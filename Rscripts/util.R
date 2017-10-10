suppressWarnings(
    suppressMessages({
        require(magrittr)
        require(ggplot2)
        require(gridExtra)
    })
)

#' Get references from BAM file
#' @param infiles Named vector of BAM files.
#' @return Data frame with reference names, lengths, and display names
get_references <- function(infiles,
                           display.fun=function(.x) gsub('^.*ref\\|(.+)\\|.*$', '\\1', .x)
                           ) {
    for (bf in infiles) {
        h <- readLines(con <- pipe(paste0('samtools view -H ', bf)))
        close(con)        
        if(sum(grep('^@SQ', h)) == 0) next
        h <- h[grep('^@SQ', h)]
        ret <- do.call(rbind, lapply(strsplit(h,'\t'), function(s){
            s <- unlist(s)
            c(substr(s[grep('^SN', s)], 4, 100000),
              as.numeric(substr(s[grep('^LN', s)], 4, 100000))
            )
        }))
        ret <- data.frame(chrom=ret[,1], len=as.numeric(ret[,2]), stringsAsFactors=F)
        ret$display <- display.fun(ret$chrom)
        ret <- ret[order(-ret$len),]
        if(length(unique(ret$display)) == length(ret$display)) {
            ret$display <- factor(ret$display, levels=ret$display)
        } else {
            ret$display <- paste(ret$display, 1:nrow(ret), sep='.')
            ret$display <- factor(ret$display, levels=ret$display)
        }
        if(nrow(ret) > 0) return(ret)
    }
    return(data.frame(chrom="", len=0, display=""))
}


#' Get contig counts from BAM files
#' 
#' Calls samtools idxstats to get the number of reads mapping to each contig.
#' 
#' @param infiles Named vector of BAM files.
#' @param refs Reference table from get_references
#' @return Data frame with number of reads mapped to all contigs
contig_counts <- function(infiles, refs=get_references(infiles)) {
    snames <- names(infiles)
    if(is.null(snames)) snames <- paste0('S', 1:length(infiles))
    names(infiles) <- snames
    
    lapply(names(infiles), function(n) {
        bf <- infiles[n]
        df <- read.table(pipe(paste0('samtools idxstats ', bf)),
                         sep='\t', stringsAsFactors=F)
        names(df) <- c('chrom', 'len', 'mapped', 'unmapped')
        df <- merge(refs, df, by=c('chrom','len'), all.x=T)
        df[is.na(df)] <- 0
        df$sample <- n
        df
    }) %>% do.call(rbind, .)
}

#' Get coverage depth from BAM files
#' 
#' Calls samtools depth to get the number of reads mapping to each position.
#' 
#' @param infiles Named vector of BAM files.
#' @return Data frame with chromosome, position, and depth in each sample.
samtools_depth <- function(infiles) {
    # samtools depth will return nothing if there are no reads in the first 
    # file. Order the input BAMs so that the first one has reads
    ccts <- contig_counts(infiles, refs <- get_references(infiles)) %>%
        dplyr::select(chrom, len, display, mapped, sample) %>%
        tidyr::spread(sample, mapped) %>%
        dplyr::select(names(infiles)) %>% colSums
    if(max(ccts) == 0) {
        warning("No samples have counts.")
        return(data.frame(chrom="", pos=1))
    }
    neworder <- infiles[order(-ccts)]
    cmd <- paste0('samtools depth ', paste(neworder, collapse=' '))
    cov.df <- read.table(con <- pipe(cmd), sep='\t', stringsAsFactors=F)
    names(cov.df) <- c('chrom', 'pos', names(neworder))
    
    # Put back in infile order
    cov.df <- cov.df[, c('chrom', 'pos', names(infiles))]
    cov.df
}


#' Share legends
marrangeGrob_sharedlegend <- function (grobs, ncol, nrow, ..., 
                                       top = quote(paste("page", g, "of", pages)),
                                       position = c("bottom", "right"),
                                       legend = 1
) 
{
    # Get the legend from the first plot
    position <- match.arg(position)
    
    if(is.numeric(legend) | (is.character(legend) & legend %in% names(grobs))) {
        g <- ggplotGrob(grobs[[legend]] + theme(legend.position = position))$grobs
        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    }
    
    switch(position,
           "bottom" = gridExtra::marrangeGrob(
               lapply(grobs, function(x) x + theme(legend.position = "none")),
               nrow = nrow, ncol = ncol, top=top, bottom=legend),
           
           "right" = gridExtra::marrangeGrob(
               lapply(grobs, function(x) x + theme(legend.position = "none")),
               nrow = nrow, ncol = ncol, top=top, right=legend)
    )
}


#' Share legends
# marrangeGrob_sharedlegend2 <- function (grobs, ncol, nrow, ..., 
#                                        top = quote(paste("page", g, "of", pages)),
#                                        position = c("bottom", "right")
# ) 
# {
#     # Get the legend from the first plot
#     position <- match.arg(position)
#     g <- ggplotGrob(grobs[[1]] + theme(legend.position = position))$grobs
#     legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#     lheight <- sum(legend$height)
#     lwidth <- sum(legend$height)    
#     # Remove the legend from all the plots
#     grobs <- lapply(grobs, function(x) x + theme(legend.position = "none"))
#     # Assign plots to pages
#     n <- length(grobs)
#     nlay <- nrow * ncol
#     pages <- n%/%nlay + as.logical(n%%nlay)
#     groups <- split(seq_along(grobs), gl(pages, nlay, n))
#     pl <- vector(mode = "list", length = pages)
#     for (page in seq_along(groups)) {
#         g <- page
#         params <- modifyList(list(...), list(top = eval(top), 
#                                              nrow = nrow, ncol = ncol))
#         
#         # Lay out the plots
#         pl[[g]] <- do.call(arrangeGrob, c(grobs[groups[[g]]], 
#                                           params))
#         # Add legend to plots
#         pl[[g]] <- arrangeGrob(pl[[g]],
#                                legend,
#                                ncol=1,
#                                heights = unit.c(unit(1, "npc") - lheight, lheight))
#     }
#     class(pl) <- c("arrangelist", class(pl))
#     pl
# }