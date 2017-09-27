suppressWarnings(
    suppressMessages({
        require(magrittr)
        require(ggplot2)
        require(gridExtra)
    })
)

#' Get references from BAM file
#' @param bf Path to BAM file
#' @return Data frame with reference names, lengths, and display names
get_references <- function(bf,
                           display.fun=function(.x) gsub('^.*ref\\|(.+)\\|.*$', '\\1', .x)
                           ) {
    h <- readLines(con <- pipe(paste0('samtools view -H ', bf)))
    close(con)
    h <- h[grep('^@SQ', h)]
    ret <- do.call(rbind, lapply(strsplit(h,'\t'), function(s){
        s <- unlist(s)
        c(substr(s[grep('^SN', s)], 4, 100000),
          as.numeric(substr(s[grep('^LN', s)], 4, 100000))
        )
    }))
    ret <- data.frame(chrom=ret[,1], len=as.numeric(ret[,2]), stringsAsFactors=F)
    ret$display <- display.fun(ret$chrom)
    ret[order(-ret$len),]
}

#' Get contig counts from BAM files
#' 
#' Calls samtools idxstats to get the number of reads mapping to each contig.
#' 
#' @param infiles Named vector of BAM files.
#' @param refs Reference table from get_references
#' @return Data frame with number of reads mapped to all contigs
contig_counts <- function(infiles, refs) {
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