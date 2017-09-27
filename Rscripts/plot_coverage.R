#! /usr/bin/env Rscript

suppressWarnings(
    suppressMessages({
        require(magrittr)
        require(dplyr)
        require(tidyr)
        require(ggplot2)
    })
)

makebins.chrom <- function(cov, reflen, nbins=min(500, reflen-1)) {
    snames <- names(cov)[3:ncol(cov)]
    bks <- as.integer(floor(seq(1, (reflen+1), length.out=nbins+1)))
    ret <- cov %>%
        mutate(bin=cut(pos, bks, 1:nbins, right=F)) %>%
        mutate(
            start=bks[1:nbins][bin],
            end=bks[2:(nbins+1)][bin]
        )

    # Include rows for missing bins
    missbins <- c(1:nbins)[! 1:nbins %in% ret$bin]
    if(length(missbins) > 0) {
        miss.df <- data.frame(chrom=ret$chrom[1],
                              bin=factor(missbins, levels=1:nbins),
                              pos=bks[1:nbins][missbins],
                              start=bks[1:nbins][missbins],
                              end=bks[2:(nbins+1)][missbins],
                              stringsAsFactors = F
        )
        for(i in 1:length(snames)) miss.df[[ snames[i] ]] <- 0
        miss.df <- miss.df[,names(ret)]
        ret <- rbind(ret, miss.df)
    }
    stopifnot(all(ret$start <= ret$pos & ret$pos < ret$end))
    stopifnot(all(c(1:nbins) %in% ret$bin))
    
    ret %>%
        group_by_at(vars(chrom,bin,start,end)) %>%
        summarise_at(snames, sum) %>%
        ungroup() %>%
        dplyr::select(chrom, start, end, snames)
}


plotbins.chrom <- function(cov.bins) {
    snames <- names(cov.bins)[4:ncol(cov.bins)]
    chromS <- min(cov.bins$start)
    chromE <- max(cov.bins$end)
    
    prop.bins <- cov.bins %>%
        mutate_at(snames, funs(./(end-start)))
    
    prop.bins %>%
        mutate(pos=(start+end)/2) %>%
        dplyr::select(chrom, pos, snames) %>%
        tidyr::gather(key='sample', value='cov', -c(chrom, pos)) %>%
        ggplot(aes(pos, cov, color=sample, fill=sample)) +
        geom_area(alpha=.4) + 
        xlim(chromS, chromE) + theme_bw()
}

makebins.genome <- function(cov, refs, nbins=500) {
    offsets <- c(0, cumsum(refs$len)[1:(nrow(refs)-1)])
    names(offsets) <- refs$chrom
    cov.merged <- cov %>%
        mutate(pos=pos+offsets[chrom]) %>%
        mutate(chrom="genome")
    
    nbins <- min(nbins, sum(refs$len))
    gbins <- makebins.chrom(cov.merged, sum(refs$len), nbins)
    gbins
}

plot_coverage_by_genome <- function(cov.df, refs, gname="Genome") {
    # Plot full genome as merged
    gbins <- makebins.genome(cov.df, refs)
    
    vlines <- cumsum(refs$len)[1:(nrow(refs)-1)]
    bounds <- c(0, vlines, sum(refs$len))
    mids <- sapply(1:(length(bounds)-1), function(i){(bounds[i] + bounds[i+1]) / 2})
    
    # Make the plots
    p <- plotbins.chrom(gbins)
    if(nrow(refs) < 15) {
        p <- p + geom_vline(xintercept=vlines, lty=3, alpha=0.7) +
                 annotate("text", x=mids, y=(max(p$data$cov)*0.9),
                          label=refs$display, size=2.5, alpha=0.7)
    }
    p + ggtitle(paste0(gname))
}

plot_coverage_by_contig <- function(cov.df, refs) {
    # Plot each contig in its own plot
    cbins.list <- lapply(1:nrow(refs), function(r) {
        makebins.chrom(cov.df[cov.df$chrom==refs[r,]$chrom, ], refs[r,]$len)
    })
    
    # Make the plots
    cbins.plots <- lapply(cbins.list, plotbins.chrom)
    # Add titles
    cbins.plots <- lapply(1:nrow(refs), function(r) {
        cbins.plots[[r]] +
            ggtitle(paste0("Contig: ", refs[r,]$display))
    })
    
    marrangeGrob_sharedlegend(cbins.plots, legend=1, nrow=2, ncol=1, top=NULL) 
    # Plot concatenated contigs
}

get_coverage <- function(infiles) {
    snames <- names(infiles)
    if(is.null(snames)) snames <- paste0('S', 1:length(infiles))
    names(infiles) <- snames
    
    cmd <- paste0('samtools depth ', paste(infiles, collapse=' '))
    cov.df <- read.table(pipe(cmd), sep='\t', stringsAsFactors=F)
    names(cov.df) <- c('chrom', 'pos', names(infiles)) 
    cov.df
}

if(!interactive()) {
    args <- commandArgs(trailingOnly=FALSE)
    script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))
    source(file.path(script.dir, 'util.R'))
    
    infiles <- args[(grep('--args', args)+1):length(args)]
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    refs <- get_references(infiles[1])
    cov.df <- get_coverage(infiles)
    
    plots.contig <- plot_coverage_by_contig(cov.df, refs)
    ggsave(file.path(outdir, 'out.coverage_contig.pdf'), plots.contig,
           width=11, height=8.5, paper='USr')
    
    gname <- basename(dirname(infiles[1]))
    plots.genome <- plot_coverage_by_genome(cov.df, refs, gname)
    ggsave(file.path(outdir, 'out.coverage_genome.pdf'), plots.genome,
           width=11, height=3.5, paper='USr')
    
} else {
    source('util.R')
    infiles <- Sys.glob("~/Projects/tmp/otu_analysis/Alcanivorax_hongdengensis_A_11_3/fixed/*.bam")
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    for(f in infiles) {
        refs <- get_references(f)
        if(nrow(refs)>0) break
    }
    cov.df <- get_coverage(infiles)

    plots.contig <- plot_coverage_by_contig(cov.df, refs)
    plots.contig
    
    plots.genome <- plot_coverage_by_genome(cov.df, refs)
    plots.genome
}


# infiles <- Sys.glob('otu_analysis/Burkholderia_cenocepacia_H111/*.pri.bam')
# infiles <- infiles[order(infiles)]
# snames <- gsub('.bam', '', basename(infiles))
# names(infiles) <- snames
# 
# outdir <- dirname(infiles[1])
# 
# # Read reference information from SAM header
# refs <- get_references(infiles[1])
# 
# 
# 
# 
# 
# 
# # pgrob <- gridExtra::marrangeGrob(cov.plots, nrow=3, ncol=2)
# pgrob <- marrangeGrob_sharedlegend(cbins.plots, nrow=2, ncol=1)
# ggsave(file.path(outdir, 'out.coverage.pdf'), width=11, height=8.5, pgrob)
# 
# 
# ### ggsave(file.path(outdir, 'out.coverage.pdf'), width=11, height=8.5, pgrob)
# 
# 
# 
# 
# 
# marrage_grob_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
#     
#     plots <- list(...)
#     position <- match.arg(position)
#     g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
#     legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#     lheight <- sum(legend$height)
#     lwidth <- sum(legend$height)
#     gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
#     gl <- c(gl, nrow = nrow, ncol = ncol)
#     
#     perpage <- nrow * ncol
#     lapply()
#     combined <- switch(position,
#                        "bottom" = marrangeGrob(do.call(marrangeGrob, gl),
#                                               legend,
#                                               ncol = 1,
#                                               heights = unit.c(unit(1, "npc") - lheight, lheight)),
#                        "right" = arrangeGrob(do.call(arrangeGrob, gl),
#                                              legend,
#                                              ncol = 2,
#                                              widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
#     grid.newpage()
#     grid.draw(combined)
#     # combined
# }
# 
# marrage_grob_shared_legend()
# 
# 
# pgrob <- gridExtra::marrangeGrob(cov.plots, nrow=3, ncol=2)
# ggsave(file.path(outdir, 'out.coverage.pdf'), width=11, height=8.5, pgrob)
# 
# 
# 
# z <- grid_arrange_shared_legend(cbins.plots[[1]], cbins.plots[[2]], cbins.plots[[3]], cbins.plots[[1]],
#                                 cbins.plots[[1]], cbins.plots[[2]], cbins.plots[[3]], cbins.plots[[1]], nrow=3, ncol=2)
#         
# 
