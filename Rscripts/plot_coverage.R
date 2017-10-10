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
        geom_area(alpha=.6) + 
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
}

if(!interactive()) {
    args <- commandArgs(trailingOnly=FALSE)
    script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))
    source(file.path(script.dir, 'util.R'))
    
    infiles <- args[(grep('--args', args)+1):length(args)]
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    # Get references information
    refs <- get_references(infiles)
    
    cov.df <- samtools_depth(infiles)
    
    plots.contig <- plot_coverage_by_contig(cov.df, refs)
    ggsave(file.path(outdir, 'out.coverage_contig.pdf'), plots.contig,
           width=11, height=8.5, paper='USr')
    
    gname <- basename(dirname(infiles[1]))
    plots.genome <- plot_coverage_by_genome(cov.df, refs, gname)
    ggsave(file.path(outdir, 'out.coverage_genome.pdf'), plots.genome,
           width=11, height=3.5, paper='USr')
    
} else {
    source('util.R')
    otu.dir <- '~/Projects/cf_hahn/otu_analysis/Pseudomonas_aeruginosa_3579'
    
    infiles <- Sys.glob(file.path(otu.dir, "*.bam"))
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    # Get references information
    refs <- get_references(infiles)

    cov.df <- samtools_depth(infiles)

}
