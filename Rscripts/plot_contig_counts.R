#! /usr/bin/env Rscript

suppressWarnings(
    suppressMessages({
        require(magrittr)
        require(dplyr)
        require(tidyr)
        require(ggplot2)
        require(MASS)
    })
)


plot_contig_counts <- function(infiles){
    snames <- names(infiles)
    if(is.null(snames)) snames <- paste0('S', 1:length(infiles))
    names(infiles) <- snames
    
    refs <- get_references(infiles[1])
    cts.all <- contig_counts(infiles, refs)
    
    mytheme <- theme_bw() + 
        theme(
            axis.text.x = element_text(angle=0, size=4),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    
    xax <- scale_x_discrete(factor(refs$display),
                            breaks=factor(refs$display)[seq(1,nrow(refs),length.out=10)],
                            labels=factor(refs$display)[seq(1,nrow(refs),length.out=10)]
    )
    
    # Plot contig lengths
    p1 <- refs %>%
        ggplot(aes(display, len / 1e6)) +
        geom_col() +
        mytheme +
        xax +
        xlab("Contig name") + ylab("Contig length (MB)") +
        ggtitle("Contig length (MB)")
    
    # Plot mapped reads
    p2 <- cts.all %>%
        ggplot(aes(display, mapped, fill=sample)) +
        geom_col(position="dodge") +
        mytheme +
        xax +
        xlab("Contig name") + ylab("Mapped reads") +
        ggtitle("Mapped reads")
    
    # Plot mapped reads per MB
    p3 <- cts.all %>%
        ggplot(aes(display, mapped / (len / 1e6), fill=sample)) +
        geom_col(position="dodge") +
        mytheme +
        xax +
        xlab("Contig name") + ylab("Mapped per MB") +
        ggtitle("Mapped per MB")
    
    list(
        contiglengths.plot=p1,
        nummapped.plot=p2,
        normmapped.plot=p3
    )
}

if(!interactive()) {
    args <- commandArgs(trailingOnly=FALSE)
    script.dir <- dirname(sub('--file=', '', args[grep('--file=', args)]))
    source(file.path(script.dir, 'util.R'))
    
    infiles <- args[(grep('--args', args)+1):length(args)]
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    ret <- plot_contig_counts(infiles)

    plot.list <- lapply(ret, function(x) {
        x + theme(axis.title.x = element_blank())
    })
    
    pg <- marrangeGrob_sharedlegend(plot.list, nrow=3, ncol=1, top=NULL, 
                                    legend=2, position="bottom")
    
    suppressMessages({    
        ggsave(file.path(outdir, 'out.contigs.pdf'), pg,
               width=8.5, height=11, paper='letter')
    })

}
