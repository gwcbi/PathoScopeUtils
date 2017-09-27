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


analyze_contig_counts <- function(infiles,
                                  weight.cutoff=1e-2,
                                  pct.cutoff=0
){
    snames <- names(infiles)
    if(is.null(snames)) snames <- paste0('S', 1:length(infiles))
    names(infiles) <- snames
    
    refs <- get_references(infiles[1])
    cts.all <- contig_counts(infiles, refs)
    
    p1 <- cts.all %>%
        ggplot2::ggplot(aes(len, mapped, color=sample, group=sample)) +
        geom_point(alpha=0.6) +
        geom_smooth(method="rlm", aes(color=sample), se=F, 
                    alpha=0.4, linetype=2, size=0.5) +
        xlab("Contig length") + ylab("# mapped") +
        ggtitle("Mapped v. length")
    
    # Make spread data frame
    cts.sp <- cts.all %>%
        dplyr::select(chrom, len, display, mapped, sample) %>%
        tidyr::spread(sample, mapped)
    
    # Detect outliers for each sample
    outliers <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        fit.summary <- summary(fit <- rlm(f , cts.sp))
        sign((fit$w < weight.cutoff) * fit$residuals)
    }) %>% do.call(cbind, .)
    
    # Possible outliers
    nsamp.cutoff <- pct.cutoff * length(snames)
    out.high <- cts.sp[rowSums(outliers==1) > nsamp.cutoff,]
    out.low <- cts.sp[rowSums(outliers==-1) > nsamp.cutoff,]
    
    list(
        out.high=out.high,
        out.low=out.low,
        scatter.plot=p1
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
    
    ret <- analyze_contig_counts(infiles)
    
    # Output scatter plot with regression lines
    suppressMessages(
        ggsave(file.path(outdir, 'out.scatter.pdf'), ret[['scatter.plot']],
               width=7, height=7, paper='USr')
    )
    
    # Output tables with outliers
    if(nrow(ret[["out.high"]]) == 0 & nrow(ret[["out.low"]]) == 0) {
        cat("Contig counts match with expectation 👍🏼\n")
    } else {
        if(nrow(ret[["out.high"]]) > 0) {
            cat("Contigs with counts greater than expected:\n")
            t <- ret[["out.high"]] %>%
                dplyr::select(display, len, names(infiles))
            print(t)
            write.table(t, file.path(outdir, 'out.high_contigs.txt'),
                        quote=F, row.names=F, sep='\t')
        }
        if(nrow(ret[["out.low"]]) > 0) {
            cat("Contigs with counts lower than expected:")
            t <- ret[["out.low"]] %>%
                dplyr::select(display, len, names(infiles))
            print(t)
            write.table(t, file.path(outdir, 'out.high_contigs.txt'),
                        quote=F, row.names=F, sep='\t')
        }
    }
}
