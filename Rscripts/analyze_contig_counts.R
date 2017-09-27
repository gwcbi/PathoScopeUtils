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
    
    # Fit robust regression
    models <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        rlm(f , cts.sp)
    })
    
    # Detect outliers for each sample
    outliers <- lapply(models, function(fit){
        sign((fit$w < weight.cutoff) * fit$residuals)
    }) %>% do.call(cbind, .)

    # Possible outliers
    nsamp.cutoff <- pct.cutoff * length(snames)
    out.high <- cts.sp[rowSums(outliers==1) > nsamp.cutoff,]
    out.low <- cts.sp[rowSums(outliers==-1) > nsamp.cutoff,]
    
    # Model params
    params <- lapply(models, function(fit){
        predicted <- predict(fit, data.frame(len=sum(refs$len)), interval='confidence')
        c(fit$coefficients, predicted)
    }) %>% do.call(rbind, .) %>%
    data.frame(., row.names=snames)
    names(params) <- c('(Intercept)', 'len', 'fit', 'lwr', 'upr')

    # Outliers removed
    cts.noout <- cts.sp %>%
        dplyr::filter(!chrom %in% out.high$chrom)
    params$out.rm <- sapply(snames, function(n) sum(cts.noout[,n]))

    list(
        out.high=out.high,
        out.low=out.low,
        scatter.plot=p1,
        params=params
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
    write.table(ret[["params"]], file.path(outdir, 'out.model_contigs.txt'),
                quote=F, row.names=T, sep='\t')
} else {
    source('util.R')
    infiles <- Sys.glob("~/Projects/tmp/otu_analysis/Alcanivorax_hongdengensis_A_11_3/fixed/*.bam")
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    refs <- get_references(infiles[1])
    
}
