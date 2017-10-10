#! /usr/bin/env Rscript

suppressWarnings(
    suppressMessages({
        require(magrittr)
        require(dplyr)
        require(tidyr)
        require(tibble)
        require(ggplot2)
        require(MASS)
        require(robust)
    })
)

plot_contig_scatter <- function(cts.all) {
    cts.all %>%
        ggplot2::ggplot(aes(len, mapped, color=sample, group=sample)) +
        geom_point(alpha=0.6) +
        geom_smooth(method="rlm", aes(color=sample), se=F, 
                    alpha=0.4, linetype=2, size=0.5) +
        xlab("Contig length") + ylab("# mapped") +
        ggtitle("Mapped v. length")
}

fit_lm <- function(cts.sp, snames) {
    models <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        suppressWarnings(lm(f, cts.sp))
    })
    names(models) <- snames
    models
}

fit_rlm <- function(cts.sp, snames) {
    models <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        suppressWarnings(MASS::rlm(f , cts.sp, maxit=100))
    })
    names(models) <- snames
    models
}

fit_glmrob <- function(cts.sp, snames) {
    models <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        robust::glmRob(f, family='poisson', data=cts.sp)
    })
    names(models) <- snames
    models
}

outliers_robust <- function(models, weight.cutoff) {
    # Detect outliers from robust regression using weights
    val <- 'weights'
    if("rlm" %in% class(models[[1]])) val <- 'w'
    
    # Detect outliers for each sample    
    outliers <- lapply(models, function(fit){
        sign((fit[[val]] < weight.cutoff) * fit$residuals)
    }) %>% do.call(cbind, .)
    outliers
}

outliers_cookD <- function(models, d.cutoff=NULL) {
    # Detect outliers from regression using Cook's distance
    if(is.null(d.cutoff)) d.cutoff <- (4 / length(models[[1]]$residuals))
    outliers <- lapply(models, function(fit){
        sign((cooks.distance(fit) > d.cutoff) * fit$residuals)
    }) %>% do.call(cbind, .)
    outliers
}

param_table <- function(cts.sp, snames, models, is_outlier, nsamp.cutoff) {
    d <- data.frame(
        adj.count=colSums(cts.sp[!is_outlier, snames]),
        exp.count=round(sapply(models, function(fit) sum(sapply(predict(fit), function(x) max(0,x))))),
        Intercept=sapply(models, function(fit) fit$coefficients[['(Intercept)']]),
        slope=sapply(models, function(fit) fit$coefficients[['len']])
    )
    row.names(d) <- snames
    d
}

model_contig_counts <- function(cts.all, refs,
                                weight.cutoff=1e-2,
                                pct.cutoff=0,
                                method=c("glmRob", "rlm", "lm")
){
    method <- match.arg(method)
    
    # Make spread data frame
    cts.sp <- cts.all %>%
        dplyr::select(chrom, len, display, mapped, sample) %>%
        tidyr::spread(sample, mapped)
    
    snames <- names(cts.sp)[4:ncol(cts.sp)]
    
    if(nrow(cts.sp)<=2) {
        return(list(
            out.high=data.frame(),
            out.low=data.frame(), 
            params=data.frame()
        ))
    }
    switch(method,
           lm = {
               #--- Fit a Linear Model ---#
               models <- fit_lm(cts.sp, snames)
               outliers <- outliers_cookD(models)
           }, rlm = {
               #--- Fit a linear model by robust regression ---#
               models <- fit_rlm(cts.sp, snames)
               outliers <- outliers_robust(models, weight.cutoff)
           }, glmRob = {
               #--- Fit a Robust Generalized Linear Model ---#
               models <- fit_glmrob(cts.sp, snames)
               outliers <- outliers_robust(models, weight.cutoff)
           }
    )
    
    # Identify outliers across samples
    nsamp.cutoff <- pct.cutoff * length(snames)
    is_outlier <- rowSums(outliers==1) > nsamp.cutoff # | rowSums(outliers==-1) > nsamp.cutoff
    out.tab <- cts.sp[is_outlier, c('display', 'len', snames)]
    params <- param_table(cts.sp, snames, models, is_outlier, nsamp.cutoff)

    list(
        out.tab=out.tab,
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
    
    # Get references information
    refs <- get_references(infiles)
    
    # Get contig counts
    cts.all <- contig_counts(infiles, refs)
    
    # Output scatter plot with regression lines
    suppressWarnings(suppressMessages({
        p <- plot_contig_scatter(cts.all)        
        ggsave(file.path(outdir, 'out.scatter.pdf'), p,
               width=7, height=7, paper='USr')
        
    }))
    
    # Model contig counts
    ret <- model_contig_counts(cts.all, refs, pct.cutoff=(1/length(infiles)))
    
    if(nrow(ret[['params']])==0) {
        cat(paste0("Model could not be fit. Found ", length(unique(cts.all$chrom)), " contigs.\n"))
        if(length(unique(cts.all$chrom)) <=2) {
            cat("Contig counts match with expectation 👍🏼\n")
        }
    } else {
        # Write parameter table
        ret[['params']] %>% tibble::rownames_to_column("sample") %>%
            write.table(file.path(outdir, 'out.countmodel.txt'),
                    quote=F, row.names=F, sep='\t')
        # Output tables with outliers
        if(nrow(ret[["out.tab"]]) == 0) {
            cat("Contig counts match with expectation 👍🏼\n")
        } else {
            # Print parameter table
            cat("Unusual contig coverage. Model parameters:\n")
            print(ret[["params"]])
            cat("Contigs with counts greater than expected:\n")
            print(ret[["out.tab"]])
            ret[["out.tab"]] %>% 
                write.table(file.path(outdir, 'out.outliers.txt'),
                            quote=F, row.names=F, sep='\t')
        }
    }
} else {
    source('util.R')
    otu.dir <- '~/Projects/cf_hahn/otu_analysis/Nocardia_brevicatena_NBRC_12119'

    infiles <- Sys.glob(file.path(otu.dir, "*.bam"))
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    # Get references information
    refs <- get_references(infiles)
    
    cts.all <- contig_counts(infiles, refs)

    # Model contig counts
    ret <- model_contig_counts(cts.all, refs, pct.cutoff=(1/length(infiles)))
}
