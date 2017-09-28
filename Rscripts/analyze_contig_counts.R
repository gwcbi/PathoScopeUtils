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

plot_contig_scatter <- function(cts.all) {
    cts.all %>%
        ggplot2::ggplot(aes(len, mapped, color=sample, group=sample)) +
        geom_point(alpha=0.6) +
        geom_smooth(method="rlm", aes(color=sample), se=F, 
                    alpha=0.4, linetype=2, size=0.5) +
        xlab("Contig length") + ylab("# mapped") +
        ggtitle("Mapped v. length")
}

model_contig_counts <- function(cts.all, refs,
                                weight.cutoff=1e-2,
                                pct.cutoff=0
){
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
    
    # Fit OLS regression
    # models.lm <- lapply(snames, function(sn){
    #     f <- as.formula(paste0('`', sn, '` ~ len'))
    #     lm(f , cts.sp)
    # })
    # 
    # Detect outliers using Cook's distance
    # cutoff.cd <- 4 / nrow(cts.sp)
    # outliers.cd <- lapply(models.lm, function(fit){
    #     sign((cooks.distance(fit) > cutoff.cd) * fit$residuals)
    # }) %>% do.call(cbind, .)

    # Fit robust regression
    models <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        rlm(f , cts.sp, maxit=100)
    })
    # models.w <- lapply(models, function(fit) fit$w) %>% do.call(cbind, .)
    
    # Detect outliers for each sample
    outliers <- lapply(models, function(fit){
        # sign((fit$w < quantile(fit$w, 0.005)) * fit$residuals)
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
    for(f in infiles) {
        refs <- get_references(f)
        if(nrow(refs)>0) break
    }
    
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
        write.table(ret[["params"]], file.path(outdir, 'out.model_contigs.txt'),
                    quote=F, row.names=T, sep='\t')
    
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
} else {
    source('util.R')
    otu.dir <- '~/Projects/tmp/otu_analysis/Georgenia_sp_SUBG003'

    
    infiles <- Sys.glob(file.path(otu.dir, "*.bam"))
    infiles <- infiles[order(infiles)]
    names(infiles) <- gsub('.bam$','',basename(infiles))
    outdir <- dirname(infiles[1])
    
    for(f in infiles) {
        refs <- get_references(f)
        if(nrow(refs)>0) break
    }
    cts.all <- contig_counts(infiles, refs)

    weight.cutoff <- 1e-2
    pct.cutoff <- 1 / length(infiles)
    
    # Make spread data frame
    cts.sp <- cts.all %>%
        dplyr::select(chrom, len, display, mapped, sample) %>%
        tidyr::spread(sample, mapped)
    
    snames <- names(cts.sp)[4:ncol(cts.sp)]
    
    # Fit OLS regression
    models.lm <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        lm(f , cts.sp)
    })
    
    # Fit robust regression
    models.rlm <- lapply(snames, function(sn){
        f <- as.formula(paste0('`', sn, '` ~ len'))
        rlm(f , cts.sp, maxit=100)
    })
    
    utoff <- 4 / (nrow(cts.sp))
    # Detect outliers for each sample
    outliers.cd <- lapply(models.lm, function(fit){
        cooks.distance(fit) > cd.cutoff
        # sign((fit$w < weight.cutoff) * fit$residuals)
    }) %>% do.call(cbind, .)
    cts.sp[rowSums(outliers.cd) > nsamp.cutoff,]

    # Detect outliers for each sample
    outliers <- lapply(models, function(fit){
        sign((fit$w < weight.cutoff) * fit$residuals)
    }) %>% do.call(cbind, .)

    models.w[rowSums(outliers==1) > nsamp.cutoff,]
    # Possible outliers
    nsamp.cutoff <- pct.cutoff * length(snames)
    out.high <- cts.sp[rowSums(outliers==1) > nsamp.cutoff,]
    out.low <- cts.sp[rowSums(outliers==-1) > nsamp.cutoff,]

    sapply(models, function(fit) fit$w[which(rowSums(outliers==1) > nsamp.cutoff)])
    outliers[rowSums(outliers==1) > nsamp.cutoff,]
}
