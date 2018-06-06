#!/usr/bin/env Rscript

if(!interactive()) {
    args <- commandArgs(trailingOnly=TRUE)
    pct <- as.numeric(args[1])
    reports <- args[2:length(args)]
    top_otus <- unique(do.call(c, lapply(reports, function(f) {
        report <- read.table(f, stringsAsFactors=F, header=T, sep='\t', skip=1)
        tot <- sum(report$Final.Best.Hit.Read.Numbers)
        top <- report[report$Final.Best.Hit.Read.Numbers > (tot * pct),]
        gsub('ti\\|(\\d+).*', '\\1', top$Genome, perl=T)
    })))
    cat(top_otus[order(top_otus)], sep='\n')
}
