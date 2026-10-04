# DESeq2 negative-binomial fit on raw counts.
# independentFiltering is off so the reported correction family is the declared BH
# over every tested feature, and size factors come from poscounts so that a feature
# with zeros cannot collapse the whole normalisation.
adapter_deseq2 <- function(entry, y, meta, opts) {
    data <- DESeq2::DESeqDataSetFromMatrix(round(y), meta, ~group)
    data <- DESeq2::DESeq(data, sfType = entry$parameters$sfType %||% "poscounts",
                          quiet = TRUE, parallel = FALSE)
    fit <- DESeq2::results(data, contrast = c("group", "Treatment", "Control"),
                           independentFiltering = isTRUE(entry$parameters$independentFiltering),
                           alpha = opts$alpha)
    list(log2fc = fit$log2FoldChange, pvalue = fit$pvalue, qvalue = fit$padj)
}