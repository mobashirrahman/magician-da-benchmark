# Wilcoxon rank test on the CLR the shared contract already produced.
# The contract owns the zero policy and the CLR geometry; this adapter only tests.
adapter_wilcoxon_clr <- function(entry, y, meta, opts) {
    control <- meta$group == "Control"
    p <- apply(y, 1, function(z) {
        if (length(unique(z)) < 2) return(1)
        suppressWarnings(stats::wilcox.test(z[!control], z[control], exact = FALSE)$p.value)
    })
    # Difference of group means in the coordinates the contract centred.
    lfc <- rowMeans(y[, !control, drop = FALSE]) - rowMeans(y[, control, drop = FALSE])
    list(log2fc = lfc, pvalue = p, qvalue = NULL)
}