# LinDA: bias-corrected log ratio to a reference feature, with adaptive zero handling.
#
# Zero handling stays inside the package: the contract adds no pseudocount. The
# adaptive branch is the declared setting and its choice against the imputation
# branch is reported, not chosen on speed. Internal prevalence filtering is
# disabled, so the multiple-testing family stays explicit.
adapter_linda <- function(entry, y, meta, opts) {
    data <- as.data.frame(y)
    fit <- MicrobiomeStat::linda(feature.dat = data, meta.dat = as.data.frame(meta),
                                 formula = "~group", feature.dat.type = "count",
                                 prev.filter = entry$parameters$prev_filter %||% 0,
                                 mean.abund.filter = 0, max.abund.filter = 0,
                                 is.winsor = isTRUE(entry$parameters$winsor),
                                 adaptive = isTRUE(entry$parameters$adaptive),
                                 zero.handling = entry$parameters$zero_handling %||% "imputation",
                                 pseudo.cnt = entry$parameters$pseudo_cnt %||% 0.5,
                                 corr.cut = entry$parameters$corr_cut %||% 0.1,
                                 p.adj.method = entry$adjustment_family, alpha = opts$alpha,
                                 n.cores = 1, verbose = FALSE)
    table <- fit$output
    if (is.null(table)) MAG$fail("LinDA returned no output table")
    # One table per model variable; the group coefficient is the Treatment contrast.
    if (!is.data.frame(table)) {
        candidates <- grep("Treatment", names(table), value = TRUE)
        if (!length(candidates)) candidates <- names(table)
        table <- table[[candidates[1]]]
    }
    if (is.null(table)) MAG$fail("LinDA returned no group coefficient table")
    index <- MAG$row_match(rownames(table), rownames(y), "LinDA results")
    message("LinDA bias (reference-constrained mode shift): ", format(fit$bias))
    list(log2fc = table$log2FoldChange[index], pvalue = table$pvalue[index],
         qvalue = table$padj[index])
}