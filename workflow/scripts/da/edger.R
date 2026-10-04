# edgeR quasi-likelihood pipeline on raw counts.
# edgeR reports no adjusted values, so the contract applies the declared BH itself.
adapter_edger <- function(entry, y, meta, opts) {
    data <- edgeR::DGEList(y, group = meta$group)
    data <- edgeR::calcNormFactors(data, method = entry$parameters$normalization %||% "TMM")
    design <- stats::model.matrix(~group, meta)
    data <- edgeR::estimateDisp(data, design)
    fit <- edgeR::glmQLFit(data, design)
    table <- edgeR::glmQLFTest(fit, coef = 2)$table
    list(log2fc = table$logFC, pvalue = table$PValue, qvalue = NULL)
}