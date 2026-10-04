# limma-voom on raw counts. limma reports no adjusted values; the contract applies BH.
adapter_limma_voom <- function(entry, y, meta, opts) {
    data <- edgeR::DGEList(y)
    data <- edgeR::calcNormFactors(data, method = entry$parameters$normalization %||% "TMM")
    design <- stats::model.matrix(~group, meta)
    voom <- limma::voom(data, design, plot = FALSE)
    fit <- limma::eBayes(limma::lmFit(voom, design))
    list(log2fc = fit$coefficients[, 2], pvalue = fit$p.value[, 2], qvalue = NULL)
}