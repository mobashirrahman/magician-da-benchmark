# Plain Wilcoxon on relative abundance: the naive rank baseline.
adapter_wilcoxon_relative <- function(entry, y, meta, opts) {
    control <- meta$group == "Control"
    p <- apply(y, 1, function(z) {
        if (length(unique(z)) < 2) return(1)
        suppressWarnings(stats::wilcox.test(z[!control], z[control], exact = FALSE)$p.value)
    })
    m_c <- rowMeans(y[, control, drop = FALSE])
    m_t <- rowMeans(y[, !control, drop = FALSE])
    lfc <- ifelse(m_c > 0 & m_t > 0, log2(m_t / m_c), NA_real_)
    list(log2fc = unname(lfc), pvalue = unname(p), qvalue = NULL)
}
