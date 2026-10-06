# Plain two-sample t-test on relative abundance: the naive baseline.
# No compositional correction, no zero model. Any method must beat this.
adapter_ttest_relative <- function(entry, y, meta, opts) {
    control <- meta$group == "Control"
    p <- apply(y, 1, function(z) {
        if (length(unique(z)) < 2) return(1)
        suppressWarnings(stats::t.test(z[!control], z[control])$p.value)
    })
    m_c <- rowMeans(y[, control, drop = FALSE])
    m_t <- rowMeans(y[, !control, drop = FALSE])
    lfc <- ifelse(m_c > 0 & m_t > 0, log2(m_t / m_c), NA_real_)
    list(log2fc = unname(lfc), pvalue = unname(p), qvalue = NULL)
}
