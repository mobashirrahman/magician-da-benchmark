# LOCOM: compositional logistic contrasts, median-centred by construction.
#
# locom() takes samples-by-features counts and a binary outcome. Filtering is owned
# by the shared contract, so filter.thresh is zero here; an additional internal
# removal would silently shrink the multiple-testing family. Taxon results and the
# community-wide test are different hypotheses: only the taxon table is used and
# p.global is recorded separately. The package does not declare a log base for its
# coefficients, so they are reported as returned and never rescaled into log2 fold
# change: fold-change error stays unavailable for this method.
adapter_locom <- function(entry, y, meta, opts) {
    data <- as.data.frame(t(y))
    rownames(data) <- meta$sample_id
    outcome <- as.integer(meta$group == "Treatment")
    perm <- as.integer(entry$parameters$perm_no %||% 999)
    if (perm < 999) MAG$fail("LOCOM n.perm.max must be at least 999 for a 0.001 resolution")
    fit <- LOCOM::locom(otu.table = data, Y = outcome,
                        filter.thresh = entry$parameters$filter_thresh %||% 0,
                        seed = opts$seed, n.perm.max = perm, n.cores = 1L,
                        permute = TRUE, verbose = FALSE)
    raw <- as.numeric(fit$p.otu[1, ])
    adjusted <- as.numeric(fit$q.otu[1, ])
    if (is.null(raw) || is.null(adjusted))
        MAG$fail("LOCOM did not return taxon p-values and q-values")
    index <- MAG$row_match(colnames(fit$p.otu), rownames(y), "LOCOM taxon results")
    effects <- as.numeric(fit$effect.size[1, ])
    # The community-wide test is a different hypothesis and is only recorded.
    if (!is.null(fit$p.global))
        message("locom p.global (community-wide, not a taxon test): ", format(fit$p.global))
    list(log2fc = effects[index], pvalue = raw[index], qvalue = adjusted[index])
}