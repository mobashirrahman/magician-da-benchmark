# ZicoSeq: reference-normalised permutation tests on counts (GUniFrac 1.9 API).
#
# 1.9 takes feature.dat as a features-by-samples matrix, not a data frame, and
# returns per-feature F statistics, R2 and raw/adjusted p-values. The F statistic is
# on transformed ratios, not a genomic log2 fold change, so the registry declares
# effect_scale=f_statistic and no fold-change error is ever computed for this
# method. Both raw and native adjusted p-values are retained. ZicoSeq's own
# reference set is reported, not hidden.
adapter_zicoseq <- function(entry, y, meta, opts) {
    data <- as.matrix(y)
    perm <- as.integer(entry$parameters$perm_no %||% 999)
    if (perm < 999) MAG$fail("ZicoSeq perm.no must be at least 999 for a 0.001 resolution")
    fit <- GUniFrac::ZicoSeq(meta.dat = as.data.frame(meta), feature.dat = data,
                             grp.name = "group", adj.name = NULL,
                             feature.dat.type = "count",
                             prev.filter = entry$parameters$prev_filter %||% 0,
                             mean.abund.filter = 0, max.abund.filter = 0,
                             is.winsor = TRUE, outlier.pct = 0.03,
                             is.post.sample = TRUE,
                             post.sample.no = as.integer(entry$parameters$post_sample_no %||% 25),
                             perm.no = perm,
                             ref.pct = entry$parameters$ref_pct %||% 0.50,
                             return.feature.dat = FALSE, verbose = FALSE)
    raw <- fit$p.raw
    adjusted <- fit$p.adj.fdr
    if (is.null(raw) || is.null(adjusted))
        MAG$fail("ZicoSeq did not expose both p.raw and p.adj.fdr")
    statistic <- fit$F0
    lfc <- rep(NA_real_, nrow(y))
    if (!is.null(statistic)) {
        statistic <- as.matrix(statistic)
        # Mirror the default max combination across link functions.
        combined <- apply(statistic, 1, function(row) suppressWarnings(max(row, na.rm = TRUE)))
        index <- MAG$row_match(rownames(statistic), rownames(y), "ZicoSeq F statistics")
        lfc <- as.numeric(combined)[index]
    }
    index <- MAG$row_match(names(raw), rownames(y), "ZicoSeq results")
    p <- as.numeric(raw)[index]
    q <- as.numeric(adjusted)[index]
    # Features ZicoSeq filtered itself are exclusions of the method, not of the truth.
    if (!is.null(fit$filter.features) && length(fit$filter.features)) {
        dropped <- match(fit$filter.features, rownames(y))
        dropped <- dropped[!is.na(dropped)]
        p[dropped] <- NA_real_
        q[dropped] <- NA_real_
        lfc[dropped] <- NA_real_
        message("ZicoSeq filtered ", length(dropped), " features internally: ",
                paste(rownames(y)[dropped], collapse = ", "))
    }
    if (!is.null(fit$ref.features))
        message("ZicoSeq reference set (", length(fit$ref.features), " features): ",
                paste(head(fit$ref.features, 10), collapse = ", "),
                if (length(fit$ref.features) > 10) ", ..." else "")
    list(log2fc = lfc, pvalue = p, qvalue = q)
}