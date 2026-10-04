# ADAPT: censored zero model with a learned reference set.
#
# adapt() returns an S4 DAresult. The full table lives in the `details` slot with
# Taxa, prevalence, log10foldchange, pval and adjusted_pval; the learned reference
# set lives in the `reference` slot. Zeros are treated as censored observations
# rather than imputed, so no pseudocount is introduced here. The learned set is
# reported, not assumed correct: the truth lane uses the typical-change counterpart
# of a learned set.
adapter_adapt <- function(entry, y, meta, opts) {
    ps <- phyloseq::phyloseq(
        phyloseq::otu_table(as.data.frame(t(y)), taxa_are_rows = FALSE),
        phyloseq::sample_data(as.data.frame(meta)))
    fit <- ADAPT::adapt(ps, cond.var = "group", base.cond = "Control",
                        prev.filter = entry$parameters$prev_filter %||% 0,
                        depth.filter = entry$parameters$depth_filter %||% 0,
                        alpha = opts$alpha)
    details <- tryCatch(methods::slot(fit, "details"), error = function(e) NULL)
    if (is.null(details)) details <- as.data.frame(summary(fit, select = "all"))
    table <- as.data.frame(details)
    index <- MAG$row_match(as.character(table$Taxa), rownames(y), "ADAPT results")
    p <- rep(NA_real_, nrow(y)); q <- rep(NA_real_, nrow(y)); lfc <- rep(NA_real_, nrow(y))
    if ("log10foldchange" %in% names(table)) lfc <- table$log10foldchange[index] * log2(10)
    if ("pval" %in% names(table)) p <- table$pval[index]
    if ("adjusted_pval" %in% names(table)) q <- table$adjusted_pval[index]
    reference <- tryCatch(methods::slot(fit, "reference"), error = function(e) character(0))
    message("ADAPT learned reference set (", length(reference), " taxa): ",
            paste(head(reference, 10), collapse = ", "),
            if (length(reference) > 10) ", ..." else "")
    list(log2fc = lfc, pvalue = p, qvalue = q)
}