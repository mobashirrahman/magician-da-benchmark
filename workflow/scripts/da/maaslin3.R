# MaAsLin 3 on length-corrected relative abundance.
#
# Parameter names follow the frozen 1.4.0 release: normalization, transform,
# correction, standardize, median_comparison_abundance, evaluate_only. The relative
# abundance lane sets median_comparison_abundance=FALSE and takes the individual
# abundance rows of the group coefficient. The composition-corrected lane sets it
# TRUE, which changes the hypothesis to change relative to the typical feature; it
# is scored against its own truth lane, never against relative-abundance truth.
# With transform LOG the coefficients are natural-log units and are converted to
# log2 before they leave this adapter.
adapter_maaslin3 <- function(entry, y, meta, opts) {
    scratch <- file.path(opts$scratch, "maaslin3")
    dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
    median_adjusted <- isTRUE(entry$parameters$median_comparison_abundance)
    fit <- tryCatch(
        maaslin3::maaslin3(input_data = t(y), input_metadata = meta, output = scratch,
                           fixed_effects = "group",
                           reference = entry$parameters$reference %||% "group,Control",
                           random_effects = character(0),
                           min_abundance = 0, min_prevalence = 0,
                           normalization = entry$parameters$normalization %||% "NONE",
                           transform = entry$parameters$transform %||% "LOG",
                           correction = entry$adjustment_family,
                           standardize = isTRUE(entry$parameters$standardize),
                           median_comparison_abundance = median_adjusted,
                           evaluate_only = entry$parameters$evaluate_only %||% "abundance",
                           # The prevalence warning is meaningless when only abundance is fitted.
                           warn_prevalence = FALSE,
                           plot_summary_plot = FALSE, plot_associations = FALSE,
                           max_pngs = 0, cores = 1,
                           save_models = FALSE, save_plots_rds = FALSE),
        finally = unlink(scratch, recursive = TRUE))
    # The function returns both hurdle halves; the abundance lane takes its own rows.
    abundance <- fit$fit_data_abundance
    results <- if (is.list(abundance) && !is.data.frame(abundance)) abundance$results else abundance
    fit <- list(results = results)
    if (is.null(fit$results) || !NROW(fit$results))
        MAG$fail("MaAsLin 3 returned no results")
    # Abundance rows only. Joint and prevalence rows belong to the hurdle output and
    # would silently widen the multiple-testing family of this lane.
    rows <- fit$results[fit$results$metadata == "group", , drop = FALSE]
    if (!nrow(rows)) MAG$fail("MaAsLin 3 returned no group abundance rows")
    message("MaAsLin 3 abundance rows: ", nrow(rows), "; metadata values seen: ",
            paste(unique(as.character(rows$value %||% rows$name %||% "group")), collapse = ", "))
    # A two-level factor yields one coefficient per feature; if the package reports a
    # row per level, keep the Treatment row so the contrast stays Treatment/Control.
    if ("value" %in% names(rows) && any(grepl("Treatment", rows$value))) {
        rows <- rows[grepl("Treatment", rows$value), , drop = FALSE]
    } else if ("name" %in% names(rows) && any(grepl("Treatment", rows$name))) {
        rows <- rows[grepl("Treatment", rows$name), , drop = FALSE]
    }
    index <- MAG$row_match(as.character(rows$feature), rownames(y), "MaAsLin 3 results")
    failed <- if ("error" %in% names(rows)) !is.na(rows$error) & nzchar(as.character(rows$error)) else rep(FALSE, nrow(rows))
    lfc <- as.numeric(rows$coef)[index] / log(2)
    p <- as.numeric(rows$pval_individual)[index]
    q <- as.numeric(rows$qval_individual)[index]
    dropped <- index[failed[index]]
    p[dropped] <- NA_real_
    q[dropped] <- NA_real_
    lfc[dropped] <- NA_real_
    if (any(failed)) message("MaAsLin 3 left ", sum(failed), " features without a fit")
    list(log2fc = lfc, pvalue = p, qvalue = q)
}