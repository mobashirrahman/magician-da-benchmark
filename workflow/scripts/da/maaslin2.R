# MaAsLin 2 on the CLR the shared contract produced.
# No normalisation or transform: the contract already applied both, and MaAsLin2's
# own corrections would change the hypothesis. Abundance rows only. Only arguments
# the installed release declares are passed.
adapter_maaslin2 <- function(entry, y, meta, opts) {
    scratch <- file.path(opts$scratch, "maaslin2")
    dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
    declared <- names(formals(Maaslin2::Maaslin2))
    arguments <- list(input_data = t(y), input_metadata = meta, output = scratch,
                      fixed_effects = "group", reference = "group,Control",
                      normalization = entry$parameters$normalization %||% "NONE",
                      transform = entry$parameters$transform %||% "NONE",
                      standardize = isTRUE(entry$parameters$standardize),
                      min_abundance = 0, min_prevalence = 0,
                      plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
    if ("save_models" %in% declared) arguments$save_models <- FALSE
    if ("correction" %in% declared) arguments$correction <- entry$adjustment_family
    arguments <- arguments[names(arguments) %in% declared]
    fit <- tryCatch(do.call(Maaslin2::Maaslin2, arguments)$results,
                    finally = unlink(scratch, recursive = TRUE))
    fit <- fit[fit$metadata == "group", , drop = FALSE]
    fit <- fit[match(rownames(y), fit$feature), , drop = FALSE]
    list(log2fc = fit$coef, pvalue = fit$pval, qvalue = fit$qval)
}