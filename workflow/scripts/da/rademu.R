# radEmu: contrasts relative to the typical taxon change, without a pseudocount.
#
# emuFit returns a coef table keyed by taxon index j and covariate index k. Only the
# treatment coefficient (k = 2 under the ~group design) is tested, declared before
# any discovery is inspected. radEmu constrains its log ratios and identifies
# contrasts only up to that constraint; it cannot recover absolute abundance in
# individual samples from sequencing alone. The constraint is declared in the
# registry and the convergence state is reported here. No adjusted values are
# returned by the package, so the contract applies the declared BH over its p-values.
adapter_rademu <- function(entry, y, meta, opts) {
    data <- as.data.frame(t(y))
    rownames(data) <- meta$sample_id
    taxa <- seq_len(ncol(data))
    fit <- radEmu::emuFit(formula = ~group, data = as.data.frame(meta), Y = data,
                          test_kj = data.frame(j = taxa, k = 2L),
                          run_score_tests = TRUE,
                          remove_zero_comparison_pvals = entry$parameters$remove_zero_comparison_pvals %||% 0.01,
                          tolerance = entry$parameters$tolerance %||% 1e-4,
                          verbose = FALSE)
    table <- fit$coef
    if (is.null(table)) MAG$fail("radEmu returned no coefficient table")
    table <- as.data.frame(table)
    # One row per taxon for the Treatment coefficient; k is the design column only in
    # the request, the answer is keyed by covariate name and taxon.
    tested <- table[table$covariate == "groupTreatment", , drop = FALSE]
    if (!nrow(tested)) MAG$fail("radEmu returned no Treatment coefficient tests")
    index <- MAG$row_match(as.character(tested$category), rownames(y), "radEmu results")
    lfc <- rep(NA_real_, nrow(y)); p <- rep(NA_real_, nrow(y)); q <- rep(NA_real_, nrow(y))
    # The documentation states the units explicitly: an estimate of 1 means exp(1).
    if ("estimate" %in% names(tested)) lfc <- as.numeric(tested$estimate) / log(2)
    if ("pval" %in% names(tested)) p <- as.numeric(tested$pval)
    if (!is.null(fit$estimation_converged) && any(!fit$estimation_converged))
        message("radEmu estimation did not converge for every feature")
    if (!is.null(fit$null_estimation_unconverged) && NROW(fit$null_estimation_unconverged))
        message("radEmu null refits did not converge for ",
                NROW(fit$null_estimation_unconverged),
                " tests; those rows keep whatever the package reported")
    message("radEmu identification constraint: pseudohuber median of log ratios is zero; ",
            "contrasts are relative to the typical taxon change")
    list(log2fc = lfc, pvalue = p, qvalue = NULL)
}