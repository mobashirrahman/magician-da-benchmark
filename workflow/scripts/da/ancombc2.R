# ANCOM-BC2 lives in the ANCOMBC package; there is no separate ANCOMBC2 package.
# Internal prevalence and library filtering are disabled: the shared contract has
# already applied the declared common filter, and a second filter inside the
# package would silently remove features from the multiple-testing family.
adapter_ancombc2 <- function(entry, y, meta, opts) {
    fit <- ANCOMBC::ancombc2(data = y, meta_data = meta, fix_formula = "group", group = "group",
                             p_adj_method = entry$adjustment_family,
                             pseudo_sens = isTRUE(entry$parameters$pseudo_sens),
                             prv_cut = 0, lib_cut = 0,
                             struc_zero = isTRUE(entry$parameters$struc_zero),
                             alpha = opts$alpha, n_cl = 1, verbose = FALSE)$res
    fit <- fit[match(rownames(y), fit$taxon), , drop = FALSE]
    q <- fit$q_groupTreatment
    # Sensitivity analysis failures are not discoveries.
    if ("passed_ss_groupTreatment" %in% names(fit) && any(!fit$passed_ss_groupTreatment))
        q[!fit$passed_ss_groupTreatment] <- 1
    list(log2fc = fit$lfc_groupTreatment / log(2),
         pvalue = fit$p_groupTreatment, qvalue = q)
}