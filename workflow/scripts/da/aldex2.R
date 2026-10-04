# ALDEx2 on raw counts. Gamma is a declared variant of one Dirichlet family, set
# before any discovery is inspected. The installed release is asked for its scale
# argument; a release without it fails visibly rather than silently using the
# previous behaviour.
adapter_aldex2 <- function(entry, y, meta, opts) {
    gamma <- as.numeric(entry$parameters$gamma)
    if (is.na(gamma)) MAG$fail("ALDEx2 gamma must be a declared number")
    mc <- as.integer(entry$parameters$mc_samples %||% opts$mc_samples)
    # The installed release decides the argument names; 1.42 calls the count matrix
    # `reads`. Ask for them instead of assuming.
    declared <- names(formals(ALDEx2::aldex.clr))
    if (!length(declared)) MAG$fail("Cannot inspect the installed ALDEx2::aldex.clr signature")
    arguments <- list()
    arguments[[declared[1]]] <- y
    if ("gamma" %in% declared) {
        if (is.na(gamma)) MAG$fail("ALDEx2 gamma must be a declared number")
        arguments$gamma <- gamma
    } else if (gamma != 0) {
        MAG$fail("Installed ALDEx2 has no gamma argument; gamma variants require it")
    }
    for (name in c("conds", "mc.samples", "denom", "verbose", "useMC")) {
        if (!name %in% declared) next
        arguments[[name]] <- switch(name,
            conds = as.character(meta$group),
            "mc.samples" = mc,
            denom = "all",
            verbose = FALSE,
            useMC = FALSE)
    }
    clr <- do.call(ALDEx2::aldex.clr, arguments)
    test <- ALDEx2::aldex.ttest(clr, paired.test = FALSE, verbose = FALSE)
    effect <- ALDEx2::aldex.effect(clr, verbose = FALSE, useMC = FALSE)
    # we.eBH is ALDEx2's own correction over the Dirichlet draws; never re-adjust it.
    list(log2fc = effect$diff.btw, pvalue = test$we.ep, qvalue = test$we.eBH)
}