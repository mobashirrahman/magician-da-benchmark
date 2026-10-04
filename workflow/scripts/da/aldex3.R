# ALDEx 3 scale-model regression on raw counts (1.3.1 API).
#
# aldex() takes the count matrix first and the formula second; gamma travels in the
# scale function's own arguments; nsample sets the Monte Carlo count; streamsize is
# the explicit chunk limit instead of the 8000 MB default. p.val and p.val.adj are
# members of the fitted object by declaration, so return.pars requests them and the
# adapter refuses to proceed if either is missing. The effect is the mean of the
# Treatment coefficient draws, an aggregate over the Monte Carlo replicates that is
# never persisted. ALDEx2 stays in the benchmark as the baseline for this family.
adapter_aldex3 <- function(entry, y, meta, opts) {
    gamma <- as.numeric(entry$parameters$gamma)
    if (is.na(gamma)) MAG$fail("ALDEx3 gamma must be a declared number")
    nsample <- as.integer(entry$parameters$mc_samples %||% opts$mc_samples)
    stream <- as.numeric(entry$parameters$stream_size_mb %||% 256)
    scale_name <- entry$parameters$scale %||% "clr.sm"
    scale_fn <- tryCatch(utils::getFromNamespace(scale_name, "ALDEx3"),
                         error = function(e) NULL)
    if (is.null(scale_fn)) MAG$fail("ALDEx3 exposes no scale function named ", scale_name)
    fit <- ALDEx3::aldex(Y = y, X = ~group, data = meta,
                         method = entry$parameters$method %||% "lm",
                         nsample = nsample, scale = scale_fn, gamma = gamma,
                         streamsize = stream, n.cores = 1L,
                         return.pars = c("estimate", "p.val", "p.val.adj"),
                         p.adjust.method = entry$adjustment_family,
                         test = entry$parameters$test %||% "t.HC3")
    if (is.null(fit$p.val) || is.null(fit$p.val.adj))
        MAG$fail("ALDEx3 fit does not expose both p.val and p.val.adj; refusing to re-adjust")
    treatment <- grep("Treatment", rownames(fit$p.val), value = TRUE)
    if (!length(treatment)) {
        candidates <- grep("Intercept", rownames(fit$p.val), value = TRUE, invert = TRUE)
        if (length(candidates) != 1)
            MAG$fail("ALDEx3 returned no identifiable Treatment coefficient")
        treatment <- candidates
    }
    message("ALDEx3 Treatment coefficient row: ", treatment,
            "; Monte Carlo draws: ", nsample, "; scale: ", scale_name,
            " with gamma ", gamma, "; streaming limit: ", stream, " MB")
    index <- MAG$row_match(colnames(fit$p.val), rownames(y), "ALDEx3 features")
    draws <- fit$estimate[treatment, index, , drop = FALSE]
    lfc <- apply(draws, 2, function(values) mean(values, na.rm = TRUE))
    list(log2fc = as.numeric(lfc),
         pvalue = as.numeric(fit$p.val[treatment, index]),
         qvalue = as.numeric(fit$p.val.adj[treatment, index]))
}