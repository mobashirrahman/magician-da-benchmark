# metaGEENOME: differential abundance from its GEE component only.
#
# The published package bundles DA with extensive downstream reporting. This adapter
# therefore requires an isolated DA entry point; if the wrapper cannot be called
# without taxonomy or diversity machinery, the job fails visibly and the method stays
# deferred rather than being scored on partial output.
adapter_metageenome <- function(entry, y, meta, opts) {
    if (!requireNamespace("metaGEENOME", quietly = TRUE))
        MAG$fail("metaGEENOME is not installed; the method is deferred until its DA component is validated")
    entry_points <- c("geeDA", "metaGEENOME", "run_DA")
    chosen <- entry_points[vapply(entry_points,
                                  function(name) exists(name, where = asNamespace("metaGEENOME"),
                                                        inherits = FALSE), logical(1))]
    if (!length(chosen)) MAG$fail("metaGEENOME exposes no isolated DA entry point")
    call <- get(chosen[1], envir = asNamespace("metaGEENOME"))
    data <- as.data.frame(t(y))
    rownames(data) <- meta$sample_id
    fit <- tryCatch(do.call(call, list(counts = data, metadata = meta, group = "group",
                                       CTF = isTRUE(entry$parameters$ctf))),
                    error = function(e) NULL)
    if (is.null(fit)) MAG$fail("metaGEENOME DA component failed on the declared inputs")
    tables <- Filter(function(x) is.data.frame(x) || is.matrix(x), fit)
    if (!length(tables)) MAG$fail("metaGEENOME returned no DA table")
    table <- as.data.frame(tables[[1]])
    index <- MAG$row_match(as.character(table[[1]]), rownames(y), "metaGEENOME results")
    lfc <- rep(NA_real_, nrow(y)); p <- rep(NA_real_, nrow(y)); q <- rep(NA_real_, nrow(y))
    for (candidate in c("beta", "estimate", "Estimate", "coef"))
        if (candidate %in% names(table)) {lfc <- table[[candidate]][index]; break}
    for (candidate in c("p.value", "pvalue", "pval"))
        if (candidate %in% names(table)) {p <- table[[candidate]][index]; break}
    list(log2fc = lfc, pvalue = p, qvalue = q)
}