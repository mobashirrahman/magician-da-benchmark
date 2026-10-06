# Random-score ranker: uniform p-values, zero effects. Without a naive
# baseline a ranking is uninterpretable; no method may lose to chance.
adapter_random_score <- function(entry, y, meta, opts) {
    n <- nrow(y)
    p <- stats::runif(n)
    list(log2fc = rep(0, n), pvalue = p, qvalue = NULL)
}
