# Shared input/output contract for every differential abundance adapter.
#
# Adapters receive one matrix, one metadata table and one feature table, and
# return effects in the units the registry declares. Everything that must hold for
# all methods lives here: strict sample and feature joins, the declared unit-aware
# zero policy, the CLR transform, native adjusted values, per-feature fit status
# and an explicit exclusion reason. Adapters never see feature names that could be
# silently rewritten by a package.

MAG <- new.env()

MAG$fail <- function(...) stop(paste0(...), call. = FALSE)

# Flag columns are written by pandas, so accept TRUE, True, true, 1 and NA alike.
MAG$as_flag <- function(values, label) {
    if (is.logical(values)) return(values)
    text <- tolower(trimws(as.character(values)))
    if (!all(text %in% c("true", "false", "1", "0", "na", "")))
        MAG$fail(label, " must contain only booleans")
    out <- text %in% c("true", "1")
    out[text %in% c("na", "")] <- NA
    out
}

# An empty catalogue is a valid outcome, not a failure. It gets its own status so
# that a method can never look ineligible for a reason the report cannot show.
MAG$signal_empty <- function() {
    stop(structure(class = c("mag_empty", "error", "condition"),
                   list(message = "No features pass common count/prevalence filtering",
                        call = NULL)))
}

MAG$read_inputs <- function(matrix_file, metadata_file, features_file) {
    x <- as.matrix(read.delim(matrix_file, row.names = 1, check.names = FALSE))
    meta <- read.delim(metadata_file, row.names = 1, check.names = FALSE)
    features <- read.delim(features_file, row.names = 1, check.names = FALSE)
    if (anyDuplicated(colnames(x)) || anyDuplicated(rownames(x)))
        MAG$fail("Matrix sample or feature IDs are duplicated")
    if (!setequal(colnames(x), rownames(meta)))
        MAG$fail("Matrix and sample metadata IDs must match exactly; refusing to intersect silently")
    # Explicit, recorded reordering. Never a silent column intersection.
    reordered <- !identical(colnames(x), rownames(meta))
    meta <- meta[colnames(x), , drop = FALSE]
    if (!identical(as.character(rownames(x)), as.character(rownames(features))))
        MAG$fail("Feature metadata order disagrees with matrix")
    if (!"retained_for_da" %in% names(features))
        MAG$fail("Feature table has no retained_for_da column")
    features$reserved_other <- if ("reserved_other" %in% names(features))
        MAG$as_flag(features$reserved_other, "reserved_other") else rep(FALSE, nrow(features))
    features$reserved_other[is.na(features$reserved_other)] <- FALSE
    if (!"read_pairs" %in% names(meta))
        MAG$fail("Sample metadata must carry read_pairs so that zero policies are unit-aware")
    if (!"group" %in% names(meta))
        MAG$fail("Sample metadata has no group column")
    if (any(!is.finite(x)) || any(x < 0))
        MAG$fail("Abundance values must be finite and nonnegative")
    list(x = x, meta = meta, features = features, reordered = reordered)
}

MAG$check_contrast <- function(meta) {
    groups <- unique(as.character(meta$group))
    if (!setequal(groups, c("Control", "Treatment")))
        MAG$fail("Exactly Control and Treatment groups are required")
    counts <- table(as.character(meta$group))
    if (any(counts < 2))
        MAG$fail("Each group needs at least two samples for a variance estimate")
    # A reversed label or a renamed sample must be rejected, never quietly accepted:
    # sample identifiers are generated in group order, so the labels have to agree.
    ordered <- meta[order(rownames(meta)), , drop = FALSE]
    labels <- as.character(ordered$group)
    blocks <- c(rep("Control", sum(labels == "Control")), rep("Treatment", sum(labels == "Treatment")))
    if (!identical(labels, blocks))
        MAG$fail("Group labels disagree with the sorted sample order; a reversed contrast is rejected")
    # Treatment against Control, always. A reversed factor would invert every sign.
    meta$group <- factor(meta$group, levels = c("Control", "Treatment"))
    meta$sample_id <- rownames(meta)
    meta
}

MAG$library <- function(meta) {
    values <- as.numeric(meta$read_pairs)
    if (any(!is.finite(values)) || any(values <= 0))
        MAG$fail("read_pairs must be positive and finite for a unit-aware zero policy")
    values
}

# Declared zero handling. `measurement` adds the pseudocount in the metric's own
# units; `count` adds it in count units, so one constant means the same thing
# whatever the input scale.
MAG$apply_zero_policy <- function(y, policy, library) {
    if (is.null(policy) || identical(policy$zero_handling, "none"))
        return(y)
    if (identical(policy$zero_handling, "pseudocount")) {
        offset <- as.numeric(policy$pseudocount)
        if (is.na(offset) || offset < 0) MAG$fail("pseudocount must be a nonnegative number")
        if (identical(policy$pseudocount_scale %||% "measurement", "measurement"))
            return(y + offset)
        return(y + matrix(offset / library, nrow = nrow(y), ncol = ncol(y), byrow = TRUE))
    }
    if (identical(policy$zero_handling, "dirichlet")) {
        alpha <- as.numeric(policy$dirichlet_alpha)
        if (is.na(alpha) || alpha <= 0) MAG$fail("dirichlet_alpha must be positive")
        totals <- colSums(y)
        posterior <- (y + alpha) / rep(totals + alpha * nrow(y), each = nrow(y))
        return(ifelse(y > 0, y, posterior))
    }
    MAG$fail("Unknown zero_handling: ", policy$zero_handling)
}

# Explicit base-2 centred log-ratio, matching the truth lane coordinates.
MAG$centre_log_ratio <- function(y) {
    values <- log2(y)
    values - rowMeans(values)
}

MAG$prepare_input <- function(entry, input) {
    keep <- MAG$as_flag(input$features$retained_for_da, "retained_for_da")
    keep[is.na(keep)] <- FALSE
    reserved <- MAG$as_flag(input$features$reserved_other, "reserved_other")
    reserved[is.na(reserved)] <- FALSE
    # A reserved `other` category carries filtered mass so totals are preserved. It
    # belongs to the composition and never to the tested set.
    y <- input$x[keep | reserved, , drop = FALSE]
    library <- MAG$library(input$meta)
    zero_policy <- entry$parameters$zero_policy
    if (!identical(entry$transform, "none")) {
        y <- MAG$apply_zero_policy(y, zero_policy, library)
    }
    if (identical(entry$transform, "clr"))
        y <- MAG$centre_log_ratio(y)
    else if (!identical(entry$transform, "none"))
        MAG$fail("Unknown transform: ", entry$transform)
    if (any(!is.finite(y))) MAG$fail("Prepared input contains non-finite values")
    # Packages may sanitize names ('-' becomes '.'); stable internal IDs keep the
    # original genome and MAG IDs in the exported result.
    rownames(y) <- sprintf("feature_%06d", seq_len(nrow(y)))
    # A reserved category joins the geometry but is never a discovery, even if a
    # contradictory table marked it retained as well.
    list(y = y, keep = keep, reserved = reserved, testable = keep & !reserved,
         library = library)
}

MAG$result_frame <- function(entry, ids) {
    n <- length(ids)
    status <- rep("no_result", n)
    reason <- rep(NA_character_, n)
    excluded <- !entry$keep & !entry$reserved
    reserved <- entry$reserved
    status[excluded] <- "filtered"; reason[excluded] <- "common count/prevalence filter"
    status[reserved] <- "reserved_category"; reason[reserved] <- "reserved other category, never a discovery"
    data.frame(feature_id = as.character(ids),
               endpoint = entry$endpoint,
               contrast = entry$contrast,
               effect_scale = entry$effect_scale,
               reference = entry$reference,
               adjustment_family = entry$adjustment_family,
               zero_policy = if (is.null(entry$parameters$zero_policy)) "none" else entry$parameters$zero_policy$name,
               log2fc = rep(NA_real_, n),
               pvalue = rep(NA_real_, n),
               qvalue = rep(NA_real_, n),
               tested = rep(FALSE, n),
               fit_status = status,
               exclusion_reason = reason,
               stringsAsFactors = FALSE)
}

# Adapters return a list of length-nrow(y) vectors, in the row order they received.
MAG$collect <- function(fit, entry, prepared, ids) {
    n <- sum(prepared$keep | prepared$reserved)
    if (is.null(fit$log2fc) || length(fit$log2fc) != n)
        MAG$fail("Adapter returned ", length(fit$log2fc), " effects for ", n, " input features")
    p <- if (is.null(fit$pvalue)) rep(NA_real_, n) else fit$pvalue
    q <- if (is.null(fit$qvalue)) rep(NA_real_, n) else fit$qvalue
    for (name in c("log2fc", "pvalue", "qvalue")) {
        v <- switch(name, log2fc = fit$log2fc, pvalue = p, qvalue = q)
        if (name != "log2fc" && any(is.finite(v) & (v < 0 | v > 1)))
            MAG$fail("Invalid ", name, " values")
    }
    if (entry$native_q && !any(is.finite(q)) && any(is.finite(p)))
        MAG$fail("Registry declares native adjusted values but the adapter returned none")
    frame <- MAG$result_frame(entry, ids)
    supplied <- prepared$keep | prepared$reserved
    frame$log2fc[supplied] <- fit$log2fc
    frame$pvalue[supplied] <- p
    frame$qvalue[supplied] <- q
    if (!entry$native_q) {
        # Never re-adjust a package's own q-values; BH over our own p-values is the
        # declared correction when the package returns none.
        adjusted <- stats::p.adjust(p, method = entry$adjustment_family)
        frame$qvalue[supplied] <- adjusted
    }
    # Only retained, non-reserved features are discoveries. A reserved category keeps
    # its value but is never tested.
    finite <- prepared$testable & is.finite(frame$pvalue) & is.finite(frame$qvalue)
    frame$tested <- finite
    frame$fit_status[finite] <- "fit"
    frame$exclusion_reason[finite] <- NA_character_
    quiet <- prepared$testable & !finite
    frame$exclusion_reason[quiet] <- "adapter returned no finite test"
    if (!any(frame$tested)) MAG$fail("Method produced no finite tests")
    frame
}

`%||%` <- function(a, b) if (is.null(a)) b else a

MAG$row_match <- function(observed, wanted, label) {
    if (is.null(observed)) MAG$fail("Adapter returned no ", label)
    idx <- match(wanted, observed)
    if (any(is.na(idx))) MAG$fail("Adapter did not report ", label, " for every tested feature")
    idx
}

MAG$report_exclusions <- function(frame) {
    reasons <- frame$exclusion_reason[!is.na(frame$exclusion_reason)]
    if (!length(reasons)) return(list())
    as.list(table(reasons))
}

MAG$write_status <- function(path, entry, extra) {
    jsonlite::write_json(c(list(method = entry$id, method_family = entry$method_family,
                                variant = entry$variant, adapter = entry$adapter,
                                package = if (is.null(entry$package)) "base R" else entry$package,
                                endpoint = entry$endpoint, effect_scale = entry$effect_scale,
                                reference = entry$reference, contrast = entry$contrast,
                                adjustment_family = entry$adjustment_family,
                                native_q = entry$native_q,
                                transform = entry$transform,
                                zero_policy = extra$zero_policy,
                                settings = entry$parameters,
                                registry_digest = extra$registry_digest),
                        extra$body),
                       path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}