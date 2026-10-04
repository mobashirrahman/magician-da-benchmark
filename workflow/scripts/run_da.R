#!/usr/bin/env Rscript
# One method per job: explicit status, version/session record, no silent skips.
#
# Usage: run_da.R <registry-entry-json> <matrix> <samples> <features>
#                  <results.tsv> <status.json> <session.txt>
#                  <alpha> <mc-samples> <seed> <timeout-seconds> <allow-failures>
#                  <case> <metric> <kind> <scratch>
#
# All shared input and output rules live in da/contract.R; a single adapter under
# da/ does the package-specific work. The registry entry is the only source of
# declared endpoint, units, correction family and settings.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 16) stop("run_da.R expects 16 arguments; see the header of this script")
REGISTRY_JSON <- args[1]
matrix_file <- args[2]; metadata_file <- args[3]; features_file <- args[4]
out <- args[5]; status_out <- args[6]; session_out <- args[7]
alpha <- as.numeric(args[8]); mc <- as.integer(args[9])
seed <- as.integer(args[10]); timeout_seconds <- as.numeric(args[11])
allow_failures <- args[12] == "true"
case <- args[13]; metric <- args[14]; kind <- args[15]; scratch <- args[16]

Sys.setenv(TZ = "UTC")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for status output")
script <- grep("^--file=", commandArgs(), value = TRUE)
here <- if (length(script)) dirname(sub("^--file=", "", script[1])) else getwd()
source(file.path(here, "da", "contract.R"))

set.seed(seed)
# One core per job, including implicit package workers.
options(mc.cores = 1, warn = 1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

start <- proc.time()[[3]]
# An environment, so that every helper sees the same declared values and the flags
# each step adds (zero policy, prepared membership) are shared rather than copied.
entry <- new.env(parent = emptyenv())
declared <- jsonlite::fromJSON(REGISTRY_JSON, simplifyVector = FALSE)
for (key in names(declared)) assign(key, declared[[key]], envir = entry)
entry$id <- entry$method_id
entry$settings <- entry$parameters
body <- list(case = case, metric = metric, kind = kind, status = "success", message = "",
             package_version = "base R", elapsed_seconds = 0, n_tested = 0L,
             n_retained = 0L, n_reserved_other = 0L, n_features = 0L, raw_p_available = FALSE,
             seed = seed, mc_samples = mc, alpha = alpha, timeout_seconds = timeout_seconds,
             exclusions = list(), message_log = character(0))
result <- NULL
all_ids <- NULL
status <- "success"
logs <- character(0)

note <- function(...) logs <<- c(logs, paste0(...))

# The adapter runs in a forked child so that an explicit deadline can be enforced and
# recorded. A timeout is a reported outcome, never a silent retry with another method.
supervise <- function(code, seconds) {
    job <- parallel::mcparallel(code, name = "magician_da", silent = TRUE)
    deadline <- Sys.time() + seconds
    repeat {
        collected <- parallel::mccollect(job, wait = FALSE)
        if (length(collected)) {
            # mccollect returns one element per finished child: our own value list.
            return(list(state = "done", value = collected[[1]]))
        }
        if (as.numeric(Sys.time() - deadline, units = "secs") >= 0) {
            # Terminate the worker outright, then reap it. A killed job is reported as a
            # timeout; it is never silently retried with a different statistical method.
            for (signal in c(tools::SIGTERM, tools::SIGKILL)) {
                try(tools::pskill(job$pid, signal), silent = TRUE)
                Sys.sleep(0.2)
            }
            suppressWarnings(parallel::mccollect(job, wait = FALSE))
            return(list(state = "timeout", value = NULL))
        }
        Sys.sleep(0.2)
    }
}

prepare <- function() {
    if (!(metric %in% unlist(entry$accepted_inputs)))
        MAG$fail("Method ", entry$id, " does not accept ", metric, " input; it accepts ",
                 paste(unlist(entry$accepted_inputs), collapse = ", "))
    entry$parameters$zero_policy <- if (identical(entry$zero_policy, "none")) NULL else entry$zero_policy
    if (identical(entry$transform, "none") && !is.null(entry$parameters$zero_policy))
        MAG$fail("A zero policy without a transform would change nothing")
    if (is.null(entry$truth_lane_policy)) entry$truth_lane_policy <- entry$zero_policy
    body$truth_lane_policy <<- entry$truth_lane_policy
    pkg <- entry$package
    if (!is.null(pkg)) {
        if (!length(find.package(pkg, quiet = TRUE))) {
            status <<- "unavailable"
            MAG$fail("Required package is missing: ", pkg)
        }
        loadNamespace(pkg)
        body$package_version <<- as.character(utils::packageVersion(pkg))
    }
    input <- MAG$read_inputs(matrix_file, metadata_file, features_file)
    all_ids <<- rownames(input$x)
    body$reserve_other <<- any(input$features$reserved_other)
    meta <- MAG$check_contrast(input$meta)
    body$n_features <<- length(all_ids)
    prepared <- MAG$prepare_input(entry, input)
    entry$keep <- prepared$keep
    entry$reserved <- prepared$reserved
    body$n_retained <<- sum(prepared$keep)
    body$n_reserved_other <<- sum(prepared$reserved)
    # A method that only accepts raw counts must receive integers.
    if (identical(unlist(entry$accepted_inputs), "counts") && metric == "counts" &&
        any(abs(prepared$y - round(prepared$y)) > 1e-8))
        MAG$fail("Counts must be integers")
    if (!sum(prepared$keep))
        MAG$signal_empty()
    adapter <- file.path(here, "da", paste0(entry$adapter, ".R"))
    if (!file.exists(adapter)) MAG$fail("Missing adapter file: ", adapter)
    list(y = prepared$y, meta = meta, library = prepared$library,
         policy = entry$parameters$zero_policy, adapter = adapter)
}

run_adapter <- function(context, meta) {
    source(context$adapter, local = TRUE)
    handler <- paste0("adapter_", entry$adapter)
    if (!exists(handler, mode = "function")) MAG$fail("Adapter does not define ", handler, "()")
    opts <- list(alpha = alpha, mc_samples = mc, seed = seed, scratch = scratch,
                 library = context$library, policy = context$policy)
    do.call(handler, list(entry = entry, y = context$y, meta = meta, opts = opts))
}

outcome <- tryCatch({
    context <- withCallingHandlers(prepare(),
        message = function(m) {note(conditionMessage(m)); invokeRestart("muffleMessage")},
        warning = function(w) {note(conditionMessage(w)); invokeRestart("muffleWarning")})
    meta <- context$meta
    supervised <- supervise({
        captured <- character(0)
        value <- withCallingHandlers(
            run_adapter(context, meta),
            message = function(m) {captured <<- c(captured, conditionMessage(m)); invokeRestart("muffleMessage")},
            warning = function(w) {captured <<- c(captured, conditionMessage(w)); invokeRestart("muffleWarning")})
        list(value = value, logs = captured)
    }, max(timeout_seconds, 0.05))
    if (identical(supervised$state, "timeout")) {
        # Plain assignment at top level: `<<-` would skip globalenv and find base::body.
        status <- "timeout"
        body$message <- paste("Exceeded the declared", timeout_seconds, "second budget")
        NULL
    } else {
        value <- supervised$value
        if (inherits(value, "try-error") || is.null(value) || is.null(value$value)) stop(value)
        note(value$logs)
        value$value
    }
}, mag_empty = function(e) {status <<- "empty"; body$message <<- conditionMessage(e); NULL},
   error = function(e) {
       if (status != "unavailable") status <<- "failed"
       body$message <<- conditionMessage(e)
       NULL
   })
body$message_log <- logs
body$elapsed_seconds <- proc.time()[[3]] - start

if (!is.null(outcome)) {
    tryCatch({
        result <- MAG$collect(outcome, entry,
                              list(keep = entry$keep, reserved = entry$reserved,
                                   testable = entry$keep & !entry$reserved),
                              all_ids)
        body$n_tested <- sum(result$tested)
        body$raw_p_available <- any(is.finite(result$pvalue))
        body$exclusions <- MAG$report_exclusions(result)
    }, error = function(e) {
        status <<- "failed"
        body$message <<- conditionMessage(e)
        result <<- NULL
    })
}

if (is.null(result)) {
    entry$keep <- rep(FALSE, body$n_features)
    entry$reserved <- rep(FALSE, body$n_features)
    ids <- all_ids %||% character(0)
    empty <- data.frame(feature_id = character(0), endpoint = character(0), contrast = character(0),
                        effect_scale = character(0), reference = character(0),
                        adjustment_family = character(0), zero_policy = character(0),
                        log2fc = numeric(0), pvalue = numeric(0), qvalue = numeric(0),
                        tested = logical(0), fit_status = character(0), exclusion_reason = character(0),
                        stringsAsFactors = FALSE)
    reason <- if (status == "unavailable") status else body$message
    reason <- gsub("[\t\r\n]+", " ", reason)
    result <- if (length(ids)) data.frame(
        feature_id = ids, endpoint = entry$endpoint, contrast = entry$contrast,
        effect_scale = entry$effect_scale, reference = entry$reference,
        adjustment_family = entry$adjustment_family,
        zero_policy = if (is.null(entry$parameters$zero_policy)) "none" else entry$parameters$zero_policy$name,
        log2fc = NA_real_, pvalue = NA_real_, qvalue = NA_real_, tested = FALSE,
        fit_status = ifelse(status == "unavailable", "not_run", "no_result"),
        exclusion_reason = reason,
        stringsAsFactors = FALSE) else empty
    body$n_features <- length(ids)
}
write.table(result, out, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
body$status <- status
MAG$write_status(status_out, entry, list(zero_policy = if (is.null(entry$parameters$zero_policy))
                                                         "none" else entry$parameters$zero_policy$name,
                                         registry_digest = entry$registry_digest, body = body))
writeLines(capture.output(sessionInfo()), session_out)
cat(entry$id, status, body$message, "\n")
if (!allow_failures && status %in% c("failed", "unavailable", "timeout")) quit(status = 1)