#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(TreeSummarizedExperiment)
  library(benchdamic)
})

# ----  Command-line args --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: run_benchdamic_ALL.R  <counts.tsv>  <metadata.tsv> ",
       "<group_column>  <out_dir>")

cts_file   <- args[1]
meta_file  <- args[2]
group_col  <- args[3]          # e.g. "group"
out_dir    <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----  Import data -------------------------------------------------------
cts  <- read.delim(cts_file,  row.names = 1, check.names = FALSE)
meta <- read.delim(meta_file, row.names = 1)

tse  <- TreeSummarizedExperiment(
          assays  = list(counts = as.matrix(cts)),
          colData = meta)

# ----  Discover every DA_*() wrapper ------------------------------------
all_wrappers <- grep("^DA_", ls("package:benchdamic"), value = TRUE)
message("Discovered ", length(all_wrappers), " DA wrappers.")

# Optional: drop wrappers you know will fail (e.g. need paired data)
# all_wrappers <- setdiff(all_wrappers, c("DA_MAST", "DA_Seurat"))

# ----  Helper: run one method safely ------------------------------------
run_one <- function(fun_name, tse_obj, grp) {
  fun <- get(fun_name, envir = asNamespace("benchdamic"))
  message("[", format(Sys.time(), "%H:%M:%S"), "]  ", fun_name)

  # Minimal arg list that works for 95 % of wrappers
  base_args <- list(
      object       = tse_obj,
      assay_name   = "counts",
      pseudo_count = FALSE
  )

  # Add *either* formula or contrast depending on the wrapper
  if ("formula"  %in% names(formals(fun)))
      base_args$formula  <- as.formula(paste("~", grp))
  if ("contrast" %in% names(formals(fun)))
      base_args$contrast <- c(grp,
                              levels(colData(tse_obj)[[grp]])[2],
                              levels(colData(tse_obj)[[grp]])[1])

  # Silently skip if the wrapper errors out
  tryCatch({
      do.call(fun, base_args)
    },
    error = function(e) {
      warning(fun_name, " failed (", e$message, "). Skipped.")
      NULL
    })
}

# ----  Run all methods ---------------------------------------------------
results <- lapply(all_wrappers, run_one, tse_obj = tse, grp = group_col)
names(results) <- gsub("^DA_", "", all_wrappers)

# Keep only successful runs
results <- Filter(Negate(is.null), results)

# ----  Save tidy tables --------------------------------------------------
export_one <- function(x, method, dir) {
  out <- file.path(dir, paste0("DA_", method, ".csv"))
  write.csv(benchdamic::exportDA(x, method = method),
            file = out, row.names = FALSE)
  out
}
out_files <- mapply(export_one,
                    x = results,
                    method = names(results),
                    MoreArgs = list(dir = out_dir),
                    SIMPLIFY = FALSE)

cat("\nFinished.\nResults written to:\n",
    paste0("  • ", unlist(out_files), "\n"), sep = "")
