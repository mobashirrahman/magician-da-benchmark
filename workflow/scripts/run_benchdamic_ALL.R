#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
#  BenchDamic – run differential‑abundance wrappers (clean version)
#  ---------------------------------------------------------------------------
#  Usage:  run_benchdamic_ALL.R  <counts.tsv>  <metadata.tsv>  <group_col>  <out_dir>


suppressPackageStartupMessages({
  library(TreeSummarizedExperiment)
  library(benchdamic)
})

# ----------  command line --------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: run_benchdamic_ALL.R  <counts.tsv>  <metadata.tsv>  <group_col>  <out_dir>")

cts_file <- args[1]; meta_file <- args[2]; grp_col <- args[3]; out_dir <- args[4]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------  read data ---------------------------------------------------------
cts  <- read.delim(cts_file, row.names = 1, check.names = FALSE, na.strings = c("", "NA"))
cts[is.na(cts)] <- 0
cts  <- data.matrix(cts)

# Convert to integers for DESeq2 compatibility (some metrics like coverage are decimals)
# Round and convert to integers - this is necessary for DESeq2 and other count-based methods
cts <- round(cts)
storage.mode(cts) <- "integer"
meta <- read.delim(meta_file, row.names = 1, check.names = FALSE)

common <- intersect(colnames(cts), rownames(meta))
if (!length(common)) stop("Counts and metadata share zero sample IDs.")
cts  <- cts[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

meta[[grp_col]] <- factor(meta[[grp_col]])
if (length(levels(meta[[grp_col]])) != 2)
  stop("group column must have exactly 2 levels – found ", length(levels(meta[[grp_col]])))
if ("Control" %in% levels(meta[[grp_col]]))
  meta[[grp_col]] <- relevel(meta[[grp_col]], ref = "Control")
levs <- levels(meta[[grp_col]])
contrast_vec <- c(grp_col, levs[2], levs[1])
coef_val     <- 2

# ----------  build TSE & add normalization factors ----------------------------
tse <- TreeSummarizedExperiment(assays = list(counts = as.matrix(cts)), colData = meta)

# add all normalisations that wrappers might look for
suppressMessages({
  tse <- norm_edgeR(tse, method = "TMM")     # adds NF.TMM
  tse <- norm_CSS(tse)                        # adds NF.CSS
  tse <- norm_DESeq2(tse, method = "ratio")  # adds NF.ratio (size factors)
  tse <- norm_DESeq2(tse, method = "poscounts")  # adds NF.poscounts (size factors)
})

# ----------  helper -----------------------------------------------------------
run_and_export <- function(fun, mname, extras) {
  cat(sprintf("\n[%s]  Running %s...\n", format(Sys.time(), "%H:%M:%S"), mname))
  res <- tryCatch(do.call(fun, c(list(object = tse, assay_name = "counts"), extras)),
                  error = function(e) {cat("   ✗ ", mname, " failed: ", conditionMessage(e), "\n", sep = ""); NULL})
  if (is.null(res)) return(NULL)
  out_csv <- file.path(out_dir, paste0("DA_", mname, ".csv"))
  
  # Enhanced export with log fold change
  export_df <- NULL
  
  # Extract p-values (always available)
  if (!is.null(res$pValMat)) {
    export_df <- res$pValMat
  } else {
    export_df <- data.frame(rawP = NA, adjP = NA)
  }
  
  # Extract log fold changes if available
  if (!is.null(res$statInfo) && "logFC" %in% colnames(res$statInfo)) {
    # Ensure same row order and add logFC
    if (nrow(res$statInfo) == nrow(export_df) && 
        all(rownames(res$statInfo) == rownames(export_df))) {
      export_df$logFC <- res$statInfo$logFC
    }
  } else if (!is.null(res$statInfo) && "log2FoldChange" %in% colnames(res$statInfo)) {
    # For DESeq2-style results
    if (nrow(res$statInfo) == nrow(export_df) && 
        all(rownames(res$statInfo) == rownames(export_df))) {
      export_df$logFC <- res$statInfo$log2FoldChange
    }
  }
  
  # Write the results
  write.csv(export_df, out_csv, row.names = TRUE)
  cat("   ✓ ", mname, " done → ", out_csv, " (columns: ", paste(colnames(export_df), collapse = ", "), ")\n", sep = "")
  invisible(res)
}

# ----------  wrapper parameter sets -------------------------------------------
wrap <- list(
  ALDEx2 = list(fun = DA_ALDEx2,
                args = list(design = grp_col, contrast = contrast_vec, test = "t")),

  ANCOM  = list(fun = DA_ANCOM,
                args = list(fix_formula = grp_col, contrast = contrast_vec, BC = TRUE)),

  basic  = list(fun = DA_basic,
                args = list(test = "t", contrast = contrast_vec)),

  corncob = list(fun = DA_corncob,
                 args = list(formula = as.formula(paste("~", grp_col)),
                              phi.formula = as.formula(paste("~", grp_col)),
                              formula_null = ~1, phi.formula_null = as.formula(paste("~", grp_col)),
                              coefficient = paste0(grp_col, levs[2]), test = "Wald")),

  dearseq = list(fun = DA_dearseq,
                 args = list(variables2test = grp_col, test = "asymptotic")),

  DESeq2 = list(fun = DA_DESeq2,
                args = list(design = as.formula(paste("~", grp_col)), contrast = contrast_vec, norm = "")),

  edgeR  = list(fun = DA_edgeR,
                args = list(design = as.formula(paste("~", grp_col)), coef = coef_val, norm = "TMM")),

  limma  = list(fun = DA_limma,
                args = list(design = as.formula(paste("~", grp_col)), coef = coef_val, norm = "TMM")),

  metagenomeSeq = list(fun = DA_metagenomeSeq,
                       args = list(design = as.formula(paste("~", grp_col)), coef = coef_val, norm = "CSS")),

  NOISeq = list(fun = DA_NOISeq,
                args = list(norm = "tmm", contrast = contrast_vec))
)

# ----------  run ----------------------------------------------------------------
results <- lapply(names(wrap), function(nm) run_and_export(wrap[[nm]]$fun, nm, wrap[[nm]]$args))

cat("\nFinished all wrappers. Results written to:\n")
for (nm in names(wrap)) cat("  • ", file.path(out_dir, paste0("DA_", nm, ".csv")), "\n", sep = "")
