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

# Optional CLI: --mode core|ext|both (default: both)
mode <- "both"
# Optional CLI: --force-counts (default: FALSE) — experimental: build pseudo-counts for non-count metrics
force_counts <- FALSE
if (length(args) > 4) {
  extra <- args[5:length(args)]
  # Support --mode value or --mode=value
  if (any(grepl("^--mode(=|$)", extra))) {
    idx <- which(grepl("^--mode(=|$)", extra))[1]
    token <- extra[idx]
    if (grepl("=", token)) {
      mode <- sub("^--mode=", "", token)
    } else if (length(extra) >= idx + 1) {
      mode <- extra[idx + 1]
    }
  }
  # --force-counts (no value or =true/1/yes)
  if (any(grepl("^--force-counts(=|$)", extra))) {
    idx <- which(grepl("^--force-counts(=|$)", extra))[1]
    token <- extra[idx]
    if (grepl("=", token)) {
      val <- tolower(sub("^--force-counts=", "", token))
      force_counts <- val %in% c("1", "true", "t", "yes", "y")
    } else {
      force_counts <- TRUE
    }
  }
  mode <- tolower(mode)
  if (!mode %in% c("core", "ext", "both")) mode <- "both"
}

# ----------  read data ---------------------------------------------------------
cts  <- read.delim(cts_file, row.names = 1, check.names = FALSE, na.strings = c("", "NA"))
cts[is.na(cts)] <- 0
cts  <- data.matrix(cts)
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

# ----------  detect metric and prepare assays ---------------------------------
# Infer metric from file name (expects e.g., counts_GC.tsv, tpm_GC.tsv, etc.)
metric_name <- sub("(_GC)?\\.tsv$", "", basename(cts_file))

# Build TreeSummarizedExperiment with appropriate assays
assays_list <- list()

# For count-based methods we need integer counts; valid for counts metric, or forced (experimental)
if (tolower(metric_name) == "counts" || force_counts) {
  if (tolower(metric_name) == "counts") {
    counts_int <- round(cts)
  } else {
    # Experimental pseudo-counts for non-count metrics:
    # Scale each sample to a target library size and round to integers.
    # This preserves relative abundances but does NOT reconstruct true counts.
    target_libsize <- 1e6
    libs <- colSums(cts, na.rm = TRUE)
    libs[!is.finite(libs) | libs <= 0] <- 1
    scale_factors <- target_libsize / libs
    scaled <- sweep(cts, 2, scale_factors, "*")
    counts_int <- round(scaled)
  }
  storage.mode(counts_int) <- "integer"
  assays_list$counts <- as.matrix(counts_int)
} else {
  # Keep a placeholder counts assay (numeric); count-based methods will be skipped
  assays_list$counts <- as.matrix(cts)
}

# Build CLR assay for compositional/continuous methods
make_clr_assay <- function(mat, pseudocount = 0.5) {
  stopifnot(all(is.finite(mat)))
  x <- mat + pseudocount
  # Per-sample CLR: center each sample (column) by its geometric mean
  # geometric mean per column
  sgm <- apply(x, 2, function(col) exp(mean(log(col))))
  clr <- sweep(x, 2, sgm, "/")
  log2(clr)
}
assays_list$clr <- make_clr_assay(cts)

# Construct TSE
tse <- TreeSummarizedExperiment(assays = assays_list, colData = meta)

# Add normalisations when counts are provided, or when pseudo-counts are forced (experimental)
if (tolower(metric_name) == "counts" || force_counts) {
  suppressMessages({
    tse <- norm_edgeR(tse, method = "TMM")     # adds NF.TMM
    tse <- norm_CSS(tse)                        # adds NF.CSS
    tse <- norm_DESeq2(tse, method = "ratio")  # adds NF.ratio (size factors)
    tse <- norm_DESeq2(tse, method = "poscounts")  # adds NF.poscounts (size factors)
  })
}

# ----------  helper -----------------------------------------------------------
run_and_export <- function(fun, mname, extras) {
  cat(sprintf("\n[%s]  Running %s...\n", format(Sys.time(), "%H:%M:%S"), mname))
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  res <- tryCatch(do.call(fun, c(list(object = tse), extras)),
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

# ----------  NEW METHOD WRAPPERS ----------------------------------------------

# ANCOM-BC (optional)
DA_ANCOMBC_wrapper <- function(object, assay_name = "counts", fix_formula, rand_formula = NULL) {
  if (!requireNamespace("ancombc", quietly = TRUE))
    stop("Package 'ancombc' not installed.")
  y <- SummarizedExperiment::assay(object, assay_name)
  md <- as.data.frame(SummarizedExperiment::colData(object))
  # Use ancombc2-style interface via ancombc::ancombc2 if available, else ancombc
  if (utils::packageVersion("ancombc") >= "2.0.0" && !is.null(getNamespace("ancombc")$ancombc2)) {
    res <- ancombc::ancombc2(data = list(abn = y, meta = md), fix_formula = fix_formula,
                             rand_formula = rand_formula, p_adj_method = "BH")
    pv <- res$pW; qv <- res$qW; tt <- res$lfc
    term <- colnames(pv)[1]
    p <- pv[, term]; q <- qv[, term]; lfc <- tt[, term]
    rn <- rownames(pv)
  } else {
    fit <- ancombc::ancombc(phyloseq = NULL, formula = fix_formula, W = y, meta_data = md,
                            p_adj_method = "BH", global = FALSE)
    nm <- fix_formula
    p <- fit$res$p_val[[nm]]; q <- fit$res$q_val[[nm]]; lfc <- fit$res$beta[[nm]]
    rn <- names(p)
  }
  list(
    pValMat  = data.frame(rawP = setNames(as.numeric(p), rn), adjP = setNames(as.numeric(q), rn)),
    statInfo = data.frame(logFC = setNames(as.numeric(lfc), rn))
  )
}

# ANCOM-BC2
DA_ANCOMBC2_wrapper <- function(object, assay_name = "counts", fix_formula, rand_formula = NULL) {
  if (!requireNamespace("ANCOMBC2", quietly = TRUE))
    stop("Package 'ANCOMBC2' not installed.")
  y  <- SummarizedExperiment::assay(object, assay_name)
  md <- as.data.frame(SummarizedExperiment::colData(object))
  res <- ANCOMBC2::ancombc2(data = list(abn = y, meta = md),
                            fix_formula = fix_formula, rand_formula = rand_formula,
                            p_adj_method = "BH")
  pv <- res$pW; qv <- res$qW; tt <- res$lfc
  term <- colnames(pv)[1]
  p <- pv[, term]; q <- qv[, term]; lfc <- tt[, term]
  rn <- rownames(pv)
  list(
    pValMat  = data.frame(rawP = setNames(as.numeric(p), rn), adjP = setNames(as.numeric(q), rn)),
    statInfo = data.frame(logFC = setNames(as.numeric(lfc), rn))
  )
}

# MaAsLin2
DA_Maaslin2_wrapper <- function(object, assay_name = "clr", fix_formula) {
  if (!requireNamespace("Maaslin2", quietly = TRUE))
    stop("Package 'Maaslin2' not installed.")
  y  <- SummarizedExperiment::assay(object, assay_name)
  md <- as.data.frame(SummarizedExperiment::colData(object))
  tmp_dir <- tempfile("maaslin2_"); dir.create(tmp_dir)
  
  # Handle fixed_effects parameter - if fix_formula is just a column name, use it directly
  fixed_effects_param <- if(grepl("^[a-zA-Z_][a-zA-Z0-9_]*$", fix_formula)) {
    fix_formula  # Simple column name
  } else {
    all.vars(as.formula(fix_formula))  # Parse formula
  }
  
  fit <- Maaslin2::Maaslin2(input_data = y, input_metadata = md,
                             output = tmp_dir, fixed_effects = fixed_effects_param,
                             normalization = "NONE", transform = "NONE")
  tab <- read.delim(file.path(tmp_dir, "all_results.tsv"), check.names = FALSE)
  eff <- if(grepl("^[a-zA-Z_][a-zA-Z0-9_]*$", fix_formula)) {
    fix_formula  # Simple column name
  } else {
    all.vars(as.formula(fix_formula))[1]  # Parse formula
  }
  sub <- tab[tab$metadata == eff, , drop = FALSE]
  rn  <- sub$feature
  p   <- sub$pval; q <- sub$qval; lfc <- sub$coef
  list(
    pValMat  = data.frame(rawP = setNames(as.numeric(p), rn), adjP = setNames(as.numeric(q), rn)),
    statInfo = data.frame(logFC = setNames(as.numeric(lfc), rn))
  )
}

# Wilcoxon on CLR (two-group)
DA_Wilcoxon_CLR_wrapper <- function(object, assay_name = "clr", group) {
  x  <- SummarizedExperiment::assay(object, assay_name)
  md <- as.data.frame(SummarizedExperiment::colData(object))
  g  <- droplevels(as.factor(md[[group]]))
  stopifnot(nlevels(g) == 2)
  g1 <- levels(g)[1]; g2 <- levels(g)[2]
  p <- vapply(rownames(x), function(feat) {
    stats::wilcox.test(x[feat, g == g1], x[feat, g == g2], exact = FALSE)$p.value
  }, numeric(1))
  lfc <- rowMeans(x[, g == g2, drop = FALSE]) - rowMeans(x[, g == g1, drop = FALSE])
  q <- p.adjust(p, method = "BH")
  list(
    pValMat  = data.frame(rawP = p, adjP = q, row.names = names(p)),
    statInfo = data.frame(logFC = lfc, row.names = names(lfc))
  )
}

# fastANCOM
DA_fastANCOM_wrapper <- function(object, assay_name = "counts", group) {
  if (!requireNamespace("fastANCOM", quietly = TRUE))
    stop("Package 'fastANCOM' not installed.")
  x  <- SummarizedExperiment::assay(object, assay_name)
  md <- as.data.frame(SummarizedExperiment::colData(object))
  g  <- droplevels(as.factor(md[[group]]))
  fit <- fastANCOM::fastANCOM(feature_table = t(x), class = g)
  rn <- rownames(x)
  p  <- fit$p_values[ rn ]
  q  <- stats::p.adjust(p, method = "BH")
  lfc <- fit$effect_size[ rn ]
  list(
    pValMat  = data.frame(rawP = p, adjP = q, row.names = rn),
    statInfo = data.frame(logFC = lfc, row.names = rn)
  )
}

# LinDA
DA_LinDA_wrapper <- function(object, assay_name = "counts", fix_formula) {
  if (!requireNamespace("LinDA", quietly = TRUE))
    stop("Package 'LinDA' not installed.")
  y  <- SummarizedExperiment::assay(object, assay_name)
  md <- as.data.frame(SummarizedExperiment::colData(object))
  fit <- LinDA::linda(W = y, meta = md, formula = fix_formula)
  tab <- fit$output
  term <- if(grepl("^[a-zA-Z_][a-zA-Z0-9_]*$", fix_formula)) {
    fix_formula  # Simple column name
  } else {
    all.vars(as.formula(fix_formula))[1]  # Parse formula
  }
  sub  <- tab[tab$term == term, , drop = FALSE]
  rn   <- sub$taxon
  p    <- sub$pvalue; q <- sub$qvalue; lfc <- sub$coef
  list(
    pValMat  = data.frame(rawP = setNames(as.numeric(p), rn), adjP = setNames(as.numeric(q), rn)),
    statInfo = data.frame(logFC = setNames(as.numeric(lfc), rn))
  )
}

# ----------  wrapper parameter sets -------------------------------------------
wrap <- list()

# ------------------ CORE methods (original 10) ------------------
# Run on real counts, or (experimental) on non-count metrics when --force-counts is set
if (mode %in% c("core", "both") && (tolower(metric_name) == "counts" || force_counts)) {
  wrap$ALDEx2 <- list(fun = DA_ALDEx2,
                      args = list(assay_name = "counts", design = grp_col, contrast = contrast_vec, test = "t"))
  wrap$ANCOM  <- list(fun = DA_ANCOM,
                      args = list(assay_name = "counts", fix_formula = grp_col, contrast = contrast_vec, BC = TRUE))
  wrap$basic  <- list(fun = DA_basic,
                      args = list(assay_name = "counts", test = "t", contrast = contrast_vec))
  wrap$corncob <- list(fun = DA_corncob,
                       args = list(assay_name = "counts",
                                   formula = as.formula(paste("~", grp_col)),
                                   phi.formula = as.formula(paste("~", grp_col)),
                                   formula_null = ~1, phi.formula_null = as.formula(paste("~", grp_col)),
                                   coefficient = paste0(grp_col, levs[2]), test = "Wald"))
  wrap$dearseq <- list(fun = DA_dearseq,
                       args = list(assay_name = "counts", variables2test = grp_col, test = "asymptotic"))
  wrap$DESeq2 <- list(fun = DA_DESeq2,
                      args = list(assay_name = "counts", design = as.formula(paste("~", grp_col)), contrast = contrast_vec))
  wrap$edgeR  <- list(fun = DA_edgeR,
                      args = list(assay_name = "counts", design = as.formula(paste("~", grp_col)), coef = coef_val, norm = "TMM"))
  wrap$limma  <- list(fun = DA_limma,
                      args = list(assay_name = "counts", design = as.formula(paste("~", grp_col)), coef = coef_val, norm = "TMM"))
  wrap$metagenomeSeq <- list(fun = DA_metagenomeSeq,
                             args = list(assay_name = "counts", design = as.formula(paste("~", grp_col)), coef = coef_val, norm = "CSS"))
  wrap$NOISeq <- list(fun = DA_NOISeq,
                      args = list(assay_name = "counts", norm = "tmm", contrast = contrast_vec))
}

# ------------------ EXTENDED methods (new) ------------------
if (mode %in% c("ext", "both")) {
  # Counts-based extended (on counts metric, or forced experimental mode)
  if (tolower(metric_name) == "counts" || force_counts) {
    if (requireNamespace("ancombc", quietly = TRUE)) {
      wrap$ANCOMBC <- list(fun = DA_ANCOMBC_wrapper,
                           args = list(assay_name = "counts", fix_formula = grp_col))
    }
    if (requireNamespace("ANCOMBC2", quietly = TRUE)) {
      wrap$ANCOMBC2 <- list(fun = DA_ANCOMBC2_wrapper,
                            args = list(assay_name = "counts", fix_formula = grp_col))
    }
    if (requireNamespace("LinDA", quietly = TRUE)) {
      wrap$LinDA <- list(fun = DA_LinDA_wrapper,
                         args = list(assay_name = "counts", fix_formula = grp_col))
    }
    if (requireNamespace("fastANCOM", quietly = TRUE)) {
      wrap$fastANCOM <- list(fun = DA_fastANCOM_wrapper,
                             args = list(assay_name = "counts", group = grp_col))
    }
  }
  # Transform-based extended (all metrics)
  wrap$Wilcoxon_CLR <- list(fun = DA_Wilcoxon_CLR_wrapper,
                            args = list(assay_name = "clr", group = grp_col))
  if (requireNamespace("Maaslin2", quietly = TRUE)) {
    wrap$MaAsLin2 <- list(fun = DA_Maaslin2_wrapper,
                          args = list(assay_name = "clr", fix_formula = grp_col))
  }
}

# ----------  run ----------------------------------------------------------------
# Ensure output directory exists even if no methods run
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

results <- lapply(names(wrap), function(nm) run_and_export(wrap[[nm]]$fun, nm, wrap[[nm]]$args))

cat("\nFinished all wrappers. Results written to:\n")
if (length(wrap) == 0) {
  cat("  (No methods ran - check metric type and --force-counts flag)\n")
} else {
  for (nm in names(wrap)) cat("  • ", file.path(out_dir, paste0("DA_", nm, ".csv")), "\n", sep = "")
}

# ----------  session info for reproducibility ---------------------------------
sink(file.path(out_dir, "SESSION_INFO.txt"))
cat(sprintf("Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
print(sessionInfo())
sink()
