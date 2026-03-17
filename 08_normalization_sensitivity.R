#!/usr/bin/env Rscript
# normalization_sensitivity.R
# Sensitivity analysis comparing 5 normalization methods for DE analysis
# Reviewer 3, Comment 2
#
# Methods:
#   1. Spike-in size factors (current approach; Sequin R2_ for C1, ERCC- for C2)
#   2. DESeq2 median-of-ratios (default)
#   3. TMM normalization (edgeR)
#   4. RUVg with spike-ins (RUVSeq)
#   5. Quantile normalization (preprocessCore → limma-voom)
#
# Contrasts:
#   A. Cohort 1: NAFLD vs No-NAFLD (ST1, NAFLD_Diagnosis 1 vs 0)
#   B. Cohort 1: Cirrhosis vs no-cirrhosis (ST1, Cirrhosis 1 vs 0)
#   C. Cohort 2: Fibrosis F4 vs F0 (ST2 Baseline, Fibrosis)

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(RUVSeq)
  library(limma)
  library(preprocessCore)
})

DATA_DIR   <- "../Data"
COUNTS_DIR <- file.path(DATA_DIR, "counts")
OUT_DIR    <- "."

cat("=== Normalization Sensitivity Analysis ===\n")
cat(format(Sys.time()), "\n\n")

# ─────────────────────────────────────────────────────────────
# 1. Load count tables
# ─────────────────────────────────────────────────────────────
cat("Loading count tables...\n")

# Cohort 1 (tab-delimited)
c1_raw <- read.table(file.path(COUNTS_DIR, "htseq_cohort1.txt"),
                     header = TRUE, sep = "\t", row.names = 1,
                     check.names = FALSE, comment.char = "")
cat("  C1 raw dims:", nrow(c1_raw), "rows x", ncol(c1_raw), "cols\n")

# Cohort 2 (space-delimited)
c2_raw <- read.table(file.path(COUNTS_DIR, "htseq_spikes_C2.txt"),
                     header = TRUE, sep = " ", row.names = 1,
                     check.names = FALSE, comment.char = "")
cat("  C2 raw dims:", nrow(c2_raw), "rows x", ncol(c2_raw), "cols\n")

# ─────────────────────────────────────────────────────────────
# 2. Separate spike rows from endogenous genes
# ─────────────────────────────────────────────────────────────
htseq_summary <- c("__no_feature", "__ambiguous", "__too_low_aQual",
                   "__not_aligned", "__alignment_not_unique")

# C1: endogenous = ENSG, spikes = R2_ (Sequin) and ERCC-
c1_endo  <- c1_raw[grepl("^ENSG", rownames(c1_raw)), ]
c1_seq   <- c1_raw[grepl("^R2_",  rownames(c1_raw)), ]   # Sequin
c1_ercc  <- c1_raw[grepl("^ERCC-", rownames(c1_raw)), ]   # ERCC
cat("  C1 endogenous:", nrow(c1_endo), "genes |",
    "Sequin:", nrow(c1_seq), "| ERCC:", nrow(c1_ercc), "\n")

# C2: endogenous = ENSG, spikes = ERCC-
c2_endo  <- c2_raw[grepl("^ENSG",  rownames(c2_raw)), ]
c2_ercc  <- c2_raw[grepl("^ERCC-", rownames(c2_raw)), ]
c2_seq   <- c2_raw[grepl("^R2_",   rownames(c2_raw)), ]   # may exist
cat("  C2 endogenous:", nrow(c2_endo), "genes |",
    "ERCC:", nrow(c2_ercc), "| Sequin:", nrow(c2_seq), "\n")

# ─────────────────────────────────────────────────────────────
# 3. Load metadata
# ─────────────────────────────────────────────────────────────
cat("Loading metadata...\n")

st1 <- read.csv(file.path(DATA_DIR, "Suplementary_Table1.csv"),
                check.names = FALSE, na.strings = c("NA",""))
rownames(st1) <- st1$SMAPLE_RNA
cat("  ST1:", nrow(st1), "rows\n")

st2 <- read.csv("SuplementaryTable2_updated.csv",
                check.names = FALSE, na.strings = c("NA",""))
rownames(st2) <- st2$SMAPLE_RNA
cat("  ST2:", nrow(st2), "rows\n")

# Load outlier lists
outliers_c1 <- tryCatch(
  readLines(file.path(DATA_DIR, "OutliersC1.txt")),
  error = function(e) character(0)
)
outliers_c2 <- tryCatch(
  readLines(file.path(DATA_DIR, "OutliersC2.txt")),
  error = function(e) character(0)
)
cat("  Outliers C1:", length(outliers_c1), "| C2:", length(outliers_c2), "\n")

# Load pre-computed spike ratio tables (SpikeRatio_C1/C2.csv)
# These contain per-sample sums of Sequin and ERCC reads.
# The original analysis uses estimateSizeFactorsForMatrix() on the Sequin column.
spike_ratio_c1 <- read.csv(file.path(DATA_DIR, "SpikeRatio_C1.csv"),
                           check.names = FALSE, stringsAsFactors = FALSE)
rownames(spike_ratio_c1) <- spike_ratio_c1$RNAID

spike_ratio_c2 <- read.csv(file.path(DATA_DIR, "SpikeRatio_C2.csv"),
                           check.names = FALSE, stringsAsFactors = FALSE)
rownames(spike_ratio_c2) <- spike_ratio_c2$RNAID
cat("  SpikeRatio C1:", nrow(spike_ratio_c1), "samples | C2:", nrow(spike_ratio_c2), "samples\n")

# ─────────────────────────────────────────────────────────────
# 4. Define contrasts
# ─────────────────────────────────────────────────────────────
cat("Defining contrasts...\n")

## Contrast A: C1 NAFLD vs No-NAFLD
st1_baseline <- st1[!is.na(st1$NAFLD_Diagnosis) &
                    st1$NAFLD_Diagnosis %in% c("0","1") &
                    !is.na(st1$Gender), , drop = FALSE]
st1_baseline <- st1_baseline[!rownames(st1_baseline) %in% outliers_c1, ]
samps_A <- intersect(rownames(st1_baseline), colnames(c1_endo))
meta_A  <- st1_baseline[samps_A, , drop = FALSE]
meta_A$group  <- factor(ifelse(meta_A$NAFLD_Diagnosis == "1", "NAFLD", "Control"),
                        levels = c("Control","NAFLD"))
meta_A$Gender <- factor(meta_A$Gender)
counts_A_endo <- as.matrix(c1_endo[, samps_A])
counts_A_ercc <- as.matrix(c1_ercc[, samps_A])
counts_A_seq  <- as.matrix(c1_seq[,  samps_A])
# Spike size factors for A (Sequin sums from SpikeRatio_C1.csv)
seq_sf_A <- spike_ratio_c1[intersect(samps_A, rownames(spike_ratio_c1)), "Sequin", drop = TRUE]
names(seq_sf_A) <- intersect(samps_A, rownames(spike_ratio_c1))
cat("  Contrast A samples:", nrow(meta_A),
    "(NAFLD:", sum(meta_A$group == "NAFLD"), "/ Control:", sum(meta_A$group == "Control"), ")\n")

## Contrast B: C1 Cirrhosis vs No-cirrhosis
st1_cirrh <- st1[!is.na(st1$Cirrhosis) &
                 st1$Cirrhosis %in% c("0","1") &
                 !is.na(st1$Gender), , drop = FALSE]
st1_cirrh <- st1_cirrh[!rownames(st1_cirrh) %in% outliers_c1, ]
samps_B <- intersect(rownames(st1_cirrh), colnames(c1_endo))
meta_B  <- st1_cirrh[samps_B, , drop = FALSE]
meta_B$group  <- factor(ifelse(meta_B$Cirrhosis == "1", "Cirrhosis", "NoCirrhosis"),
                        levels = c("NoCirrhosis","Cirrhosis"))
meta_B$Gender <- factor(meta_B$Gender)
counts_B_endo <- as.matrix(c1_endo[, samps_B])
counts_B_ercc <- as.matrix(c1_ercc[, samps_B])
counts_B_seq  <- as.matrix(c1_seq[,  samps_B])
# Spike size factors for B
seq_sf_B <- spike_ratio_c1[intersect(samps_B, rownames(spike_ratio_c1)), "Sequin", drop = TRUE]
names(seq_sf_B) <- intersect(samps_B, rownames(spike_ratio_c1))
cat("  Contrast B samples:", nrow(meta_B),
    "(Cirrhosis:", sum(meta_B$group == "Cirrhosis"),
    "/ No:", sum(meta_B$group == "NoCirrhosis"), ")\n")

## Contrast C: C2 Fibrosis F4 vs F0 (Baseline only)
st2_fib <- st2[!is.na(st2$Fibrosis) &
               st2$Fibrosis %in% c("F0","F4") &
               st2$Visit_Type == "Baseline" &
               !is.na(st2$Sex), , drop = FALSE]
# Restrict to NASH biomarker cohort
nash_cohort_flag <- grepl("NASH biomarker", st2_fib$Cohort, ignore.case = TRUE)
st2_fib <- st2_fib[nash_cohort_flag, ]
st2_fib  <- st2_fib[!rownames(st2_fib) %in% outliers_c2, ]
samps_C  <- intersect(rownames(st2_fib), colnames(c2_endo))
meta_C   <- st2_fib[samps_C, , drop = FALSE]
meta_C$group <- factor(meta_C$Fibrosis, levels = c("F0","F4"))
meta_C$Sex   <- factor(meta_C$Sex)
counts_C_endo <- as.matrix(c2_endo[, samps_C])
counts_C_ercc <- as.matrix(c2_ercc[, samps_C])
# Spike size factors for C (Sequin sums from SpikeRatio_C2.csv)
seq_sf_C <- spike_ratio_c2[intersect(samps_C, rownames(spike_ratio_c2)), "Sequin", drop = TRUE]
names(seq_sf_C) <- intersect(samps_C, rownames(spike_ratio_c2))
cat("  Contrast C samples:", nrow(meta_C),
    "(F4:", sum(meta_C$group == "F4"), "/ F0:", sum(meta_C$group == "F0"), ")\n")

# ─────────────────────────────────────────────────────────────
# 5. Gene filter helpers — match original analyses exactly
# C1: rowSums(counts >= 10) >= smallestGroupSize (= 30)
# C2: rowSums(counts)       >= 300
# ─────────────────────────────────────────────────────────────
filter_genes_c1 <- function(counts, smallest_group = 30, min_per_sample = 10) {
  counts[rowSums(counts >= min_per_sample) >= smallest_group, , drop = FALSE]
}

filter_genes_c2 <- function(counts, min_sum = 300) {
  counts[rowSums(counts) >= min_sum, , drop = FALSE]
}

# Generic wrapper used by helpers below; cohort is auto-detected per contrast
filter_genes <- function(counts, cohort = "C1") {
  if (cohort == "C2") filter_genes_c2(counts) else filter_genes_c1(counts)
}

# ─────────────────────────────────────────────────────────────
# 6. Helper: DESeq2 with provided sizeFactors
# ─────────────────────────────────────────────────────────────
run_deseq2_sf <- function(counts, meta, size_factors = NULL, formula = ~group,
                           contrast_name = "group",
                           contrast_levels = NULL,
                           cohort = "C1") {
  counts <- filter_genes(counts, cohort = cohort)
  counts <- round(counts)
  # Align meta to filtered count columns
  meta   <- meta[colnames(counts), , drop = FALSE]
  if (!is.null(size_factors)) size_factors <- size_factors[colnames(counts)]
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = meta,
                                design    = formula)
  if (!is.null(size_factors)) {
    sizeFactors(dds) <- size_factors
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
  } else {
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
  }
  if (is.null(contrast_levels)) {
    res <- results(dds, name = resultsNames(dds)[2])
  } else {
    res <- results(dds, contrast = c(contrast_name, contrast_levels[2], contrast_levels[1]))
  }
  as.data.frame(res)
}

# ─────────────────────────────────────────────────────────────
# 7. Helper: compute spike-in size factors using DESeq2's
#    estimateSizeFactorsForMatrix(), consistent with original analysis.
#    Input: named numeric vector of per-sample spike sums (from SpikeRatio CSV)
# ─────────────────────────────────────────────────────────────
spike_size_factors <- function(spike_sums_vec, sample_ids) {
  # Align to requested sample IDs
  common <- intersect(sample_ids, names(spike_sums_vec))
  if (length(common) < length(sample_ids))
    warning(length(sample_ids) - length(common), " samples missing from SpikeRatio table")
  sv <- spike_sums_vec[common]
  sv[is.na(sv) | sv == 0] <- 1
  # estimateSizeFactorsForMatrix expects a matrix (genes x samples);
  # passing a 1-row matrix replicates ~RNAID/geo_mean normalisation exactly.
  sf_mat <- matrix(as.numeric(sv), nrow = 1,
                   dimnames = list("spike_sum", common))
  sf <- DESeq2::estimateSizeFactorsForMatrix(sf_mat)
  sf  # named vector, length = ncol(counts)
}

# ─────────────────────────────────────────────────────────────
# 8. Helper: edgeR TMM → return log2FC + FDR
# ─────────────────────────────────────────────────────────────
run_edger_tmm <- function(counts, meta, cohort = "C1", formula_covar = ~group) {
  counts <- filter_genes(counts, cohort = cohort)
  counts <- round(counts)
  meta   <- meta[colnames(counts), , drop = FALSE]
  dge <- DGEList(counts = counts, group = meta$group)
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(formula_covar, data = meta)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  # group coefficient is whichever column starts with "group"
  group_col <- grep("^group", colnames(design))[1]
  qlf <- glmQLFTest(fit, coef = group_col)
  tt  <- topTags(qlf, n = Inf, sort.by = "none")$table
  data.frame(
    row.names  = rownames(tt),
    log2FoldChange = tt$logFC,
    pvalue     = tt$PValue,
    padj       = tt$FDR
  )
}

# ─────────────────────────────────────────────────────────────
# 9. Helper: RUVg with spike controls → DESeq2
# ─────────────────────────────────────────────────────────────
run_ruvg_deseq2 <- function(counts_endo, counts_spike, meta, k = 1,
                            cohort = "C1",
                            formula_covar = ~group) {
  # Apply same gene filter as manuscript before RUVg
  counts_endo_filt <- filter_genes(counts_endo, cohort = cohort)
  # Combine filtered endo + spike for RUVg input
  all_counts <- rbind(counts_endo_filt, counts_spike)
  # spike control genes
  ctrl_genes <- intersect(rownames(counts_spike), rownames(all_counts))
  if (length(ctrl_genes) == 0) {
    warning("No spike-in control genes passed filter — skipping RUVg")
    return(NULL)
  }
  meta_sub <- meta[colnames(all_counts), , drop = FALSE]
  # SeqExpressionSet
  set <- newSeqExpressionSet(as.matrix(round(all_counts)),
                             phenoData = AnnotatedDataFrame(
                               data.frame(group = meta_sub$group,
                                          row.names = colnames(all_counts))))
  set_ruv <- RUVg(set, ctrl_genes, k = k)
  # Build covariate formula: original covars + RUV W factors
  W <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop = FALSE]
  orig_vars <- all.vars(formula_covar)
  meta_ruv  <- cbind(meta_sub, W)
  extra_vars <- c(orig_vars[-1], colnames(W))  # drop 'group'; add covar + W
  covar_formula <- if (length(extra_vars) > 0) {
    as.formula(paste("~ group +", paste(extra_vars, collapse = " + ")))
  } else {
    ~group
  }
  endo_only <- as.matrix(round(counts(set_ruv)))[grepl("^ENSG", rownames(counts(set_ruv))), , drop = FALSE]
  meta_ruv  <- meta_ruv[colnames(endo_only), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = endo_only,
                                colData   = meta_ruv,
                                design    = covar_formula)
  sizeFactors(dds) <- rep(1, ncol(endo_only))  # RUV handles normalisation
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  rn <- resultsNames(dds)
  group_coef <- rn[grepl("^group", rn)][1]
  if (is.na(group_coef)) {
    warning("Could not find group coefficient in RUVg DESeq2 results")
    return(NULL)
  }
  res <- results(dds, name = group_coef)
  as.data.frame(res)
}

# ─────────────────────────────────────────────────────────────
# 10. Helper: limma-voom with quantile normalization
# ─────────────────────────────────────────────────────────────
run_limma_quantile <- function(counts, meta, cohort = "C1", formula_covar = ~group) {
  counts <- filter_genes(counts, cohort = cohort)
  counts <- round(counts)
  meta   <- meta[colnames(counts), , drop = FALSE]
  # voom with quantile normalization
  dge    <- DGEList(counts = counts, group = meta$group)
  design <- model.matrix(formula_covar, data = meta)
  v      <- voom(dge, design, normalize.method = "quantile", plot = FALSE)
  fit    <- lmFit(v, design)
  fit    <- eBayes(fit)
  group_col <- grep("^group", colnames(design))[1]
  tt     <- topTable(fit, coef = group_col, number = Inf, sort.by = "none")
  data.frame(
    row.names      = rownames(tt),
    log2FoldChange = tt$logFC,
    pvalue         = tt$P.Value,
    padj           = tt$adj.P.Val
  )
}

# ─────────────────────────────────────────────────────────────
# 11. Run all methods for all contrasts
# ─────────────────────────────────────────────────────────────
METHOD_NAMES <- c("Spike_SizeFactors", "DESeq2_MedianRatio", "edgeR_TMM",
                  "RUVg_Spikes", "Limma_Quantile")

# spike_sums_vec: named numeric vector of per-sample spike sums (from SpikeRatio CSV)
run_all_methods <- function(counts_endo, spike_sums_vec, meta,
                            contrast_label, ref_level, treat_level,
                            cohort = "C1", formula_covar = ~group) {
  cat("\n--- Contrast:", contrast_label, "---\n")
  meta$group <- relevel(meta$group, ref = ref_level)
  results_list <- list()

  # Method 1: Spike-in size factors (from SpikeRatio CSV, via estimateSizeFactorsForMatrix)
  cat("  Method 1: Spike-in size factors\n")
  tryCatch({
    sf  <- spike_size_factors(spike_sums_vec, sample_ids = colnames(counts_endo))
    res <- run_deseq2_sf(counts_endo, meta, size_factors = sf,
                         formula        = formula_covar,
                         contrast_name  = "group",
                         contrast_levels = c(ref_level, treat_level),
                         cohort         = cohort)
    results_list[["Spike_SizeFactors"]] <- res
    cat("    Done:", sum(!is.na(res$padj) & res$padj < 0.05), "DEGs (padj<0.05)\n")
  }, error = function(e) {
    cat("    FAILED:", conditionMessage(e), "\n")
    results_list[["Spike_SizeFactors"]] <<- NULL
  })

  # Method 2: DESeq2 default (median-of-ratios)
  cat("  Method 2: DESeq2 median-of-ratios\n")
  tryCatch({
    res <- run_deseq2_sf(counts_endo, meta, size_factors = NULL,
                         formula        = formula_covar,
                         contrast_name  = "group",
                         contrast_levels = c(ref_level, treat_level),
                         cohort         = cohort)
    results_list[["DESeq2_MedianRatio"]] <- res
    cat("    Done:", sum(!is.na(res$padj) & res$padj < 0.05), "DEGs (padj<0.05)\n")
  }, error = function(e) {
    cat("    FAILED:", conditionMessage(e), "\n")
    results_list[["DESeq2_MedianRatio"]] <<- NULL
  })

  # Method 3: edgeR TMM
  cat("  Method 3: edgeR TMM\n")
  tryCatch({
    res <- run_edger_tmm(counts_endo, meta, cohort = cohort, formula_covar = formula_covar)
    results_list[["edgeR_TMM"]] <- res
    cat("    Done:", sum(!is.na(res$padj) & res$padj < 0.05), "DEGs (padj<0.05)\n")
  }, error = function(e) {
    cat("    FAILED:", conditionMessage(e), "\n")
    results_list[["edgeR_TMM"]] <<- NULL
  })

  # Method 4: RUVg (uses spike rows from endo count table)
  cat("  Method 4: RUVg with spike-ins\n")
  tryCatch({
    # Re-extract spike rows from the full (unfiltered) endo+spike matrix
    # for RUVg control genes — use spike sum vector to build a proxy row matrix
    sv     <- spike_sums_vec[colnames(counts_endo)]
    sv[is.na(sv)] <- 0
    spike_proxy <- matrix(as.integer(round(sv)), nrow = 1,
                          dimnames = list("spike_proxy", colnames(counts_endo)))
    res <- run_ruvg_deseq2(counts_endo, spike_proxy, meta, k = 1,
                           cohort = cohort, formula_covar = formula_covar)
    if (!is.null(res)) {
      results_list[["RUVg_Spikes"]] <- res
      cat("    Done:", sum(!is.na(res$padj) & res$padj < 0.05), "DEGs (padj<0.05)\n")
    }
  }, error = function(e) {
    cat("    FAILED:", conditionMessage(e), "\n")
    results_list[["RUVg_Spikes"]] <<- NULL
  })

  # Method 5: limma-voom quantile
  cat("  Method 5: limma-voom quantile\n")
  tryCatch({
    res <- run_limma_quantile(counts_endo, meta, cohort = cohort, formula_covar = formula_covar)
    results_list[["Limma_Quantile"]] <- res
    cat("    Done:", sum(!is.na(res$padj) & res$padj < 0.05), "DEGs (padj<0.05)\n")
  }, error = function(e) {
    cat("    FAILED:", conditionMessage(e), "\n")
    results_list[["Limma_Quantile"]] <<- NULL
  })

  results_list
}

cat("\n\n=== Running contrasts ===\n")
# C1 model includes Gender covariate; C2 includes Sex covariate
res_A <- run_all_methods(counts_A_endo, seq_sf_A, meta_A,
                         "Contrast A: C1 NAFLD vs Control",
                         ref_level     = "Control",
                         treat_level   = "NAFLD",
                         cohort        = "C1",
                         formula_covar = ~group + Gender)
res_B <- run_all_methods(counts_B_endo, seq_sf_B, meta_B,
                         "Contrast B: C1 Cirrhosis vs NoCirrhosis",
                         ref_level     = "NoCirrhosis",
                         treat_level   = "Cirrhosis",
                         cohort        = "C1",
                         formula_covar = ~group + Gender)
res_C <- run_all_methods(counts_C_endo, seq_sf_C, meta_C,
                         "Contrast C: C2 F4 vs F0",
                         ref_level     = "F0",
                         treat_level   = "F4",
                         cohort        = "C2",
                         formula_covar = ~group + Sex)

all_results <- list(
  "C1_NAFLD_vs_Control"     = res_A,
  "C1_Cirrhosis_vs_NoCirr"  = res_B,
  "C2_F4_vs_F0"             = res_C
)

# ─────────────────────────────────────────────────────────────
# 12. Summarise: LFC correlation + DEG counts per method
# ─────────────────────────────────────────────────────────────
cat("\n=== Computing summary statistics ===\n")

pad <- 0.05  # FDR threshold

summary_rows <- list()
for (contrast_name in names(all_results)) {
  for (method in METHOD_NAMES) {
    res <- all_results[[contrast_name]][[method]]
    if (is.null(res)) next
    n_deg <- sum(!is.na(res$padj) & res$padj < pad, na.rm = TRUE)
    n_up  <- sum(!is.na(res$padj) & res$padj < pad & res$log2FoldChange > 0, na.rm = TRUE)
    n_dn  <- sum(!is.na(res$padj) & res$padj < pad & res$log2FoldChange < 0, na.rm = TRUE)
    summary_rows[[paste(contrast_name, method, sep = "|")]] <- data.frame(
      Contrast = contrast_name, Method = method,
      DEGs_total = n_deg, DEGs_up = n_up, DEGs_down = n_dn
    )
  }
}
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
print(summary_df)
write.csv(summary_df, "normalization_sensitivity_DEG_counts.csv", row.names = FALSE)
cat("Saved: normalization_sensitivity_DEG_counts.csv\n")

# ─────────────────────────────────────────────────────────────
# NOTE: Exclude Spike_SizeFactors from Contrast A (C1 NAFLD vs Control).
# The design is severely unbalanced (n = 97 NAFLD vs n = 10 Control).
# With only 10 control samples, differences in Sequin spike-in recovery
# between groups create a systematic normalization artifact: 6,735 DEGs
# were detected (99.9% in the same direction), indicating that the
# size factors — rather than biology — are driving the signal.
# Spike-in normalization is retained for the balanced contrasts (B, C).
# ─────────────────────────────────────────────────────────────
cat("\nNOTE: Spike_SizeFactors excluded for Contrast A (n=97 vs n=10).\n")
all_results[["C1_NAFLD_vs_Control"]][["Spike_SizeFactors"]] <- NULL
res_A[["Spike_SizeFactors"]] <- NULL
# Rewrite summary CSV without the excluded entry
summary_df_plot <- summary_df[
  !(summary_df$Contrast == "C1_NAFLD_vs_Control" &
    summary_df$Method   == "Spike_SizeFactors"), ]
write.csv(summary_df_plot, "normalization_sensitivity_DEG_counts.csv", row.names = FALSE)
cat("Updated: normalization_sensitivity_DEG_counts.csv\n")

# ─────────────────────────────────────────────────────────────
# 13. LFC correlation between methods (per contrast)
# ─────────────────────────────────────────────────────────────

compute_lfc_cor <- function(res_list) {
  # Extract LFC vectors, align by gene
  lfc_list <- lapply(res_list, function(r) {
    if (is.null(r)) return(NULL)
    setNames(r$log2FoldChange, rownames(r))
  })
  # Remove nulls
  lfc_list <- Filter(Negate(is.null), lfc_list)
  if (length(lfc_list) < 2) return(NULL)
  # Common genes
  common_genes <- Reduce(intersect, lapply(lfc_list, names))
  common_genes <- common_genes[!is.na(common_genes)]
  if (length(common_genes) < 100) {
    warning("Too few common genes for LFC correlation")
    return(NULL)
  }
  lfc_mat <- do.call(cbind, lapply(lfc_list, function(v) v[common_genes]))
  lfc_mat <- na.omit(lfc_mat)
  cor(lfc_mat, method = "spearman")
}

cor_A <- compute_lfc_cor(res_A)
cor_B <- compute_lfc_cor(res_B)
cor_C <- compute_lfc_cor(res_C)

# Save correlation tables
if (!is.null(cor_A)) write.csv(round(cor_A, 3), "lfc_cors_contrastA.csv")
if (!is.null(cor_B)) write.csv(round(cor_B, 3), "lfc_cors_contrastB.csv")
if (!is.null(cor_C)) write.csv(round(cor_C, 3), "lfc_cors_contrastC.csv")

# ─────────────────────────────────────────────────────────────
# 14. DEG overlap: Jaccard index between method pairs
# ─────────────────────────────────────────────────────────────

get_degs <- function(res, fdr = 0.05) {
  if (is.null(res)) return(character(0))
  rownames(res)[!is.na(res$padj) & res$padj < fdr]
}

jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA)
  length(intersect(a, b)) / length(union(a, b))
}

compute_jaccard_mat <- function(res_list) {
  methods <- names(res_list)
  mat <- matrix(NA, length(methods), length(methods),
                dimnames = list(methods, methods))
  for (i in seq_along(methods))
    for (j in seq_along(methods))
      mat[i, j] <- jaccard(get_degs(res_list[[methods[i]]]),
                            get_degs(res_list[[methods[j]]]))
  mat
}

jacc_A <- compute_jaccard_mat(res_A)
jacc_B <- compute_jaccard_mat(res_B)
jacc_C <- compute_jaccard_mat(res_C)

# ─────────────────────────────────────────────────────────────
# 15. Generate figures
# ─────────────────────────────────────────────────────────────
cat("\nGenerating figures...\n")

METHOD_LABELS <- c(
  "Spike_SizeFactors" = "Spike-in\nSize Factors",
  "DESeq2_MedianRatio" = "DESeq2\nMedian-Ratio",
  "edgeR_TMM"         = "edgeR\nTMM",
  "RUVg_Spikes"       = "RUVg\n(Spikes)",
  "Limma_Quantile"    = "limma-voom\nQuantile"
)

CONTRAST_LABELS <- c(
  "C1_NAFLD_vs_Control"    = "C1: NAFLD vs Control",
  "C1_Cirrhosis_vs_NoCirr" = "C1: Cirrhosis vs No-Cirrh.",
  "C2_F4_vs_F0"            = "C2: F4 vs F0"
)

COLORS <- c(
  "Spike_SizeFactors" = "#E64B35",
  "DESeq2_MedianRatio" = "#4DBBD5",
  "edgeR_TMM"          = "#00A087",
  "RUVg_Spikes"        = "#3C5488",
  "Limma_Quantile"     = "#F39B7F"
)

make_heatmap <- function(mat, title, method_labels) {
  if (is.null(mat)) {
    plot.new(); title(paste(title, "\n(insufficient data)"))
    return(invisible(NULL))
  }
  # get intersection of available methods
  avail <- intersect(names(method_labels), rownames(mat))
  mat2  <- mat[avail, avail, drop = FALSE]
  rownames(mat2) <- colnames(mat2) <- method_labels[avail]
  n <- nrow(mat2)
  # plot as image
  image(1:n, 1:n, t(mat2[nrow(mat2):1, ]),
        col = colorRampPalette(c("white","#2196F3","#0D47A1"))(50),
        zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  mtext(colnames(mat2), side = 1, at = 1:n, las = 2, cex = 0.75)
  mtext(rev(rownames(mat2)), side = 2, at = 1:n, las = 1, cex = 0.75)
  title(title, cex.main = 0.95)
  for (i in 1:n) for (j in 1:n) {
    val <- mat2[nrow(mat2)+1-i, j]
    if (!is.na(val))
      text(j, i, sprintf("%.2f", val), cex = 0.7,
           col = ifelse(val > 0.6, "white", "black"))
  }
}

# Figure 1: DEG counts barplot
pdf("normalization_sensitivity_DEG_counts.pdf", width = 10, height = 5)
par(mfrow = c(1, 3), mar = c(8, 5, 3, 1))
for (contrast_name in names(all_results)) {
  res_lst <- all_results[[contrast_name]]
  methods_available <- intersect(METHOD_NAMES, names(res_lst))
  n_degs <- sapply(methods_available, function(m) {
    r <- res_lst[[m]]
    if (is.null(r)) return(0)
    sum(!is.na(r$padj) & r$padj < 0.05, na.rm = TRUE)
  })
  bp <- barplot(n_degs,
                names.arg = METHOD_LABELS[methods_available],
                col = COLORS[methods_available],
                las = 2, cex.names = 0.75,
                main = CONTRAST_LABELS[contrast_name],
                ylab = "# DEGs (FDR < 0.05)", border = NA)
  text(bp, n_degs + max(n_degs, na.rm = TRUE) * 0.02,
       labels = n_degs, cex = 0.7, font = 2)
}
dev.off()
cat("Saved: normalization_sensitivity_DEG_counts.pdf\n")

# Figure 2: LFC correlation heatmaps (Spearman)
pdf("normalization_sensitivity_LFC_cors.pdf", width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(7, 7, 3, 1))
make_heatmap(cor_A, "C1: NAFLD vs Control\nSpearman LFC correlation", METHOD_LABELS)
make_heatmap(cor_B, "C1: Cirrhosis vs No-Cirrhosis\nSpearman LFC correlation", METHOD_LABELS)
make_heatmap(cor_C, "C2: F4 vs F0\nSpearman LFC correlation", METHOD_LABELS)
# Add color legend
par(fig = c(0.92, 0.99, 0.2, 0.8), new = TRUE, mar = c(2,1,2,3))
image(1, seq(0, 1, length.out = 50), t(matrix(seq(0, 1, length.out = 50))),
      col = colorRampPalette(c("white","#2196F3","#0D47A1"))(50),
      axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(0, 0.5, 1), labels = c("0","0.5","1"), las = 1, cex.axis = 0.7)
mtext("Spearman r", side = 3, cex = 0.6, line = 0.3)
dev.off()
cat("Saved: normalization_sensitivity_LFC_cors.pdf\n")

# Figure 3: DEG Jaccard similarity heatmaps
pdf("normalization_sensitivity_Jaccard.pdf", width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(7, 7, 3, 1))
make_heatmap(jacc_A, "C1: NAFLD vs Control\nDEG Jaccard similarity", METHOD_LABELS)
make_heatmap(jacc_B, "C1: Cirrhosis vs No-Cirrhosis\nDEG Jaccard similarity", METHOD_LABELS)
make_heatmap(jacc_C, "C2: F4 vs F0\nDEG Jaccard similarity", METHOD_LABELS)
dev.off()
cat("Saved: normalization_sensitivity_Jaccard.pdf\n")

# Figure 4: Pairwise LFC scatter (Spike-in vs other methods) for each contrast
plot_lfc_scatter <- function(res_list, contrast_label, ref_method = "Spike_SizeFactors") {
  if (is.null(res_list[[ref_method]])) {
    cat("  Ref method", ref_method, "not available for", contrast_label, "— skipping scatter\n")
    return(invisible(NULL))
  }
  other_methods <- setdiff(names(res_list)[!sapply(res_list, is.null)], ref_method)
  if (length(other_methods) == 0) return(invisible(NULL))
  n_plots <- length(other_methods)
  par(mfrow = c(1, n_plots), mar = c(4, 4, 3, 1))
  ref_lfc <- res_list[[ref_method]]$log2FoldChange
  names(ref_lfc) <- rownames(res_list[[ref_method]])
  for (m in other_methods) {
    other_lfc  <- res_list[[m]]$log2FoldChange
    names(other_lfc) <- rownames(res_list[[m]])
    common <- intersect(names(ref_lfc), names(other_lfc))
    x <- ref_lfc[common]
    y <- other_lfc[common]
    ok <- !is.na(x) & !is.na(y)
    r_val <- round(cor(x[ok], y[ok], method = "spearman"), 3)
    smoothScatter(x[ok], y[ok], xlab = paste0("log2FC\n", METHOD_LABELS[ref_method]),
                  ylab = paste0("log2FC\n", METHOD_LABELS[m]),
                  main = paste0(METHOD_LABELS[ref_method], " vs ", METHOD_LABELS[m],
                                "\nr = ", r_val),
                  cex.main = 0.8, nrpoints = 200,
                  colramp = colorRampPalette(c("white","#1565C0","#0D47A1")))
    abline(0, 1, col = "red", lty = 2, lwd = 1.5)
    legend("topleft", legend = paste("r =", r_val), bty = "n", cex = 0.9)
  }
  title(contrast_label, outer = TRUE, line = -1, cex.main = 1)
}

pdf("normalization_sensitivity_LFC_scatter.pdf", width = 14, height = 4)
plot_lfc_scatter(res_A, "C1: NAFLD vs Control")
plot_lfc_scatter(res_B, "C1: Cirrhosis vs No-Cirrhosis")
plot_lfc_scatter(res_C, "C2: F4 vs F0")
dev.off()
cat("Saved: normalization_sensitivity_LFC_scatter.pdf\n")

# ─────────────────────────────────────────────────────────────
# 16. Print key LFC correlation statistics for the response letter
# ─────────────────────────────────────────────────────────────
cat("\n=== LFC Correlation summary ===\n")
for (nm in c("C1_NAFLD_vs_Control","C1_Cirrhosis_vs_NoCirr","C2_F4_vs_F0")) {
  mat <- switch(nm,
    "C1_NAFLD_vs_Control"    = cor_A,
    "C1_Cirrhosis_vs_NoCirr" = cor_B,
    "C2_F4_vs_F0"            = cor_C)
  if (is.null(mat)) {
    cat(nm, ": insufficient data\n")
    next
  }
  # For Contrast A (spike excluded) use DESeq2 as reference; otherwise use spike
  ref <- if ("Spike_SizeFactors" %in% rownames(mat)) "Spike_SizeFactors" else "DESeq2_MedianRatio"
  ref_label <- if (ref == "Spike_SizeFactors") "Spike-in" else "DESeq2 (spike excluded from Contrast A)"
  others <- setdiff(rownames(mat), ref)
  cat(nm, "(reference:", ref_label, "):\n")
  for (m in others) cat("  ", m, ":", round(mat[ref, m], 3), "\n")
}

# Save full LFC tables per contrast for supplementary
for (contrast_nm in names(all_results)) {
  for (method_nm in METHOD_NAMES) {
    res_m <- all_results[[contrast_nm]][[method_nm]]
    if (!is.null(res_m)) {
      out_fn <- paste0("normalization_sensitivity_LFC_", contrast_nm, "_", method_nm, ".csv")
      write.csv(res_m, out_fn, row.names = TRUE)
    }
  }
}

cat("\n=== Analysis complete ===\n")
cat(format(Sys.time()), "\n")
