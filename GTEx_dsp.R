###############################################################################
# Helper functions for GTEx poison-exon differential splicing analysis.
#
# Sourced by GTEx_differential_splicing_analysis.R.
#
# Dependencies:
#   - matrixStats
#   - auto_deSV.R  (provides buildDEmodel(); must be in the working directory)
###############################################################################

library(matrixStats)
source("./auto_deSV.R")


# Run a tissue-vs-tissue differential PSI test for one pair of tissues.
# Writes beta / stderr / fdr / pval tables to output_dir.
test_one_pair <- function(x, meta, tissue_col, cutoff, output_dir, folder, pe_ids = NULL){
  tis1 <- x[1]
  tis2 <- x[2]
  s1_1 <- read.table(paste0(folder, tis1, ".poison_exon.numerator.out.tab.gz"),   header = TRUE, row.names = 1, check.names = FALSE)
  s1_2 <- read.table(paste0(folder, tis1, ".poison_exon.denominator.out.tab.gz"), header = TRUE, row.names = 1, check.names = FALSE)
  s2_1 <- read.table(paste0(folder, tis2, ".poison_exon.numerator.out.tab.gz"),   header = TRUE, row.names = 1, check.names = FALSE)
  s2_2 <- read.table(paste0(folder, tis2, ".poison_exon.denominator.out.tab.gz"), header = TRUE, row.names = 1, check.names = FALSE)
  grp <- meta[which(meta[[tissue_col]] %in% c(tis1, tis2)), ]
  grp[[tissue_col]] <- factor(as.vector(grp[[tissue_col]]), levels = c(tis1, tis2))

  # Restrict to filtered PE set if provided
  if (!is.null(pe_ids)) {
    common_pe <- intersect(pe_ids, rownames(s1_1))
    s1_1 <- s1_1[common_pe, , drop = FALSE]
    s1_2 <- s1_2[common_pe, , drop = FALSE]
    s2_1 <- s2_1[common_pe, , drop = FALSE]
    s2_2 <- s2_2[common_pe, , drop = FALSE]
  }

  s1_1 <- as.matrix(s1_1)
  s1_2 <- as.matrix(s1_2)
  s2_1 <- as.matrix(s2_1)
  s2_2 <- as.matrix(s2_2)

  rmed1_1 <- rowMedians(s1_1)
  rmed1_2 <- rowMedians(s1_2)
  rmed2_1 <- rowMedians(s2_1)
  rmed2_2 <- rowMedians(s2_2)

  # Keep PEs that are expressed (median numerator OR denominator >= cutoff) in BOTH tissues
  exp <- rep(FALSE, length(rmed1_1))
  for (i in 1:length(exp)){
    if ((rmed1_1[i] >= cutoff || rmed1_2[i] >= cutoff) && (rmed2_1[i] >= cutoff || rmed2_2[i] >= cutoff)){
      exp[i] <- TRUE
    }
  }

  mat1 <- cbind(s1_1, s2_1)
  mat2 <- cbind(s1_2, s2_2)
  # Ensure metadata and PSI matrix have the same samples in the same order
  common_samples <- intersect(rownames(grp), colnames(mat1))
  grp  <- grp[common_samples, ]
  mat1 <- mat1[, common_samples]
  mat2 <- mat2[, common_samples]
  stopifnot(identical(colnames(mat1), rownames(grp)))
  mat1 <- mat1[which(exp), ]
  mat2 <- mat2[which(exp), ]

  # PSI with pseudocount to avoid 0/0; clip to (0,1) so qlogis is finite
  psi <- mat1 / (mat1 + mat2 + 0.01)
  psi <- pmax(pmin(psi, 1 - 1e-6), 1e-6)
  mat <- qlogis(psi)  # logit transform of PSI

  method1 <- "glm"
  ms  <- tissue_col
  mf1 <- paste0(tis2, "_VS_", tis1)
  test <- buildDEmodel(feature_matrix = mat, sampleTable = grp,
                       method = method1, estimate_offset = FALSE, offset_str = NULL,
                       apply_sva = FALSE, test_terms = ms,
                       remove_terms = NULL, random_terms = NULL,
                       distro_family = "gaussian")

  write.table(test@beta,   paste0(output_dir, "/model_params_beta_",   method1, "_", ms, "_", mf1, ".txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(test@stderr, paste0(output_dir, "/model_params_stderr_", method1, "_", ms, "_", mf1, ".txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(test@fdr,    paste0(output_dir, "/model_params_fdr_",    method1, "_", ms, "_", mf1, ".txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(test@pval,   paste0(output_dir, "/model_params_pval_",   method1, "_", ms, "_", mf1, ".txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
}


# For a single tissue, return a 0/1 expression call per PE id:
# 1 if median numerator OR denominator >= cutoff, 0 otherwise.
expression_one_tissue <- function(x, ids, cutoff, folder){
  s1_1 <- read.table(paste0(folder, x, ".poison_exon.numerator.out.tab.gz"),   header = TRUE, row.names = 1, check.names = FALSE)
  s1_2 <- read.table(paste0(folder, x, ".poison_exon.denominator.out.tab.gz"), header = TRUE, row.names = 1, check.names = FALSE)

  s1_1 <- as.matrix(s1_1)
  s1_2 <- as.matrix(s1_2)
  rmed1_1 <- rowMedians(s1_1)
  rmed1_2 <- rowMedians(s1_2)
  names(rmed1_1) <- rownames(s1_1)
  names(rmed1_2) <- rownames(s1_2)

  tmp1_1 <- rep(0, length(ids))
  tmp1_2 <- rep(0, length(ids))
  names(tmp1_1) <- ids
  names(tmp1_2) <- ids
  common1 <- intersect(names(rmed1_1), ids)
  common2 <- intersect(names(rmed1_2), ids)
  tmp1_1[common1] <- rmed1_1[common1]
  tmp1_2[common2] <- rmed1_2[common2]

  res <- (pmax(tmp1_1, tmp1_2, na.rm = TRUE) >= cutoff)
  return(as.numeric(res))
}


# Return PE ids that are expressed in at least one target_tissue sample
# (e.g. any "Brain" sub-tissue) and in zero non-target tissues.
find_tissue_uniq <- function(tissues, ids, target_tissue, folder){
  mat <- sapply(tissues, FUN = expression_one_tissue, cutoff = 3, ids = ids, folder = folder)
  rownames(mat) <- ids
  idx <- grep(target_tissue, tissues)
  rsums1 <- rowSums(mat[, idx,  drop = FALSE])
  rsums2 <- rowSums(mat[, -idx, drop = FALSE])

  return(ids[which((rsums1 > 0) & (rsums2 == 0))])
}
