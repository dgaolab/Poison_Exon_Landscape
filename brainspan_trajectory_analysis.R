###############################################################################
# PSI GAM Analysis — BrainSpan developmental trajectory
###############################################################################

# ---- 0. Configuration (edit here only) ----------------------
CONFIG <- list(
  # Paths — edit these to point at your local copies
  path_bs    = "path/to/brainspan_poison_exon_output",
  pe_ref     = "spliceIDs_ref.poison_exon.no_ir.GRCh38.GencodeV26.filtered.txt",
  output_dir = "./3_output_brainSpan",
  csv_dir    = "./3_output_brainSpan/2_brain_span_files",

  # Sample filters
  min_region_n    = 30, # keep brain regions with > this many samples total
  min_age_group_n = 3,  # drop age groups that are singleton (n=1) in more than this many tissues

  # Feature filters
  coverage_threshold = 5,    # min read count per sample
  coverage_min_pct   = 0.5,  # min fraction of samples passing coverage per age group
  var_cutoff         = 1e-5, # min row variance to keep feature

  # GAM settings
  gam_k         = 9,       
  psi_transform = "logit", 
  fdr_method    = "fdr"
)

# ---- 1. Libraries -------------------------------------------
library(tidyverse)
library(data.table)
library(R.utils)
library(mgcv)
library(matrixStats)

# ---- 2. Load & filter metadata ------------------------------
load_metadata <- function(path_bs, min_region_n, min_age_group_n) {
  meta <- read.table(file.path(path_bs, "meta.tsv"),
                     header = TRUE, sep = "\t", check.names = FALSE)

  # Filter 1: remove age groups that have only 1 sample in more than min_age_group_n tissues
  age_tissue_counts <- meta %>%
    dplyr::count(age_group, brain_region, name = "n_samples") %>%
    group_by(age_group) %>%
    summarise(n_sparse_tissues = sum(n_samples == 1), .groups = "drop") %>%
    arrange(age_group) %>%
    mutate(kept = n_sparse_tissues <= min_age_group_n)

  for (i in seq_len(nrow(age_tissue_counts))) {
    status <- if (age_tissue_counts$kept[i]) "KEEP" else "EXCLUDE"
    message(sprintf("  %-45s  singleton_tissues = %d  [%s]",
                    age_tissue_counts$age_group[i],
                    age_tissue_counts$n_sparse_tissues[i],
                    status))
  }

  keep_ages <- age_tissue_counts$age_group[age_tissue_counts$kept]
  excluded  <- age_tissue_counts$age_group[!age_tissue_counts$kept]
  if (length(excluded) > 0) {
    message(sprintf("=> Excluding %d age group(s): %s\n",
                    length(excluded), paste(excluded, collapse = ", ")))
  } else {
    message("=> All age groups retained.\n")
  }

  meta_f <- meta %>% filter(age_group %in% keep_ages)

  # Filter 2: remove brain regions with too few samples
  region_counts <- meta_f %>% dplyr::count(brain_region)
  small_regions <- region_counts$brain_region[region_counts$n < min_region_n]
  meta_f <- meta_f %>% filter(!(brain_region %in% small_regions))

  return(meta_f)
}

# ---- 3. Load PSI count matrices -----------------------------
load_psi_counts <- function(path_bs, keep_samples) {
  files_path <- file.path(path_bs, "poison_exon_spliceIDcounts_perTissue")
  files      <- list.files(files_path, full.names = TRUE)

  read_and_bind <- function(pattern) {
    flist <- grep(pattern, files, value = TRUE)
    mat   <- lapply(flist, function(f)
      read.table(gzfile(f), header = TRUE, row.names = 1, check.names = FALSE))
    dplyr::bind_cols(mat)[, keep_samples]
  }

  list(
    numerator   = read_and_bind(".*Numerator"),
    R3          = read_and_bind(".*R3")
  )
}

# ---- 4. Feature filtering helpers ---------------------------
# Keep features with >= threshold reads in >= sample_pct of samples
low_coverage_filtering <- function(input_matrix, threshold, sample_pct) {
  m <- as.matrix(input_matrix)
  keep <- apply(m, 1, function(x) mean(x >= threshold, na.rm = TRUE) >= sample_pct)
  rownames(input_matrix)[keep]
}

# ---- 5. GAM fitting -----------------------------------------
get_gam <- function(x, PSI_tissue, tissue_meta, gam_k = 9,
                    psi_transform = "none") {
  PE_name <- rownames(PSI_tissue)[x]

  df <- t(PSI_tissue[x, , drop = FALSE]) %>%
    as.data.frame() %>%
    setNames("PSI") %>%
    filter(!is.na(PSI)) %>%
    mutate(`entity:sample_id` = rownames(.)) %>%
    left_join(tissue_meta, by = "entity:sample_id")

  if (psi_transform == "logit") {
    df$PSI <- qlogis(pmax(pmin(df$PSI, 0.999), 0.001))
  }

  fit <- tryCatch(
    gam(PSI ~ s(age_group_int, k = gam_k), data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # Model with sex and ancestry as parametric covariates (for GCV comparison only)
  fit_cov <- tryCatch(
    gam(PSI ~ s(age_group_int, k = gam_k) + gender + ethn, data = df),
    error = function(e) NULL
  )

  s <- summary(fit)
  fitted_vals           <- t(as.data.frame(fitted(fit)))
  rownames(fitted_vals) <- PE_name
  colnames(fitted_vals) <- df$`entity:sample_id`

  list(
    summary = data.frame(PE_name      = PE_name,
                         EDF          = s$edf,
                         s_pv         = s$s.pv,
                         dev_expl     = s$dev.expl,
                         r_sq         = s$r.sq,
                         GCV_base     = fit$gcv.ubre.dev,
                         GCV_cov      = if (!is.null(fit_cov)) fit_cov$gcv.ubre.dev else NA_real_),
    fitValues = fitted_vals
  )
}

# ---- 6. Per-region pipeline ---------------------------------
run_region <- function(region, metadata, PSI_n, PSI_R3, cfg) {

  tmeta <- metadata %>%
    filter(brain_region == region) %>%
    select(`entity:sample_id`, age_group, brain_region, gender, ethn) %>%
    mutate(age_group_int = as.integer(factor(age_group, levels = sort(unique(age_group))))) %>%
    arrange(age_group_int)

  sids <- tmeta$`entity:sample_id`

  # Compute PSI = numerator / (numerator + R3)
  n_mat <- PSI_n[,  sids]
  d_mat <- PSI_R3[, sids]
  den   <- n_mat + d_mat

  # Filter 1: low coverage — per age group, keep union of passing features
  age_levels    <- sort(unique(tmeta$age_group))
  keep_coverage <- unique(unlist(lapply(age_levels, function(ag) {
    s <- tmeta$`entity:sample_id`[tmeta$age_group == ag]
    low_coverage_filtering(den[, s, drop = FALSE],
                           cfg$coverage_threshold, cfg$coverage_min_pct)
  })))

  psi <- (n_mat / den)[keep_coverage, ]

  psi <- as.matrix(psi)
  psi[is.nan(psi)] <- 0  # NaN = 0/0 (no reads); treat as PSI = 0
  psi <- as.data.frame(psi)

  # Filter 2: low variance
  rv  <- rowVars(as.matrix(psi), na.rm = TRUE)
  psi <- psi[!is.na(rv) & rv >= cfg$var_cutoff, ]
  message("  Features after filtering: ", nrow(psi))

  psi <- psi[, sids]  # ensure sample order matches metadata
  stopifnot(all(colnames(psi) == sids))

  # Fit GAMs across all features
  results <- lapply(seq_len(nrow(psi)), function(i)
    get_gam(i, psi, tmeta, cfg$gam_k, cfg$psi_transform))
  results <- Filter(Negate(is.null), results)

  # Collect summaries and apply FDR correction
  rdf <- do.call(rbind, lapply(results, `[[`, "summary")) %>%
    filter(!is.na(s_pv)) %>%
    mutate(p.adjust = p.adjust(s_pv, method = cfg$fdr_method))

  # Collect fitted values (back-transformed to PSI scale if logit was used)
  fit_df <- bind_rows(lapply(lapply(results, `[[`, "fitValues"), as.data.frame))

  write.csv(rdf,    file.path(cfg$csv_dir,
            paste0("PE_gam_results_",       region, ".csv")), row.names = FALSE)
  write.csv(fit_df, file.path(cfg$csv_dir,
            paste0("PE_gam_fitted_values_", region, ".csv")), row.names = TRUE)

  invisible(list(summary = rdf, fitted = fit_df))
}

# ---- 7. Main ------------------------------------------------
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$csv_dir,    showWarnings = FALSE, recursive = TRUE)

bs_meta_f     <- load_metadata(CONFIG$path_bs, CONFIG$min_region_n, CONFIG$min_age_group_n)
save(bs_meta_f, file = file.path(CONFIG$output_dir, "bs_meta_f.RData"))
counts        <- load_psi_counts(CONFIG$path_bs, bs_meta_f$`entity:sample_id`)

# Restrict to reference PE list
pe_ref <- read.table(CONFIG$pe_ref, header = TRUE, sep = "\t", check.names = FALSE)
pe_ids <- pe_ref$spliceID
counts$numerator   <- counts$numerator[rownames(counts$numerator)   %in% pe_ids, ]
counts$R3          <- counts$R3[rownames(counts$R3) %in% pe_ids, ]

brain_regions <- sort(unique(bs_meta_f$brain_region))

# -- GAM model fitting across all brain regions ----
for (region in brain_regions) {
  run_region(region, bs_meta_f, counts$numerator, counts$R3, CONFIG)
}
