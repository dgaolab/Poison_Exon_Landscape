###############################################################################
# GTEx differential splicing analysis for poison exons
#
# Inputs:
#   - gtex_dir : directory containing
#       * poison_exon_spliceIDcounts_perTissue/  (per-tissue numerator/denominator
#         count tables, named <tissue>.poison_exon.numerator.out.tab.gz and
#         <tissue>.poison_exon.denominator.out.tab.gz)
#       * GTEx_v8_RNASeq_meta.txt                (sample metadata with SMTSD column)
#   - pe_ref  : filtered poison-exon reference table
#                 (spliceIDs_ref.poison_exon.no_ir.GRCh38.GencodeV26.filtered.txt)
#   - GTEx_dsp.R : helper script defining test_one_pair() and find_tissue_uniq()
#
###############################################################################

library(ggplot2)
library(combinat)
library(data.table)
library(gplots)
library(hues)
library(matrixStats)
library(rwantshue)
library(dplyr)
library(tibble)

source("./GTEx_dsp.R")

# ---- CONFIG -----------------------------------------------------------------
gtex_dir <- "path/to/gtex_poison_exon_output"  
pe_ref   <- "spliceIDs_ref.poison_exon.no_ir.GRCh38.GencodeV26.filtered.txt"
outdir   <- "./4_output_GTEx"
# -----------------------------------------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "RData"), showWarnings = FALSE)
dir.create(file.path(outdir, "pdf"),   showWarnings = FALSE)

dd  <- read.table(pe_ref, header = TRUE, sep = "\t", check.names = FALSE)
pe  <- as.vector(dd$spliceID)
dd1 <- dd
rownames(dd1) <- as.vector(dd$spliceID)
dd1 <- dd1[, -1]

folder <- file.path(gtex_dir, "poison_exon_spliceIDcounts_perTissue/")
fl <- list.files(folder)
tl <- unique(sapply(strsplit(fl, "\\."), "[[", 1))

ps <- t(combn(tl, 2))

meta <- fread(file.path(gtex_dir, "GTEx_v8_RNASeq_meta.txt"), sep = "\t")
meta <- data.frame(meta, row.names = 1)
meta$tissue <- as.vector(meta$SMTSD)
meta$tissue <- gsub(" - ", "_", meta$tissue)
meta$tissue <- gsub(" ",   "_", meta$tissue)
meta$tissue <- gsub("\\(", "",  meta$tissue)
meta$tissue <- gsub("\\)", "",  meta$tissue)

diff_dir <- file.path(outdir, "differential_splicing_test")
dir.create(diff_dir, showWarnings = FALSE)
apply(ps, 1, FUN = test_one_pair, meta = meta, tissue_col = "tissue",
      cutoff = 3, output_dir = diff_dir, folder = folder, pe_ids = pe)


##################################################################################
# ------  extract values from models ------
extract_values <- function(x, ref, keyword, col_index){
  fld <- file.path(outdir, "differential_splicing_test")
  pre <- paste0("model_params_", keyword, "_glm_tissue")
  f   <- paste0(fld, "/", pre, "_", x[2], "_VS_", x[1], ".txt")
  d   <- read.table(f, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
  vals <- rep(NA, length(ref))
  names(vals) <- ref
  vals[rownames(d)] <- d[rownames(d), col_index]
  return(vals)
}

has_file <- apply(ps, 1, function(x){
  file.exists(paste0(diff_dir, "/model_params_beta_glm_tissue_", x[2], "_VS_", x[1], ".txt"))
})
ps <- ps[has_file, ]
tissues_kept <- unique(as.vector(ps))
cat(sprintf("Pairs with model output: %d / %d  |  Tissues represented: %d\n",
            sum(has_file), sum(has_file) + sum(!has_file), length(tissues_kept)))

df_beta <- c()
df_fdr  <- c()
df_base <- c()

for (i in 1:nrow(ps)){
  tmp <- extract_values(ps[i, ], ref = pe, keyword = "beta", col_index = 1)
  df_base <- cbind(df_base, tmp)
  tmp <- extract_values(ps[i, ], ref = pe, keyword = "beta", col_index = 2)
  df_beta <- cbind(df_beta, tmp)
  tmp <- extract_values(ps[i, ], ref = pe, keyword = "fdr",  col_index = 2)
  df_fdr  <- cbind(df_fdr, tmp)
}

df_target <- df_base + df_beta

coln <- apply(ps, 1, FUN = function(x){ paste0(x[2], "_VS_", x[1]) })

colnames(df_beta)   <- coln
colnames(df_base)   <- coln
colnames(df_target) <- coln
colnames(df_fdr)    <- coln

saveRDS(df_beta,   file = file.path(outdir, "RData/df_beta.rds"))
saveRDS(df_base,   file = file.path(outdir, "RData/df_base.rds"))
saveRDS(df_target, file = file.path(outdir, "RData/df_target.rds"))
saveRDS(df_fdr,    file = file.path(outdir, "RData/df_fdr.rds"))

res <- matrix(0, ncol = ncol(df_beta), nrow = nrow(df_beta))

# Back-transform from logit scale to PSI scale, then threshold on |delta PSI| >= 0.2
psi_base   <- plogis(df_base)
psi_target <- plogis(df_target)
delta_psi  <- psi_target - psi_base

res[Reduce(intersect, list(which(delta_psi >=  0.2), which(df_fdr < 0.1), which(psi_base < 0.5 & psi_target < 0.5)))] <-  1
res[Reduce(intersect, list(which(delta_psi <= -0.2), which(df_fdr < 0.1), which(psi_base < 0.5 & psi_target < 0.5)))] <- -1

rownames(res) <- rownames(df_beta)
colnames(res) <- colnames(df_beta)

rmaxs <- rowMaxs(abs(res))
res1  <- res[which(rmaxs != 0), ]

##################################################################################
# Brain-vs-other-tissue poison exons
brain_comp       <- grep("Brain", colnames(res))
brain_inter_comp <- grep("^Brain_.+_VS_Brain", colnames(res))
brain_comp_ok    <- setdiff(brain_comp, brain_inter_comp)
res2   <- res1[, brain_comp_ok]
rmaxs2 <- rowMaxs(abs(res2))
res2   <- res2[which(rmaxs2 != 0), ]

delta_psi_brain <- as.matrix(delta_psi[rownames(res2), colnames(res2)])
delta_psi_brain[which(is.na(delta_psi_brain))] <- 0

# Flip direction so every column reads <non-brain>_VS_Brain
brain_first          <- grep("^Brain", colnames(res2))
brain_first_coln     <- colnames(res2)[brain_first]
brain_first_coln_inv <- unname(sapply(brain_first_coln, FUN = function(x){
  ll <- strsplit(x, "_VS_"); paste0(ll[[1]][2], "_VS_", ll[[1]][1])
}))

tmp <- res2[, brain_first] * (-1)
colnames(tmp) <- brain_first_coln_inv
res3 <- cbind(tmp, res2[, -brain_first])

tmp <- delta_psi_brain[, brain_first] * (-1)
colnames(tmp) <- brain_first_coln_inv
delta_psi_brain <- cbind(tmp, delta_psi_brain[, -brain_first])
delta_psi_brain <- delta_psi_brain[, sort(colnames(delta_psi_brain))]

coln <- sapply(strsplit(colnames(delta_psi_brain), "_"), "[[", 1)
coln[which(coln == "Minor")] <- "Salivary_Gland"
coln[which(coln == "Small")] <- "Intestine"
coln[which(coln == "Whole")] <- "Blood"
non_brain_tis <- unique(coln)

scheme   <- iwanthue(seed = 122, force_init = TRUE)
myColors <- scheme$hex(length(non_brain_tis))
names(myColors) <- non_brain_tis
color_vec <- myColors[coln]

pdf(file.path(outdir, "pdf/test_Brain_GTEx_poison_exon_spliceID_heatmap.pdf"))
heatmap.2(delta_psi_brain, Colv = FALSE, dendrogram = "row",
          scale = "none", col = "bluered", trace = "none",
          ColSideColors = color_vec, density.info = "none",
          labRow = FALSE, labCol = FALSE)
dev.off()


##################################################################################
# brain-unique PEs
brain_uniq <- find_tissue_uniq(tl, pe, "Brain", folder = folder)
brain_diff <- rownames(delta_psi_brain)

brain_poison  <- rbind(data.frame(ID = brain_uniq, Class = "brain_uniq"),
                       data.frame(ID = brain_diff, Class = "brain_diff"))
brain_poison1 <- dd1[as.vector(brain_poison$ID), ]
brain_poison1$Tissue <- brain_poison$Class

