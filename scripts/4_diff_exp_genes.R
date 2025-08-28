#!/usr/bin/env Rscript
# Pseudobulk + limma-voom DE: Post vs Pre for pt3, pt7, and combined (paired)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ------------------- Paths -------------------
rds_path <- "../Michigan_data/singlecell/seurat.RDS"
out_dir  <- "../analyses"
fig_dir  <- "../figures/volcanoes"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------- Load Seurat -------------------
obj <- readRDS(rds_path)
stopifnot(inherits(obj, "Seurat"))

assay <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cts   <- GetAssayData(obj, assay = assay, slot = "counts")
stopifnot(inherits(cts, "dgCMatrix"))

md <- obj@meta.data
stopifnot(all(c("time", "orig.ident") %in% colnames(md)))

# ------------------- Parse metadata: time & participant -------------------
clean <- function(x) trimws(tolower(as.character(x)))
time0 <- clean(md$time)
time  <- ifelse(time0 %in% c("pre","post"), time0, NA_character_)

oi <- as.character(md$orig.ident)
pid_num <- ifelse(grepl("[Pp][0-9]+", oi), sub(".*[Pp]([0-9]+).*", "\\1", oi), NA_character_)
pid <- ifelse(is.na(pid_num), NA_character_, paste0("pt", pid_num))

keep <- !is.na(time) & pid %in% c("pt3","pt7")
if (!any(keep)) stop("No cells with time in {pre, post} and pid in {pt3, pt7}.")

cts  <- cts[, keep, drop = FALSE]
time <- time[keep]
pid  <- pid[keep]

# ------------------- Pseudobulk by (pid x time) -------------------
# sample_id like "pt3_Pre", "pt3_Post", etc.
sample_id <- factor(paste0(pid, "_", ifelse(time=="pre","Pre","Post")))
S <- sparse.model.matrix(~ sample_id - 1)    # cells x samples
agg <- cts %*% S                             # genes x samples (sparse)

colnames(agg) <- sub("^sample_id", "", colnames(S))
# convert to ordinary integer matrix for downstream tools
agg <- as.matrix(agg)
mode(agg) <- "numeric"
agg <- round(agg)
storage.mode(agg) <- "integer"

# sample metadata
col_pid  <- sub("^(pt[0-9]+)_.*$", "\\1", colnames(agg))
col_time <- tolower(sub("^.*_(Pre|Post)$", "\\1", colnames(agg)))
coldata <- data.frame(
  row.names = colnames(agg),
  pid  = factor(col_pid,  levels = c("pt3","pt7")),
  time = factor(col_time, levels = c("pre","post"))
)

# ------------------- edgeR DGEList + TMM + voom -------------------
dge <- DGEList(counts = agg)
keep_genes <- filterByExpr(dge, design = model.matrix(~ pid + time, data = coldata))
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ pid + time, data = coldata)   # paired model
v <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Combined Post vs Pre (blocking by pid)
res_comb <- topTable(fit, coef = "timepost", number = Inf, sort.by = "P") %>%
  rownames_to_column("gene")
write_tsv(res_comb, file.path(out_dir, "limma_pts_combined_Post_vs_Pre.tsv"))

# ------------------- Per-patient logFCs & direction -------------------
# Compute logCPM (voom) differences Post-Pre per gene per participant
logcpm <- v$E  # genes x samples on log2 scale after TMM+voom
# sample order sanity: names match coldata rows
stopifnot(identical(colnames(logcpm), rownames(coldata)))

# helper to fetch columns for each pid/time
get_col <- function(who_pid, who_time) {
  which(coldata$pid == who_pid & coldata$time == who_time)
}

i_pt3_pre  <- get_col("pt3","pre")
i_pt3_post <- get_col("pt3","post")
i_pt7_pre  <- get_col("pt7","pre")
i_pt7_post <- get_col("pt7","post")

# compute deltas where both timepoints exist
has_pt3 <- length(i_pt3_pre)  == 1 && length(i_pt3_post) == 1
has_pt7 <- length(i_pt7_pre)  == 1 && length(i_pt7_post) == 1

delta_pt3 <- if (has_pt3) logcpm[, i_pt3_post] - logcpm[, i_pt3_pre] else rep(NA_real_, nrow(logcpm))
delta_pt7 <- if (has_pt7) logcpm[, i_pt7_post] - logcpm[, i_pt7_pre] else rep(NA_real_, nrow(logcpm))

per_patient <- tibble(
  gene = rownames(logcpm),
  logFC_pt3 = as.numeric(delta_pt3),
  logFC_pt7 = as.numeric(delta_pt7),
  same_direction = sign(logFC_pt3) == sign(logFC_pt7)
)

# Merge with combined limma stats
res_with_patients <- res_comb %>%
  select(gene, logFC_combined = logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  left_join(per_patient, by = "gene")

write_tsv(res_with_patients, file.path(out_dir, "limma_with_per_patient_logFCs.tsv"))

# ------------------- Optional: naive paired t-test on two deltas (df = 1) -------------------
# This is for reference only; with n=2 pairs the p-values are not reliable.
paired_t <- function(x, y) {
  d <- cbind(x, y)
  ok <- stats::complete.cases(d)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 1L || length(y) < 1L) return(c(t = NA, p = NA))
  diff <- x - y
  if (length(diff) != 2L) return(c(t = NA, p = NA))  # we expect exactly two patients
  tstat <- mean(diff) / (stats::sd(diff) / sqrt(2))
  pval  <- 2 * stats::pt(-abs(tstat), df = 1)
  c(t = tstat, p = pval)
}

pt_test <- mapply(function(g) {
  paired_t(per_patient$logFC_pt7[per_patient$gene == g],
           per_patient$logFC_pt3[per_patient$gene == g])
}, per_patient$gene)

# pack results (keeping order of genes)
pt_tbl <- tibble(
  gene = per_patient$gene,
  t_paired_df1 = as.numeric(pt_test["t", ]),
  p_paired_df1 = as.numeric(pt_test["p", ])
)

res_with_t <- res_with_patients %>% left_join(pt_tbl, by = "gene")
write_tsv(res_with_t, file.path(out_dir, "limma_plus_naivePairedT.tsv"))


# ------------------- Volcano plot: Post vs Pre (participants combined) -------------------
stopifnot(all(c("gene","logFC","P.Value") %in% colnames(res_comb)))

# thresholds
lfc_thresh <- 1
p_thresh   <- 0.05

volcano_df <- res_comb %>%
  dplyr::mutate(
    neglog10P = -log10(pmax(P.Value, .Machine$double.xmin)),
    sig = (abs(logFC) >= lfc_thresh & P.Value < p_thresh),
    direction_label = dplyr::case_when(
      sig & logFC > 0  ~ "↑Post-iso",
      sig & logFC < 0  ~ "↑Pre-iso",
      TRUE             ~ "NS"
    )
  )

# pick top N significant genes per side
lab_pos <- volcano_df %>%
  dplyr::filter(sig, logFC > 0) %>%
  dplyr::arrange(P.Value) %>%
  dplyr::slice_head(n = 20)

lab_neg <- volcano_df %>%
  dplyr::filter(sig, logFC < 0) %>%
  dplyr::arrange(P.Value) %>%
  dplyr::slice_head(n = 20)

lab_df <- dplyr::bind_rows(lab_pos, lab_neg)

p_volcano <- ggplot2::ggplot(
  volcano_df,
  ggplot2::aes(x = logFC, y = neglog10P, color = direction_label)
) +
  ggplot2::geom_point(size = 1, alpha = 0.85) +
  ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = "dotted", color = "black") +
  ggplot2::scale_color_manual(
    values = c("↑Pre-iso" = "#f8766d",
               "↑Post-iso" = "#00bfc4",
               "NS"    = "grey70"),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Differentially expressed genes (participants combined)",
    x = "log2 Fold Change (Post / Pre)",
    y = "-log10(p-value)"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    plot.title   = ggplot2::element_text(hjust = 0.5),
    legend.position = c(0.9, 0.1)  # x=80% across, y=20% up from bottom
  )

# add labels for the top 20 per side
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(x = logFC, y = neglog10P, label = gene, color = direction_label),
      inherit.aes = FALSE,
      size = 2,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.alpha = 0.6,
      seed = 42
    )
} else {
  p_volcano <- p_volcano +
    ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(x = logFC, y = neglog10P, label = gene, color = direction_label),
      size = 2, vjust = -0.4
    )
}

# save
if (requireNamespace("ragg", quietly = TRUE)) {
  ggplot2::ggsave(
    filename = file.path(fig_dir, "volcano_pts_combined_Post_vs_Pre_labeled.png"),
    plot = p_volcano, width = 7, height = 6, dpi = 600,
    device = ragg::agg_png, limitsize = FALSE
  )
} else {
  ggplot2::ggsave(
    filename = file.path(fig_dir, "volcano_pts_combined_Post_vs_Pre_labeled.png"),
    plot = p_volcano, width = 8, height = 6, dpi = 600, limitsize = FALSE
  )
}

# ------------------- Export significant DE genes -------------------
sig_genes <- volcano_df %>%
  dplyr::filter(sig) %>%
  dplyr::transmute(
    gene,
    direction = ifelse(logFC > 0, "Post-associated", "Pre-associated"),
    logFC,
    P.Value,
    neglog10P
  )

out_sig <- file.path(out_dir, "significant_genes_combined.tsv")
readr::write_tsv(sig_genes, out_sig)

message("Significant gene list saved to: ", normalizePath(out_sig))


message("Done. Results in: ", normalizePath(out_dir))
