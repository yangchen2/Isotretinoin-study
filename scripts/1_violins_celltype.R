#!/usr/bin/env Rscript
setwd("/Users/yangchen/PhD/Gallo_lab/Isotretinoin-study/scripts")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# -------- Load Seurat object --------
rds_path <- "../Michigan_data/singlecell/seurat.RDS"
obj <- readRDS(rds_path)
stopifnot(inherits(obj, "Seurat"))

# -------- Metadata & setup (memory-light) --------
md <- obj@meta.data
stopifnot(all(c("celltype","time","orig.ident") %in% names(md)))

# choose metric
value_col <- dplyr::case_when(
  "nFeature_RNA" %in% names(md) ~ "nFeature_RNA",
  "nCount_RNA"   %in% names(md) ~ "nCount_RNA",
  TRUE ~ NA_character_
)
if (is.na(value_col)) stop("Set `value_col` to a numeric metadata column (e.g., 'nFeature_RNA').")
stopifnot(value_col %in% colnames(md))

# clean time and parse participant from orig.ident (…_P7_… → pt7)
clean <- function(x) trimws(tolower(as.character(x)))
time_clean <- clean(md$time)
pid_num <- ifelse(grepl("[Pp][0-9]+", md$orig.ident),
                  sub(".*[Pp]([0-9]+).*", "\\1", md$orig.ident),
                  NA_character_)
pid_vec <- ifelse(is.na(pid_num), NA_character_, paste0("pt", pid_num))

# restrict global stats to the two participants to save memory
keep_pid <- pid_vec %in% c("pt3","pt7")
keep_time <- time_clean %in% c("pre","post")
ok <- keep_pid & keep_time & !is.na(md$celltype) & !is.na(md[[value_col]])
global_y <- range(md[[value_col]][ok], na.rm = TRUE)

# consistent celltype order (pt3+pt7 only)
ct_order <- md[ok, ] |>
  dplyr::mutate(celltype = factor(celltype)) |>
  dplyr::group_by(celltype) |>
  dplyr::summarize(med = stats::median(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(med)) |>
  dplyr::pull(celltype)

# output dir
out_dir <- "../figures/violins/celltype"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# jitter dot settings (smaller to lower memory)
max_points_per_group <- 1000L
point_size  <- 0.2
point_alpha <- 0.25
jitter_w    <- 0.15

# prefer ragg device if available
save_png <- function(path, plot, width_in = 4.2, height_in = 4.8, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(path, plot, width = width_in, height = height_in, dpi = dpi,
                    device = ragg::agg_png, limitsize = FALSE)
  } else {
    ggplot2::ggsave(path, plot, width = width_in, height = height_in, dpi = dpi, limitsize = FALSE)
  }
}

set.seed(1)

# -------- Per-participant plotting (pt3, pt7) --------
for (pid_target in c("pt3","pt7")) {
  idx <- which(pid_vec == pid_target & keep_time)
  if (!length(idx)) { message("No data for ", pid_target, "; skipping."); next }
  
  # small per-participant df
  dat_pid <- data.frame(
    celltype = factor(md$celltype[idx], levels = ct_order),
    time     = factor(time_clean[idx], levels = c("pre","post")),
    value    = md[[value_col]][idx],
    stringsAsFactors = FALSE
  )
  dat_pid <- dat_pid[!is.na(dat_pid$celltype) & !is.na(dat_pid$time) & !is.na(dat_pid$value), , drop = FALSE]
  if (!nrow(dat_pid)) { message("No usable rows for ", pid_target); next }
  
  # loop over this participant's cell types
  for (ct in levels(dat_pid$celltype)) {
    subdat <- dat_pid[dat_pid$celltype == ct, , drop = FALSE]
    if (!nrow(subdat)) next
    
    # t-test (needs ≥2 per group)
    x_pre  <- subdat$value[subdat$time == "pre"]
    x_post <- subdat$value[subdat$time == "post"]
    if (sum(!is.na(x_pre)) >= 2 && sum(!is.na(x_post)) >= 2) {
      tt <- t.test(x_pre, x_post, var.equal = FALSE)
      p_text <- paste0("t-test p = ", formatC(tt$p.value, format = "e", digits = 2))
    } else {
      p_text <- "t-test: NA (need ≥2 per group)"
    }
    y_pos <- global_y[2] - 0.02 * diff(global_y)
    
    # dot sampling per group (lower to reduce memory)
    subdat_pts <- subdat |>
      dplyr::group_by(time) |>
      dplyr::group_modify(function(.x, .y) {
        k <- min(nrow(.x), max_points_per_group)
        if (k < nrow(.x)) .x[sample.int(nrow(.x), k), , drop = FALSE] else .x
      }) |>
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(subdat, ggplot2::aes(x = time, y = value, fill = time)) +
      ggplot2::geom_violin(trim = TRUE, scale = "width") +
      ggplot2::geom_jitter(data = subdat_pts, width = jitter_w, height = 0,
                           size = point_size, alpha = point_alpha) +
      ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.25, alpha = 0.9) +
      ggplot2::annotate("text", x = 1.5, y = y_pos, label = p_text, size = 4) +
      ggplot2::labs(title = paste0(as.character(ct), " — ", toupper(pid_target)),
                    x = NULL, y = "Gene Expression") +
      ggplot2::coord_cartesian(ylim = global_y) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    outfile <- file.path(out_dir, pid_target, paste0("violin_", pid_target, "_",
                                         gsub("[^A-Za-z0-9_\\-]+","_", as.character(ct)), ".png"))
    save_png(outfile, p)
    message("Saved: ", normalizePath(outfile))
    
    # aggressively free per-plot temporaries
    rm(subdat, subdat_pts, p, x_pre, x_post, tt); gc(verbose = FALSE)
  }
  
  rm(dat_pid); gc(verbose = FALSE)
}
message("Done. Files in: ", normalizePath(out_dir))

# -------- Combined pt3+pt7 plotting --------
out_dir_comb <- file.path("../figures/violins/celltype", "pts_combined")
dir.create(out_dir_comb, showWarnings = FALSE, recursive = TRUE)

idx_comb <- which(pid_vec %in% c("pt3","pt7") & keep_time)
if (length(idx_comb)) {
  dat_comb <- data.frame(
    celltype = factor(md$celltype[idx_comb], levels = ct_order),
    time     = factor(time_clean[idx_comb], levels = c("pre","post")),
    value    = md[[value_col]][idx_comb],
    stringsAsFactors = FALSE
  )
  dat_comb <- dat_comb[!is.na(dat_comb$celltype) & !is.na(dat_comb$time) & !is.na(dat_comb$value), , drop = FALSE]
  
  for (ct in levels(dat_comb$celltype)) {
    subdat <- dat_comb[dat_comb$celltype == ct, , drop = FALSE]
    if (!nrow(subdat)) next
    
    # t-test
    x_pre  <- subdat$value[subdat$time == "pre"]
    x_post <- subdat$value[subdat$time == "post"]
    if (sum(!is.na(x_pre)) >= 2 && sum(!is.na(x_post)) >= 2) {
      tt <- t.test(x_pre, x_post, var.equal = FALSE)
      p_text <- paste0("t-test p = ", formatC(tt$p.value, format = "e", digits = 2))
    } else {
      p_text <- "t-test: NA (need ≥2 per group)"
    }
    y_pos <- global_y[2] - 0.02 * diff(global_y)
    
    # dot sampling
    subdat_pts <- subdat |>
      dplyr::group_by(time) |>
      dplyr::group_modify(function(.x, .y) {
        k <- min(nrow(.x), max_points_per_group)
        if (k < nrow(.x)) .x[sample.int(nrow(.x), k), , drop = FALSE] else .x
      }) |>
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(subdat, ggplot2::aes(x = time, y = value, fill = time)) +
      ggplot2::geom_violin(trim = TRUE, scale = "width") +
      ggplot2::geom_jitter(data = subdat_pts, width = jitter_w, height = 0,
                           size = point_size, alpha = point_alpha) +
      ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.25, alpha = 0.9) +
      ggplot2::annotate("text", x = 1.5, y = y_pos, label = p_text, size = 4) +
      ggplot2::labs(title = paste0(as.character(ct), " — PT3+PT7"),
                    x = NULL, y = "Gene Expression") +
      ggplot2::coord_cartesian(ylim = global_y) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    outfile <- file.path(out_dir_comb, paste0("violin_pts_combined_",
                                              gsub("[^A-Za-z0-9_\\-]+","_", as.character(ct)), ".png"))
    save_png(outfile, p)
    message("Saved combined: ", normalizePath(outfile))
    
    rm(subdat, subdat_pts, p, x_pre, x_post, tt); gc(verbose = FALSE)
  }
}


#################### GENE INFO #########################
# pick the assay you care about
assay <- DefaultAssay(obj)           # or "RNA", "SCT", etc.
feat  <- rownames(obj[[assay]])      # all genes/features in that assay

# all genes for the assay
genes <- rownames(obj[["RNA"]])
length(genes); head(genes)

# Check if genes of interest are present and annotated
genes_to_check <- c("PLA2G7", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3")
present <- genes_to_check %in% feat
data.frame(gene = genes_to_check, present)
