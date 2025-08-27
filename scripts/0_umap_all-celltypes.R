#!/usr/bin/env Rscript
setwd("/Users/yangchen/PhD/Gallo_lab/Isotretinoin-study/scripts")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

# ---------- Inputs ----------
rds_path <- "../Michigan_data/singlecell/seurat.RDS"
out_dir  <- "../figures/umaps"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# UMAP settings
npcs     <- 30
dims_use <- 1:npcs
pt_size  <- 0.2
seed     <- 1234

# ---------- Load ----------
obj <- readRDS(rds_path)
stopifnot(inherits(obj, "Seurat"))
stopifnot("celltype" %in% colnames(obj@meta.data))

# Prefer SCT if present; otherwise RNA
DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"

# ---------- Ensure embeddings ----------
if (!"umap" %in% Reductions(obj)) {
  # Make sure variable features exist
  if (length(VariableFeatures(obj)) == 0) {
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  }
  # PCA if missing
  if (!"pca" %in% Reductions(obj)) {
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
  }
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE, seed.use = seed)
}

# ***** ---------- Plot: all samples colored by celltype ---------- *****
Idents(obj) <- "celltype"
p_celltype <- DimPlot(
  obj, reduction = "umap",
  group.by = "celltype",
  label = TRUE, repel = TRUE, label.size = 3,
  pt.size = pt_size
) + ggplot2::labs(title = "UMAP: by celltype of all samples") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(out_dir, "umap_celltype_allsamples.png"),
  plot = p_celltype,
  width = 8, height = 6, dpi = 600, limitsize = FALSE
)

# ***** ---------- Plot: all samples split pre- and post- colored by celltype ---------- *****
# Clean time to Pre/Post and keep order
time_clean_raw <- trimws(tolower(as.character(obj$time)))
obj$time_clean <- factor(
  ifelse(time_clean_raw %in% c("pre","post"), time_clean_raw, NA_character_),
  levels = c("pre","post"),
  labels = c("Pre","Post")
)

# Base UMAP split by time, colored by celltype
Idents(obj) <- "celltype"
p2 <- DimPlot(
  obj, reduction = "umap",
  group.by = "celltype",
  split.by = "time_clean",
  ncol = 2, label = FALSE, pt.size = pt_size
) + ggplot2::labs(title = "UMAP: by cell type of all samples") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# ----- Build label centroids for selected celltypes -----
emb <- as.data.frame(Embeddings(obj, "umap"))
emb$cell <- rownames(emb)
md  <- obj@meta.data

emb$celltype   <- md[emb$cell, "celltype", drop = TRUE]
emb$time_clean <- md[emb$cell, "time_clean", drop = TRUE]

targets <- with(emb, tolower(celltype) %in% c("sebaceous gland cells"))
label_df <- emb[targets & !is.na(emb$time_clean), , drop = FALSE] |>
  dplyr::group_by(celltype, time_clean) |>
  dplyr::summarise(
    UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
    UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

# Overlay labels (will appear in the correct Pre/Post facet)
p2 <- p2 +
  ggrepel::geom_label_repel(
    data = label_df,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, label = celltype),
    inherit.aes = FALSE,
    size = 3,
    label.padding = grid::unit(0.15, "lines"),
    label.r = grid::unit(0.1, "lines"),
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 42
  )

# Save @ 600 dpi
ggsave(file.path(out_dir, "umap_celltype_pre_vs_post.png"),
       p2, width = 12, height = 6, dpi = 600, limitsize = FALSE)
message("Saved: ", normalizePath(file.path(out_dir, "umap_celltype_pre_vs_post.png")))

# ***** ---------- Build df_all for per-participant Pre/Post ggplots ---------- *****
stopifnot("orig.ident" %in% colnames(obj@meta.data))
emb_all <- as.data.frame(Embeddings(obj, "umap"))[, 1:2, drop = FALSE]

# Use the cleaned time we already constructed (Pre/Post)
time_fac <- obj$time_clean

# Parse participant from orig.ident like ..._P7_...
oi <- as.character(obj@meta.data$orig.ident)
pid <- ifelse(grepl("[Pp][0-9]+", oi),
              paste0("pt", sub(".*[Pp]([0-9]+).*", "\\1", oi)),
              NA_character_)

df_all <- data.frame(
  UMAP_1   = emb_all$UMAP_1,
  UMAP_2   = emb_all$UMAP_2,
  celltype = obj@meta.data$celltype,
  time     = time_fac,
  pid      = pid,
  check.names = FALSE
)
df_all <- df_all[!is.na(df_all$time) & !is.na(df_all$pid) & !is.na(df_all$celltype), , drop = FALSE]
df_all$celltype <- factor(df_all$celltype)
df_all$time     <- factor(df_all$time, levels = c("Pre","Post"))
df_all$pid      <- factor(df_all$pid)

# Global limits & palette for consistent look
xlim_glob <- range(df_all$UMAP_1, na.rm = TRUE)
ylim_glob <- range(df_all$UMAP_2, na.rm = TRUE)
ct_levels <- levels(df_all$celltype)
pal <- setNames(scales::hue_pal(h = c(0, 360) + 15, c = 110, l = 55)(length(ct_levels)), ct_levels)

# ***** ---------- Plot: split by patient (pt3 / pt7), Pre vs Post, colored by celltype ---------- *****
# helper: prefer ragg (no Cairo/XQuartz), else base device
save_png <- function(path, plot, width_in = 12, height_in = 6, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(path, plot, width = width_in, height = height_in, dpi = dpi,
           device = ragg::agg_png, limitsize = FALSE)
  } else {
    ggsave(path, plot, width = width_in, height = height_in, dpi = dpi, limitsize = FALSE)
  }
}

draw_umap_prepost_gg <- function(data_df, pid, outfile,
                                 xlim_glob = range(data_df$UMAP_1, na.rm = TRUE),
                                 ylim_glob = range(data_df$UMAP_2, na.rm = TRUE),
                                 pt_size = 1.1,
                                 pt_alpha = 1,
                                 legend_pt_size = 4,
                                 base_size = 11,
                                 title_size = 20,
                                 axis_title_size = 14,
                                 axis_text_size  = 12,
                                 legend_title_size = 12,
                                 legend_text_size  = 10,
                                 strip_text_size  = 14,
                                 max_points_per_panel = 100000L,
                                 width_in = 12, height_in = 6, dpi = 600,
                                 label_celltype = "Sebaceous Gland Cells",
                                 label_size = 3.5,
                                 label_alpha = 0.9,
                                 label_color = "black") {
  
  stopifnot(is.data.frame(data_df))
  need <- c("UMAP_1","UMAP_2","celltype","time","pid")
  stopifnot(all(need %in% names(data_df)))
  
  d <- data_df %>% dplyr::filter(pid == !!pid, !is.na(time), !is.na(celltype))
  if (nrow(d) == 0L) { message("No cells for ", pid, "; skipping."); return(invisible(NULL)) }
  
  set.seed(1)
  d_pts <- d %>%
    dplyr::group_by(time) %>%
    dplyr::group_modify(function(.x, .y) {
      k <- min(nrow(.x), max_points_per_panel)
      if (k < nrow(.x)) .x[sample.int(nrow(.x), k), , drop = FALSE] else .x
    }) %>% dplyr::ungroup()
  
  # centroid label per facet for target celltype
  label_df <- d %>%
    dplyr::filter(tolower(celltype) == tolower(label_celltype)) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
                     UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(label = label_celltype)
  
  p <- ggplot2::ggplot(d_pts, ggplot2::aes(UMAP_1, UMAP_2, color = celltype)) +
    ggplot2::geom_point(size = pt_size, alpha = pt_alpha) +
    ggplot2::facet_wrap(~ time, ncol = 2, scales = "fixed") +
    ggplot2::coord_cartesian(xlim = xlim_glob, ylim = ylim_glob) +
    ggplot2::scale_color_manual(
      values = pal, limits = ct_levels, drop = FALSE, name = "Cell type",
      guide = ggplot2::guide_legend(override.aes = list(size = legend_pt_size, alpha = 1))
    ) +
    ggplot2::labs(
      title = paste0("UMAP: by cell type — ", toupper(pid)),
      x = "UMAP_1", y = "UMAP_2"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text  = ggplot2::element_text(size = axis_text_size),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text  = ggplot2::element_text(size = legend_text_size),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = strip_text_size, face = "bold"),
      panel.border = ggplot2::element_blank()
    )
  
  if (nrow(label_df)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        data = label_df,
        ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE,
        size = label_size, color = label_color,
        label.size = 0.2, box.padding = 0.25, point.padding = 0.3,
        min.segment.length = 0, segment.alpha = 0.6, max.overlaps = Inf
      )
    } else {
      p <- p + ggplot2::geom_label(
        data = label_df,
        ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE,
        size = label_size, color = label_color,
        alpha = label_alpha, label.size = 0.2
      )
    }
  }
  
  save_png(outfile, p, width_in = width_in, height_in = height_in, dpi = dpi)
  message("Saved: ", normalizePath(outfile))
}

# --- Make your two plots (per participant) ---
draw_umap_prepost_gg(
  df_all, "pt3",
  file.path(out_dir, "umap_celltype_pre_vs_post_pt3.png"),
  xlim_glob = xlim_glob, ylim_glob = ylim_glob,
  pt_size = 1.1, pt_alpha = 1, legend_pt_size = 4
)

draw_umap_prepost_gg(
  df_all, "pt7",
  file.path(out_dir, "umap_celltype_pre_vs_post_pt7.png"),
  xlim_glob = xlim_glob, ylim_glob = ylim_glob,
  pt_size = 1.1, pt_alpha = 1, legend_pt_size = 4
)

message("Saved per-participant UMAPs to: ", normalizePath(out_dir))
