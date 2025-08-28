#!/usr/bin/env Rscript
setwd("/Users/yangchen/PhD/Gallo_lab/Isotretinoin-study/scripts")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# Inputs
# -----------------------------
rds_path <- "../Michigan_data/singlecell/seurat.RDS"
target_celltypes <- c("Differentiated Keratinocytes")
# "Myeloid Cells", "Fibroblasts", "Sebaceous Gland Cells", "Differentiated Keratinocytes", "Basal Keratinocytes"

#genes <- c("PLA2G7", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3")
genes <- c("PLA2G2A")

# -----------------------------
# Load object and metadata
# -----------------------------
obj <- readRDS(rds_path)
stopifnot(inherits(obj, "Seurat"))
md <- obj@meta.data
stopifnot(all(c("celltype","time","orig.ident") %in% colnames(md)))

if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"

# Clean time to Pre/Post
clean <- function(x) trimws(tolower(as.character(x)))
md$time_clean <- clean(md$time)
keep_levels <- c("pre","post")
md <- md %>%
  mutate(time_fac = factor(ifelse(time_clean %in% keep_levels, time_clean, NA_character_),
                           levels = keep_levels,
                           labels = c("Pre","Post")))

# Parse participant IDs
pid_num <- ifelse(grepl("[Pp][0-9]+", md$orig.ident),
                  sub(".*[Pp]([0-9]+).*", "\\1", md$orig.ident),
                  NA_character_)
pid_vec <- ifelse(is.na(pid_num), NA_character_, paste0("pt", pid_num))

# -----------------------------
# Helper / settings
# -----------------------------
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)
max_points_per_group <- 5000L
point_size  <- 0.2
point_alpha <- 0.25
jitter_w    <- 0.15
set.seed(1)

# -----------------------------
# Loop over cell types
# -----------------------------
for (target_celltype in target_celltypes) {
  message("Processing celltype: ", target_celltype)
  
  for (pid_target in c("pt3","pt7","pts_combined")) {
    message("  Subsetting: ", pid_target)

    # Make directory for this celltype + participant
    out_dir <- file.path("../figures/violins/prostaglandin_genes",
                         safe_name(target_celltype),
                         pid_target)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Subset cells
    if (pid_target == "pts_combined") {
      sel_cells <- rownames(md)[!is.na(md$time_fac) &
                                tolower(trimws(md$celltype)) == tolower(trimws(target_celltype)) &
                                pid_vec %in% c("pt3","pt7")]
    } else {
      sel_cells <- rownames(md)[!is.na(md$time_fac) &
                                tolower(trimws(md$celltype)) == tolower(trimws(target_celltype)) &
                                pid_vec == pid_target]
    }
    if (length(sel_cells) == 0) {
      message("    No cells for ", target_celltype, " / ", pid_target, "; skipping.")
      next
    }
    
    obj_sub <- subset(obj, cells = sel_cells)
    obj_sub@meta.data$time_fac <- md[colnames(obj_sub), "time_fac", drop = TRUE]
    
    # Filter genes present
    genes_present <- intersect(genes, rownames(obj_sub))
    genes_missing <- setdiff(genes, rownames(obj_sub))
    if (length(genes_missing)) {
      message("    Skipping missing genes: ", paste(genes_missing, collapse = ", "))
    }
    if (!length(genes_present)) {
      message("    None of the requested genes present; skipping.")
      next
    }
    
    # Loop over genes
    for (g in genes_present) {
      df <- FetchData(obj_sub, vars = c(g, "time_fac"))
      colnames(df) <- c("expr", "time_fac")
      df <- df %>% filter(!is.na(time_fac))
      if (nrow(df) == 0) next
      
      # t-test
      x_pre  <- df$expr[df$time_fac == "Pre"]
      x_post <- df$expr[df$time_fac == "Post"]
      if (sum(!is.na(x_pre)) >= 2 && sum(!is.na(x_post)) >= 2) {
        tt <- t.test(x_pre, x_post, var.equal = FALSE)
        p_text <- paste0("t-test p = ", formatC(tt$p.value, format = "e", digits = 2))
      } else {
        p_text <- "t-test: NA (need ≥2 per group)"
      }
      
      # Sample points
      df_pts <- df %>%
        group_by(time_fac) %>%
        group_modify(function(.x, .y) {
          k <- min(nrow(.x), max_points_per_group)
          if (k < nrow(.x)) .x[sample.int(nrow(.x), size = k), , drop = FALSE] else .x
        }) %>%
        ungroup()
      
      # Plot
      y_lim <- range(df$expr, na.rm = TRUE)
      y_pos <- y_lim[2] - 0.02 * diff(y_lim)
      
      p <- ggplot(df, aes(x = time_fac, y = expr, fill = time_fac)) +
        geom_violin(trim = TRUE, scale = "width") +
        geom_jitter(data = df_pts, width = jitter_w, height = 0,
                    size = point_size, alpha = point_alpha) +
        geom_boxplot(width = 0.12, outlier.size = 0.25, alpha = 0.9) +
        annotate("text", x = 1.5, y = y_pos, label = p_text, size = 3.5) +
        labs(
          title = paste0(target_celltype, " — ", g, " — ", toupper(pid_target)),
          x = NULL, y = "Gene Expression"
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))
      
      out_file <- file.path(out_dir, paste0("violin_", safe_name(target_celltype),
                                            "_", safe_name(pid_target),
                                            "_", safe_name(g), ".png"))
      ggsave(out_file, p, width = 4.2, height = 4.8, dpi = 600, limitsize = FALSE)
      message("    Saved: ", normalizePath(out_file))
    }
  }
}
message("All done.")
