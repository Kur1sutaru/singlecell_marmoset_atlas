## ====== Packages ======
pkgs <- c("dplyr","tidyr","ggplot2","scales","forcats","ggalluvial","patchwork","Seurat")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(forcats); library(ggalluvial); library(patchwork); library(Seurat)
})

## ====== Assumptions ======
# predtable and obj already exist in your environment

## ====== Normalize / assert columns ======
names(predtable) <- trimws(names(predtable))
need_cols <- c("Cluster","Cells","Major label","Subtype label")
missing <- setdiff(need_cols, names(predtable))
if (length(missing)) stop("Missing columns in predtable: ", paste(missing, collapse=", "))

cluster_col <- "Cluster"
count_col   <- "Cells"
label_major <- "Major label"
label_sub   <- "Subtype label"

## ====== Output folder ======
base_out <- "figs_no_sample"
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

## ====== Helper: order factor levels by weighted/global abundance ======
order_by_abundance <- function(x, w = NULL) {
  if (is.null(w)) {
    lev <- names(sort(table(x), decreasing = TRUE))
  } else {
    agg <- tapply(w, x, sum)
    lev <- names(sort(agg, decreasing = TRUE))
  }
  factor(x, levels = lev)
}

## ====== Consistent palette across ALL plots ======
all_labels <- sort(unique(c(predtable[[label_major]], predtable[[label_sub]])))
base_cols <- c(
  "#66C2A5","#1F78B4","#FC8D62","#8DA0CB","#A6D854",
  "#E31A1C","#FFD92F","#E78AC3","#33A02C","#B15928",
  "#B3B3B3","#CAB2D6","#FDBF6F","#A6CEE3"
)
if (length(all_labels) > length(base_cols)) {
  extra <- scales::hue_pal()(length(all_labels) - length(base_cols))
  all_cols <- c(base_cols, extra)
} else {
  all_cols <- base_cols[seq_along(all_labels)]
}
celltype_cols <- setNames(all_cols, all_labels)
locks <- c(
  "Fibroblast"="#66C2A5", "Endothelial"="#1F78B4", "Keratinocyte"="#8DA0CB",
  "NK cell"="#FB9A99", "NKT cell"="#FDBF6F", "Monocyte"="#B15928",
  "Neutrophil"="#E31A1C", "Mast cell"="#E78AC3",
  "Vascular smooth muscle cell"="#33A02C"
)
for (nm in names(locks)) if (nm %in% names(celltype_cols)) celltype_cols[nm] <- locks[nm]

## ====== Core plotting (NO Sample used) ======
make_plots_global <- function(df_in,
                              label_choice = c("major","subtype"),
                              out_dir = file.path(base_out, "major")) {
  label_choice <- match.arg(label_choice)
  label_col <- if (label_choice == "major") label_major else label_sub
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  df <- df_in %>%
    dplyr::transmute(
      Cluster = as.character(.data[[cluster_col]]),
      Label   = as.character(.data[[label_col]]),
      Cells   = as.numeric(.data[[count_col]])
    )
  
  ## ===== Global composition (pooled) =====
  comp <- df %>%
    dplyr::group_by(.data$Label) %>%
    dplyr::summarise(Cells = sum(.data$Cells, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Proportion = .data$Cells / sum(.data$Cells)) %>%
    dplyr::mutate(Label = order_by_abundance(.data$Label, w = .data$Cells))
  
  # Bar (global)
  p_bar <- ggplot(comp, aes(x = .data$Label, y = .data$Proportion, fill = .data$Label)) +
    geom_col(color = "white") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = celltype_cols, drop = FALSE) +
    coord_flip() +
    labs(x = NULL, y = "Proportion (global)", fill = "Cell type",
         title = paste0("Global composition — ", label_col)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "global_bar.png"), p_bar, width = 7, height = 6, dpi = 300)
  
  # Donut (global)
  comp_donut <- comp %>%
    dplyr::arrange(dplyr::desc(.data$Proportion)) %>%
    dplyr::mutate(ymax = cumsum(.data$Proportion),
                  ymin = dplyr::lag(.data$ymax, default = 0))
  p_donut <- ggplot(comp_donut, aes(ymin = .data$ymin, ymax = .data$ymax, xmin = 3, xmax = 4, fill = .data$Label)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y") + xlim(c(0,4)) +
    scale_fill_manual(values = celltype_cols, drop = FALSE) +
    theme_void(base_size = 12) +
    theme(legend.position = "right") +
    labs(title = paste("Global donut —", label_col), fill = "Cell type")
  ggsave(file.path(out_dir, "global_donut.png"), p_donut, width = 8, height = 6, dpi = 300)
  
  ## ===== Alluvial: Cluster → Label =====
  flow_cluster <- df %>%
    dplyr::group_by(.data$Cluster, .data$Label) %>%
    dplyr::summarise(Cells = sum(.data$Cells, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      Label   = order_by_abundance(.data$Label, w = .data$Cells),
      Cluster = order_by_abundance(.data$Cluster, w = .data$Cells)
    )
  
  p_alluvial <- ggplot(flow_cluster,
                       aes(axis1 = .data$Cluster, axis2 = .data$Label, y = .data$Cells)) +
    geom_alluvium(aes(fill = .data$Label), knot.pos = 0.4, alpha = 0.9) +
    geom_stratum(width = 0.15, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = celltype_cols, drop = FALSE) +
    labs(x = NULL, y = "Cells", fill = "Cell type",
         title = paste0("Alluvial: Cluster → ", label_col)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  ggsave(file.path(out_dir, "alluvial_cluster_to_label.png"), p_alluvial, width = 12, height = 6, dpi = 300)
  
  message("Saved figures to: ", normalizePath(out_dir))
}

## ====== Run (GLOBAL, no Sample) for Major and Subtype ======
make_plots_global(predtable, "major",   file.path(base_out, "major"))
make_plots_global(predtable, "subtype", file.path(base_out, "subtype"))

## ====== Merge labels into obj (by seurat_clusters) and plot UMAPs ======
add_labels_to_obj <- function(obj, pred_df) {
  map_major <- pred_df %>%
    dplyr::distinct(Cluster = .data[[cluster_col]], Major = .data[[label_major]]) %>%
    dplyr::mutate(Cluster = as.character(.data$Cluster))
  map_sub <- pred_df %>%
    dplyr::distinct(Cluster = .data[[cluster_col]], Subtype = .data[[label_sub]]) %>%
    dplyr::mutate(Cluster = as.character(.data$Cluster))
  
  md <- obj@meta.data
  md$Cluster <- as.character(obj$seurat_clusters)
  md <- md %>% dplyr::left_join(map_major, by = "Cluster") %>% dplyr::left_join(map_sub, by = "Cluster")
  names(md)[names(md) == "Major"]   <- "Major_label"
  names(md)[names(md) == "Subtype"] <- "Subtype_label"
  
  if (anyNA(md$Major_label))   md$Major_label   <- ifelse(is.na(md$Major_label),   paste0("Cluster ", md$Cluster), md$Major_label)
  if (anyNA(md$Subtype_label)) md$Subtype_label <- ifelse(is.na(md$Subtype_label), paste0("Cluster ", md$Cluster), md$Subtype_label)
  
  md$Major_label   <- order_by_abundance(md$Major_label)
  md$Subtype_label <- order_by_abundance(md$Subtype_label)
  obj@meta.data <- md
  obj
}

# Ensure UMAP exists
if (!"umap" %in% names(obj@reductions)) {
  if (!"pca" %in% names(obj@reductions)) {
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
  }
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
}

# Add labels into obj
obj <- add_labels_to_obj(obj, predtable)

# UMAP colored by Major (no Sample)
Idents(obj) <- obj$Major_label
p_umap_major <- DimPlot(obj, label = TRUE, repel = TRUE, cols = celltype_cols) +
  ggtitle("UMAP — Major labels")
ggsave(file.path(base_out, "UMAP_major.png"), p_umap_major, width = 8, height = 6, dpi = 300)

# UMAP colored by Subtype (no Sample)
Idents(obj) <- obj$Subtype_label
p_umap_sub <- DimPlot(obj, label = TRUE, repel = TRUE, cols = celltype_cols) +
  ggtitle("UMAP — Subtype labels")
ggsave(file.path(base_out, "UMAP_subtype.png"), p_umap_sub, width = 8, height = 6, dpi = 300)

message("All figures saved under: ", normalizePath(base_out))



## ===== Major cell types =====
# mapping: cluster IDs → major labels
major_map <- c(
  "0" = "Fibroblast",
  "1" = "Fibroblast",
  "2" = "Fibroblast",
  "3" = "Fibroblast",
  "4" = "Endothelial",
  "5" = "NKT cell",
  "6" = "Vascular smooth muscle cell",
  "7" = "Vascular smooth muscle cell",
  "8" = "Fibroblast",
  "9" = "Neutrophil",
  "10" = "Mast cell",
  "11" = "Keratinocyte",
  "12" = "Vascular smooth muscle cell",
  "13" = "Fibroblast",
  "14" = "Fibroblast",
  "15" = "Fibroblast",
  "16" = "Monocyte",
  "17" = "NK cell"
 
)

# set Idents and rename
Idents(obj) <- "seurat_clusters"   # make sure cluster identities are active
obj <- RenameIdents(obj, major_map)

# store into metadata
obj$MajorType <- Idents(obj)


## ===== Subtypes =====
# mapping: cluster IDs → subtype labels
subtype_map <- c(
  "0" = "",
  "1" = "",
  "2" = "",
  "3" = "",
  "4" = "",
  "5" = "",
  "6" = "",
  "7" = "",
  "8" = "",
  "9" = "",
  "10" = "",
  "11" = "",
  "12" = "",
  "13" = "",
  "14" = "",
  "15" = "",
  "16" = "",
  "17" = ""

)

Idents(obj) <- "seurat_clusters"
obj <- RenameIdents(obj, subtype_map)

# store into metadata
obj$Subtype <- Idents(obj)

library(dplyr)

# 1) Build mapping vectors from your predtable
cluster_col  <- "Cluster"
label_major  <- "Major label"
label_sub    <- "Subtype label"

map_major <- predtable %>%
  distinct(Cluster = .data[[cluster_col]], Major = .data[[label_major]]) %>%
  mutate(Cluster = as.character(Cluster))

map_sub <- predtable %>%
  distinct(Cluster = .data[[cluster_col]], Subtype = .data[[label_sub]]) %>%
  mutate(Cluster = as.character(Cluster))

major_vec <- setNames(map_major$Major, map_major$Cluster)
sub_vec   <- setNames(map_sub$Subtype, map_sub$Cluster)

# 2) Per-cell cluster as character (keys for indexing)
cluster_chr <- as.character(obj$seurat_clusters)

# 3) WRITE metadata columns WITHOUT joins (preserves row order/rownames)
obj$Major_label   <- unname(major_vec[cluster_chr])
obj$Subtype_label <- unname(sub_vec[cluster_chr])

# Optional: sanity checks
setdiff(unique(cluster_chr), names(major_vec))   # clusters missing a major label
setdiff(unique(cluster_chr), names(sub_vec))     # clusters missing a subtype label

# 4) Plot
DimPlot(obj, group.by = "Major_label", label = TRUE, repel = TRUE)
DimPlot(obj, group.by = "Subtype_label", label = TRUE, repel = TRUE)


## Ensure per-cell labels are assigned (vectorized, no joins)
cluster_col  <- "Cluster"
label_major  <- "Major label"
label_sub    <- "Subtype label"

map_major <- predtable |>
  dplyr::distinct(Cluster = .data[[cluster_col]], Major = .data[[label_major]]) |>
  dplyr::mutate(Cluster = as.character(.data$Cluster))

major_vec <- setNames(map_major$Major, map_major$Cluster)

cluster_chr <- as.character(obj$seurat_clusters)
obj$Major_label <- unname(major_vec[cluster_chr])

# quick sanity checks
cat("NAs in Major_label:", sum(is.na(obj$Major_label)), "\n")
print(head(table(obj$Major_label)))

## Make it a factor with informative order (by abundance)
order_by_abundance <- function(x) {
  lev <- names(sort(table(x), decreasing = TRUE))
  factor(x, levels = lev)
}
obj$Major_label <- order_by_abundance(obj$Major_label)

## Build a palette for ONLY the labels present
present_labels <- na.omit(unique(as.character(obj$Major_label)))
# If you already built a big 'celltype_cols' map, subset it; otherwise make one:
if (exists("celltype_cols")) {
  pal_major <- celltype_cols[present_labels]
  # fill any missing names with a hue palette
  if (anyNA(pal_major)) {
    miss <- is.na(pal_major)
    pal_major[miss] <- scales::hue_pal()(sum(miss))
    names(pal_major) <- present_labels
  }
} else {
  pal_major <- setNames(scales::hue_pal()(length(present_labels)), present_labels)
}

## Plot: use group.by = "Major_label", and pass the palette
# (labels = TRUE is optional; it just tries to label cluster/group centers)
p_umap_major <- DimPlot(
  obj,
  group.by = "Major_label",
  label = TRUE, repel = TRUE,
  cols = pal_major
) + ggtitle("UMAP — Major labels")

ggsave(
  filename = file.path(getwd(), "UMAP_major.png"),
  plot = p_umap_major,
  width = 10, height = 8, dpi = 600
)



suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tibble) })
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(rlang) })

# Safe helper: returns a named vector  Cluster -> Label
build_map <- function(df, cluster_col, label_col) {
  cc <- rlang::sym(cluster_col)
  lc <- rlang::sym(label_col)
  df %>%
    transmute(
      Cluster = trimws(as.character(!!cc)),
      Label   = trimws(as.character(!!lc))
    ) %>%
    distinct(Cluster, .keep_all = TRUE) %>%   # one label per cluster
    arrange(Cluster) %>%
    tibble::deframe()                         # names = Cluster, values = Label
}
.build_map <- build_map
cluster_col <- "Cluster"
major_col   <- "Major label"
sub_col     <- "Subtype label"

major_map <- build_map(predtable, cluster_col, major_col)
sub_map   <- build_map(predtable, cluster_col, sub_col)
library(Seurat)

# identities = clusters as character (levels "0","1",...)
Idents(obj) <- factor(as.character(obj$seurat_clusters))

# ---- Major ----
# keep only keys that exist in current levels (avoids warnings)
major_use <- major_map[names(major_map) %in% levels(Idents(obj))]

# optional: stop if any cluster is unmapped
missing_major <- setdiff(levels(Idents(obj)), names(major_use))
if (length(missing_major)) stop("Missing Major for clusters: ", paste(missing_major, collapse=", "))

obj <- RenameIdents(obj, major_use)
obj$Major_label <- as.character(Idents(obj))   # write to metadata

# ---- Subtype ----
Idents(obj) <- factor(as.character(obj$seurat_clusters))  # reset to clusters
sub_use <- sub_map[names(sub_map) %in% levels(Idents(obj))]
missing_sub <- setdiff(levels(Idents(obj)), names(sub_use))
if (length(missing_sub)) stop("Missing Subtype for clusters: ", paste(missing_sub, collapse=", "))

obj <- RenameIdents(obj, sub_use)
obj$Subtype_label <- as.character(Idents(obj)) # write to metadata
Idents(obj) <- factor(obj$Major_label)

# Simple plots (no custom palette needed to test)
p_major <- DimPlot(obj, label = TRUE, repel = TRUE) + ggtitle("UMAP — Major")
p_sub   <- DimPlot(obj, group.by = "Subtype_label", label = TRUE, repel = TRUE) + ggtitle("UMAP — Subtype")

ggsave("UMAP_major.png", p_major, width = 12, height = 9, dpi = 600)
ggsave("UMAP_subtype.png", p_sub,   width = 12, height = 9, dpi = 600)

new.cluster.ids <- c("Fibroblast", "Fibroblast","Fibroblast", "Fibroblast","Endothelial", "NKT cell","Vascular smooth muscle cell", 
                     "Vascular smooth muscle cell","Fibroblast", "Neutrophil","Mast cell", "Keratinocyte","Vascular smooth muscle cell", 
                     "Fibroblast","Fibroblast", "Fibroblast","Monocyte", "NK cell")
names(new.cluster.ids) <- levels(rat)
rat <- RenameIdents(rat, new.cluster.ids)
DimPlot(rat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(rat, "rat_CRA017434.rds")


library(dplyr)
library(dplyr)

## --- pick the right column names from your table ---
pred_df <- predtable  # (use the object you actually have)
cluster_col <- "Cluster"
sub_col <- if ("Subtype label" %in% names(pred_df)) "Subtype label" else "Subtype.label"

## --- build mapping: cluster -> subtype (as character) ---
sub_map <- pred_df %>%
  transmute(
    Cluster = as.character(.data[[cluster_col]]),
    Subtype = trimws(as.character(.data[[sub_col]]))
  ) %>%
  distinct(Cluster, .keep_all = TRUE) %>%
  arrange(Cluster)

sub_vec <- setNames(sub_map$Subtype, sub_map$Cluster)

## --- assign per-cell (make vector UNNAMED so Seurat doesn't try to match names) ---
# use your Seurat object name; you used 'rat' in your command
clust_chr <- as.character(rat$seurat_clusters)

# safest: write straight into @meta.data (preserves order, no name matching)
rat@meta.data$Subtype_label <- unname(sub_vec[clust_chr])

## sanity checks
cat("NAs in Subtype_label:", sum(is.na(rat@meta.data$Subtype_label)), "\n")
print(head(table(rat@meta.data$Subtype_label, useNA = "ifany")))

## plot
DimPlot(rat, reduction = "umap", group.by = "Subtype_label", label = TRUE, repel = TRUE) +
  ggtitle("UMAP — Subtypes from predtable")


## plot with subtype annotation
DimPlot(rat, reduction = "umap", group.by = "Subtype_label", label = TRUE, repel = TRUE) +
  ggtitle("UMAP — Subtypes from rat tongue")


