# ============================================================
# Fibroblast subsetting, reclustering, validation & plotting
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)        # v5+
  library(Seuratmouseect)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(plyr)
})

# ----------------- paths and setup -----------------
in_rds   <- mouse        # <- your input file
out_dir  <- "fibro_mouse_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------- 1) load & subset fibro -----------------
mouse <- readRDS(in_rds)
DefaultAssay(mouse) <- if ("RNA" %in% Assays(mouse)) "RNA" else DefaultAssay(mouse)

# columns to search for fibro labels
candidates <- c(
  "celltype_merged","label_coarse","label_fine","celltype","cell_type","CellType",
  "annotation","annotations","predicted.celltype.l2","predicted.celltype.l1",
  "seurat_annotations","broad_label","majorType","fine_label"
)
candidates <- intersect(candidates, colnames(mouse@meta.data))
if (length(candidates) == 0) stop("No known label columns found in metadata.")

# match fibro in any candidate column
fibro_mask <- Reduce(`|`, lapply(candidates, function(cn)
  grepl("fibro", mouse@meta.data[[cn]], ignore.case = TRUE)
))
if (!any(fibro_mask)) stop("No fibroblast-like labels matched in metadata.")

# (optional) exclude common contaminants if mislabeled under fibro
drop_mask <- Reduce(`|`, lapply(candidates, function(cn)
  grepl("(peri|pericyte|smc|smooth|myocyte|endo|vascular|krt|keratin|epithel)",
        mouse@meta.data[[cn]], ignore.case = TRUE)
))
fibro_mask <- fibro_mask & !drop_mask

fib <- subset(mouse, cells = colnames(mouse)[fibro_mask])
if (ncol(fib) < 50) warning("Fewer than 50 fibro cells kept — proceed carefully.")
fib$parent_label <- mouse@meta.data[colnames(fib), intersect("celltype_merged", candidates), drop = TRUE]

saveRDS(fib, file.path(out_dir, "mouse_fibro_only_initial.rds"))

# ----------------- 2) normalize & reduce -----------------
# Option A: LogNormalize (simple & fast)
fib <- NormalizeData(fib)
fib <- FindVariableFeatures(fib, nfeatures = 26207)
fib <- ScaleData(fib, features = VariableFeatures(fib), verbose = FALSE)
fib <- RunPCA(fib, features = VariableFeatures(fib), npcs = 50, verbose = FALSE)

# # Option B: SCTransform (uncomment if preferred)
# fib <- SCTransform(fib, verbose = FALSE)
# fib <- RunPCA(fib, npcs = 50, verbose = FALSE)

# ----------------- 3) neighbors, UMAP, clustering --------
fib <- FindNeighbors(fib, dims = 1:30)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (r in resolutions) {
  fib <- FindClusters(fib, resolution = r, algorithm = 1)  # Louvain
}
fib <- RunUMAP(fib, dims = 1:30, n.neighbors = 30, min.dist = 0.25)

# Preview UMAPs across resolutions
p_list <- lapply(resolutions, function(r)
  DimPlot(fib, reduction = "umap",
          group.by = paste0("RNA_snn_res.", r),
          label = TRUE, repel = TRUE) + ggtitle(paste("res", r))
)
ggsave(file.path(out_dir, "UMAP_all_resolutions.png"),
       wrap_plots(p_list, ncol = 2), width = 12, height = 10, dpi = 300)

# choose default resolution for downstream
default_res <- 0.6
cluster_col <- paste0("RNA_snn_res.", default_res)
Idents(fib) <- fib@meta.data[[cluster_col]]

# ----------------- 4) DEGs per subcluster -----------------
min_cells <- 30
keep_clusters <- names(which(table(Idents(fib)) >= min_cells))
fib_deg <- subset(fib, idents = keep_clusters)

markers <- FindAllMarkers(
  fib_deg,
  assay = DefaultAssay(fib_deg),
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.10,
  return.thresh = 0.05
) %>% arrange(cluster, desc(avg_log2FC))
write.csv(markers, file.path(out_dir, "markers_per_cluster.csv"), row.names = FALSE)

# ----------------- 5) validation panels -----------------
# (A) fibro & subtypes (adjust as needed)
panels <- list(
  Fibro_core          = c("Col1a1","Col1a2","Col3a1","Dcn","Lum","Dpt","Fn1","Col15a1","Sparc","Col18a1","Cdh11","Pdgfra"),
  Pi16_fibro          = c("Pi16","Dpp4","Pdgfra"),
  Myofibro            = c("Acta2","Tagln","Myl9","Myh11","Cnn1"),
  Inflammatory_fibro  = c("Cxcl14","Cxcl12","Il6","Ccl2"),
  Chondro_fibro       = c("Acan","Sox9","Col2a1","Comp")    # keep if relevant
)
# (B) contaminants to check
panels$Perivascular   <- c("Rgs5","Kcnj8","Abcc9","Pdgfrb","Mcam")              # pericyte/SMC-like
panels$Keratinocyte   <- c("Krt14","Krt5","Krt17","Krt6a","Krt6b","Krt1","Krt10","Krt8","Krt18","Krt15","Epcam")
# (C) proliferation
panels$Proliferation  <- c("Mki67","Top2a","Pcna","Cdk1","Ube2c","Birc5","Ccnb1","Cenpf")

# keep genes present
panels <- lapply(panels, function(v) intersect(v, rownames(fib)))
panels <- panels[sapply(panels, length) > 0]

# per-cell module scores
for (nm in names(panels)) {
  fib <- AddModuleScore(fib, features = list(panels[[nm]]), name = paste0(nm,"_Score"))
  sc_col <- paste0(nm,"_Score1")
  names(fib@meta.data)[names(fib@meta.data) == sc_col] <- paste0(nm,"_Score")
}
score_cols <- paste0(names(panels), "_Score")
score_cols <- score_cols[score_cols %in% colnames(fib@meta.data)]

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---------- prep ----------
if (!exists("out_dir")) out_dir <- "fibro_mouse_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
mod_dir <- file.path(out_dir, "module_scores"); dir.create(mod_dir, FALSE, TRUE)

# Detect cluster column if not set
if (!exists("cluster_col")) {
  cand <- grep("_snn_res\\.", colnames(fib@meta.data), value = TRUE)
  if (length(cand) == 0) stop("No clustering column found (.*_snn_res.*).")
  cluster_col <- cand[which.max(nchar(cand))]  # pick the most specific
}
Idents(fib) <- fib@meta.data[[cluster_col]]

# Define panels if missing
if (!exists("panels")) {
  panels <- list(
    Fibro_core          = c("Col1a1","Col1a2","Col3a1","Dcn","Lum","Dpt","Fn1","Col15a1","Sparc","Col18a1","Cdh11","Pdgfra"),
    Pi16_fibro          = c("Pi16","Dpp4","Pdgfra"),
    Myofibro            = c("Acta2","Tagln","Myl9","Myh11","Cnn1"),
    Inflammatory_fibro  = c("Cxcl14","Cxcl12","Il6","Ccl2"),
    Perivascular        = c("Rgs5","Kcnj8","Abcc9","Pdgfrb","Mcam"),
    Keratinocyte        = c("Krt14","Krt5","Krt17","Krt6a","Krt6b","Krt1","Krt10","Krt8","Krt18","Krt15","Epcam"),
    Proliferation       = c("Mki67","Top2a","Pcna","Cdk1","Ube2c","Birc5","Ccnb1","Cenpf")
  )
}

# Keep only genes present
panels <- lapply(panels, function(v) intersect(v, rownames(fib)))
panels <- panels[sapply(panels, length) > 0]

# Ensure module scores exist; if missing, compute
for (nm in names(panels)) {
  sc <- paste0(nm, "_Score")
  if (!sc %in% colnames(fib@meta.data)) {
    fib <- AddModuleScore(fib, features = list(panels[[nm]]), name = paste0(nm, "_Score"))
    # AddModuleScore appends a "1"
    colnames(fib@meta.data)[colnames(fib@meta.data) == paste0(nm, "_Score1")] <- sc
  }
}

# Helper: per-module plots
plot_module <- function(seu, module_name, genes, cluster_col, out_dir) {
  score_col <- paste0(module_name, "_Score")
  genes <- intersect(genes, rownames(seu))
  if (!score_col %in% colnames(seu@meta.data)) return(invisible(NULL))
  
  # order clusters by mean module score (nice for Dot/Violin)
  ord <- seu@meta.data |>
    mutate(cluster = seu@meta.data[[cluster_col]]) |>
    group_by(cluster) |>
    summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") |>
    arrange(desc(mean_score)) |>
    pull(cluster) |> as.character()
  
  # UMAP colored by module score (quantile cutoffs to stabilize palette)
  p_umap <- FeaturePlot(seu, features = score_col, reduction = "umap",
                        min.cutoff = "q05", max.cutoff = "q95") +
    ggtitle(paste0(module_name, " score (UMAP)"))
  
  # Violin by cluster (ordered)
  seu@meta.data[[cluster_col]] <- factor(seu@meta.data[[cluster_col]], levels = ord)
  p_vln <- VlnPlot(seu, features = score_col, group.by = cluster_col,
                   pt.size = 0, assay = DefaultAssay(seu)) +
    ggtitle(paste0(module_name, " score (Violin)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # DotPlot of the **score** by cluster (size = % cells with score > 0)
  p_dot_score <- DotPlot(seu, features = score_col, group.by = cluster_col) +
    RotatedAxis() + ggtitle(paste0(module_name, " score (DotPlot)"))
  
  # DotPlot of the module **genes** by cluster
  p_dot_genes <- NULL
  if (length(genes) > 0) {
    p_dot_genes <- DotPlot(seu, features = genes, group.by = cluster_col) +
      RotatedAxis() + ggtitle(paste0(module_name, " genes (DotPlot)"))
  }
  
  # Save
  ggsave(file.path(out_dir, paste0("UMAP_", module_name, "_score.png")),
         p_umap, width = 6.5, height = 5.5, dpi = 320)
  ggsave(file.path(out_dir, paste0("Violin_", module_name, "_score.png")),
         p_vln, width = 8, height = 5.5, dpi = 320)
  ggsave(file.path(out_dir, paste0("DotPlot_", module_name, "_score.png")),
         p_dot_score, width = 6.5, height = 5.5, dpi = 320)
  if (!is.null(p_dot_genes)) {
    ggsave(file.path(out_dir, paste0("DotPlot_", module_name, "_genes.png")),
           p_dot_genes, width = 10, height = 6.5, dpi = 320)
  }
  
  # combined sheet (UMAP + Violin + ScoreDot + GenesDot)
  comb <- if (is.null(p_dot_genes)) (p_umap | p_vln) / p_dot_score else (p_umap | p_vln) / (p_dot_score | p_dot_genes)
  ggsave(file.path(out_dir, paste0("COMBO_", module_name, ".png")),
         comb, width = 14, height = if (is.null(p_dot_genes)) 7.5 else 11, dpi = 320)
  
  invisible(list(umap = p_umap, violin = p_vln, dot_score = p_dot_score, dot_genes = p_dot_genes))
}

# ---------- run per module ----------
all_mod_plots <- list()
for (nm in names(panels)) {
  all_mod_plots[[nm]] <- plot_module(fib, nm, panels[[nm]], cluster_col, mod_dir)
}

# ---------- optional: one big combined figure across modules ----------
# (UMAPs of all scores side-by-side; adjust ncol if many)
umaps <- lapply(names(panels), function(nm) all_mod_plots[[nm]]$umap + ggtitle(nm))
umaps <- umaps[!sapply(umaps, is.null)]
if (length(umaps) > 0) {
  ggsave(file.path(mod_dir, "ALL_modules_UMAP_scores.png"),
         wrap_plots(umaps, ncol = 3), width = 16, height = ceiling(length(umaps)/3)*5.5, dpi = 300)
}

# Export cluster means of each module score (handy for tables/heatmaps)
score_cols <- paste0(names(panels), "_Score")
avg_mod <- fib@meta.data |>
  mutate(cluster = fib@meta.data[[cluster_col]]) |>
  group_by(cluster) |>
  summarise(across(all_of(score_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") |>
  arrange(as.numeric(as.character(cluster)))
write.csv(avg_mod, file.path(mod_dir, "module_scores_cluster_means.csv"), row.names = FALSE)

# --- ensure output dir for modules ---
if (!exists("out_dir")) out_dir <- "fibro_mouse_outputs"
mod_dir <- file.path(out_dir, "module_scores"); dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)

# --- detect cluster column if not set ---
if (!exists("cluster_col")) {
  cand <- grep("_snn_res\\.", colnames(fib@meta.data), value = TRUE)
  stopifnot(length(cand) > 0)
  cluster_col <- cand[which.max(nchar(cand))]
}
Idents(fib) <- fib@meta.data[[cluster_col]]

# ---------- safer helper (no .data / across) ----------
cluster_order_by_score <- function(seu, score_col, cluster_col) {
  cl <- seu@meta.data[[cluster_col]]
  sc <- seu@meta.data[[score_col]]
  m <- tapply(sc, cl, function(v) mean(v, na.rm = TRUE))
  names(sort(m, decreasing = TRUE))
}

plot_module2 <- function(seu, module_name, genes, cluster_col, out_dir) {
  score_col <- paste0(module_name, "_Score")
  if (!score_col %in% colnames(seu@meta.data)) return(invisible(NULL))
  genes <- intersect(genes, rownames(seu))
  
  # order clusters by mean score
  ord <- cluster_order_by_score(seu, score_col, cluster_col)
  
  # UMAP colored by module score
  p_umap <- FeaturePlot(seu, features = score_col, reduction = "umap",
                        min.cutoff = "q05", max.cutoff = "q95", order = TRUE) +
    ggtitle(paste0(module_name, " score (UMAP)"))
  
  # Violin (ordered)
  seu@meta.data[[cluster_col]] <- factor(seu@meta.data[[cluster_col]], levels = ord)
  p_vln <- VlnPlot(seu, features = score_col, group.by = cluster_col, pt.size = 0) +
    ggtitle(paste0(module_name, " score (Violin)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # DotPlot of the score (DotPlot can use metadata columns via FetchData)
  p_dot_score <- DotPlot(seu, features = score_col, group.by = cluster_col) +
    RotatedAxis() + ggtitle(paste0(module_name, " score (DotPlot)"))
  
  # DotPlot of module genes
  p_dot_genes <- NULL
  if (length(genes) > 0) {
    p_dot_genes <- DotPlot(seu, features = genes, group.by = cluster_col) +
      RotatedAxis() + ggtitle(paste0(module_name, " genes (DotPlot)"))
  }
  
  # Save all pieces
  ggsave(file.path(out_dir, paste0("UMAP_", module_name, "_score.png")),
         p_umap, width = 6.5, height = 5.5, dpi = 320)
  ggsave(file.path(out_dir, paste0("Violin_", module_name, "_score.png")),
         p_vln, width = 8, height = 5.5, dpi = 320)
  ggsave(file.path(out_dir, paste0("DotPlot_", module_name, "_score.png")),
         p_dot_score, width = 6.5, height = 5.5, dpi = 320)
  if (!is.null(p_dot_genes)) {
    ggsave(file.path(out_dir, paste0("DotPlot_", module_name, "_genes.png")),
           p_dot_genes, width = 10, height = 6.5, dpi = 320)
  }
  
  comb <- if (is.null(p_dot_genes)) (p_umap | p_vln) / p_dot_score else (p_umap | p_vln) / (p_dot_score | p_dot_genes)
  ggsave(file.path(out_dir, paste0("COMBO_", module_name, ".png")),
         comb, width = 14, height = if (is.null(p_dot_genes)) 7.5 else 11, dpi = 320)
  
  invisible(list(umap = p_umap, violin = p_vln, dot_score = p_dot_score, dot_genes = p_dot_genes))
}

# ---------- run per module with safety (won’t abort on one failure) ----------
all_mod_plots <- setNames(vector("list", length(panels)), names(panels))
for (nm in names(panels)) {
  all_mod_plots[[nm]] <- tryCatch(
    plot_module2(fib, nm, panels[[nm]], cluster_col, mod_dir),
    error = function(e) { message("Module ", nm, " failed: ", e$message); NULL }
  )
}

# ---------- combined UMAP sheet (only modules that succeeded) ----------
ok_umaps <- lapply(names(all_mod_plots), function(nm) {
  x <- all_mod_plots[[nm]]
  if (is.null(x)) return(NULL)
  x$umap + ggtitle(nm)
})
ok_umaps <- ok_umaps[!sapply(ok_umaps, is.null)]
if (length(ok_umaps) > 0) {
  ggsave(file.path(mod_dir, "ALL_modules_UMAP_scores.png"),
         wrap_plots(ok_umaps, ncol = 3), width = 16,
         height = ceiling(length(ok_umaps)/3) * 5.5, dpi = 300)
}

# ---------- export cluster means for each module score (base R aggregate) ----------
score_cols <- paste0(names(panels), "_Score")
score_cols <- intersect(score_cols, colnames(fib@meta.data))
df <- fib@meta.data[, c(cluster_col, score_cols), drop = FALSE]

avg_mod <- aggregate(df[, score_cols, drop = FALSE],
                     by = list(cluster = df[[cluster_col]]),
                     FUN = function(x) mean(as.numeric(x), na.rm = TRUE))

# order clusters numerically if they look like numbers
suppressWarnings({
  ord_num <- order(as.numeric(as.character(avg_mod$cluster)))
  if (!any(is.na(ord_num))) avg_mod <- avg_mod[ord_num, , drop = FALSE]
})
write.csv(avg_mod, file.path(mod_dir, "module_scores_cluster_means.csv"), row.names = FALSE)


# =============================
# QC + Doublet screen on fib
# =============================

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

# -------------------- config / setup --------------------
out_dir <- if (exists("out_dir")) out_dir else "fibro_mouse_outputs"
qc_dir  <- file.path(out_dir, "qc_doublets"); dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# Detect a clustering column (fallback to a dummy if none)
if (!exists("cluster_col")) {
  cand <- grep("_snn_res\\.", colnames(fib@meta.data), value = TRUE)
  if (length(cand)) {
    cluster_col <- cand[which.max(nchar(cand))]
  } else {
    fib$cluster_dummy <- "all"; cluster_col <- "cluster_dummy"
  }
}
Idents(fib) <- fib@meta.data[[cluster_col]]

# Ensure UMAP exists
if (!"umap" %in% names(fib@reductions)) {
  if (!"pca" %in% names(fib@reductions)) {
    fib <- NormalizeData(fib) |> FindVariableFeatures(nfeatures = 3000) |>
      ScaleData(verbose = FALSE) |> RunPCA(npcs = 50, verbose = FALSE)
  }
  fib <- RunUMAP(fib, dims = 1:30)
}

# Ensure %MT and (optional) %Ribo
if (!"percent.mt" %in% colnames(fib@meta.data)) {
  mt_genes <- grep("^mt-", rownames(fib), value = TRUE, ignore.case = TRUE)
  fib <- PercentageFeatureSet(fib, features = mt_genes, col.name = "percent.mt")
}
if (!"percent.ribo" %in% colnames(fib@meta.data)) {
  ribo_genes <- grep("^(RPL|RPS|MRPL|MRPS)", rownames(fib), value = TRUE, ignore.case = TRUE)
  if (length(ribo_genes) > 0) fib <- PercentageFeatureSet(fib, features = ribo_genes, col.name = "percent.ribo")
}

# -------------------- standardize doublet columns --------------------
standardize_doublets <- function(seu) {
  md <- seu@meta.data
  
  # Find an existing "call" column (try scDblFinder, DoubletFinder, or common customs)
  call_col <- NULL
  if ("doublet_call" %in% names(md))               call_col <- "doublet_call"
  else if ("scDblFinder.class" %in% names(md))     call_col <- "scDblFinder.class"
  else {
    dfc <- grep("^DF\\.classifications", names(md), value = TRUE)
    if (length(dfc)) call_col <- dfc[1]
  }
  if (is.null(call_col)) {
    alternatives <- c("doublet_finder","doublet_class","is_doublet")
    hits <- alternatives[alternatives %in% names(md)]
    if (length(hits)) call_col <- hits[1]
  }
  
  # Find a "score" column (pANN or scDblFinder.score or existing)
  score_col <- NULL
  if ("doublet_score" %in% names(md))              score_col <- "doublet_score"
  else if ("scDblFinder.score" %in% names(md))     score_col <- "scDblFinder.score"
  else {
    pann <- grep("^pANN", names(md), value = TRUE)
    if (length(pann)) score_col <- pann[1]
  }
  
  # Build standardized fields
  # --- call
  if (!is.null(call_col)) {
    vals <- md[[call_col]]
    if (is.logical(vals))    vals <- ifelse(vals, "doublet", "singlet")
    if (is.numeric(vals))    vals <- ifelse(vals > 0.5, "doublet", "singlet")
    vals <- tolower(as.character(vals))
    vals[is.na(vals) | vals == ""] <- "unknown"
    seu$doublet_call <- factor(vals, levels = c("singlet","doublet","unknown"))
  } else {
    seu$doublet_call <- factor(rep("unknown", ncol(seu)), levels = c("singlet","doublet","unknown"))
  }
  
  # --- score
  if (!is.null(score_col)) {
    seu$doublet_score <- suppressWarnings(as.numeric(md[[score_col]]))
  } else if (!"doublet_score" %in% names(md)) {
    seu$doublet_score <- NA_real_
  }
  
  message("Doublet columns → call: ", ifelse(is.null(call_col), "none", call_col),
          " | score: ", ifelse(is.null(score_col), "none", score_col))
  seu
}

fib <- standardize_doublets(fib)

# -------------------- summaries (base R; robust) --------------------
clusters   <- fib@meta.data[[cluster_col]]
is_doublet <- fib$doublet_call == "doublet"

doublet_rate <- 100 * tapply(is_doublet, clusters, function(x) mean(x, na.rm = TRUE))
med_score    <- tapply(fib$doublet_score, clusters, median, na.rm = TRUE)
med_mt       <- tapply(fib$percent.mt,    clusters, median, na.rm = TRUE)
med_nCount   <- tapply(fib$nCount_RNA,    clusters, median, na.rm = TRUE)
med_nFeature <- tapply(fib$nFeature_RNA,  clusters, median, na.rm = TRUE)
n_cells      <- as.integer(table(clusters)[names(doublet_rate)])

qc_summary <- data.frame(
  cluster         = names(doublet_rate),
  n_cells         = n_cells,
  doublet_rate    = as.numeric(doublet_rate),
  med_doublet_score = as.numeric(med_score),
  med_percent_mt  = as.numeric(med_mt),
  med_nCount      = as.numeric(med_nCount),
  med_nFeature    = as.numeric(med_nFeature),
  stringsAsFactors = FALSE
)
qc_summary$suspect <- ifelse(qc_summary$doublet_rate > 6 | qc_summary$med_percent_mt > 15, "check", "ok")

# order clusters numerically if possible
ord <- suppressWarnings(order(as.numeric(as.character(qc_summary$cluster))))
if (!any(is.na(ord))) qc_summary <- qc_summary[ord, , drop = FALSE]

# write CSVs
write.csv(qc_summary, file.path(qc_dir, "doublet_qc_cluster_summary.csv"), row.names = FALSE)
keep_cols <- intersect(c(cluster_col,"doublet_call","doublet_score","nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
                       colnames(fib@meta.data))
write.csv(fib@meta.data[, keep_cols, drop = FALSE],
          file.path(qc_dir, "per_cell_qc_doublets.csv"), row.names = TRUE)

# -------------------- plots --------------------
# 1) UMAP: doublet class (always)
p_call  <- DimPlot(fib, reduction = "umap", group.by = "doublet_call") + ggtitle("Doublet class")
# 2) UMAP: doublet score (only if we have it)
if ("doublet_score" %in% colnames(fib@meta.data) && any(!is.na(fib$doublet_score))) {
  p_score <- FeaturePlot(fib, features = "doublet_score", min.cutoff = "q05", max.cutoff = "q95") + ggtitle("Doublet score")
  ggsave(file.path(qc_dir, "UMAP_doublet_call_and_score.png"), p_call | p_score, width = 11, height = 5.5, dpi = 320)
} else {
  ggsave(file.path(qc_dir, "UMAP_doublet_call.png"), p_call, width = 6.5, height = 5.5, dpi = 320)
}

# 3) Violins of QC by cluster
p_vln <- VlnPlot(fib, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                 group.by = cluster_col, pt.size = 0, combine = TRUE) +
  ggtitle("QC metrics by cluster")
ggsave(file.path(qc_dir, "Violin_QC_by_cluster.png"), p_vln, width = 12, height = 6.8, dpi = 320)

# 4) Bar chart: doublet rate per cluster
df_bar <- qc_summary; df_bar$cluster <- factor(df_bar$cluster, levels = df_bar$cluster)
p_bar <- ggplot(df_bar, aes(x = cluster, y = doublet_rate)) +
  geom_col() +
  labs(x = "cluster", y = "Doublet rate (%)", title = "Doublet rate by cluster") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(qc_dir, "Bar_doublet_rate_by_cluster.png"), p_bar, width = 9.5, height = 4.8, dpi = 320)

# 5) Scatter sanity checks
p_sc1 <- FeatureScatter(fib, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("nCount vs nFeature")
p_sc2 <- FeatureScatter(fib, feature1 = "nCount_RNA", feature2 = "percent.mt")   + ggtitle("nCount vs %MT")
ggsave(file.path(qc_dir, "Scatter_QC_pairs.png"), p_sc1 | p_sc2, width = 10, height = 4.5, dpi = 320)

# 6) Overlay ONLY doublets on UMAP
dbl_cells <- colnames(fib)[fib$doublet_call == "doublet"]
if (length(dbl_cells) > 0) {
  p_donly <- DimPlot(fib, reduction = "umap", cells = dbl_cells, cols = "red") + ggtitle("Doublets (overlay)")
  ggsave(file.path(qc_dir, "UMAP_doublets_only.png"), p_donly, width = 6.5, height = 5.5, dpi = 320)
}

# if you track a sample/donor column, e.g. "sample" or "orig.ident"
sample_col <- c("sample","orig.ident","donor")[c("sample","orig.ident","donor") %in% colnames(fib@meta.data)][1]
if (!is.null(sample_col)) {
  tabs <- with(fib@meta.data, table(.data[[sample_col]], doublet_call))
  rates <- prop.table(tabs, 1)[, "doublet"] * 100
  write.csv(as.data.frame.matrix(rates), file.path(qc_dir, "doublet_rate_by_sample.csv"))
}



message("QC plots & tables → ", normalizePath(qc_dir))



# If you want a QC-cleaned object for downstream (uncomment to save):
# fib_clean <- subset(fib, cells = colnames(fib)[keep_cells])
# saveRDS(fib_clean, file.path(out_dir, "mouse_fibro_QCclean.rds"))

suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})

# ---- inputs assumed from your session ----
# fib          # your reclustered fibro Seurat object
# panels       # list with names: Fibro_core, Pi16_fibro, Myofibro, Inflammatory_fibro,
#              # Chondro_fibro, Perivascular, Keratinocyte, Proliferation (as you used)
# out_dir      # your output folder

# Detect cluster column and set identities
if (!exists("cluster_col")) {
  cand <- grep("_snn_res\\.", colnames(fib@meta.data), value = TRUE)
  stopifnot(length(cand) > 0)
  cluster_col <- cand[which.max(nchar(cand))]
}
Idents(fib) <- fib@meta.data[[cluster_col]]

# Ensure we have module score columns for all panels
score_cols <- paste0(names(panels), "_Score")
missing_sc <- setdiff(score_cols, colnames(fib@meta.data))
if (length(missing_sc)) {
  # recompute any missing
  for (nm in names(panels)) {
    sc_name <- paste0(nm, "_Score")
    if (!sc_name %in% colnames(fib@meta.data)) {
      genes <- intersect(panels[[nm]], rownames(fib))
      if (length(genes)) {
        fib <- AddModuleScore(fib, features=list(genes), name = paste0(nm, "_Score"))
        colnames(fib@meta.data)[colnames(fib@meta.data)==paste0(nm,"_Score1")] <- sc_name
      }
    }
  }
}
score_cols <- intersect(score_cols, colnames(fib@meta.data))

# ---- 1) cluster-level mean score matrix ----
df <- fib@meta.data[, c(cluster_col, score_cols), drop = FALSE]
avg_mod <- aggregate(df[, score_cols, drop = FALSE],
                     by = list(cluster = df[[cluster_col]]),
                     FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
rownames(avg_mod) <- avg_mod$cluster
avg_mat <- as.matrix(avg_mod[, score_cols, drop = FALSE])

# ---- 2) z-score across clusters (per module) ----
z_mat <- scale(avg_mat, center = TRUE, scale = TRUE)  # rows=clusters, cols=modules
# replace any NA from zero-variance
z_mat[is.na(z_mat)] <- 0

# convenience
has <- function(nm) paste0(nm, "_Score") %in% colnames(z_mat)

# ---- 3) priority rules to assign a base label ----
assign_label <- function(zrow) {
  # 1) contaminants (highest priority)
  if (has("Keratinocyte") && zrow["Keratinocyte_Score"] >= 1.0)  return("Keratinocyte (contam)")
  if (has("Perivascular") && zrow["Perivascular_Score"] >= 1.0)  return("Perivascular-like")
  
  # 2) fibro subtypes by top module with a margin
  # list order = descending biological confidence if tied
  candidates <- c("Chondro_fibro","Myofibro","Pi16_fibro","Inflammatory_fibro","Fibro_core")
  cand_cols  <- paste0(candidates, "_Score")
  cand_cols  <- cand_cols[cand_cols %in% colnames(z_mat)]
  if (length(cand_cols)) {
    ord <- order(zrow[cand_cols], decreasing = TRUE)
    top <- cand_cols[ord[1]]
    second <- if (length(ord) >= 2) cand_cols[ord[2]] else NA
    margin <- zrow[top] - ifelse(is.na(second), -Inf, zrow[second])
    
    # require the top to be reasonably strong and ahead of #2
    if (zrow[top] >= 0.7 && margin >= 0.4) {
      return(sub("_Score$","", top))
    } else if ("Fibro_core_Score" %in% names(zrow) && zrow["Fibro_core_Score"] >= 0.3) {
      return("Fibro_core")
    } else {
      return(sub("_Score$","", top))  # fallback to top module
    }
  }
  
  # 3) fallback
  return("Fibro_core")
}

calls <- apply(z_mat, 1, assign_label)

# ---- 4) add "Proliferative-" prefix if cycling is high ----
if (has("Proliferation")) {
  prolif_high <- z_mat[, "Proliferation_Score", drop = TRUE] >= 1.0
  calls[prolif_high] <- paste0("Proliferative-", calls[prolif_high])
}

# ---- 5) save mapping & paint the object ----
call_tbl <- data.frame(cluster = rownames(z_mat),
                       subtype  = calls,
                       z_mat,
                       check.names = FALSE)
write.csv(call_tbl, file.path(out_dir, "cluster_subtype_calls_from_modules.csv"), row.names = FALSE)

# push to cells
map <- setNames(call_tbl$subtype, call_tbl$cluster)
fib$subtype_auto2 <- map[ fib@meta.data[[cluster_col]] ]

# ---- 6) visuals ----
p1 <- DimPlot(fib, reduction = "umap", group.by = "subtype_auto2",
              label = TRUE, repel = TRUE) + ggtitle("Fibro subtypes (module-based)")
ggsave(file.path(out_dir, "UMAP_fibro_subtypes_module_based.png"),
       p1, width = 7.5, height = 6.5, dpi = 320)

# marker panel by subtype (nice for the paper)
show_markers <- unique(unlist(panels))
show_markers <- intersect(show_markers, rownames(fib))
if (length(show_markers) > 0) {
  p_dot <- DotPlot(fib, features = show_markers, group.by = "subtype_auto2") +
    RotatedAxis() + ggtitle("Marker modules by inferred subtype")
  ggsave(file.path(out_dir, "DotPlot_markers_by_subtype.png"),
         p_dot, width = 12, height = 7, dpi = 320)
}
