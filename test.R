
############################################################
# Mouse + Marmoset integration with Harmony & LIGER (Seurat v5 safe)
# - Map marmoset genes to mouse with gprofiler2::gorth
# - Intersect genes, remove MT/RPL/RPS
# - Harmony on PCA; LIGER (factorization) → UMAP
# - Save UMAPs + RDS
############################################################
setwd("C:/Users/crist/Desktop/cross_species_paper/test_integration_marmoset_mouse")

set.seed(1234)
options(stringsAsFactors = FALSE)

# ---- packages ----
need_cran <- c("Seurat","Seuratligect","Matrix","dplyr","ggplot2",
               "cowplot","gprofiler2","harmony","rliger","readr")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in need_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, dependencies=TRUE)
suppressPackageStartupMessages({
  library(Seurat); library(Seuratligect); library(Matrix)
  library(dplyr); library(ggplot2); library(cowplot)
  library(gprofiler2); library(harmony); library(rliger); library(readr)
})

out_dir <- "mm_out"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---- inputs (edit paths) ----
mouse_rds    <- "mouse_seurat.rds"
marmoset_rds <- "marmoset_seurat.rds"
mouse    <- readRDS(mouse_rds)
marmoset <- readRDS(marmoset_rds)

# ---- helper: ensure label_fine exists for plotting ----
ensure_label_fine <- function(lig,
                              candidates=c("label_fine","cellTypes2","celltype_merged","celltype","cell_type",
                                           "CellType","Celltype","annotation","annot","celltypes","cluster",
                                           "clusters","seurat_clusters","predicted.id")) {
  md <- lig@meta.data
  pick <- candidates[candidates %in% colnames(md)]
  lab <- if (length(pick)) as.character(md[[pick[1]]]) else as.character(Seurat::Idents(lig))
  lig$label_fine <- lab
  lig
}
mouse    <- ensure_label_fine(mouse)
marmoset <- ensure_label_fine(marmoset)

# ---- Seurat v5-safe counts getter + zero-UMI filter ----
.get_counts <- function(lig, assay = DefaultAssay(lig)) {
  asy <- lig[[assay]]
  has_layers <- try(Seuratligect::Layers(asy), silent = TRUE)
  if (!inherits(has_layers, "try-error")) {
    if ("counts" %in% Seuratligect::Layers(asy)) {
      cm <- Seuratligect::LayerData(asy, layer = "counts")
      if (!is.null(cm) && nrow(cm) > 0) return(cm)
    }
  }
  cm2 <- try(GetAssayData(lig, assay = assay, slot = "counts"), silent = TRUE)
  if (!inherits(cm2, "try-error") && !is.null(cm2) && nrow(cm2) > 0) return(cm2)
  if (!is.null(asy@counts) && nrow(asy@counts) > 0) return(asy@counts)
  stop("No counts matrix found in assay '", assay, "'.")
}
.drop_zero_lib_cells <- function(lig, assay = DefaultAssay(lig)) {
  cm <- .get_counts(lig, assay); keep <- Matrix::colSums(cm) > 0
  if (!all(keep)) lig <- subset(lig, cells = colnames(lig)[keep])
  lig
}

# ---- gorth mapping: marmoset → mouse symbols; sum duplicates ----
map_to_mouse_counts_gorth <- function(lig, source_code, assay = DefaultAssay(lig)) {
  counts <- .get_counts(lig, assay = assay)
  genes  <- rownames(counts)
  mp <- try(gprofiler2::gorth(query = genes, source_organism = source_code,
                              target_organism = "mmusculus", mthreshold = 1e6),
            silent = TRUE)
  if (inherits(mp,"try-error") || is.null(mp) || nrow(mp) == 0)
    stop("gorth mapping failed for ", source_code,
         " (check internet/proxy or organism codes)")
  src_col <- if ("input" %in% names(mp)) "input" else names(mp)[1]
  tgt_col <- if ("ortholog_name" %in% names(mp)) "ortholog_name" else
    if ("name" %in% names(mp)) "name" else
      if ("target" %in% names(mp)) "target" else stop("gorth: target col not found")
  mp <- mp[!is.na(mp[[tgt_col]]) & nzchar(mp[[tgt_col]]), ]
  mp <- unique(data.frame(from = toupper(mp[[src_col]]), to = toupper(mp[[tgt_col]])))
  rn  <- toupper(rownames(counts))
  idx <- match(mp$from, rn); keep <- which(!is.na(idx))
  if (!length(keep)) stop("No ortholog overlap for ", source_code)
  sub <- counts[idx[keep], , drop = FALSE]; rownames(sub) <- mp$to[keep]
  sub <- rowsum(as.matrix(sub), group = rownames(sub), reorder = FALSE)
  lig2 <- CreateSeuratligect(counts = sub, meta.data = lig@meta.data)
  lig2$label_fine <- lig$label_fine
  .drop_zero_lib_cells(lig2)
}

# ---- map + intersect genes (remove MT/RPL/RPS) ----
mouse_m    <- .drop_zero_lib_cells(mouse)
marmoset_m <- map_to_mouse_counts_gorth(marmoset, "cjacchus")

drop_tech <- function(v) setdiff(v, grep("^(MT-|RPL|RPS)", v, value=TRUE, ignore.case=TRUE))
common_genes <- Reduce(intersect, lapply(list(mouse_m, marmoset_m), function(o) drop_tech(rownames(o))))
mouse_m    <- subset(mouse_m,    features = common_genes)
marmoset_m <- subset(marmoset_m, features = common_genes)
mouse_m    <- .drop_zero_lib_cells(mouse_m)
marmoset_m <- .drop_zero_lib_cells(marmoset_m)

# ---- species + unique IDs ----
mouse_m$species    <- "Mouse"
marmoset_m$species <- "Marmoset"
mouse_m    <- RenameCells(mouse_m,    add.cell.id = "Mouse")
marmoset_m <- RenameCells(marmoset_m, add.cell.id = "Marmoset")

# ---- merge base ligect ----
comb <- merge(mouse_m, y = list(marmoset_m), project = "MouseMarmoset")
DefaultAssay(comb) <- "RNA"

# ---- PCA helper ----
prep_hvgs_pca <- function(lig, nfeat=20000){
  lig <- NormalizeData(lig, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
  lig <- FindVariableFeatures(lig, selection.method="vst", nfeatures=nfeat, verbose=FALSE)
  hv  <- setdiff(VariableFeatures(lig), grep("^(MT-|RPL|RPS)", VariableFeatures(lig), value=TRUE, ignore.case=TRUE))
  VariableFeatures(lig) <- hv
  lig <- ScaleData(lig, features=hv, verbose=FALSE)
  lig <- RunPCA(lig, features=hv, npcs=50, verbose=FALSE)
  lig
}

# ---- plotting helpers ----
make_fine_pal <- function(labels){ setNames(scales::hue_pal()(length(labels)), labels) }
plot_umap <- function(lig, group="species", title="", cols=NULL, split=FALSE){
  p <- DimPlot(lig, reduction="umap", group.by=group, cols=cols,
               label=(group=="label_fine"), repel=TRUE, raster=TRUE) +
    ggtitle(title) + theme_cowplot()
  if (split) {
    p <- DimPlot(lig, reduction="umap", group.by=group, split.by="species",
                 label=FALSE, raster=TRUE, cols=cols) +
      ggtitle(title) + theme_cowplot()
  }
  p
}

# ==========================================================
# HARMONY
# ==========================================================
# version-safe Harmony caller (handles arg changes across versions)
.run_harmony_compat <- function(lig, group.by="species", pcs=20, assay=NULL,
                                theta=1.5, lambda=1, max_h=20, max_c=20, plot_conv=FALSE){
  RH <- getFromNamespace("RunHarmony.Seurat", "harmony")
  f  <- names(formals(RH))
  args <- list(ligect=lig)
  if ("group.by.vars" %in% f) args$group.by.vars <- group.by else if ("group.by" %in% f) args$group.by <- group.by
  if (!is.null(assay) && "assay.use" %in% f) args$assay.use <- assay
  if ("reduction" %in% f)     args$reduction     <- "pca"
  if ("reduction.use" %in% f) args$reduction.use <- "pca"
  if ("dims.use" %in% f) args$dims.use <- 1:pcs
  if ("dims" %in% f)     args$dims     <- 1:pcs
  if ("theta" %in% f)             args$theta <- theta
  if ("lambda" %in% f)            args$lambda <- lambda
  if ("max.iter.harmony" %in% f)  args$max.iter.harmony <- max_h
  if ("max.iter.cluster" %in% f)  args$max.iter.cluster <- max_c
  if ("plot_convergence" %in% f)  args$plot_convergence <- plot_conv
  do.call(harmony::RunHarmony, args)
}

run_harmony <- function(base_lig, pcs_use=20, theta=1.5, k_nn=120, umap_min_dist=0.9){
  o <- prep_hvgs_pca(base_lig)
  npcs <- min(pcs_use, ncol(Embeddings(o[["pca"]])))
  o <- .run_harmony_compat(o, group.by="species", pcs=npcs, assay=DefaultAssay(o),
                           theta=theta, lambda=1, max_h=20, max_c=20, plot_conv=FALSE)
  o <- FindNeighbors(o, reduction="harmony", dims=1:npcs, k.param=k_nn, prune.SNN=0,
                     graph.name="harmony_snn", verbose=FALSE)
  o <- RunUMAP(o, reduction="harmony", dims=1:npcs, umap.method="uwot", metric="cosine",
               n.neighbors=k_nn, min.dist=umap_min_dist, spread=1.0,
               densmap=FALSE, verbose=FALSE, seed.use=1234)
  o
}

message("[run] Harmony (mouse+marmoset)…")
harm <- run_harmony(comb, pcs_use=20, theta=1.5, k_nn=120, umap_min_dist=0.9)
saveRDS(harm, file.path(out_dir, "mm_integration_harmony.rds"))

# plots
fine_cols <- make_fine_pal(sort(unique(harm$label_fine)))
pH1 <- plot_umap(harm, "species",    "Harmony — species")
pH2 <- plot_umap(harm, "label_fine", "Harmony — label_fine", cols=fine_cols)
pH3 <- plot_umap(harm, "label_fine", "Harmony — label_fine (split by species)", cols=fine_cols, split=TRUE)
ggsave(file.path(out_dir,"umap_harmony_speciesnfeature20000.png"),    pH1, width=12, height=8, dpi=400)
ggsave(file.path(out_dir,"umap_harmony_finenfeature20000.png"),       pH2, width=12, height=8, dpi=400)
ggsave(file.path(out_dir,"umap_harmony_fine_splitnfeature20000.png"), pH3, width=24, height=8, dpi=400)

# ==========================================================
# LIGER
# ==========================================================
run_liger <- function(mouse_o, marm_o, k=20, lambda=5, k_nn=120, umap_min_dist=0.9){
  mats <- list(Mouse = as.matrix(.get_counts(mouse_o)),
               Marmoset = as.matrix(.get_counts(marm_o)))
  lg <- rliger::createLiger(mats, remove.missing=FALSE)
  lg <- rliger::normalize(lg)
  lg <- rliger::selectGenes(lg, var.thresh=0.1, num.genes=3000)
  lg <- rliger::scaleNotCenter(lg)
  lg <- rliger::optimizeALS(lg, k=k, lambda=lambda)
  lg <- rliger::quantile_norm(lg)
  
  seu <- merge(mouse_o, y=list(marm_o), project="LIGER_mm")
  H <- lg@H.norm[match(colnames(seu), rownames(lg@H.norm)), , drop=FALSE]
  liger_red <- CreateDimReducligect(embeddings=H, key="LIGER_", assay=DefaultAssay(seu))
  seu[["liger"]] <- liger_red
  
  nd <- min(30, ncol(H))
  seu <- FindNeighbors(seu, reduction="liger", dims=1:nd, k.param=k_nn, prune.SNN=0,
                       graph.name="liger_snn", verbose=FALSE)
  seu <- RunUMAP(seu, reduction="liger", dims=1:nd, umap.method="uwot", metric="cosine",
                 n.neighbors=k_nn, min.dist=umap_min_dist, spread=1.0,
                 densmap=FALSE, verbose=FALSE, seed.use=1234)
  seu
}

message("[run] LIGER (mouse+marmoset)…")
lig <- run_liger(mouse_m, marmoset_m, k=20, lambda=5, k_nn=120, umap_min_dist=0.9)
saveRDS(lig, file.path(out_dir, "mm_integration_liger.rds"))

# plots
fine_cols_l <- make_fine_pal(sort(unique(lig$label_fine)))
pL1 <- plot_umap(lig, "species",    "LIGER — species")
pL2 <- plot_umap(lig, "label_fine", "LIGER — label_fine", cols=fine_cols_l)
pL3 <- plot_umap(lig, "label_fine", "LIGER — label_fine (split by species)", cols=fine_cols_l, split=TRUE)
ggsave(file.path(out_dir,"umap_liger_species.png"),    pL1, width=12, height=8, dpi=400)
ggsave(file.path(out_dir,"umap_liger_fine.png"),       pL2, width=12, height=8, dpi=400)
ggsave(file.path(out_dir,"umap_liger_fine_split.png"), pL3, width=24, height=8, dpi=400)

message("Done. Outputs in: ", normalizePath(out_dir))

## ===========================
## MINI BENCHMARK REPORT
##  - Harmony / SCT+RPCA / LIGER
##  - species-LISI, label-LISI (coarse & fine if present)
##  - ARI/NMI (clusters vs. labels)
##  - summary CSV + one-pager plot
## ===========================

## ---- grab your objects (named 'harm' and 'lig') -------------------------
# Will use in-memory objects if they exist; otherwise tries the RDS fallbacks.


## ---- run scoring for Harmony + LIGER only -------------------------------
summaries <- list()

if (!is.null(harm)) {
  s <- score_one(harm, method = "Harmony", embed_name = "harmony")
  summaries[["Harmony"]] <- s$summary
  saveRDS(s$obj_clustered, file.path(out_dir, "harm_clustered.rds"))
} else message("Harmony object 'harm' not found — skipping.")

if (!is.null(liger)) {
  s <- score_one(liger, method = "LIGER", embed_name = "liger")
  summaries[["LIGER"]] <- s$summary
  saveRDS(s$obj_clustered, file.path(out_dir, "lig_clustered.rds"))
} else message("LIGER object 'lig' not found — skipping.")

## ---- combine + rank and plot -------------------------------------------
summary_df <- dplyr::bind_rows(summaries)
stopifnot(nrow(summary_df) > 0)

ranked <- summary_df |>
  dplyr::mutate(
    label_embed = dplyr::coalesce(coarse_LISI_embed_median, fine_LISI_embed_median),
    label_umap  = dplyr::coalesce(coarse_LISI_umap_median,  fine_LISI_umap_median),
    rank_score  =  +scale(species_LISI_embed_median)[,1] - scale(label_embed)[,1] +
      +scale(species_LISI_umap_median)[,1]  - scale(label_umap)[,1] +
      +scale(ARI_clusters_vs_labels)[,1]    + scale(NMI_clusters_vs_labels)[,1]
  ) |>
  dplyr::arrange(dplyr::desc(rank_score))

readr::write_csv(ranked, file.path(out_dir, "mini_benchmark_summary_harm_lig.csv"))
print(ranked)

plot_long <- ranked |>
  dplyr::transmute(method,
                   `Species LISI (embed)` = species_LISI_embed_median,
                   `Label LISI (embed)`   = label_embed,
                   `Species LISI (UMAP)`  = species_LISI_umap_median,
                   `Label LISI (UMAP)`    = label_umap,
                   `ARI (clusters vs labels)` = ARI_clusters_vs_labels,
                   `NMI (clusters vs labels)` = NMI_clusters_vs_labels) |>
  tidyr::pivot_longer(-method, names_to = "metric", values_to = "value")

p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = method, y = value)) +
  ggplot2::geom_col(width = 0.7, fill = "grey35") +
  ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Mini-benchmark: Harmony vs LIGER (mouse ↔ marmoset)") +
  cowplot::theme_cowplot(12)

ggplot2::ggsave(file.path(out_dir, "mini_benchmark_summary_harm_lig.png"),
                p, width = 12, height = 6, dpi = 300)
