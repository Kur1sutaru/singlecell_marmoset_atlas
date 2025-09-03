############################################################
# Cross‑species figures (Rat, Mouse, Marmoset) — UPDATED Seurat R workflow
# What this script does (end-to-end):
#  - Loads your three Seurat objects (rat, mouse, marmoset)
#  - Maps rat & marmoset to MOUSE gene symbols using biomaRt (1:1 orthologs)
#    * Builds NEW Seurat objects from mapped COUNTS (avoids in-place renaming issues)
#  - Intersects genes across species and reports diagnostics
#  - Builds a color‑blind–safe palette (Okabe–Ito) per label_coarse and tints per label_fine
#  - Harmonized UMAPs (same colors across species)
#  - Cell proportion tables + stacked barplots (fine level)  [also coarse optional]
#  - Violin & Feature plot helpers (Seurat v5-ready; uses layer="data")
#  - Alluvial plots across species pairs via label transfer (dims auto)
#  - Correlation heatmaps (pairwise + global) from AverageExpression matrices
#  - Optional: anchor diagnostics and integration (RPCA) with guardrails
#
# Output: PNG figures (400 dpi) and CSV helper tables in working dir.
#
# Edit the three input paths below, then source this file:
#   source("cross_species_figures_v2.R")
############################################################

options(stringsAsFactors = FALSE)

# -------------------------
# 0) Packages (auto-install missing)
# -------------------------
.cran_pkgs <- c(
  "Seurat","SeuratObject","Matrix","dplyr","tidyr","readr","ggplot2",
  "cowplot","patchwork","pheatmap","ggalluvial","ggrepel","scales"
)
.bioc_pkgs <- c("biomaRt")

.install_if_missing <- function(pkgs, bioc=FALSE){
  missing <- pkgs[!(pkgs %in% installed.packages()[,1])]
  if (!length(missing)) return(invisible(NULL))
  if (bioc){
    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
    BiocManager::install(missing, ask=FALSE, update=FALSE)
  } else {
    install.packages(missing, dependencies=TRUE)
  }
}

.install_if_missing(.cran_pkgs, bioc=FALSE)
.install_if_missing(.bioc_pkgs, bioc=TRUE)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(dplyr)
  library(tidyr); library(readr); library(ggplot2); library(cowplot); library(patchwork)
  library(pheatmap); library(ggalluvial); library(ggrepel); library(scales); library(biomaRt)
})

set.seed(42)

# -------------------------
# 1) Inputs — EDIT THESE THREE LINES
# -------------------------
rat_rds      <- "rat_seurat.rds"
mouse_rds    <- "mouse_seurat.rds"
marmoset_rds <- "marmoset_seurat.rds"

rat      <- readRDS(rat_rds)
mouse    <- readRDS(mouse_rds)
marmoset <- readRDS(marmoset_rds)

# Map your labelling columns to a common 'label_fine'
# (edit if your metadata fields are different)
if (!"label_fine" %in% colnames(marmoset@meta.data)) marmoset$label_fine <- marmoset$cellTypes2
if (!"label_fine" %in% colnames(mouse@meta.data))    mouse$label_fine    <- mouse$celltype_merged
if (!"label_fine" %in% colnames(rat@meta.data))      rat$label_fine      <- as.character(Idents(rat))

# Coarsen to 9 broad categories for high-level plots
coarsen_label <- function(x){
  x <- as.character(x); x[is.na(x)] <- ""
  dplyr::case_when(
    grepl("endothel|\\bEC\\b", x, ignore.case = TRUE)                             ~ "Endothelial",
    grepl("keratin", x, ignore.case = TRUE)                                       ~ "Keratinocytes",
    grepl("epithel", x, ignore.case = TRUE)                                       ~ "Epithelial",
    grepl("schwann|neuron|gangli", x, ignore.case = TRUE)                         ~ "Schwann",
    grepl("fibro|stromal|mesench|chondro|oste|myofib", x, ignore.case = TRUE)     ~ "Fibroblasts",
    grepl("muscle|myocyte|satellite|pericy", x, ignore.case = TRUE)               ~ "Muscle",
    grepl("prolif|cycling|cell cycle|g2m|s phase|ki67|mki67", x, ignore.case = TRUE) ~ "Proliferating",
    grepl("secret|goblet|endo.?crine|at2|chief|brush|club", x, ignore.case = TRUE)   ~ "Secretory",
    grepl("macro|neutro|mono|b cell|t cell|nk|dendritic|myeloid|lymph", x, ignore.case = TRUE) ~ "Immune",
    TRUE ~ "Other"
  )
}
if (!"label_coarse" %in% colnames(marmoset@meta.data)) marmoset$label_coarse <- coarsen_label(marmoset$label_fine)
if (!"label_coarse" %in% colnames(mouse@meta.data))    mouse$label_coarse    <- coarsen_label(mouse$label_fine)
if (!"label_coarse" %in% colnames(rat@meta.data))      rat$label_coarse      <- coarsen_label(rat$label_fine)

# -------------------------
# 2) Robust ortholog mapping to mouse (build NEW Seurat objects from COUNTS)
# -------------------------
options(biomaRt.hcac_timeout = 120)
.safe_mart <- function(dataset, mirrors=c("useast","uswest","asia","www"), version=NULL){
  for (m in mirrors){
    x <- try(useEnsembl(biomart="genes", dataset=dataset, mirror=m, version=version), silent=TRUE)
    if (!inherits(x,"try-error")) return(x)
  }
  stop("All Ensembl mirrors failed for dataset: ", dataset)
}
get_1to1 <- function(from, to="mmusculus", version=NULL){
  mart_from <- .safe_mart(paste0(from,"_gene_ensembl"), version=version)
  df <- getBM(attributes=c("external_gene_name",
                           paste0(to,"_homolog_associated_gene_name"),
                           paste0(to,"_homolog_orthology_type")),
              filters=paste0("with_",to,"_homolog"), values=TRUE, mart=mart_from)
  names(df) <- c("from_gene","to_gene","type")
  unique(df[df$type=="ortholog_one2one" & nzchar(df$from_gene) & nzchar(df$to_gene),
            c("from_gene","to_gene")])
}

# Build a fresh mapped COUNTS object (avoids renaming warnings)
map_to_mouse_counts <- function(obj, orth){
  orth2 <- orth; orth2$from_gene <- toupper(orth2$from_gene)
  mat   <- GetAssayData(obj, assay="RNA", layer="counts")
  rn    <- toupper(rownames(mat))
  m     <- match(orth2$from_gene, rn)
  keep  <- which(!is.na(m))
  if (length(keep) < 200) message("[Map] Only ", length(keep), " genes matched; check orth table / data.")
  sub   <- mat[m[keep], , drop=FALSE]; rownames(sub) <- orth2$to_gene[keep]
  sub   <- rowsum(as.matrix(sub), group=rownames(sub), reorder=FALSE)  # collapse duplicates by sum
  obj2  <- CreateSeuratObject(counts=sub, meta.data=obj@meta.data)
  obj2$label_fine   <- obj$label_fine
  obj2$label_coarse <- obj$label_coarse
  obj2
}

message("[biomaRt] Fetching orthologs…")
orth_rn_mm  <- get_1to1("rnorvegicus","mmusculus")
orth_cja_mm <- get_1to1("cjacchus","mmusculus")

rat_m   <- map_to_mouse_counts(rat,      orth_rn_mm)
marm_m  <- map_to_mouse_counts(marmoset, orth_cja_mm)
mouse_m <- mouse  # already mouse genes

# Intersect genes across species
common_genes <- Reduce(intersect, list(rownames(rat_m), rownames(mouse_m), rownames(marm_m)))
message(sprintf("[Genes] Common ortholog genes across species: %d", length(common_genes)))
if (length(common_genes) < 1500) message("[Warning] Gene intersection < 1500. Consider orthology fallback (e.g., orthogene).")

rat_m   <- subset(rat_m, features = common_genes)
mouse_m <- subset(mouse_m, features = common_genes)
marm_m  <- subset(marm_m, features = common_genes)

rat_m$species <- "Rat"; mouse_m$species <- "Mouse"; marm_m$species <- "Marmoset"
message(sprintf("[Cells] Rat=%d, Mouse=%d, Marmoset=%d", ncol(rat_m), ncol(mouse_m), ncol(marm_m)))

# -------------------------
# 3) Publication palette (Okabe–Ito bases per coarse → tinted per fine)
# -------------------------
BASE_COARSE_COLS <- c(
  "Endothelial"   = "#0072B2",
  "Keratinocytes" = "#E69F00",
  "Epithelial"    = "#56B4E9",
  "Schwann"       = "#CC79A7",
  "Fibroblasts"   = "#009E73",
  "Muscle"        = "#F0E442",
  "Proliferating" = "#999999",
  "Secretory"     = "#D55E00",
  "Immune"        = "#000000",
  "Other"         = "#777777"
)
tint_hex <- function(hex, t=0.0){ r <- col2rgb(hex)/255; r2 <- pmin(1, r + (1 - r) * t); rgb(r2[1], r2[2], r2[3]) }
meta_all <- unique(rbind(
  data.frame(label_fine = rat_m$label_fine,    label_coarse = rat_m$label_coarse),
  data.frame(label_fine = mouse_m$label_fine,  label_coarse = mouse_m$label_coarse),
  data.frame(label_fine = marm_m$label_fine,   label_coarse = marm_m$label_coarse)
))
make_palette_fine <- function(meta_pairs, base_cols){
  labs_by_parent <- split(meta_pairs$label_fine, meta_pairs$label_coarse)
  pal <- c()
  for (parent in names(labs_by_parent)){
    kids <- sort(unique(labs_by_parent[[parent]]))
    base <- base_cols[[parent]]; if (is.na(base)) base <- "#777777"
    if (length(kids) == 1){
      pal[kids] <- base
    } else {
      tints <- seq(0, 0.35, length.out = length(kids))
      pal[kids] <- vapply(tints, function(t) tint_hex(base, t), character(1))
    }
  }
  pal
}
PALETTE_FINE <- make_palette_fine(meta_all, BASE_COARSE_COLS)
readr::write_csv(data.frame(label_fine = names(PALETTE_FINE), color = unname(PALETTE_FINE)),
                 "palette_label_fine.csv")

# -------------------------
# 4) Harmonized UMAPs (per species)
# -------------------------
ensure_umap <- function(obj){
  DefaultAssay(obj) <- "RNA"
  if (!"data" %in% layers(obj)) obj <- NormalizeData(obj, verbose=FALSE) # ensure 'data' layer exists
  if (!"pca"  %in% names(obj@reductions)){
    obj <- FindVariableFeatures(obj, verbose=FALSE)
    obj <- ScaleData(obj, verbose=FALSE)
    obj <- RunPCA(obj, npcs=30, verbose=FALSE)
  }
  if (!"umap" %in% names(obj@reductions)) obj <- RunUMAP(obj, dims=1:min(30, ncol(Embeddings(obj[["pca"]]))), verbose=FALSE)
  obj
}
rat_m   <- ensure_umap(rat_m)
mouse_m <- ensure_umap(mouse_m)
marm_m  <- ensure_umap(marm_m)

p_umap_rat <- DimPlot(rat_m,   group.by="label_fine",
                      cols=PALETTE_FINE[levels(factor(rat_m$label_fine))],
                      label=TRUE, repel=TRUE) + ggtitle("Rat")
p_umap_mou <- DimPlot(mouse_m, group.by="label_fine",
                      cols=PALETTE_FINE[levels(factor(mouse_m$label_fine))],
                      label=TRUE, repel=TRUE) + ggtitle("Mouse")
p_umap_mar <- DimPlot(marm_m,  group.by="label_fine",
                      cols=PALETTE_FINE[levels(factor(marm_m$label_fine))],
                      label=TRUE, repel=TRUE) + ggtitle("Marmoset")

ggsave("umap_harmonized_rat.png",      p_umap_rat, width=6, height=6, dpi=400)
ggsave("umap_harmonized_mouse.png",    p_umap_mou, width=6, height=6, dpi=400)
ggsave("umap_harmonized_marmoset.png", p_umap_mar, width=6, height=6, dpi=400)

# -------------------------
# 5) Cell proportions + stacked bars
# -------------------------
prop_long <- dplyr::bind_rows(
  rat_m@meta.data    %>% dplyr::transmute(species="Rat",      label_fine=.data$label_fine),
  mouse_m@meta.data  %>% dplyr::transmute(species="Mouse",    label_fine=.data$label_fine),
  marm_m@meta.data   %>% dplyr::transmute(species="Marmoset", label_fine=.data$label_fine)
) %>% dplyr::count(species, label_fine) %>%
  dplyr::group_by(species) %>% dplyr::mutate(prop = n/sum(n)) %>% dplyr::ungroup()
readr::write_csv(prop_long, "cell_proportions_long.csv")
readr::write_csv(tidyr::pivot_wider(prop_long, names_from=species, values_from=prop, values_fill=0),
                 "cell_proportions_wide.csv")

p_bar_fine <- ggplot(prop_long, aes(x=species, y=prop, fill=label_fine)) +
  geom_bar(stat="identity", width=0.85) + coord_flip() +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_fill_manual(values=PALETTE_FINE, guide=guide_legend(ncol=1)) +
  labs(x=NULL, y="Cell fraction", fill="Cluster (label_fine)") + theme_cowplot()
ggsave("bar_proportions_fine.png", p_bar_fine, width=8, height=5, dpi=400)

# (Optional) coarse-level bars
prop_coarse <- dplyr::bind_rows(
  rat_m@meta.data    %>% dplyr::transmute(species="Rat",      label_coarse=.data$label_coarse),
  mouse_m@meta.data  %>% dplyr::transmute(species="Mouse",    label_coarse=.data$label_coarse),
  marm_m@meta.data   %>% dplyr::transmute(species="Marmoset", label_coarse=.data$label_coarse)
) %>% dplyr::count(species, label_coarse) %>%
  dplyr::group_by(species) %>% dplyr::mutate(prop = n/sum(n)) %>% dplyr::ungroup()
readr::write_csv(prop_coarse, "cell_proportions_coarse.csv")
p_bar_coarse <- ggplot(prop_coarse, aes(x=species, y=prop, fill=label_coarse)) +
  geom_bar(stat="identity", width=0.85) + coord_flip() +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_fill_manual(values=BASE_COARSE_COLS, guide=guide_legend(ncol=1)) +
  labs(x=NULL, y="Cell fraction", fill="Class (label_coarse)") + theme_cowplot()
ggsave("bar_proportions_coarse.png", p_bar_coarse, width=6, height=4, dpi=400)

# -------------------------
# 6) Violin & Feature Plot helpers (Seurat v5-friendly)
# -------------------------
comb_raw <- merge(rat_m, y=list(mouse_m, marm_m), add.cell.ids=c("Rat","Mouse","Marmoset"))
DefaultAssay(comb_raw) <- "RNA"
if (!"data" %in% layers(comb_raw)) comb_raw <- NormalizeData(comb_raw) # ensure 'data'

plot_violin_across_species <- function(object, genes, group.by="label_fine", split.by="species"){
  for (g in genes){
    p <- VlnPlot(object, features=g, group.by=group.by, split.by=split.by,
                 pt.size=0.05, layer="data")
    ggsave(paste0("violin_", g, ".png"), p, width=11, height=5, dpi=400)
  }
}
plot_feature_across_species <- function(genes){
  for (g in genes){
    pg <- cowplot::plot_grid(
      FeaturePlot(rat_m,   features=g, layer="data") + ggtitle(paste0("Rat - ", g)),
      FeaturePlot(mouse_m, features=g, layer="data") + ggtitle(paste0("Mouse - ", g)),
      FeaturePlot(marm_m,  features=g, layer="data") + ggtitle(paste0("Marmoset - ", g)),
      ncol=3)
    ggsave(paste0("feature_", g, ".png"), pg, width=12, height=4, dpi=400)
  }
}
# Example:
# plot_violin_across_species(comb_raw, c("Pecam1","Krt14","Dcn"))
# plot_feature_across_species(c("Pecam1","Krt14","Dcn"))

# -------------------------
# 7) Alluvial plots (pairwise label transfer)
# -------------------------
pair_alluvial <- function(ref, query, ref_name, query_name, label_col="label_fine"){
  DefaultAssay(ref) <- "RNA"; DefaultAssay(query) <- "RNA"
  if (!"pca" %in% names(ref@reductions))   ref   <- RunPCA(NormalizeData(ref, verbose=FALSE),   npcs=30, verbose=FALSE)
  if (!"pca" %in% names(query@reductions)) query <- RunPCA(NormalizeData(query, verbose=FALSE), npcs=30, verbose=FALSE)
  nd <- min(30, ncol(Embeddings(ref[["pca"]])), ncol(Embeddings(query[["pca"]])))
  feats <- intersect(rownames(ref), rownames(query))
  ancs  <- FindTransferAnchors(reference=ref, query=query, features=feats, dims=1:nd)
  pred  <- TransferData(anchorset=ancs, refdata=ref[[label_col]][,1])
  query$pred_from_ref <- pred$predicted.id

  flows <- dplyr::count(query@meta.data, pred_from_ref, !!rlang::sym(label_col))
  names(flows) <- c("ref","query","n")
  df <- dplyr::transmute(flows, axis1=paste(ref_name, ref, sep=":"), axis2=paste(query_name, query, sep=":"), n=n, ref=ref)

  p <- ggplot(df, aes(y=n, axis1=axis1, axis2=axis2)) +
    ggalluvial::geom_alluvium(aes(fill=ref), width=1/12, alpha=0.9) +
    ggalluvial::geom_stratum(width=1/12) +
    ggalluvial::geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    theme_cowplot() + labs(title=paste0(ref_name," → ",query_name," (",label_col,")"), y="Cells", x=NULL)
  ggsave(paste0("alluvial_", tolower(ref_name), "_to_", tolower(query_name), ".png"), p, width=11, height=6, dpi=400)
  readr::write_csv(flows, paste0("alluvial_table_", tolower(ref_name), "_to_", tolower(query_name), ".csv"))
}
# Example:
# pair_alluvial(mouse_m, rat_m,      "Mouse","Rat")
# pair_alluvial(mouse_m, marm_m,     "Mouse","Marmoset")
# pair_alluvial(marm_m,  rat_m,      "Marmoset","Rat")

# -------------------------
# 8) Correlation heatmaps (AverageExpression by fine labels; Spearman)
# -------------------------
avg_by_fine <- function(obj){
  Idents(obj) <- obj$label_fine
  as.matrix(AverageExpression(obj, assays="RNA", layer="data")$RNA)
}
A_rat_f   <- avg_by_fine(rat_m)
A_mouse_f <- avg_by_fine(mouse_m)
A_marm_f  <- avg_by_fine(marm_m)

C_mouse_rat  <- cor(A_mouse_f, A_rat_f,  method="spearman")
C_mouse_marm <- cor(A_mouse_f, A_marm_f, method="spearman")
C_marm_rat   <- cor(A_marm_f,  A_rat_f,  method="spearman")
pheatmap::pheatmap(C_mouse_rat,  filename="corr_mouse_vs_rat.png",        width=8, height=8)
pheatmap::pheatmap(C_mouse_marm, filename="corr_mouse_vs_marmoset.png",   width=8, height=8)
pheatmap::pheatmap(C_marm_rat,   filename="corr_marmoset_vs_rat.png",     width=8, height=8)

colnames(A_rat_f)   <- paste0("Rat:",      colnames(A_rat_f))
colnames(A_mouse_f) <- paste0("Mouse:",    colnames(A_mouse_f))
colnames(A_marm_f)  <- paste0("Marmoset:", colnames(A_marm_f))
A_all_lbls <- cbind(A_mouse_f, A_rat_f, A_marm_f)
pheatmap::pheatmap(cor(A_all_lbls, method="spearman"),
                   filename="corr_all_species_all_clusters.png", width=12, height=10)

# -------------------------
# 9) Optional: Integration diagnostics + guarded RPCA integration
# -------------------------
diag_pair_anchors <- function(A, B){
  prep <- function(o){
    DefaultAssay(o) <- "RNA"
    if (!"pca" %in% names(o@reductions)){
      o <- NormalizeData(o, verbose=FALSE) %>% FindVariableFeatures(nfeatures=3000, verbose=FALSE) %>%
        ScaleData(features=VariableFeatures(o), verbose=FALSE) %>% RunPCA(features=VariableFeatures(o), npcs=30, verbose=FALSE)
    }
    o
  }
  A2 <- prep(A); B2 <- prep(B)
  feats <- SelectIntegrationFeatures(object.list=list(A2,B2), nfeatures=3000)
  anc  <- tryCatch(FindIntegrationAnchors(object.list=list(A2,B2),
                                          anchor.features=feats, reduction="rpca",
                                          dims=1:min(30, ncol(Embeddings(A2[["pca"]])), ncol(Embeddings(B2[["pca"]]))),
                                          k.filter=NA), error=function(e) NULL)
  n <- if (is.null(anc)) 0 else nrow(anc@anchors)
  data.frame(pair = paste(unique(A$species)[1], unique(B$species)[1], sep="-"), anchors=n)
}
anc_tbl <- rbind(diag_pair_anchors(mouse_m, rat_m),
                 diag_pair_anchors(mouse_m, marm_m),
                 diag_pair_anchors(marm_m,  rat_m))
print(anc_tbl)
readr::write_csv(anc_tbl, "diagnostics_anchor_counts.csv")

if (sum(anc_tbl$anchors) >= 20){
  message("[Integration] Sufficient anchors detected; running IntegrateData (RPCA).")
  obj.list <- list(rat_m, mouse_m, marm_m)
  obj.list <- lapply(obj.list, function(o){
    DefaultAssay(o) <- "RNA"
    o <- NormalizeData(o, verbose=FALSE) %>% FindVariableFeatures(selection.method="vst", nfeatures=3000, verbose=FALSE) %>%
      ScaleData(features=VariableFeatures(o), verbose=FALSE) %>% RunPCA(features=VariableFeatures(o), npcs=30, verbose=FALSE)
  })
  features <- SelectIntegrationFeatures(object.list=obj.list, nfeatures=3000)
  nd <- min(30, sapply(obj.list, function(o) ncol(Embeddings(o[["pca"]]))))
  anchors  <- FindIntegrationAnchors(object.list=obj.list, anchor.features=features,
                                     reduction="rpca", dims=1:nd, k.filter=NA)
  comb <- IntegrateData(anchorset=anchors, dims=1:nd)
  DefaultAssay(comb) <- "integrated"
  comb <- ScaleData(comb, verbose=FALSE) %>% RunPCA(npcs=30, verbose=FALSE) %>% RunUMAP(dims=1:min(30, ncol(Embeddings(comb[["pca"]]))), verbose=FALSE)
  saveRDS(comb, "integrated_combined.rds")
} else {
  message("[Integration] Too few anchors — skipping integration.")
}

message("Done. Figures and tables written to the working directory.")
