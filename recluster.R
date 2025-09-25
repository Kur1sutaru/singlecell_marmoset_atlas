library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(patchwork)

set.seed(1234)

# Assume 'fibro' is a Seurat obj with counts in assay = "RNA"
DefaultAssay(fibro) <- "RNA"

# Mitochondrial fraction
fibro[["percent.mt"]] <- PercentageFeatureSet(fibro, pattern = "^mt-")  # fibro
# If human: pattern = "^MT-"

VlnPlot(fibro, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)



## ----- 1) QC metrics -----
DefaultAssay(fibro) <- "RNA"
fibro[["percent.mt"]] <- PercentageFeatureSet(fibro, pattern = "^mt-")  # fibro; use "^MT-" for human

VlnPlot(fibro, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

## ----- 2) QC filters (tuned to your plots) -----
min_feats  <- 200
max_feats  <- 5000      # >5k likely doublets/multiplets
max_counts <- 30000     # trim extreme UMI outliers; relax if you expect very large cells
max_mt     <- 15        # most cells are below this

fibro <- subset(
  fibro,
  subset = nFeature_RNA >= min_feats &
    nFeature_RNA <= max_feats &
    nCount_RNA   <= max_counts &
    percent.mt   <= max_mt
)

# Quick check after filtering
VlnPlot(fibro, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

## ----- 3) SCTransform (do NOT remove MT genes; just model percent.mt) -----
fibro <- SCTransform(
  fibro,
  assay = "RNA",
  vst.flavor = "v2",
  method = "glmGamPoi",
  vars.to.regress = "percent.mt",
  verbose = TRUE
)

## ----- 4) DR + clustering -----
fibro <- RunPCA(fibro, assay = "SCT")
dims_use <- 1:30
fibro <- RunUMAP(fibro, dims = dims_use)
fibro <- FindNeighbors(fibro, dims = dims_use)
fibro <- FindClusters(fibro, resolution = 0.3)

DimPlot(fibro, label = TRUE, repel = TRUE) + NoLegend()

## ----- 5) Markers (use SCT assay) -----
DefaultAssay(fibro) <- "SCT"
markers <- FindAllMarkers(
  fibro, only.pos = TRUE, test.use = "wilcox",
  logfc.threshold = 0.25, min.pct = 0.1
)

# make sure the gene column exists (older Seurat sometimes puts genes in rownames)
if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)

# save everything
write.csv(markers, "markers_all_clusters_wilcox_fibrores03.csv", row.names = FALSE)

## ----- 6) Feature plots on SCT
DefaultAssay(fibro) <- "SCT"
FeaturePlot(fibro, features = c("Acta2","Des","Col1a1"), ncol = 3, order = TRUE)

# Dot plot for a marker panel
genes <- c("Pecam1","Kdr","Col1a1","Pdgfra","Aif1","Cx3cr1")
DotPlot(fibro, features = genes, assay = "SCT", cols = c("lightgrey","red")) +
  RotatedAxis() + ggtitle("Marker panel")

# Ridge plot (expression distributions)
RidgePlot(fibro, features = c("Pecam1","Col1a1"), ncol = 2)

## ------ 6.5) markers
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

## ===== 0) Assumptions =====
# - Your fibroect is 'fibro' (assay RNA/SCT already set up)
# - Identities are meaningful (e.g., seurat_clusters or your annotations)
DefaultAssay(fibro) <- if ("SCT" %in% names(fibro@assays)) "SCT" else "RNA"

## ===== 1) Marker panels (from your list) =====
panels_raw <- list(
  Fibroblast_general   = c("Col1a1","Col1a2","Pdgfra","Dcn","Lum"),
  Fibroblast_perineur  = c("Slc2a1","Itgb4","Cldn1"),
  Fibroblast_endoneur  = c("Pdgfra","Dcn","Col6a1"),
  Fibroblast_myofibro  = c("Acta2","Tagln","Myh11"),
  Fibroblast_perivasc  = c("Col15a1","Fn1","Vtn"),
  
  Keratinocyte_basal   = c("Krt14","Krt5","Itga6","Trp63"),
  Keratinocyte_spinous = c("Krt1","Krt10","Krt16"),
  Keratinocyte_granular= c("Flg","Lor","Ivl"),
  Keratinocyte_cornif  = c("Sprr1a","Lce3d"),
  
  Endo_pan             = c("Pecam1","Cdh5","Vwf","Tie1"),
  Endo_arterial        = c("Efnb2","Gja5","Hey1"),
  Endo_venous          = c("Nr2f2","Ephb4","Flt4"),
  Endo_capillary       = c("Mfsd2a","Car4","Slc2a1"),
  Endo_lymphatic       = c("Prox1","Lyve1","Pdpn","Flt4"),
  
  Imm_Macrophage       = c("Adgre1","Cd68","C1qa","Mrc1"),
  Imm_Dendritic        = c("Itgax","H2-Ab1","Cd209a"),
  Imm_Monocyte         = c("Ly6c2","Ccr2"),
  Imm_Tcell            = c("Cd3d","Cd3e","Cd4","Cd8a"),
  Imm_Bcell            = c("Cd79a","Cd19","Ms4a1"),
  Imm_NK               = c("Nkg7","Klrb1c","Gzmb"),
  Imm_Mast             = c("Kit","Cma1","Mcpt4"),
  
  Schwann_mye          = c("Mpz","Mbp","Plp1","Prx"),
  Schwann_nonmye       = c("Ngfr","L1cam","Gap43"),
  Schwann_repair       = c("Shh","Gdnf","Sox2"),
  Satellite_glia       = c("Fabp7","Glul","Apoe","Gfap"),
  
  Pericyte_cap         = c("Pdgfrb","Cspg4","Rgs5","Kcnj8"),
  Pericyte_venular     = c("Des","Acta2","Abcc9"),
  
  SMC_contractile      = c("Acta2","Myh11","Tagln","Cnn1"),
  SMC_synthetic        = c("Fn1","Mgp","Ctgf","Thbs1")
)

## ===== 2) Map to fibroect features (case-insensitive, drop missing) =====
gene_map <- setNames(rownames(fibro), toupper(rownames(fibro)))
map_genes <- function(v) unique(na.omit(gene_map[toupper(v)]))
panels <- lapply(panels_raw, map_genes)
panels <- panels[sapply(panels, length) > 0]

## ===== 3) Plot helpers =====
outdir <- "plots_markers"; dir.create(outdir, showWarnings = FALSE)
group.by <- if ("label_coarse" %in% colnames(fibro@meta.data)) "label_coarse" else "seurat_clusters"

plot_panel <- function(fibro, panel_name, feats, pt.size=0.2) {
  message(">> Panel: ", panel_name, " (", length(feats), " genes)")
  feats <- intersect(feats, rownames(fibro))
  if (length(feats) == 0) return(invisible(NULL))
  
  # FeaturePlot grid (on UMAP if available)
  fp <- tryCatch({
    FeaturePlot(fibro, features = feats, order = TRUE, label = FALSE)
  }, error=function(e) NULL)
  
  # Violin (expression by identity)
  vp <- VlnPlot(fibro, features = feats, group.by = group.by, pt.size = 0.1, combine = TRUE, stack = FALSE)
  
  # Ridge (distributions)
  rp <- RidgePlot(fibro, features = feats, group.by = group.by, ncol = 2)
  
  # DotPlot (scaled avg + pct)
  dp <- DotPlot(fibro, features = feats, group.by = group.by, cols = c("lightgrey","red")) +
    RotatedAxis() + ggtitle(panel_name)
  
  # Save
  if (!is.null(fp)) {
    ggsave(file.path(outdir, paste0(panel_name,"__FeaturePlot.png")), fp, width=10, height=8, dpi=300)
    ggsave(file.path(outdir, paste0(panel_name,"__FeaturePlot.pdf")), fp, width=10, height=8)
  }
  ggsave(file.path(outdir, paste0(panel_name,"__Violin.png")), vp, width=12, height=8, dpi=300)
  ggsave(file.path(outdir, paste0(panel_name,"__Violin.pdf")), vp, width=12, height=8)
  
  ggsave(file.path(outdir, paste0(panel_name,"__Ridge.png")),  rp, width=12, height=8, dpi=300)
  ggsave(file.path(outdir, paste0(panel_name,"__Ridge.pdf")),  rp, width=12, height=8)
  
  ggsave(file.path(outdir, paste0(panel_name,"__Dot.png")),    dp, width=10, height=6, dpi=300)
  ggsave(file.path(outdir, paste0(panel_name,"__Dot.pdf")),    dp, width=10, height=6)
  
  invisible(list(FeaturePlot=fp, Violin=vp, Ridge=rp, Dot=dp))
}

## ===== 4) Run for all panels =====
plots <- lapply(names(panels), function(nm) plot_panel(fibro, nm, panels[[nm]]))

## ===== 5) (Optional) One big DotPlot across all genes with panel facets =====
all_feats <- unique(unlist(panels))
dp_all <- {
  dp <- DotPlot(fibro, features = all_feats, group.by = group.by, cols = c("lightgrey","red")) + theme_bw(11)
  df <- dp$data
  # Map feature -> panel
  f2p <- rep(names(panels), lengths(panels)); names(f2p) <- unlist(panels)
  df$Panel <- f2p[as.character(df$features.plot)]
  df$features.plot <- factor(df$features.plot, levels = rev(all_feats))
  library(ggplot2)
  ggplot(df, aes(id, features.plot, size=pct.exp, color=avg.exp.scaled)) +
    geom_point() +
    scale_size(range=c(0,6)) +
    scale_color_gradient(low="lightgrey", high="red") +
    facet_grid(Panel ~ ., scales="free_y", space="free_y") +
    labs(x=group.by, y=NULL, title="Marker panels (DotPlot)") +
    theme_bw(11)
}
ggsave(file.path(outdir, "ALL_PANELS__DotPlot.png"), dp_all, width=10, height=14, dpi=300)
ggsave(file.path(outdir, "ALL_PANELS__DotPlot.pdf"), dp_all, width=10, height=14)

saveRDS(fibro, "fibro_clustering_sct_res03.rds")

DefaultAssay(fibro) <- if ("SCT" %in% names(fibro@assays)) "SCT" else "RNA"
grp <- if ("label_coarse" %in% colnames(fibro@meta.data)) "label_coarse" else "seurat_clusters"

panels <- list(
  Macrophage = c("Adgre1","Cd68","Mrc1","C1qa","Apoe","H2-Ab1"),
  Monocyte   = c("Ly6c2","Ccr2","Plac8","Lyz2"),
  cDC1       = c("Xcr1","Clec9a","Batf3","H2-Ab1"),
  cDC2       = c("Itgam","Sirpa","Klf4","H2-Ab1"),
  Langerhans = c("Cd207","Itgax","H2-Ab1","Epcam"),
  T_CD4      = c("Cd3d","Cd3e","Cd4","Il7r","Foxp3"),
  T_CD8      = c("Cd3d","Cd3e","Cd8a","Nkg7","Gzmb","Prf1"),
  B_Plasma   = c("Ms4a1","Cd19","Cd79a","Jchain","Sdc1"),
  Mast       = c("Kit","Tpsb2","Cma1","Cpa3","Mcpt4"),
  NK_ILC     = c("Ncr1","Klrc1","Klrd1","Prf1","Gzmb","Gata3","Il1rl1"),
  Neutro     = c("Ly6g","S100a8","S100a9","Lcn2"),
  Eosino     = c("Epx","Prg2","Siglecf")
)

genes <- unique(unlist(panels))
DotPlot(fibro, features = genes, group.by = grp, cols = c("lightgrey","red")) + RotatedAxis()

## check some redundant clusters
# wrapper that uses your current Seurat object `fibro`
map_genes <- function(genes) map_genes_verbose(fibro, genes)

# now your original calls work:
mfib  <- map_genes(c("Acta2","Tagln","Pdgfra","Col1a1","Col1a2","Dcn","Fn1","Postn","Itga11"))
smc   <- map_genes(c("Myh11","Cnn1","Tagln","Acta2","Myocd","Actg2","Smtn","Tpm2","Myl9"))
peric <- map_genes(c("Pdgfrb","Rgs5","Cspg4","Mcam","Des","Kcnj8","Abcc9","Notch3"))


fibro <- AddModuleScore(fibro, features=list(mfib),  name="MFib")
fibro <- AddModuleScore(fibro, features=list(smc),   name="SMC")
fibro <- AddModuleScore(fibro, features=list(peric), name="Peri")

# inspect per cluster
scores <- aggregate(cbind(MFib1,SMC1,Peri1) ~ fibro$seurat_clusters, fibro@meta.data, mean)
scores[order(-scores$MFib1), ]

# rename cluster 11
fibro$celltype <- as.character(fibro$celltype %||% fibro$seurat_clusters)
fibro$celltype[fibro$seurat_clusters == "11"] <- "Myofibroblast"
fibro$celltype <- factor(fibro$celltype)



mural_panels <- list(Myofibroblast = mfib, SMC = smc, Pericyte = peric)
ms <- add_scores(obj, mural_panels, prefix="MS_mural_"); obj <- ms$obj
mural_cols <- ms$cols   # e.g., "Myofibroblast_Score", "SMC_Score", "Pericyte_Score"
p_umap_scores <- FeaturePlot(obj, features = mural_cols, reduction = "umap",
                             cols = c("grey90","blue"), order = TRUE)
ggsave("umap_mural_programs.png", p_umap_scores, width=10, height=4, dpi=300)
p_vln <- VlnPlot(obj, features = mural_cols, group.by = "seurat_clusters",
                 pt.size = 0, stack = TRUE, flip = TRUE, same.y.lims = TRUE)
ggsave("vln_mural_programs_by_cluster.png", p_vln, width=10, height=6, dpi=300)
calls <- means |>
  group_by(cluster) |>
  slice_max(order_by = score, n = 1, with_ties = FALSE) |>
  transmute(cluster = as.character(cluster), best_signature = signature, best_score = score)

print(calls)           # shows which program wins per cluster
write.csv(calls, "mural_program_best_by_cluster.csv", row.names = FALSE)

#----- 7) Quick label transfer (manual mapping example)
# make sure Idents(fibro) are the numeric cluster IDs "0","1",...
Idents(fibro) <- fibro$seurat_clusters

# rename identities using your map
label_map <- c("0"=" PI16/DPP4⁺ Sca1⁺ adventitial fibroblasts","1"="ECM-remodeling fibroblasts",
               "2"="Inflammatory / stress-activated fibroblasts","3"="Inflammatory / stress-activated fibroblasts",
               "4"="Doublet","5"="Cxcl12⁺ perivascular fibroblasts",
               "6"="Perivascular fibroblasts","7"="Tenocyte fibroblasts",
               "8"="PI16/DPP4⁺ Sca1⁺ adventitial fibroblasts","9"="Perivascular fibroblasts",
               "10"="WNT-responsive fibroblasts","11"="PI16/DPP4⁺ Sca1⁺ adventitial fibroblasts",
               "12"="Perivascular fibroblasts","13"="Tenocyte fibroblasts")

fibro <- RenameIdents(fibro, label_map)

# if you also want a metadata column:
fibro$labels <- Idents(fibro)
# remove all cells whose current identity is "Doublet"
fibro_noDoublet <- subset(fibro, subset = labels != "Doublets")

saveRDS(fibro_noDoublet, file = "fibro_sct_annotated.rds")

###########################################################################
### Shivani annot
## --- Setup ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
DefaultAssay(fibro) <- "RNA"

## --- 1) Define marker groups (mouse TitleCase) ----
gene_groups <- list(
  `FIBROBLASTS CORE MARKERS` = c("Col1a1","Col3a1","Dcn","Lum","Pdgfra"),
  
  `Cluster 0 • Classic Perineural` = c("Cldn3","Ildr2","Opcml","Nell2","Lgi2",
                                       # helper to rule in/out PI16+ adventitial
                                       "Pi16","Dpp4","Cd34"),
  `Cluster 1 • Perivascular`       = c("Hmcn2","Col8a1","Fbln7","Meox1","Bmper","Eln"),
  `Cluster 2 • Inflammatory`       = c("Ccl2","Ccl7","Cxcl1","Cxcl2","Ptgs2","Il6"),
  `Cluster 3 • Stress/IER`         = c("Fos","Jun","Atf3","Egr1","Zfp36","Socs3","Has1","Serpine1"),
  
  `Cluster 4 • Contaminant check`  = c("Col1a1","Dcn","Pdgfra","Col3a1",  # fibro check
                                       "Pecam1","Prtn3","Muc13","Nrcam","Il31ra","Trpv2"), # non-fibro
  `Cluster 5 • Neurotrophic Perivasc` = c("Cygb","Abcc9","Kcnj8","Apoe","Ntf3","Penk",
                                          "Fgf13","Gabra3","Dner","Foxp2","Trpm3"),
  `Cluster 6 • Neuron-support Perineural (male)` = c("Lgi1","Adam23","Lrrtm3","Ctnna3","Vsnl1",
                                                     "Kdm5d","Uty","Eif2s3y"),
  `Cluster 7 • Wnt-Responsive Perivasc` = c("Myl10","Wif1","Wnt10a","Dkk3","Nkd1","Ndnf"),
  `Cluster 8 • Stress-responsive Perivasc` = c("Cygb","Kcnj8","Ccn2","Mgp","Gfra2","Hspa1a"),
  `Cluster 9 • Progenitor Perivasc` = c("Col4a6","Fbn2","Parvb","Lgr5","Grem1","Synpr"),
  `Cluster 10 • Tenocyte-like` = c("Scx","Tnmd","Fmod","Col11a1","Comp","Prg4"),
  `Cluster 11 • Endothelial check` = c("Pecam1","Kdr","Cdh5","Vwf","Emcn",
                                       "Col1a1","Dcn","Pdgfra","Col3a1"), # fibro check
  `Cluster 12 • Mesenchymal Perineural` = c("Grid2","Cntn1","Cadps","Ndnf","Slitrk6","Pou3f3","Hand1"),
  `Cluster 13 • Mechanosensory Perineural` = c("Piezo2","Ngfr","Gfra1","Grin2b","Reln")
)

## Normalize case and keep only genes present
fix_case <- function(x) stringr::str_to_title(x)
gene_groups <- lapply(gene_groups, fix_case)
present    <- rownames(fibro)
gene_groups <- lapply(gene_groups, function(g) intersect(g, present))

## Optional: set plot identities to your subtype labels if available
if ("labels" %in% colnames(fibro@meta.data)) {
  Idents(fibro) <- fibro$labels
}

## --- 2) Build DotPlot data and keep ≥35% expressing (no avg.exp filter) ----
dot_genes <- unique(unlist(gene_groups))
dp <- Seurat::DotPlot(fibro, features = dot_genes, scale = FALSE)
df <- dp$data

# harmonize column names from Seurat versions
if ("features.plot" %in% names(df)) df <- dplyr::rename(df, feature = features.plot)
if ("features"      %in% names(df)) df <- dplyr::rename(df, feature = features)
if ("ident.1"       %in% names(df)) df <- dplyr::rename(df, ident   = ident.1)
if ("id"            %in% names(df)) df <- dplyr::rename(df, ident   = id)

# make pct.exp 0–100 if needed
if (max(df$pct.exp, na.rm = TRUE) <= 1) df$pct.exp <- 100 * df$pct.exp

df_filt <- df %>%
  dplyr::filter(pct.exp >= 35) %>%            # keep only ≥35% expressing
  dplyr::mutate(
    feature = factor(feature, levels = dot_genes),
    ident   = factor(ident)
  )

## --- 3) DotPlot (cleaner) ----
p_dot <- ggplot(df_filt, aes(x = feature, y = ident,
                             size = pct.exp, color = avg.exp)) +
  geom_point() +
  scale_size_continuous(limits = c(35, 100), range = c(1.8, 7),
                        breaks = c(35, 50, 75, 100),
                        name = "Percent Expressed") +
  scale_color_gradient(low = "grey90", high = "blue4", name = "Avg expression") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(title = "Fibroblast subtype markers (≥35% expressing)",
       x = "Genes", y = "Group / cluster")

print(p_dot)
ggsave("dotplot_fibro_subtypes_ge35pct.png", p_dot, width = 14, height = 7, dpi = 300)

## --- 4) Quick FeaturePlot panels (pick a few per group) ----
feat_highlights <- c(
  "Col1a1","Col3a1","Dcn","Lum","Pdgfra",           # core
  "Cldn3","Ildr2","Opcml","Nell2","Pi16","Dpp4",    # Cl0 + check
  "Hmcn2","Col8a1","Fbln7","Meox1","Bmper",         # Cl1
  "Ccl2","Ccl7","Cxcl1","Ptgs2","Il6",              # Cl2
  "Fos","Jun","Atf3","Egr1","Has1",                 # Cl3
  "Cygb","Abcc9","Kcnj8","Apoe","Ntf3","Penk",      # Cl5
  "Lgi1","Adam23","Lrrtm3","Ctnna3","Vsnl1",        # Cl6
  "Wif1","Wnt10a","Dkk3","Nkd1","Ndnf",             # Cl7
  "Ccn2","Mgp","Gfra2","Hspa1a",                    # Cl8
  "Col4a6","Fbn2","Parvb","Lgr5","Grem1","Synpr",   # Cl9
  "Scx","Tnmd","Fmod","Col11a1","Comp","Prg4",      # Cl10
  "Grid2","Cntn1","Cadps","Ndnf","Slitrk6","Pou3f3","Hand1",  # Cl12
  "Piezo2","Ngfr","Gfra1","Grin2b","Reln"           # Cl13
)
feat_highlights <- intersect(stringr::str_to_title(feat_highlights), present)

if (length(feat_highlights) > 0) {
  p_feat <- FeaturePlot(fibro, features = feat_highlights, cols = c("grey90","red3"),
                        reduction = "umap", ncol = 5, pt.size = 0.3, order = TRUE)
  print(p_feat)
  ggsave("featureplot_fibro_subtypes_umap.png", p_feat, width = 16, height = 12, dpi = 300)
}




# Rename cluster identities
new_cluster_names <- c(
  "Cluster 0 Classic Perineural",
  "Cluster 1 Perivascular",
  "Cluster 2 Inflammatory",
  "Cluster 3 Stress/IER",
  "Cluster 4 Endothelial",
  "Cluster 5 Neurotrophic Perivasc",
  "Cluster 6 Neuron-support Perineural",
  "Cluster 7 Wnt-Responsive Perivasc",
  "Cluster 8 Stress-responsive Perivasc",
  "Cluster 9 Progenitor Perivasc",
  "Cluster 10 Tenocyte-like",
  "Cluster 11 Endothelial",
  "Cluster 12 Mesenchymal Perineural",
  "Cluster 13 Mechanosensory Perineural"
)

# Apply new names to Seurat object
names(new_cluster_names) <- levels(fibro) 
fibro <- RenameIdents(fibro, new_cluster_names)

# Generate UMAP plot
DimPlot(fibro, reduction = "umap", label = TRUE, label.size = 4) + 
  ggtitle("UMAP of Renamed Fibroblast Clusters") +
  theme_minimal()


# Core fibroblast markers
core_fibroblast <- c("Col1a1", "Col3a1", "Dcn", "Lum", "Pdgfra")

# Subtype-specific markers
classic_perineural <- c("Cldn3", "Ildr2", "Opcml", "Nell2", "Lgi2")
classic_perineural <- c("Hmcn2", "Col8a1", "Fbln7", "Meox1", "Bmper", "Eln")
inflammatory <- c("Ccl2",  "Cxcl1", "Cxcl2", "Ptgs2", "Il6")
stress_ier <- c("Fos", "Jun", "Atf3", "Egr1", "Zfp36", "Socs3", "Has1", "Serpine1")
Endothelial <- c("Pecam1", "Vwf", "Prom1", "Nrcam", "Cd34", "Fabp4")
neurotrophic_perivasc <- c("Cygb", "Abcc9", "Kcnj8", "Apoe", "Ntf3", "Penk", "Fgf13", "Gabra3", "Dner", "Foxp2", "Trpm3")
wnt_responsive <- c("Myl10", "Wif1", "Wnt10a", "Dkk3", "Nkd1", "Ndnf")
stress_perivasc <- c("Cygb", "Kcnj8", "Ccn2", "Mgp", "Gfra2")
progenitor_perivasc <- c("Col4a6", "Fbn2", "Parvb", "Lgr5", "Grem1", "Synpr")
tenocyte_like <- c("Scx", "Tnmd", "Fmod",  "Comp", "Prg4")
mesenchymal_perineural <- c("Grid2", "Cntn1", "Cadps", "Ndnf", "Slitrk6", "Pou3f3", "Hand1")
mechanosensory_perineural <- c("Piezo2", "Ngfr", "Grin2b", "Reln")




fibro <- NormalizeData(fibro)
fibro <- ScaleData(fibro, features = unique(c(
  core_fibroblast, classic_perineural, perivascular, inflammatory, stress_ier,
  Endothelial, neurotrophic_perivasc, neuron_support_male, wnt_responsive,
  stress_perivasc, progenitor_perivasc, tenocyte_like, mesenchymal_perineural,
  mechanosensory_perineural
)))


  
# Combine all gene sets into one vector
all_genes <- unique(c(
  core_fibroblast,
  classic_perineural,
  perivascular,
  inflammatory,
  stress_ier,
  Endothelial,
  neurotrophic_perivasc,
  wnt_responsive,
  stress_perivasc,
  progenitor_perivasc,
  tenocyte_like,
  mesenchymal_perineural,
  mechanosensory_perineural
))



# Create a named list with filtered unique genes per group
gene_groups <- list(
  "Core Fibroblast" = intersect(core_fibroblast, all_genes),
  "Classic Perineural" = intersect(classic_perineural, all_genes),
  "Perivascular" = intersect(perivascular, all_genes),
  "Inflammatory" = intersect(inflammatory, all_genes),
  "Stress/IER" = intersect(stress_ier, all_genes),
  "Endothelial" = intersect(Endothelial, all_genes),
  "Neurotrophic Perivasc" = intersect(neurotrophic_perivasc, all_genes),
  "Wnt-Responsive Perivasc" = intersect(wnt_responsive, all_genes),
  "Stress-Responsive Perivasc" = intersect(stress_perivasc, all_genes),
  "Progenitor Perivasc" = intersect(progenitor_perivasc, all_genes),
  "Tenocyte-Like" = intersect(tenocyte_like, all_genes),
  "Mesenchymal Perineural" = intersect(mesenchymal_perineural, all_genes),
  "Mechanosensory Perineural" = intersect(mechanosensory_perineural, all_genes)
)


DotPlot(fibro, features = gene_groups) +
  RotatedAxis() +
  theme_minimal()

# Flatten all genes from the grouped list
flat_genes <- unique(unlist(gene_groups))
DotPlot(fibro, features = flat_genes) +
  RotatedAxis() +
  theme_minimal()

DotPlot(fibro, features = flat_genes) +
  RotatedAxis() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("DotPlot of Fibroblast Marker Genes")


# Example dot plot object
p <- DotPlot(fibro, features = all_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(group ~ ., scales = "free_y", space = "free") +
  theme(strip.text.y = element_text(angle = 0))

# Create the dot plot first
p <- DotPlot(fibro, features = all_genes)

# Add group info manually
gene_df <- data.frame(
  gene = unlist(gene_groups),
  group = rep(names(gene_groups), times = lengths(gene_groups))
)

# Match group labels to the plot data
p$data$group <- gene_df$group[match(p$data$features.plot, gene_df$gene)]

# Now facet safely
p <- p +
  facet_grid(group ~ ., scales = "free_y", space = "free") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("dotplot_fibroblast_subtypes.tiff", plot = p,
       width = 10, height = 12, dpi = 600, units = "in", compression = "lzw")




# Use your flat gene list
vln<- VlnPlot(fibro, features = flat_genes, pt.size = 0.1, group.by = "seurat_clusters") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold")
  )
ggsave("violin_plot.png", plot = vln, width = 12, height = 8, dpi = 300)




# Loop through each gene and save one violin plot per page
for (gene in flat_genes) {
  vln <- VlnPlot(fibro,
                 features = gene,
                 pt.size = 0.1,
                 group.by = "cluster_labels") +  # Uses renamed cluster labels
    theme_minimal() +
    ggtitle(paste("Violin Plot:", gene)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  ggsave(filename = paste0("violin_", gene, ".png"),
         plot = vln,
         width = 7, height = 5, dpi = 300)
}


# Rename cluster identities
new_cluster_names <- c(
  "Pi16⁺ Adventitial Fibroblasts",
  "ECM-rich Myofibroblast-like Fibroblasts (Matrifibrocytes)",
  "Activated Inflammatory Fibroblasts",
  "Stress-Responsive / Wound-Activated Fibroblasts",
  "Perivascular Fibroblasts (Endothelial-like)",
  "Lipofibroblast-like AY036118⁺/Plpp3⁺",
  "Mesenchymal Signaling Fibroblasts (Cilp⁺/Thbs4⁺/Rspo1⁺)",
  "Chondrogenic/Tendon-like Fibroblasts (Fmod⁺/Prg4⁺/Tnmd⁺)",
  "Perivascular/Angiogenic Fibroblasts (Esm1⁺/Mgp⁺/Ccn2⁺)",
  "Endothelial Cells (Cdh5⁺/Kdr⁺/Flt1⁺)",
  "Perineurial Fibroblasts (Tenm2⁺/Reln⁺/Sox9⁺)"
)

# Apply new names to Seurat object
names(new_cluster_names) <- levels(fibro)
fibro <- RenameIdents(fibro, new_cluster_names)

# Generate UMAP with new labels
DimPlot(fibro, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("Fibroblast and Endothelial Subtypes in Mouse Tongue") +
  theme_minimal()

# Rename cluster identities
new_cluster_names <- c(
  "Pi16⁺ Adventitial Fibroblasts",
  "ECM-rich Myofibroblast-like Fibroblasts (Matrifibrocytes)",
  "Activated Inflammatory Fibroblasts",
  "Stress-Responsive / Wound-Activated Fibroblasts",
  "Perivascular Fibroblasts (Endothelial-like)",
  "Lipofibroblast-like AY036118⁺/Plpp3⁺",
  "Mesenchymal Signaling Fibroblasts (Cilp⁺/Thbs4⁺/Rspo1⁺)",
  "Chondrogenic/Tendon-like Fibroblasts (Fmod⁺/Prg4⁺/Tnmd⁺)",
  "Perivascular/Angiogenic Fibroblasts (Esm1⁺/Mgp⁺/Ccn2⁺)",
  "Endothelial Cells (Cdh5⁺/Kdr⁺/Flt1⁺)",
  "Perineurial Fibroblasts (Tenm2⁺/Reln⁺/Sox9⁺)"
)

# Apply new names to Seurat object
names(new_cluster_names) <- levels(fibroect)
fibroect <- RenameIdents(fibroect, new_cluster_names)

png("UMAP_Fibroblast_Subtypes.png", width = 1200, height = 900, res = 150)
DimPlot(fibro, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("Fibroblast subtypes") +
  theme_minimal()
dev.off()

saveRDS(fibro, file = "fibro_correct_annotatedres03.rds")


# Define your gene sets
marker_genes <- list(
  "Core Fibroblast" = c("Col1a1", "Col3a1", "Dcn", "Lum", "Pdgfra"),
  "Classic Perineural" = c("Cldn3", "Ildr2", "Opcml", "Nell2", "Lgi2"),
  "Perivascular" = c("Hmcn2", "Col8a1", "Fbln7", "Meox1", "Bmper", "Eln"),
  "Inflammatory" = c("Ccl2", "Cxcl1", "Cxcl2", "Ptgs2", "Il6"),
  "Stress IER" = c("Fos", "Jun", "Atf3", "Egr1", "Zfp36", "Socs3", "Has1", "Serpine1"),
  "Endothelial" = c("Pecam1", "Vwf", "Prom1", "Nrcam", "Cd34", "Fabp4"),
  "Neurotrophic Perivasc" = c("Cygb", "Abcc9", "Kcnj8", "Apoe", "Ntf3", "Penk", "Fgf13", "Gabra3", "Dner", "Foxp2", "Trpm3"),
  "Wnt Responsive" = c("Myl10", "Wif1", "Wnt10a", "Dkk3", "Nkd1", "Ndnf"),
  "Stress Perivasc" = c("Cygb", "Kcnj8", "Ccn2", "Mgp", "Gfra2"),
  "Progenitor Perivasc" = c("Col4a6", "Fbn2", "Parvb", "Lgr5", "Grem1", "Synpr"),
  "Tenocyte-like" = c("Scx", "Tnmd", "Fmod", "Comp", "Prg4"),
  "Mesenchymal Perineural" = c("Grid2", "Cntn1", "Cadps", "Ndnf", "Slitrk6", "Pou3f3", "Hand1"),
  "Mechanosensory Perineural" = c("Piezo2", "Ngfr", "Grin2b", "Reln")
)


library(dplyr)

# Create a summary table of counts and proportions
composition_table <- fibro@meta.data %>%
  group_by(annotres03, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(annotres03) %>%
  mutate(freq = n / sum(n))  # Proportion within each cluster

library(ggplot2)

ggplot(composition_table, aes(x = annotres03, y = freq, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster (annotres03)", y = "Composition (%)", fill = "Cell Type",
       title = "Cluster Composition by Cell Type") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cluster_composition_by_celltype.png", width = 10, height = 6, dpi = 600)





library(dplyr)

# Summarize counts and percentages
composition_table <- fibro@meta.data %>%
  group_by(annotres03, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(annotres03) %>%
  mutate(freq = n / sum(n),
         label = paste0(n, " (", round(freq * 100, 1), "%)"))  # Label with count and %

write.csv(composition_table, "cellpropfibroanotres03.csv")
library(ggplot2)

ggplot(composition_table, aes(x = annotres03, y = freq, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster (annotres03)", y = "Composition (%)", fill = "Cell Type",
       title = "Cluster Composition by Cell Type with Counts and Percentages") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cluster_composition_labeled.png", width = 20, height = 10, dpi = 600)


library(dplyr)
library(ggplot2)

# Step 1: Summarize cell counts
composition_counts <- fibro@meta.data %>%
  group_by(annotres03, celltype) %>%
  summarise(count = n(), .groups = "drop")

# Step 2: Plot with count labels
ggplot(composition_counts, aes(x = annotres03, y = count, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  labs(x = "Cluster (annotres03)", y = "Cell Count", fill = "Cell Type",
       title = "Cluster Composition by Cell Type (Cell Counts Only)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# Define your gene sets
marker_genes <- list(
  "Core Fibroblast" = c("Col1a1", "Col3a1", "Dcn", "Lum", "Pdgfra"),
  "Classic Perineural" = c("Cldn3", "Ildr2", "Opcml", "Nell2", "Lgi2"),
  "Perivascular" = c("Hmcn2", "Col8a1", "Fbln7", "Meox1", "Bmper", "Eln"),
  "Inflammatory" = c("Ccl2", "Cxcl1", "Cxcl2", "Ptgs2", "Il6"),
  "Stress IER" = c("Fos", "Jun", "Atf3", "Egr1", "Zfp36", "Socs3", "Has1", "Serpine1"),
  "Endothelial" = c("Pecam1", "Vwf", "Prom1", "Nrcam", "Cd34", "Fabp4"),
  "Neurotrophic Perivasc" = c("Cygb", "Abcc9", "Kcnj8", "Apoe", "Ntf3", "Penk", "Fgf13", "Gabra3", "Dner", "Foxp2", "Trpm3"),
  "Wnt Responsive" = c("Myl10", "Wif1", "Wnt10a", "Dkk3", "Nkd1", "Ndnf"),
  "Stress Perivasc" = c("Cygb", "Ccn2", "Mgp", "Gfra2"),
  "Progenitor Perivasc" = c("Col4a6", "Fbn2", "Parvb", "Lgr5", "Grem1", "Synpr"),
  "Tenocyte-like" = c("Scx", "Tnmd", "Fmod", "Comp", "Prg4"),
  "Mesenchymal Perineural" = c("Grid2", "Cntn1", "Cadps", "Ndnf", "Slitrk6", "Pou3f3", "Hand1"),
  "Mechanosensory Perineural" = c("Piezo2", "Ngfr", "Grin2b", "Reln")
)

# If you built marker_genes from a list like marker_list
marker_genes <- unique(unlist(marker_genes))  # Flatten and deduplicate

fibro$annotres03 <- as.factor(fibro$annotres03)
dot_plot <- DotPlot(fibro, features = marker_genes, group.by = "annotres03") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("A. Marker Expression by Cluster")

violin_plot <- VlnPlot(fibro, features = marker_genes[1:4], group.by = "annotres03", pt.size = 0) +
  ggtitle("B. Violin Plot of Key Markers")

feature_plot <- FeaturePlot(fibro, features = marker_genes[1:4], reduction = "umap", cols = c("lightgrey", "blue")) &
  theme(legend.position = "none") &
  ggtitle("C. Feature Plot on UMAP")

umap_plot <- DimPlot(fibro, reduction = "umap", group.by = "annotres03", label = TRUE) +
  ggtitle("D. UMAP of Annotated Clusters")


library(patchwork)

# Combine into a 2x2 layout
final_figure <- (dot_plot | violin_plot) / (feature_plot | umap_plot) +
  plot_annotation(title = "Fibroblast and Endothelial Subtype Landscape",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("multi_panel_subtype_figure.png", plot = final_figure,
       width = 16, height = 12, dpi = 600)



library(dplyr)
library(tidyr)

# Step 1: Count cells per subtype per sample
composition_table <- fibro@meta.data %>%
  group_by(celltype, annotres03) %>%
  summarise(cell_count = n(), .groups = "drop")

# Step 2: Pivot to wide format (optional, for easier viewing)
composition_wide <- composition_table %>%
  pivot_wider(names_from = annotres03, values_from = cell_count, values_fill = 0)

# View the table
print(composition_wide)
write.csv(composition_wide, "compositionwidefibro.csv")


composition_percent <- composition_table %>%
  group_by(celltype) %>%
  mutate(percent = round(cell_count / sum(cell_count) * 100, 2))
print(composition_percent)
write.csv(composition_percent, "composition_percent_fibro.csv")


install.packages("gridExtra")
install.packages("grid")
install.packages("formattable")
install.packages("webshot")
library(formattable)

# Style the percent column with a color gradient
styled_table <- formattable(composition_df, list(
  percent = color_tile("white", "steelblue")
))

library(webshot)
library(formattable)
# Save as HTML first
export_path <- "subtype_composition_table.html"
export_image <- "subtype_composition_table.png"

export_formattable <- function(f, file) {
  htmlwidgets::saveWidget(as.htmlwidget(f), file, selfcontained = TRUE)
}
styled_table <- formattable(composition_percent, list(
  percent = color_tile("white", "steelblue")
))
export_formattable(styled_table, export_path)

# Convert HTML to PNG
webshot(export_path, export_image, vwidth = 1200, vheight = 800)
