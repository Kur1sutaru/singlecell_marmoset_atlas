setwd("C:/Users/crist/Desktop/cross_species_paper/subset_all_major_types/Fibroblast/fibro_predictions")

###### prediction of fibroblast subtypes with packages #####
## Named marker list: each element is a character vector of genes
fibro_markers <- list(
  "Pi16/dpp4+ papillary adventitial" = c("pi16","dpp4","ly6a","ly6c1","pdgfra","col15a1"),
  "ECM-myofibroblast"              = c("tagln","acta2","myh11","col1a1","col1a2","fn1"),
  "Inflammatory fibroblasts"         = c("ccl2","ccl7","il6","il1b","cxcl1","nfkbia"),
  "Stress-activated infl. fibro"     = c("fos","jun","atf3","hspa1a","hspa1b","dusp1"),
  "Matrix ecm-remodeling"           = c("mmp2","mmp14","timp2","sparc","col3a1","lox"),
  "Perivascular fibroblasts"        = c("col15a1","vit","vwa1","rgs5","pdgfrb","kcnj8","abcc9"),
  "Car-like adventitial"            = c("cxcl12","lepr","rgs5","kdr","vcam1"),
  "WNT-responsive fibroblasts"      = c("axin2","wif1","rspo1","rspo3","dkk3","fgf10"),
  "Tenocyte fibroblasts"            = c("tnmd","scx","thbs4","col11a1","col12a1","prg4"),
  "Adventitial pi16+/col15a1+ perivascular" = c("pi16","col15a1","dpp4","pdgfra"),
  "Adventitial perivascular"        = c("col15a1","vit","vwa1","pdgfrb")
)


## Make sure gene case matches your object (Seurat is usually uppercase for mouse)
for(nm in names(fibro_markers)) fibro_markers[[nm]] <- intersect(toupper(fibro_markers[[nm]]), rownames(fibro))


####### SCINA
# install.packages("SCINA")
library(SCINA)
library(Seurat)
library(Matrix)
library(preprocessCore)


# use the default assay you normalized (e.g., "RNA" after SCTransform/NormalizeData)
expr <- GetAssayData(fibro, slot = "data")  # genes x cells
# example: take your list and coerce to the same case as expr
to_mouse_case <- function(v) {
  sapply(v, function(x) {
    x <- tolower(x)
    paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))
  }, USE.NAMES = FALSE)
}

# your original (edit as needed); convert case:
fibro_markers <- list(
  `PI16/DPP4 papillary adventitial` = c("Pi16","Dpp4","Ly6a","Ly6c1","Pdgfra","Col15a1"),
  `ECM-myofibroblast`               = c("Tagln","Acta2","Myh11","Col1a1","Col1a2","Fn1"),
  `Inflammatory fibroblasts`        = c("Ccl2","Ccl7","Il6","Il1b","Cxcl1","Nfkbia"),
  `Stress-activated infl. fibro`    = c("Fos","Jun","Atf3","Hspa1a","Hspa1b","Dusp1"),
  `Matrix ECM-remodeling`           = c("Mmp2","Mmp14","Timp2","Sparc","Col3a1","Lox"),
  `Perivascular fibroblasts`        = c("Col15a1","Vit","Vwa1","Rgs5","Pdgfrb","Kcnj8","Abcc9"),
  `WNT-responsive fibroblasts`      = c("Axin2","Wif1","Rspo1","Rspo3","Dkk3","Fgf10"),
  `Tenocyte fibroblasts`            = c("Tnmd","Scx","Thbs4","Col11a1","Col12a1","Prg4")
)

# make absolutely sure the case matches the matrix:
rown_low <- tolower(rownames(expr))
norm_markers <- lapply(fibro_markers, function(g) {
  # map case-insensitively to rownames(expr)
  hit <- match(tolower(g), rown_low, nomatch = 0)
  rownames(expr)[hit[hit > 0]]
})

# drop signatures with too few present genes; keep ≥3 (adjust if needed)
norm_markers <- lapply(norm_markers, unique)
norm_markers <- norm_markers[lengths(norm_markers) >= 3]

# sanity check: you MUST have at least 2 signatures after cleaning
lengths(norm_markers)
if (length(norm_markers) < 2) stop("After cleaning, fewer than 2 signatures remain.")

sig_union <- unique(unlist(norm_markers))
expr_sub  <- expr[intersect(sig_union, rownames(expr)), , drop = FALSE]

# drop zero-variance genes to help covariance estimation
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
expr_sub <- expr_sub[matrixStats::rowSds(expr_sub) > 0, , drop = FALSE]

# double-check we still have ≥2 signatures
norm_markers <- lapply(norm_markers, function(g) intersect(g, rownames(expr_sub)))
norm_markers <- norm_markers[lengths(norm_markers) >= 3]
if (length(norm_markers) < 2) stop("After cleaning, fewer than 2 signatures remain.")

scina_res <- SCINA(
  exp          = as.matrix(expr_sub),
  signatures   = norm_markers,
  max_iter     = 200,
  convergence_n= 10,
  rm_overlap   = FALSE,   # or FALSE if you prefer to keep shared markers
  allow_unknown= TRUE
)

# labels per cell
count<-table(scina_res$cell_labels)
write.csv(count, "scinaresultsfibro.csv")
fibro$SCINA <- scina_res$cell_labels  # add to Seurat object metadata

DimPlot(fibro, group.by = "SCINA", label = TRUE, repel = TRUE)




###ScType (fast, simple scoring)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)


## Build positive & negative gene sets (we’ll just use positives here)
gs_pos <- fibro_markers
gs_neg <- lapply(gs_pos, function(x) character(0))

scores <- sctype_score(scRNAseqData = GetAssayData(fibro, slot="data"),
                       scaled = FALSE, # TRUE if you pass scaled data
                       gs = gs_pos, gs2 = gs_neg)

## Assign top type per cell
fibro$ScType_label <- colnames(scores)[max.col(t(scores), ties.method = "first")]
DimPlot(fibro, group.by = "ScType_label", label=TRUE)


### scSorter (semi-supervised)
# install.packages("scSorter")
library(scSorter)

all_genes <- unique(unlist(fibro_markers))
sig <- matrix(0, nrow = length(all_genes), ncol = length(fibro_markers),
              dimnames = list(all_genes, names(fibro_markers)))
for(nm in names(fibro_markers)) sig[fibro_markers[[nm]], nm] <- 1

pred_scSorter <- scSorter(expdata = as.matrix(GetAssayData(fibro, slot="data")[rownames(sig), ]),
                          signature = sig)
fibro$scSorter_label <- pred_scSorter$Subtype

##Cell-ID (signature enrichment / reference-free)
# remotes::install_github("carmonalab/CellID")
library(CelliD)
library(Seurat)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

## =========================
## 0) Packages
## =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(CelliD)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

DefaultAssay(fibro) <- "RNA"

## =========================
## 1) Marker gene sets (~10 each)
##    (mouse TitleCase symbols)
## =========================
fibro_markers_extended <- list(
  `PI16/DPP4 papillary adventitial` = c(
    "Pi16","Dpp4","Ly6a","Ly6c1","Pdgfra","Col15a1",
    "Igfbp7","Mfap5","Gsn","Anxa3"
  ),
  `ECM-myofibroblast` = c(
    "Tagln","Acta2","Myh11","Col1a1","Col1a2","Fn1",
    "Cnn1","Tpm2","Postn","Lmod1"
  ),
  `Inflammatory fibroblasts` = c(
    "Ccl2","Ccl7","Il6","Il1b","Cxcl1","Nfkbia",
    "Cxcl2","Ccl8","Saa3","Ptgs2"
  ),
  `Stress-activated infl. fibro` = c(
    "Fos","Jun","Atf3","Hspa1a","Hspa1b","Dusp1",
    "Egr1","Junb","Fosl1","Ier3"
  ),
  `Matrix ECM-remodeling` = c(
    "Mmp2","Mmp14","Timp2","Sparc","Col3a1","Lox",
    "Mmp3","Mmp9","Col5a1","Thbs1"
  ),
  `Perivascular fibroblasts` = c(
    "Col15a1","Vit","Vwa1","Rgs5","Pdgfrb","Kcnj8","Abcc9",
    "Cspg4","Notch3","Mcam","Anpep"
  ),
  `WNT-responsive fibroblasts` = c(
    "Axin2","Wif1","Rspo1","Rspo3","Dkk3","Fgf10",
    "Lef1","Tcf7","Wnt2","Wnt5a"
  ),
  `Tenocyte fibroblasts` = c(
    "Tnmd","Scx","Thbs4","Col11a1","Col12a1","Prg4",
    "Fmod","Comp","Mkx","Col14a1"
  )
)

## =========================
## 2) Intersect with data + keep sets with >=10 present genes
## =========================
library(purrr)

all_genes <- rownames(fibro)
genesets <- fibro_markers_extended |>
  purrr::imap(function(g, nm){
    gg <- unique(intersect(g, all_genes))
    if (length(gg) < 10) message(sprintf("[Warn] %s has only %d genes present", nm, length(gg)))
    gg
  }) |>
  purrr::keep(~ length(.x) >= 10)


## =========================
## 3) Run MCA (stored as reduction "MCA")
## =========================
fibro <- RunMCA(fibro, assay = "RNA", slot = "data")
# Note: the reduction name is uppercase "MCA" in Seurat after RunMCA()

## Optional: run UMAP on MCA coords (nice for visualization)

fibro <- RunUMAP(fibro, reduction = "mca", dims = 1:50)


## =========================
## 4) Per-cell enrichment with CelliD (Hypergeometric test)
## =========================
# HGT: pathways x cells
enr_mat <- as.matrix(HGT)

# make it cells x pathways
enr_cellxpath <- t(enr_mat)

# best pathway per cell
best_ix <- max.col(enr_cellxpath, ties.method = "first")
labels  <- colnames(enr_cellxpath)[best_ix]   # pathway (subtype) names

# store & use as identities
fibro$CellID_label <- labels
Idents(fibro) <- "CellID_label"

# plots
DimPlot(fibro, reduction = "umap", group.by = "CellID_label",
                  label = TRUE, repel = TRUE) + ggtitle("CelliD labels (UMAP)")
 DimPlot(fibro, reduction = "mca",  group.by = "CellID_label",
                  label = TRUE, repel = TRUE) + ggtitle("CelliD labels (MCA)")
print(p_umap); print(p_mca)

marker_union <- unique(unlist(genesets))
DotPlot(fibro, features = head(marker_union, 60)) +
  RotatedAxis() + ggtitle("Marker DotPlot by CelliD subtype")


## =========================
## 5) DimPlots (UMAP & MCA) colored by CelliD labels
## =========================
p_umap <- DimPlot(fibro, reduction = "umap", group.by = "CellID_label", label = TRUE, repel = TRUE) +
  ggtitle("CelliD labels (UMAP)")
p_mca  <- DimPlot(fibro, reduction = "mca",  group.by = "CellID_label", label = TRUE, repel = TRUE) +
  ggtitle("CelliD labels (MCA)")

print(p_umap)
print(p_mca)

## =========================
## 6) DotPlot of markers across CelliD labels
## =========================
marker_union <- unique(unlist(genesets))
# Show up to, say, 60 markers if you want to keep the plot readable
marker_plot <- head(marker_union, 60)

# Reorder identities by CelliD label for the plot
Idents(fibro) <- fibro$CellID_label

p_dot <- DotPlot(fibro, features = marker_plot) +
  RotatedAxis() +
  ggtitle("Marker DotPlot by CelliD subtype")

print(p_dot)

## =========================
## (Optional) Save results
## =========================
# Save labels
# write.csv(data.frame(cell = colnames(fibro), label = fibro$CellID_label),
#           "cellid_labels.csv", row.names = FALSE)

# Save plots
# ggsave("cellid_umap.png", p_umap, width = 8, height = 6, dpi = 300)
# ggsave("cellid_mca.png",  p_mca,  width = 8, height = 6, dpi = 300)
# ggsave("cellid_dotplot.png", p_dot, width = 10, height = 7, dpi = 300)
