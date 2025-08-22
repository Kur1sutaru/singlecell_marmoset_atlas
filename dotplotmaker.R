## ---- Setup ----
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

setwd("C:/Users/crist/Downloads/work/meta_atlas/rat/rat_annotation")
## ===========================
## EDIT ME (1): your object & identities
## ===========================
rat <- your_seurat_object              # <— replace with your Seurat object
Idents(rat) <- rat$subtype             # <— column that holds your subtypes

## (Optional) Order subtypes the way you want on the Y axis:
# Idents(rat) <- factor(Idents(rat), levels = c("Endothelial", "Lymphatic endothelial", ...))

## ===========================
## Marker panels (same layout as your figure)
## Tip: If your data uses Title-Case (e.g., Pecam1) keep that casing.
## ===========================
features_list <- list(
  "Endothelial" = c("Pecam1","Kdr","Klf2","Klf4","Tek","Cdh5","Icam1","Vwf"),
  "Lymphatic endothelial" = c("Prox1","Pdpn","Reln","Mmrn1"),
  "Pericytes" = c("Pdgfrb","Rgs5","Abcc9","Kcnj8","Mcam","Cspg4"),
  "Smooth Muscle" = c("Acta2","Tagln","Myh11","Cnn1","Notch3"),
  "Skeletal Muscle" = c("Des","Ckm","Ckmt2","Ttn","Myh1"),
  "skeletal\nmuscle\nsatellite cells" = c("Pax7","Six1","Itga7","Myf5"),
  "Myogenic\nprogenitor\ncells" = c("Myod1","Myog","Tcf3","Tcf12"),
  "Fibroblasts" = c("Col1a1","Col1a2","Dcn","Lum","Pdgfra","S100a4","Ctgf"),
  "Epithelial Cells" = c("Krt5","Krt14","Krt8","Krt18","Krt13","Krt10","Cldn1","Epcam","Ocln"),
  "T cell" = c("Cd3d","Cd3e","Trac","Lck","Cd2","Cd7","Cxcr4","Il7r"),
  "Macrophages" = c("Lyz2","Lcp1","Ptprc","Adgre1","Csf1r","Itgam"),
  "Schwann\nCells" = c("Mpz","Pmp22","Prx","Plp1","Mbp"),
  "Mast" = c("Kit","Cma1","Tpsab1","Gata2")
)

## ===========================
## (Optional) Map markers to rat symbols
## Skip if your markers already match rownames(rat)
## ===========================
map_to_rat <- function(genes){
  # try to match case-insensitively to your matrix first
  rn <- rownames(rat)
  found <- sapply(genes, function(g){
    ix <- which(tolower(rn) == tolower(g))
    if (length(ix)) rn[ix[1]] else g
  })
  unname(found)
}
features_list <- lapply(features_list, map_to_rat)

## Report missing genes (helps spot symbol differences)
all_feats <- unique(unlist(features_list))
missing <- setdiff(all_feats, rownames(rat))
if (length(missing)) {
  message("Markers not found in object (check casing/aliases): ",
          paste(missing, collapse=", "))
}

## ---------------------------
## Plot
## ---------------------------
p <- DotPlot(
  rat,
  features = features_list,
  assay = DefaultAssay(rat),
  cols = c("#eaefff","#1d2f99"),     # light→dark blue, similar to your legend
  col.min = 0,
  col.max = 2.5,                     # matches “Average Expression up to 2.5”
  dot.scale = 4
) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.major = element_line(size = 0.2, colour = "grey90")
  ) +
  labs(title = "Rat tongue", x = NULL, y = NULL)

## Show it
print(p)

## Save high‑res (adjust size to your subtype count)
ggsave("rat_dotplot.png", p, width = 16, height = 7, dpi = 400, bg = "white")
