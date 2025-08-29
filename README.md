scrna_mucositis_tutorial

# 🧬 Single-Cell RNA-seq Tutorial in R  
*Using mucositis dataset (GSE131630) as an example*  

---

## 📌 1. Introduction  
This tutorial introduces **R basics** and then walks through a **single-cell RNA-seq analysis** workflow using data from the paper:  

> Zhao et al. *Single-cell RNA sequencing analysis reveals alginate oligosaccharides preventing chemotherapy-induced mucositis*. Mucosal Immunology (2020)  

We will:  
1. Learn R fundamentals  
2. Explore single-cell analysis with **Seurat**  
3. Identify intestinal cell types and treatment effects  
4. Perform simple differential expression analysis  

---

## 📘 2. R Basics  

### Install and load packages
```r
# Install (only once)
install.packages(c("tidyverse", "Seurat"))

# Load
library(tidyverse)
library(Seurat)

## Data types & structures

x <- c(1, 2, 3, 4)         # numeric vector
y <- c("Busulfan", "AOS")  # character vector
df <- data.frame(sample = y, value = c(10, 20)) # data frame

head(df)

## Simple plotting

ggplot(df, aes(x = sample, y = value, fill = sample)) +
  geom_bar(stat = "identity") +
  theme_minimal()


## Load scRNA-seq Data

# Raw data is available at GEO: GSE131630. After downloading and processing with cellranger, load it into R:

# Load 10x counts (folder must contain barcodes.tsv, genes.tsv, matrix.mtx)
counts <- Read10X(data.dir = "path/to/GSE131630/")

# Create Seurat object
seu <- CreateSeuratObject(counts = counts, project = "Mucositis", min.cells = 3, min.features = 200)

# Add mitochondrial content
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

# Quick QC
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalization and Clustering

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu)

# Clustering
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:20)

DimPlot(seu, label = TRUE, reduction = "umap")

## Cell Type Identification

Use known marker genes from the paper:

Stem/TA cells: Olfm4, Slc12a2

Enterocytes: Mttp, Clca1

EEC: Fxyd3

Tuft: Rac2

Goblet: Ang4

Paneth: Defa31

```
FeaturePlot(seu, features = c("Olfm4", "Slc12a2", "Mttp", "Clca1", "Fxyd3", "Rac2", "Ang4", "Defa31"))
```

Assign cluster identities:

```
new.cluster.ids <- c("Stem", "TA", "Enterocyte", "Goblet", "Paneth", "EEC", "Tuft")
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
DimPlot(seu, reduction = "umap", label = TRUE)

```


Differential Expression (Treatment Effects)

Assuming metadata has treatment labels (AOS vs Busulfan):

Idents(seu) <- "treatment"
degs <- FindMarkers(seu, ident.1 = "Busulfan", ident.2 = "AOS", min.pct = 0.25)
head(degs)


Plot top DEGs:

VlnPlot(seu, features = c("Vil1","Defa31","Clca1"), group.by = "treatment")

🚀 7. Extensions (Optional Advanced)

Trajectory inference: use monocle3 for pseudotime

Regulatory networks: use SCENIC

Pathway enrichment: use clusterProfiler

✅ 8. Summary

Learned R basics (data frames, plotting)

Performed Seurat workflow (QC → clustering → UMAP)

Identified intestinal cell types

Found treatment differences between Busulfan vs AOS

📂 References

Zhao et al., Mucosal Immunology 2020

Seurat tutorials: https://satijalab.org/seurat/

Monocle: http://cole-trapnell-lab.github.io/monocle3/
