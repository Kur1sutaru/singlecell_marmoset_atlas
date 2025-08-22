## ================ CONFIG =================
marker_csv      <- "C:/Users/crist/Downloads/work/meta_atlas/rat/rat_annotation/rat_tongue_markers_FULL_PLUS.csv"
out_html        <- "C:/Users/crist/Downloads/work/meta_atlas/rat/rat_annotation/celltype_subtype_predictions_WITH_WIDE.html"
res_use         <- 0.05   # set to the resolution that gives you ~18 clusters
force_recluster <- FALSE  # TRUE = always recluster to res_use; FALSE = use existing obj$seurat_clusters
## =========================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(readr); library(stringr)
  library(DT); library(htmlwidgets); library(htmltools)
})
if (!requireNamespace("plyr", quietly = TRUE)) install.packages("plyr")  # for plyr::mapvalues

stopifnot(exists("obj"), inherits(obj, "Seurat"))
DefaultAssay(obj) <- "RNA"

## 0) Repair duplicate meta column names (prevents dplyr errors)
if (anyDuplicated(colnames(obj@meta.data))) {
  colnames(obj@meta.data) <- make.unique(colnames(obj@meta.data), sep = "__")
  message("Repaired duplicate meta column names (pre-scoring).")
}

## 1) Ensure clustering (force to res_use if requested)
if (force_recluster || !"seurat_clusters" %in% colnames(obj@meta.data)) {
  if (!"pca" %in% names(obj@reductions)) {
    obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData(verbose = FALSE) |> RunPCA(verbose = FALSE)
  }
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = res_use, verbose = FALSE)
}
cat("Clusters now:", length(unique(obj$seurat_clusters)), "at res =", res_use, "\n")

## 2) Read marker sets and build MAJOR vs SUBTYPE panels
mk <- readr::read_csv(marker_csv, show_col_types = FALSE) |>
  dplyr::mutate(
    gene     = as.character(gene),
    celltype = as.character(celltype),
    major    = dplyr::if_else(stringr::str_detect(celltype, ":"), sub(":.*","",celltype), celltype)
  )

present     <- function(v) unique(v[v %in% rownames(obj)])
sets_major  <- lapply(split(mk$gene, mk$major), present);    sets_major <- sets_major[lengths(sets_major) > 0]
sets_sub    <- lapply(split(mk$gene, mk$celltype), present); sets_sub   <- sets_sub[lengths(sets_sub)   > 0]

## 3) Score with AddModuleScore, rename columns as "MAJ|xxx" / "SUB|yyy"
score_sets <- function(object, sets, tag) {
  if (!length(sets)) return(object)
  name_prefix <- paste0(tag, "_")
  object <- AddModuleScore(object, features = sets, name = name_prefix, search = FALSE)
  added <- paste0(name_prefix, seq_along(sets))
  newnames <- paste0(tag, "|", names(sets))
  colnames(object@meta.data)[match(added, colnames(object@meta.data))] <- newnames
  object
}
obj <- score_sets(obj, sets_major, "MAJ")
obj <- score_sets(obj, sets_sub,   "SUB")

## (defensive) Fix any new duplicate names
if (anyDuplicated(colnames(obj@meta.data))) {
  colnames(obj@meta.data) <- make.unique(colnames(obj@meta.data), sep = "__")
  message("Repaired duplicate meta column names (post-scoring).")
}

## 4) Summarize scores by cluster and pick best labels with margins (join by cluster)
clus <- "seurat_clusters"
summarize_calls <- function(meta, pick_prefix, cluster_col = "seurat_clusters") {
  cols <- grep(paste0("^", pick_prefix, "\\|"), colnames(meta), value = TRUE)
  stopifnot(length(cols) > 0)
  meta |>
    dplyr::select(dplyr::all_of(cluster_col), dplyr::all_of(cols)) |>
    dplyr::group_by(.data[[cluster_col]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(cols), ~mean(.x, na.rm=TRUE)), .groups="drop") |>
    tidyr::pivot_longer(cols = -dplyr::all_of(cluster_col), names_to = "label_pref", values_to = "score") |>
    dplyr::mutate(label = sub("^.*\\|", "", label_pref)) |>
    dplyr::group_by(.data[[cluster_col]]) |>
    dplyr::arrange(dplyr::desc(score), .by_group = TRUE) |>
    dplyr::summarise(
      seurat_clusters = as.character(dplyr::first(.data[[cluster_col]])),
      top_label   = dplyr::first(label),
      top_score   = dplyr::first(score),
      second_score = dplyr::nth(score, 2, default = NA_real_),
      margin      = top_score - dplyr::coalesce(second_score, 0),
      .groups     = "drop"
    ) |>
    dplyr::mutate(conf = dplyr::case_when(
      margin >= 0.10 ~ "high",
      margin >= 0.05 ~ "medium",
      TRUE ~ "low"
    ))
}
major_calls <- summarize_calls(obj@meta.data, "MAJ") |>
  dplyr::rename(major_label = top_label, major_score = top_score,
                major_margin = margin, major_conf = conf)
sub_calls <- summarize_calls(obj@meta.data, "SUB") |>
  dplyr::rename(sub_label = top_label, sub_score = top_score,
                sub_margin = margin, sub_conf = conf)

calls <- dplyr::left_join(major_calls, sub_calls, by = "seurat_clusters")

## 5) Map labels back to the object (use subtype when confident)
obj$auto_major     <- plyr::mapvalues(obj$seurat_clusters, from = calls$seurat_clusters, to = calls$major_label)
obj$auto_subtype   <- plyr::mapvalues(obj$seurat_clusters, from = calls$seurat_clusters, to = calls$sub_label)
obj$final_celltype <- ifelse(calls$sub_margin[match(obj$seurat_clusters, calls$seurat_clusters)] >= 0.05,
                             obj$auto_subtype, obj$auto_major)

## 6) Cluster-level table: include zero-rows so all clusters appear per sample
md <- obj@meta.data
sample_col <- c("sample","orig.ident","batch","library","dataset")
sample_col <- sample_col[sample_col %in% colnames(md)][1]
if (is.na(sample_col)) { md$sample <- "sample1"; sample_col <- "sample"; obj@meta.data <- md }

all_clusters <- as.character(sort(unique(obj$seurat_clusters)))

base_counts <- md |>
  dplyr::transmute(Sample = .data[[sample_col]], Cluster = as.character(.data[[clus]])) |>
  dplyr::count(Sample, Cluster, name = "Cells") |>
  tidyr::complete(Sample, Cluster = all_clusters, fill = list(Cells = 0)) |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(Proportion = dplyr::if_else(sum(Cells) > 0, Cells / sum(Cells), 0)) |>
  dplyr::ungroup()

pred_table <- base_counts |>
  dplyr::left_join(
    calls |>
      dplyr::mutate(Cluster = seurat_clusters) |>
      dplyr::select(Cluster, major_label, major_score, major_margin, major_conf,
                    sub_label, sub_score, sub_margin, sub_conf),
    by = "Cluster"
  ) |>
  dplyr::rename(`Major label` = major_label, `Major score` = major_score,
                `Major margin` = major_margin, `Major conf` = major_conf,
                `Subtype label` = sub_label, `Subtype score` = sub_score,
                `Subtype margin` = sub_margin, `Subtype conf` = sub_conf) |>
  dplyr::mutate(Cluster = as.integer(Cluster)) |>
  dplyr::arrange(Sample, Cluster)

## 7) Wide table: Sample × Final cell type proportions
wide <- obj@meta.data |>
  dplyr::transmute(Sample = .data[[sample_col]], Final = obj$final_celltype) |>
  dplyr::count(Sample, Final, name = "Cells") |>
  dplyr::group_by(Sample) |> dplyr::mutate(Proportion = Cells/sum(Cells)) |> dplyr::ungroup() |>
  dplyr::select(-Cells) |>
  tidyr::pivot_wider(names_from = Final, values_from = Proportion, values_fill = 0) |>
  dplyr::arrange(Sample)

## 8) Build DT widgets
dt1 <- DT::datatable(
  pred_table,
  rownames = FALSE, filter = "top",
  class = "stripe hover row-border order-column",
  extensions = c("Buttons","FixedHeader","ColReorder","Scroller"),
  options = list(
    dom="Bfrtip", buttons=c("copy","csv","excel","pdf","print"),
    pageLength=25, lengthMenu=c(10,25,50,100),
    scrollX=TRUE, deferRender=TRUE, scroller=TRUE, fixedHeader=TRUE, autoWidth=TRUE,
    order = list(
      list(which(names(pred_table)=="Sample")-1, "asc"),
      list(which(names(pred_table)=="Cluster")-1, "asc")
    )
  )
) |>
  DT::formatPercentage("Proportion", 1) |>
  DT::formatRound(c("Major score","Major margin","Subtype score","Subtype margin"), 3)

dt2 <- DT::datatable(
  wide, rownames = FALSE, filter = "top",
  class = "stripe hover row-border order-column",
  extensions = c("Buttons","FixedHeader","ColReorder","Scroller"),
  options = list(
    dom="Bfrtip", buttons=c("copy","csv","excel","pdf","print"),
    pageLength=25, lengthMenu=c(10,25,50,100),
    scrollX=TRUE, deferRender=TRUE, scroller=TRUE, fixedHeader=TRUE, autoWidth=TRUE
  )
)
for (c in setdiff(colnames(wide), "Sample")) dt2 <- DT::formatPercentage(dt2, c, 1)

## 9) Robust save: single page → non-self-contained → wrapper with iframes (Windows-proof)
page <- htmltools::browsable(htmltools::tagList(
  htmltools::tags$h2("Cluster-level predictions (counts, proportions, and calls)"),
  dt1,
  htmltools::tags$hr(),
  htmltools::tags$h2(htmltools::HTML("Sample &times; Final cell type proportions (wide)")),
  dt2
))

libdir <- file.path(dirname(out_html), paste0(tools::file_path_sans_ext(basename(out_html)), "_files"))

ok <- try({
  htmlwidgets::saveWidget(page, out_html, selfcontained = TRUE, title = "Cell type predictions")
}, silent = TRUE)

if (!inherits(ok, "try-error")) {
  message("Saved: ", normalizePath(out_html))
} else {
  message("Self-contained save failed; trying non-self-contained...")
  ok2 <- try({
    htmlwidgets::saveWidget(page, out_html, selfcontained = FALSE, libdir = libdir,
                            title = "Cell type predictions")
  }, silent = TRUE)
  
  if (!inherits(ok2, "try-error")) {
    message("Saved (non-self-contained): ", normalizePath(out_html))
  } else {
    message("Combined page still failed. Saving tables separately and creating an index wrapper...")
    
    base <- tools::file_path_sans_ext(basename(out_html))
    dir  <- normalizePath(dirname(out_html), winslash = "/", mustWork = TRUE)
    
    html1 <- file.path(dir, paste0(base, "_table1.html"))
    html2 <- file.path(dir, paste0(base, "_table2.html"))
    lib1  <- file.path(dir, paste0(base, "_table1_files"))
    lib2  <- file.path(dir, paste0(base, "_table2_files"))
    
    htmlwidgets::saveWidget(dt1, html1, selfcontained = FALSE, libdir = lib1, title = "Cluster-level predictions")
    htmlwidgets::saveWidget(dt2, html2, selfcontained = FALSE, libdir = lib2, title = "Sample x Final proportions")
    
    wrapper <- sprintf(
      '<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Cell type predictions</title>
<style>
  body { font-family: system-ui, Segoe UI, Arial, sans-serif; margin: 20px; }
  h2 { margin: 16px 0 8px; }
  .card { margin: 12px 0; border: 1px solid #e5e7eb; border-radius: 8px; padding: 12px; }
  iframe { width: 100%%; height: 70vh; border: 0; }
</style>
</head>
<body>
  <div class="card">
    <h2>Cluster-level predictions (counts, proportions, and calls)</h2>
    <iframe src="%s"></iframe>
  </div>
  <div class="card">
    <h2>Sample &times; Final cell type proportions (wide)</h2>
    <iframe src="%s"></iframe>
  </div>
</body>
</html>', basename(html1), basename(html2))
    
    writeLines(wrapper, out_html)
    message("Saved wrapper index: ", normalizePath(out_html))
    message("Keep with the wrapper:\n  - ", html1, "\n  - ", html2, "\n  - ", lib1, "\n  - ", lib2)
  }
}
