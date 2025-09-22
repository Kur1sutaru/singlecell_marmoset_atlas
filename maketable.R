# Install if needed:
# install.packages("gt")
library(gt)

# Data frame summarizing your 21 genes
fibro_genes_info <- data.frame(
  Gene = c("Fmod","Chad","Cilp2","Col28a1","Prg4","Spp1","Grem1","Wif1","Socs3","Tnfrsf12a",
           "Pik3r5","Far2","Tfap2b","Zbtb7c","Tenm2","Amph","Angpt4","Hhatl","Kank1","Smim41","Prss12"),
  Function = c(
    "Collagen fibrillogenesis regulator",
    "Cartilage adhesion ECM protein",
    "Cartilage ECM structural protein",
    "Collagen XXVIII, basement/ECM",
    "Lubricin, cartilage/tendon lubrication",
    "Osteopontin, matricellular, inflammation",
    "BMP antagonist, stromal progenitors",
    "Wnt inhibitor, signaling modulation",
    "JAK/STAT negative feedback regulator",
    "TWEAK receptor, fibroblast activation",
    "PI3Kγ regulatory subunit, signaling",
    "Lipid metabolism enzyme",
    "Transcription factor AP-2β, mesenchyme",
    "BTB-zinc finger transcription repressor",
    "Adhesion receptor (teneurin family)",
    "Endocytosis adaptor, neuronal",
    "Angiopoietin-4, vascular factor",
    "Hedgehog acyltransferase-like",
    "Cytoskeletal/tumor suppressor",
    "Poorly characterized SMIM",
    "Neurotrypsin protease"
  ),
  scRNAseq_evidence = c(
    "Strong — tendon, dermal, synovial fibroblasts",
    "Strong — tendon/ligament fibroblasts",
    "Strong — joint/tendon fibroblasts",
    "Moderate — adventitial fibroblasts",
    "Strong — synovial/tendon fibroblasts",
    "Strong — activated/inflammatory fibroblasts",
    "Strong — stromal fibroblasts (intestine, skin)",
    "Strong — WNT-responsive fibroblasts",
    "Strong — inflammatory/stress-activated fibroblasts",
    "Strong — injury/fibrosis fibroblasts",
    "Moderate — inflammatory fibroblast subsets",
    "Moderate — dermal fibroblast subsets",
    "Moderate — mesenchymal/fibroblast subsets",
    "Moderate — stromal subsets",
    "Moderate — fibroblasts in tongue/joint datasets",
    "Weak — mostly neuronal, rare fibroblast detection",
    "Moderate — perivascular/adventitial fibroblasts",
    "Weak — rare expression",
    "Moderate — scattered fibroblasts",
    "Weak — very limited data",
    "Weak — not canonical fibroblast"
  ),
  Notes = c(
    "Well-established ECM regulator",
    "ECM adhesion protein",
    "Cartilage matrix protein",
    "Specialized collagen in stromal fibroblasts",
    "Synovial/tendon marker",
    "Inflammatory fibrosis marker",
    "Stromal stem-like fibroblast marker",
    "Canonical WNT antagonist",
    "Stress response regulator",
    "Fibroblast activation receptor",
    "Signaling regulator in fibroblasts",
    "Metabolic stromal gene",
    "Developmental regulator",
    "Transcriptional control",
    "Adhesion, less specific",
    "Possible neuronal contamination",
    "Vascular signaling, peri-fibroblasts",
    "Unclear role",
    "Cytoskeleton regulation",
    "Uncharacterized",
    "Neuronal protease"
  )
)

# Make a nice gt table
fibro_table <- fibro_genes_info |>
  gt() |>
  tab_header(title = md("**Fibroblast-relevant genes and single-cell evidence**")) |>
  cols_label(
    Gene = "Gene",
    Function = "Known Function",
    scRNAseq_evidence = "Single-cell fibroblast evidence",
    Notes = "Notes"
  ) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = Gene)) |>
  tab_options(table.font.size = 11)

# Save to PDF/PNG
gtsave(fibro_table, "fibro_genes_summary_table.png")
gtsave(fibro_table, "fibro_genes_summary_table.pdf")
