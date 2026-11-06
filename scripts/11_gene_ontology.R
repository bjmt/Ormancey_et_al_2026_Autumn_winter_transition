library(ComplexHeatmap)
library(simplifyEnrichment)

source("scripts/lib/common.R")


# Read the GO enrichment tables.
# --refresh reruns the same g:Profiler queries against the current database.

args <- commandArgs(trailingOnly = TRUE)

cls_d <- read_tsv(ref("rna_clusters.tsv"))
cls_d <- cls_d[grepl("cold", cls_d$Class), ]

levels <- c(
  "Inc_cold_early", "Inc_cold_linear", "Inc_cold_late",
  "Dec_cold_early", "Dec_cold_linear", "Dec_cold_late"
)

cls_d$Class <- factor(cls_d$Class, levels = levels)
cls_d <- cls_d[order(cls_d$Class), ]

goBP <- structure(vector("list", length(levels)), names = levels)

p <- outdir("gene_ontology")

for (cl in levels) {
  if ("--refresh" %in% args) {
    go <- gprofiler2::gost(
      cls_d$Genes[cls_d$Class == cl],
      "athaliana",
      correction_method = "fdr",
      user_threshold = 0.1
    )$result

    goBP[[cl]] <- go[go$source == "GO:BP", ]
  } else {
    goBP[[cl]] <- read_tsv(ref(paste0("go/", cl, ".tsv")))
  }

  write_tsv(goBP[[cl]], file.path(p, paste0(cl, ".tsv")))
}


# Add GO word clouds to the expression heatmap.

gowc <- anno_word_cloud_from_GO(
  as.character(cls_d$Class),
  lapply(goBP, function(x) x$term_id),
  fontsize_range = c(3, 10),
  bg_gp = grid::gpar(col = "#000000", fill = "#ffffff"),
  word_cloud_grob_param = list(
    max_width = grid::unit(30, "mm"),
    word_space = grid::unit(2, "pt"),
    line_space = grid::unit(2, "pt")
  ),
  stat = "pvalue",
  max_words = 15,
  min_stat = 1,
  exclude_words = c(
    "process", "metabolic", "response", "cellular", "cell", "compound",
    "biological", "macromolecule", "levels", "decreased", "organic"
  )
)

gcmr <- rna_expression()

z <- t(scale(t(log2(1 + gcmr))))
z[z > 2] <- 2
z[z < -2] <- -2

pdf(file.path(p, "RNAseq_clusters_GO.pdf"), width = 8.4, height = 5)

draw(Heatmap(z[cls_d$Genes, c(gS, gA[-(1:3)], gF[-(1:3)])],
  name = "Z-score",
  row_split = cls_d$Class,
  column_split = c(rep(1, 3), rep(2, 12), rep(3, 12), rep(4, 12)),
  row_title_rot = 0,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  use_raster = TRUE,
  raster_quality = 3,
  right_annotation = rowAnnotation(wc = gowc)
))

dev.off()

