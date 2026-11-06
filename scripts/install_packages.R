if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}


cran <- c(
  "readr", "dplyr", "ggplot2", "ggrepel", "cowplot", "igraph", "matrixStats",
  "gprofiler2", "circlize", "viridisLite"
)

bioc <- c(
  "edgeR", "RUVSeq", "GenomicFeatures", "rtracklayer", "Biostrings", "ComplexHeatmap",
  "simplifyEnrichment", "universalmotif", "variancePartition"
)


missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}


missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
