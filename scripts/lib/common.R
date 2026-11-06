# Run scripts from the repository root.
stopifnot(file.exists("metadata/rnaseq_samples.tsv"))


ref <- function(name) {
  file.path("data/reference", name)
}


outdir <- function(name) {
  path <- file.path("results", name)

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  path
}


read_tsv <- function(path, ...) {
  read.delim(path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    ...
  )
}


write_tsv <- function(x, path, row.names = FALSE) {
  write.table(x, path,
    sep = "\t",
    quote = FALSE,
    row.names = row.names,
    col.names = TRUE,
    na = "NA"
  )
}


read_matrix <- function(path, ...) {
  as.matrix(read_tsv(path, row.names = 1, ...))
}


samples <- function(assay) {
  read_tsv(paste0("metadata/", assay, "_samples.tsv"))
}


rna_expression <- function() {
  p <- Sys.getenv("RNA_RPKM", "data/processed/GeneRPKM.tsv")

  if (!file.exists(p)) {
    stop("Obtain processed matrices first: bash scripts/download_processed.sh")
  }

  m <- read_matrix(p)
  s <- samples("rnaseq")

  if (all(s$sample %in% colnames(m))) {
    m <- m[, s$sample, drop = FALSE]
  } else if (all(s$analysis %in% colnames(m))) {
    m <- m[, s$analysis, drop = FALSE]
  } else {
    stop("RNA-seq matrix columns do not match metadata/rnaseq_samples.tsv")
  }

  colnames(m) <- s$analysis

  m
}


atac_expression <- function() {
  p <- Sys.getenv("ATAC_RPKM", "data/processed/atac_peak_rpkm.tsv")

  if (!file.exists(p)) {
    stop("Obtain processed matrices first: bash scripts/download_processed.sh")
  }

  m <- read_matrix(p)
  s <- samples("atac")

  if (all(s$sample %in% colnames(m))) {
    m <- m[, s$sample, drop = FALSE]
  } else if (all(s$analysis %in% colnames(m))) {
    m <- m[, s$analysis, drop = FALSE]
  } else {
    stop("ATAC-seq matrix columns do not match metadata/atac_samples.tsv")
  }

  colnames(m) <- s$analysis

  m
}


gS <- c(1:3, 10:12, 19:21, 28:30, 37:39)
gA <- c(1:3, 4:6, 13:15, 22:24, 31:33)
gF <- c(1:3, 7:9, 16:18, 25:27, 34:36)

colorBlindBlack8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
