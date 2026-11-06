source("scripts/lib/common.R")

p <- outdir("motifs/figures")


# Enrichment of the discovered motifs in the three accessibility clusters.

sea <- list()
cen <- list()

read_motif_results <- function(path) {
  header <- readLines(path, n = 1)
  columns <- strsplit(header, "\t", fixed = TRUE)[[1]]

  tab <- read.delim(path,
    skip = 1,
    header = FALSE,
    col.names = c(columns, "Unused"),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "#",
    fill = TRUE
  )

  tab[, columns]
}

for (i in 1:3) {
  sea[[i]] <- read_motif_results(
    ref(paste0("motif_results/sea2/C", i, ".tsv"))
  )

  cen[[i]] <- read_motif_results(
    ref(paste0("motif_results/cen2/C", i, ".tsv"))
  )
}

tfs <- unique(c(
  sea[[1]]$ID, sea[[2]]$ID, sea[[3]]$ID,
  cen[[1]]$motif_id, cen[[2]]$motif_id, cen[[3]]$motif_id
))

tfs <- tfs[!is.na(tfs)]

nt1 <- function(x) {
  x[is.na(x)] <- 1

  x
}

mres2 <- data.frame(row.names = tfs)

for (i in 1:3) {
  sea_q <- structure(sea[[i]]$QVALUE, names = sea[[i]]$ID)
  cen_q <- structure(cen[[i]]$`adj_p-value`, names = cen[[i]]$motif_id)

  mres2[[paste0("sea", i)]] <- nt1(sea_q[tfs])
  mres2[[paste0("cen", i)]] <- nt1(cen_q[tfs])
}

# Order the enrichment columns and select motifs.

mres2 <- mres2[, c("sea1", "sea2", "sea3", "cen1", "cen2", "cen3")]
mres2L <- -log10(as.matrix(mres2))

mres2L4 <- mres2L[matrixStats::rowMaxs(mres2L) > 20, ]
mres2L4 <- mres2L4[c(1, 3, 4, 5, 6, 8, 13, 15), ]

write_tsv(
  data.frame(Motif = paste0("M", 1:8), Source = rownames(mres2L4), mres2L4),
  file.path(p, "motif_enrichment.tsv")
)

mres2L42 <- mres2L4
mres2L42[mres2L42 > 20] <- 20

if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)

  pdf(file.path(p, "motif_enrichment.pdf"), width = 4.5, height = 3.5)

  draw(Heatmap(mres2L42[1:7, 1:3],
    name = "-log10(Adj. P-value)",
    row_labels = paste0("M", 1:7),
    column_labels = c("A1", "A2", "A3"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = viridisLite::magma(10)
  ))

  dev.off()
}

if (requireNamespace("universalmotif", quietly = TRUE)) {
  library(universalmotif)

  motifs <- read_meme(ref("motifs.meme"))

  plot <- view_motifs(motifs[1:7],
    names.pos = "right",
    tryRC = FALSE,
    min.overlap = 99
  )

  ggplot2::ggsave(file.path(p, "motif_logos.pdf"),
    plot = plot,
    width = 7, height = 5
  )
}

