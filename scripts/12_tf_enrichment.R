source("scripts/lib/common.R")


# Average expression and accessibility across replicates.

rna <- rna_expression()
atac <- atac_expression()

rnaM2 <- as.matrix(data.frame(row.names = rownames(rna),
  Con_3W = rowMeans(rna[, c("3WS1", "3WS3", "3WS5")]),
  Con_4W = rowMeans(rna[, c("4WS1", "4WS3", "4WS5")]),
  Con_5W = rowMeans(rna[, c("5WS1", "5WS3", "5WS5")]),
  Con_6W = rowMeans(rna[, c("6WS1", "6WS3", "6WS5")]),
  Flu_3W = rowMeans(rna[, c("3WF1", "3WF3", "3WF5")]),
  Flu_4W = rowMeans(rna[, c("4WF1", "4WF3", "4WF5")]),
  Flu_5W = rowMeans(rna[, c("5WF1", "5WF3", "5WF5")]),
  Flu_6W = rowMeans(rna[, c("6WF1", "6WF3", "6WF5")])
))

atacM2 <- as.matrix(data.frame(row.names = rownames(atac),
  Con_3W = rowMeans(atac[, c("Con_3W_1", "Con_3W_2")]),
  Con_4W = rowMeans(atac[, c("Con_4W_1", "Con_4W_2")]),
  Con_5W = rowMeans(atac[, c("Con_5W_1", "Con_5W_2")]),
  Con_6W = rowMeans(atac[, c("Con_6W_1", "Con_6W_2")]),
  Flu_3W = rowMeans(atac[, c("Flu_3W_1", "Flu_3W_2")]),
  Flu_4W = rowMeans(atac[, c("Flu_4W_1", "Flu_4W_2")]),
  Flu_5W = rowMeans(atac[, c("Flu_5W_1", "Flu_5W_2")]),
  Flu_6W = rowMeans(atac[, c("Flu_6W_1", "Flu_6W_2")])
))

rnaM2z <- t(scale(t(rnaM2)))
atacM2z <- t(scale(t(atacM2)))


# Associate TF motif enrichments with accessibility and expression changes.

clusters <- read_tsv(ref("atac_clusters.tsv"))
gene_classes <- read_tsv(ref("rna_clusters.tsv"))
symbols <- read_tsv(ref("gene_symbols.tsv"))
de <- read_tsv(ref("de/SvF_Wk6.tsv.gz"), row.names = 1)

directory <- Sys.getenv("TF_ENRICHMENT_DIR", ref("tf_enrichment"))
sea <- list()

for (cluster in c("C1", "C2", "C3")) {
  path <- file.path(directory, paste0(cluster, ".tsv.gz"))

  if (!file.exists(path)) {
    path <- file.path(directory, cluster, "sea.tsv")
  }

  table <- read_tsv(path, comment.char = "#")
  table <- table[table$ALT_ID %in% rownames(rna), ]

  peaks <- clusters$Peaks[clusters$Cluster == cluster]
  accessibility <- colMeans(atacM2z[peaks, ])

  table$Cluster <- sub("C", "A", cluster)

  table$Pearson <- vapply(table$ALT_ID, function(gene) {
    cor(rnaM2z[gene, ], accessibility)
  }, numeric(1))

  sea[[cluster]] <- table[, setdiff(names(table), "DB")]
}

sea <- do.call(rbind, sea)
rownames(sea) <- NULL

sea$Symbol <- structure(symbols$Symbol, names = symbols$Gene)[sea$ALT_ID]
sea$Gene_Cluster <- structure(gene_classes$Class, names = gene_classes$Genes)[sea$ALT_ID]
sea$Flu6W_logFC <- de[sea$ALT_ID, "logFC"]

sea$Symbol[sea$ALT_ID == "AT5G11260"] <- "HY5"
sea$Symbol[sea$ALT_ID == sea$Symbol] <- NA

sea <- sea[!is.na(sea$Flu6W_logFC), ]
sea <- sea[sea$QVALUE < .05 & abs(sea$Flu6W_logFC) > 0, ]

p <- outdir("motifs/tf_enrichment")

write_tsv(sea, file.path(p, "TF_enrichment.tsv"))


# Plot motif enrichment against TF expression changes.

packages <- c("ggplot2", "ggrepel", "cowplot")

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  library(ggplot2)
  library(ggrepel)

  plots <- list()

  for (cluster in c("A1", "A2", "A3")) {
    keep <- which(
      sea$Cluster == cluster & sea$QVALUE < .001 & sea$Pearson > .5
    )

    selected <- sea[keep, ]

    xlimits <- if (cluster == "A3") c(NA, 0) else c(0, NA)

    plots[[cluster]] <- ggplot(selected,
      aes(Flu6W_logFC, -LOG_QVALUE, label = Symbol)) +
      geom_point() +
      geom_text_repel() +
      scale_x_continuous(limits = xlimits, name = "log FC in FLU week 4") +
      scale_y_continuous(limits = c(0, NA), name = "-log Q-value of enrichment") +
      ggtitle(cluster) +
      theme_bw() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        plot.title = element_text(hjust = .5)
      )
  }

  cowplot::save_plot(file.path(p, "TF_enrichment.pdf"),
    cowplot::plot_grid(plotlist = plots, nrow = 1),
    base_height = 2.75,
    base_width = 10
  )
} else {
  message("Install ggplot2, ggrepel and cowplot to draw the TF enrichment plots.")
}
