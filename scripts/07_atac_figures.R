source("scripts/lib/common.R")


# Scale accessibility within peaks and group by cluster.

rp2 <- atac_expression()
rnaRuvM <- rna_expression()

p2cls <- read_tsv(ref("atac_clusters.tsv"))
prom2peak <- read_tsv(ref("prom2peak.tsv"))

p <- outdir("atac/figures")

rp2zO <- t(scale(t(log2(1 + rp2))))

rp2f <- rp2[(matrixStats::rowMaxs(rp2) - matrixStats::rowMins(rp2)) > 28, ]
rp2f <- rp2f[!matrixStats::rowAnyNAs(rp2f), ]

rp2z <- t(scale(t(log2(1 + rp2f))))
rp2z[rp2z > 2] <- 2
rp2z[rp2z < -2] <- -2

write_tsv(
  data.frame(Peak = p2cls$Peaks, rp2z[p2cls$Peaks, ]),
  file.path(p, "cluster_accessibility_zscores.tsv")
)

if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)

  pdf(file.path(p, "ATACseq_clusters.pdf"), width = 4, height = 5)

  draw(Heatmap(rp2z[p2cls$Peaks, ],
    name = "Z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    row_split = factor(p2cls$Cluster, levels = c("C1", "C2", "C3")),
    column_split = factor(substr(colnames(rp2), 1, 3), levels = c("Con", "Flu")),
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    use_raster = TRUE,
    raster_quality = 3
  ))

  dev.off()
} else {
  message("Install ComplexHeatmap to draw the cluster heatmap.")
}


# Promoter accessibility at the top 500 up- and down-regulated genes.

prom2peak3 <- prom2peak[
  prom2peak$name %in% rownames(rp2) &
    prom2peak$Gene %in% rownames(rnaRuvM),
]

con_flu_6w_de <- read_tsv(ref("de/SvF_Wk6.tsv.gz"), row.names = 1)

ord1 <- colnames(rp2)[order(substr(colnames(rp2), 1, 3))]

rp2zO2 <- rp2zO
rp2zO2[rp2zO > 2] <- 2
rp2zO2[rp2zO < -2] <- -2

n_peaks <- c()

for (direction in c("up", "down")) {
  con_flu_6w_de <- con_flu_6w_de[
    order(con_flu_6w_de$logFC, decreasing = direction == "up"),
  ]

  peak_ids <- prom2peak3$name[
    prom2peak3$Gene %in% rownames(con_flu_6w_de)[1:500]
  ]

  accessibility <- rp2zO2[peak_ids, ord1]

  pdf(file.path(p, paste0("top500_", direction, ".pdf")), width = 4.5, height = 5)

  par(mar = c(7, 4, 1, 1))

  boxplot(accessibility,
    ylab = "Z-score",
    xlab = NA,
    las = 2,
    cex.axis = .7
  )

  abline(h = 0, lty = 2, lwd = 2, col = "red")

  dev.off()

  write_tsv(
    data.frame(Peak = peak_ids, accessibility),
    file.path(p, paste0("top500_", direction, ".tsv"))
  )

  n_peaks[direction] <- length(peak_ids)
}

write_tsv(
  data.frame(Direction = names(n_peaks), Peaks = unname(n_peaks)),
  file.path(p, "top500_peak_counts.tsv")
)

