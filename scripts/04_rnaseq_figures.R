source("scripts/lib/common.R")


# Differentially expressed genes against the matched constant-temperature controls.

p <- outdir("rnaseq/figures")
de_counts <- list()

for (condition in c("A", "F")) {
  for (week in 3:6) {
    de <- read_tsv(
      ref(paste0("de/Sv", condition, "_Wk", week, ".tsv.gz")),
      row.names = 1
    )

    label <- if (condition == "A") "Average" else "Fluctuating"

    de_counts[[paste(condition, week)]] <- data.frame(
      Week = week - 2,
      Condition = label,
      Direction = c("Up", "Down"),
      Count = c(
        sum(de$logFC >= 1 & de$FDR < .01),
        sum(de$logFC <= -1 & de$FDR < .01)
      )
    )
  }
}

de_counts <- do.call(rbind, de_counts)
rownames(de_counts) <- NULL

write_tsv(de_counts, file.path(p, "differential_gene_counts.tsv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  for (direction in c("Up", "Down")) {
    plot <- ggplot(de_counts[de_counts$Direction == direction, ]) +
      geom_bar(aes(Week, Count, fill = Condition),
        position = position_dodge(.9), stat = "identity") +
      scale_y_continuous(limits = c(0, 2500),
        name = paste0(direction, "-regulated genes")) +
      scale_fill_manual(values = colorBlindBlack8[c(1, 6)]) +
      xlab("Week since start of cold") +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "bottom")

    ggsave(file.path(p, paste0(direction, "_regulated_genes.pdf")),
      plot = plot,
      width = 3.5, height = 4
    )
  }
}


# Plot gene expression profiles.

gcmr <- rna_expression()
p <- outdir("rnaseq/figures")

make_plot <- function(g, name = NULL, z = FALSE, isExpr = FALSE, xlab = NA) {
  if (is.null(name)) {
    if (!isExpr)
      name <- g
    else
      name <- ""
  }
  if (isExpr) {
    e <- g
  } else {
    e <- gcmr[g, ]
  }
  if (z) e <- (e - mean(e)) / sd(e)
  gA <- c(1:3, 4:6, 13:15, 22:24, 31:33)
  gS <- c(1:3, 10:12, 19:21, 28:30, 37:39)
  gF <- c(1:3, 7:9, 16:18, 25:27, 34:36)
  y_max <- max(e) * 1.1
  y_min <- 0
  if (z) y_min <- abs(min(e)) * -1.1
  plot(rep(1:5, each = 3) + rnorm(15, 0, .1), e[gS], xlim = c(0.7, 5.3),
    pch = 19, ylim = c(y_min, y_max), xaxt = "n", xlab = xlab,
    col = "#afababff", bty = "n",
    xaxs="i", yaxs="i",
    ylab = if (z) "Z-score" else "RPKM", main = name,
  )
  axis(1, 1:5, paste0(0:4, "W"))
  points(rep(2:5, each = 3) + rnorm(12, 0, .1), e[gF[-(1:3)]], pch = 19, col = "#0070c0ff")
  points(rep(2:5, each = 3) + rnorm(12, 0, .1), e[gA[-(1:3)]], pch = 19, col = "#ffc000ff")

  a1 <- tapply(e[c(1:3, gF[-(1:3)])], rep(1:5, each = 3), mean)
  s1 <- tapply(e[c(1:3, gF[-(1:3)])], rep(1:5, each = 3), sd) / sqrt(3)
  arrows(1:5, a1 - s1, 1:5, a1 + s1, length=0.05, angle=90, code=3)
  a2 <- tapply(e[c(1:3, gA[-(1:3)])], rep(1:5, each = 3), mean)
  s2 <- tapply(e[c(1:3, gA[-(1:3)])], rep(1:5, each = 3), sd) / sqrt(3)
  arrows(1:5, a2 - s2, 1:5, a2 + s2, length=0.05, angle=90, code=3)
  a3 <- tapply(e[gS], rep(1:5, each = 3), mean)
  s3 <- tapply(e[gS], rep(1:5, each = 3), sd) / sqrt(3)
  arrows(1:5, a3 - s3, 1:5, a3 + s3, length=0.05, angle=90, code=3)

  lines(1:5, tapply(e[c(1:3, gF[-(1:3)])], rep(1:5, each = 3), mean), lwd = 2, col = "#0070c0ff")
  lines(1:5, tapply(e[c(1:3, gA[-(1:3)])], rep(1:5, each = 3), mean), lwd = 2, col = "#ffc000ff")
  lines(1:5, tapply(e[gS], rep(1:5, each = 3), mean), lwd = 2, col = "#afababff")

  abline(h=0, lwd = 2)

  # legend("topleft", legend = c("Constant", "Average cold", "Fluctuating cold"),
  #   pch = 19, col = c("#afababff", "#ffc000ff", "#0070c0ff"), bty = "n")
}



genes <- c(
  EMB2746 = "AT5G63420",
  KIN1 = "AT5G15960",
  RD29B = "AT5G52300",
  HCAR = "AT1G04620",
  SCPL49 = "AT3G10410",
  RD29A = "AT5G52310",
  TPR = "AT1G80130",
  EARLI1 = "AT4G12480",
  AZI1 = "AT4G12470",
  VIN3 = "AT5G57380",
  CBF1 = "AT4G25490",
  CBF2 = "AT4G25470",
  CBF3 = "AT4G25480",
  CCA1 = "AT2G46830",
  HY5 = "AT5G11260",
  ELF4 = "AT2G40080",
  COR28 = "AT4G33980",
  HYH = "AT3G17609",
  COR27 = "AT5G42900",
  CO = "AT5G15840",
  FKF1 = "AT1G68050",
  PRR5 = "AT5G24470",
  PRR9 = "AT2G46790"
)

for (name in names(genes)) {
  pdf(file.path(p, paste0(name, ".pdf")), width = 5.8, height = 4)

  make_plot(genes[name], name, xlab = "Week since start of cold")

  dev.off()
}


# Scale expression within genes and group by cluster.

cls_d <- read_tsv(ref("rna_clusters.tsv"))
cls_d <- cls_d[grepl("cold", cls_d$Class), ]

levels <- c(
  "Inc_cold_early", "Inc_cold_linear", "Inc_cold_late",
  "Dec_cold_early", "Dec_cold_linear", "Dec_cold_late"
)

cls_d$Class <- factor(cls_d$Class, levels = levels)
cls_d <- cls_d[order(cls_d$Class), ]


# Plot the first 20 genes in the early up-regulated cluster.

early_genes <- cls_d$Genes[cls_d$Class == "Inc_cold_early"][1:20]

pdf(file.path(p, "U1_gene_profiles.pdf"), width = 12, height = 12)

par(mfrow = c(5, 4), mar = c(3, 3, 2, 1), mgp = c(2, .6, 0))

for (gene in early_genes) {
  make_plot(gene, gene)
}

dev.off()


z <- t(scale(t(log2(1 + gcmr))))
z[z > 2] <- 2
z[z < -2] <- -2

write_tsv(
  data.frame(Gene = cls_d$Genes, z[cls_d$Genes, ]),
  file.path(p, "cluster_expression_zscores.tsv")
)


# Heatmap and cluster profiles.

if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)

  pdf(file.path(p, "RNAseq_clusters.pdf"), width = 8.4, height = 5)

  draw(Heatmap(z[cls_d$Genes, c(gS, gA[-(1:3)], gF[-(1:3)])],
    name = "Z-score",
    row_split = cls_d$Class,
    column_split = c(rep(1, 15), rep(2, 12), rep(3, 12)),
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    use_raster = TRUE,
    raster_quality = 3
  ))

  dev.off()
} else {
  message("Install ComplexHeatmap to draw the cluster heatmap.")
}

for (cl in levels) {
  pdf(file.path(p, paste0(cl, ".pdf")), width = 8, height = 4)

  boxplot(z[cls_d$Genes[cls_d$Class == cl], c(gS, gA[-(1:3)], gF[-(1:3)])],
    outline = FALSE,
    ylab = "Z-score",
    xlab = "Sample",
    main = cl,
    las = 2
  )

  dev.off()
}


# Average expression within each cluster and condition.

z_mean <- t(scale(t(log2(1 + gcmr))))
profiles <- list()

for (cl in levels) {
  genes <- cls_d$Genes[cls_d$Class == cl]

  for (condition in c("22C", "AVG", "FLU")) {
    columns <- switch(condition, "22C" = gS, "AVG" = gA, "FLU" = gF)
    means <- colMeans(z_mean[genes, columns])

    profiles[[paste(cl, condition)]] <- data.frame(
      Class = cl,
      Condition = condition,
      Week = 0:4,
      Zscore = as.numeric(tapply(means, rep(0:4, each = 3), mean))
    )
  }
}

profiles <- do.call(rbind, profiles)
rownames(profiles) <- NULL

write_tsv(profiles, file.path(p, "cluster_mean_profiles.tsv"))

pdf(file.path(p, "cluster_mean_profiles.pdf"), width = 10, height = 6)

par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

for (cl in levels) {
  profile <- profiles[profiles$Class == cl, ]
  colours <- c("22C" = "#afababff", "AVG" = "#ffc000ff", "FLU" = "#0070c0ff")

  plot(NA,
    xlim = c(0, 4), ylim = range(profile$Zscore),
    xlab = "Week since start of cold", ylab = "Mean Z-score", main = cl
  )

  for (condition in names(colours)) {
    values <- profile[profile$Condition == condition, ]

    lines(values$Week, values$Zscore,
      type = "b", pch = 19, col = colours[condition], lwd = 2
    )
  }
}

dev.off()


# Compare the up-regulated clusters with the Liu et al. expression time course.

liu <- read_tsv(ref("Liu2022_fpkm.tsv.gz"))
liuM <- as.matrix(liu[, -1])
rownames(liuM) <- liu[[1]]

liuMz <- t(scale(t(liuM)))
liuMz[liuMz > 2] <- 2
liuMz[liuMz < -2] <- -2

for (cl in levels[1:3]) {
  g <- cls_d$Genes[cls_d$Class == cl]
  g <- g[g %in% rownames(liuMz)]

  pdf(file.path(p, paste0("Liu2022_", cl, ".pdf")), width = 7, height = 4)

  boxplot(liuMz[g, 1:12],
    outline = FALSE,
    ylab = "Z-score",
    main = cl,
    las = 2
  )

  dev.off()
}


# Plot VIN3 and genes with expression correlated with VIN3.

vin3 <- "AT5G57380"
liu2 <- liuM[, 1:12]
liu2z <- t(scale(t(log2(1 + liu2))))

pdf(file.path(p, "Liu2022_VIN3.pdf"), width = 5.8, height = 4)

par(mar = c(6, 4, 2, 1))

plot(1:12, liu2z[vin3, ],
  pch = 19, ylim = c(-2, 2), xaxt = "n",
  xlab = NA, ylab = "Z-score", main = "VIN3"
)

axis(1, 1:12, colnames(liu2), las = 2, cex.axis = .8)

dev.off()

for (cl in levels[1:3]) {
  genes <- cls_d$Genes[cls_d$Class == cl]
  genes <- genes[genes %in% rownames(liuM)]

  pearson <- cor(liuM[vin3, ], t(liuM[genes, ]))[1, ]
  selected <- names(pearson)[which(pearson > .9)]

  write_tsv(
    data.frame(Gene = selected, Pearson = unname(pearson[selected])),
    file.path(p, paste0("VIN3_correlated_", cl, ".tsv"))
  )

  pdf(file.path(p, paste0("VIN3_correlated_", cl, ".pdf")), width = 7, height = 4)

  par(mar = c(6, 4, 2, 1))

  boxplot(liu2z[selected, , drop = FALSE],
    outline = FALSE, ylim = c(-2, 2),
    ylab = "Z-score", main = cl, las = 2, cex.axis = .8
  )

  dev.off()
}

