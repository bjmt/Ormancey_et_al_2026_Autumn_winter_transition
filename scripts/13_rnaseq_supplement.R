source("scripts/lib/common.R")

p <- outdir("rnaseq/supplement")


# Calculate sample correlations using the 5,000 most variable genes.

expression <- rna_expression()
variable_genes <- order(matrixStats::rowVars(expression), decreasing = TRUE)[1:5000]

pcc <- cor(log2(1 + expression[variable_genes, ]))
labels <- gsub("_", "-", samples("rnaseq")$sample)

write_tsv(data.frame(Sample = rownames(pcc), pcc),
  file.path(p, "sample_correlations.tsv"))

pdf(file.path(p, "sample_correlations.pdf"), width = 9, height = 9)

heatmap(pcc,
  symm = TRUE,
  scale = "none",
  labRow = labels,
  labCol = labels,
  col = colorRampPalette(c("white", "#0070c0"))(100),
  margins = c(8, 8),
  cexRow = .6,
  cexCol = .6
)

dev.off()


# Count expression changes relative to the initial time point.

counts <- list()

for (condition in c("S", "A", "F")) {
  for (week in 3:6) {
    de <- read_tsv(
      ref(paste0("de/W0v", condition, week, ".tsv.gz")),
      row.names = 1
    )

    counts[[paste(condition, week)]] <- data.frame(
      Condition = unname(c(S = "22C", A = "AVG", F = "FLU")[condition]),
      Week = week - 2,
      Direction = c("Up", "Down"),
      Count = c(
        sum(de$logFC >= 1 & de$FDR < .01),
        sum(de$logFC <= -1 & de$FDR < .01)
      )
    )
  }
}

counts <- do.call(rbind, counts)
rownames(counts) <- NULL

write_tsv(counts, file.path(p, "initial_timepoint_comparisons.tsv"))

pdf(file.path(p, "initial_timepoint_comparisons.pdf"), width = 8, height = 4)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
colours <- c("22C" = "#afababff", "AVG" = "#ffc000ff", "FLU" = "#0070c0ff")

for (direction in c("Up", "Down")) {
  selected <- counts[counts$Direction == direction, ]
  values <- xtabs(Count ~ Condition + Week, selected)
  values <- values[names(colours), ]

  barplot(values,
    beside = TRUE,
    col = colours,
    names.arg = paste0("W", 1:4),
    xlab = "Week since start of cold",
    ylab = "Genes",
    main = paste0(direction, "-regulated"),
    legend.text = names(colours),
    args.legend = list(bty = "n", cex = .7)
  )
}

dev.off()


# Count expression differences between AVG and FLU at each week.

counts <- list()

for (week in 3:6) {
  de <- read_tsv(ref(paste0("de/AvF_Wk", week, ".tsv.gz")), row.names = 1)

  counts[[as.character(week)]] <- data.frame(
    Week = week - 2,
    Direction = c("Up", "Down"),
    Count = c(
      sum(de$logFC >= 1 & de$FDR < .01),
      sum(de$logFC <= -1 & de$FDR < .01)
    )
  )
}

counts <- do.call(rbind, counts)
rownames(counts) <- NULL

write_tsv(counts, file.path(p, "AVG_FLU_weekly_comparisons.tsv"))

pdf(file.path(p, "AVG_FLU_weekly_comparisons.pdf"), width = 5, height = 4)

barplot(xtabs(Count ~ Direction + Week, counts),
  beside = TRUE,
  col = c("#999999", "#009E73"),
  names.arg = paste0("W", 1:4),
  xlab = "Week since start of cold",
  ylab = "Genes",
  legend.text = c("Down", "Up"),
  args.legend = list(bty = "n")
)

dev.off()
