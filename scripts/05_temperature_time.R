library(variancePartition)
library(ggplot2)

source("scripts/lib/common.R")


# Sampling temperatures and week number for the 39 RNA-seq libraries.

m <- as.data.frame(rna_expression())
colnames(m) <- samples("rnaseq")$sample

p <- outdir("temperature_time")

gS <- c(1:3, 10:12, 19:21, 28:30, 37:39)
gA <- c(1:3, 4:6, 13:15, 22:24, 31:33)
gF <- c(1:3, 7:9, 16:18, 25:27, 34:36)

temps <- read_tsv(ref("temperatures.tsv"))

tempsA <- temps[temps$Condition == "Average", ]
tempsF <- temps[temps$Condition == "Fluctuating", ]

tempsA <- tempsA[rep(1:5, each = 3), ]
tempsF <- tempsF[rep(1:5, each = 3), ]
tempsC <- tempsA
tempsC$Condition <- "Constant"
tempsC$TempSampling <- 22
tempsC$TempMax <- 22
tempsC$TempMin <- 22
tempsC$TempAvg <- 22

tempsA <- as.data.frame(tempsA)
tempsF <- as.data.frame(tempsF)
tempsC <- as.data.frame(tempsC)

rownames(tempsC) <- colnames(m)[grepl("22C", colnames(m))]
tempsAA <- tempsA[-(1:3), ]
tempsFF <- tempsF[-(1:3), ]
rownames(tempsAA) <- colnames(m)[grepl("AVG", colnames(m))]
rownames(tempsFF) <- colnames(m)[grepl("FLU", colnames(m))]

tempsAll <- rbind(tempsC, tempsAA, tempsFF)



# Partition expression variance across temperature and time.

tempsAllm5 <- ~ TempSampling + TempMax + TempMin + Week

tempsAllv5 <- fitExtractVarPartModel(
  m[, rownames(tempsAll)], tempsAllm5, tempsAll
)

write_tsv(tempsAllv5, file.path(p, "AllGenes.tsv"), row.names = TRUE)

ggsave(file.path(p, "AllGenes.pdf"),
  plot = plotVarPart(tempsAllv5, label.angle = 45),
  width = 7, height = 5
)


# Reverse the order of the variables, then repeat the model for each gene class.

tempsAllm6 <- ~ Week + TempMin + TempMax + TempSampling

tempsAllv6 <- fitExtractVarPartModel(
  m[, rownames(tempsAll)], tempsAllm6, tempsAll
)

write_tsv(tempsAllv6, file.path(p, "AllGenes_rev.tsv"), row.names = TRUE)

ggsave(file.path(p, "AllGenes_rev.pdf"),
  plot = plotVarPart(tempsAllv6, label.angle = 65),
  width = 7, height = 5
)

cls <- read_tsv(ref("rna_clusters.tsv"))

for (cl in c("early", "linear", "late")) {
  genes <- cls$Genes[grepl(paste0("cold_", cl), cls$Class)]

  v <- fitExtractVarPartModel(
    m[genes, rownames(tempsAll)], tempsAllm5, tempsAll
  )

  write_tsv(v, file.path(p, paste0(cl, ".tsv")), row.names = TRUE)

  ggsave(file.path(p, paste0(cl, ".pdf")),
    plot = plotVarPart(v, label.angle = 65),
    width = 7, height = 5
  )
}

