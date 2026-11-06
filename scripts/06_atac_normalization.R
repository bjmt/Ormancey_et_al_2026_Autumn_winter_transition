library(edgeR)
library(RUVSeq)
library(rtracklayer)
library(matrixStats)

source("scripts/lib/common.R")


# Read peak counts; an optional argument supplies a counts file.

args <- commandArgs(trailingOnly = TRUE)
count_file <- if (length(args)) args[1] else ref("atac_counts.tsv.gz")

counts <- read_matrix(count_file)

ACRs <- import(ref("atac_summits.bed"))
ACRsW <- import(ref("atac_peaks.bed"))
BL <- import(ref("blacklist.bed"))

ACRsWF <- ACRsW[!overlapsAny(ACRsW, BL)]
ACRsF <- ACRs[ACRs$name %in% ACRsWF$name]
countsW <- counts[ACRsF$name, ]

grps <- gsub("_\\d$", "", colnames(counts))

p <- outdir("atac")


# Estimate dispersion and correct batch effects.

countsW2 <- countsW[, !grepl("Avg_", colnames(countsW))]
grps2 <- grps[!grepl("Avg_", grps)]

gdge <- DGEList(countsW2, group = grps2)
gdge <- gdge[rowMedians(gdge$counts) >= 120, , keep.lib.sizes = FALSE]
gdge <- calcNormFactors(gdge, method = "upperquartile")
gdge <- estimateDisp(gdge, robust = TRUE)

peak_lengths <- structure(width(ACRsWF), names = ACRsF$name)

gdge$genes <- data.frame(
  Length = peak_lengths[rownames(gdge$counts)]
)

rp1 <- rpkm(gdge)


# Fit condition-specific time trends and obtain deviance residuals.

grps_con1 <- ifelse(grepl("Con_", grps2), 1, 0)
grps_flu1 <- ifelse(grepl("Flu_", grps2), 1, 0)
grps_allt <- rep(1:4, each = 4)

grps_cont <- grps_allt
grps_flut <- grps_allt

grps_cont[grps_con1 != 1] <- 0
grps_flut[grps_flu1 != 1] <- 0

grps_cont2 <- grps_cont
grps_flut2 <- grps_flut

grps_cont2[grps_cont2 == 0] <- 1
grps_flut2[grps_flut2 == 0] <- 1

gFit <- glmFit(gdge, model.matrix(~grps_flu1 + grps_cont2 + grps_flut2))
res <- residuals(gFit, type = "deviance")


# Correct batch effects, then fit the sample groups.

ruv_g <- 1:nrow(gdge$counts)
nrgc <- RUVr(gdge$counts, ruv_g, 2, res)

gdge2 <- DGEList(nrgc$normalizedCounts, group = grps2)
gdge2 <- calcNormFactors(gdge2)
gdge2 <- estimateDisp(gdge2)

gdge2$genes <- data.frame(
  Length = peak_lengths[rownames(gdge2$counts)]
)

fit2 <- glmQLFit(gdge2, model.matrix(~grps2))


# Accessibility changes at weeks 2 and 4 of cold.

qlf2 <- glmQLFTest(fit2,
  contrast = makeContrasts(
    "grps2Flu_6W-grps2Con_6W",
    levels = model.matrix(~grps2)
  )
)

cVf6 <- as.data.frame(topTags(qlf2, n = Inf, p.value = 0.05))

qlf3 <- glmQLFTest(fit2,
  contrast = makeContrasts(
    "grps2Flu_4W-grps2Con_4W",
    levels = model.matrix(~grps2)
  )
)

cVf4 <- as.data.frame(topTags(qlf3, n = Inf, p.value = 0.05))


rp2 <- rpkm(gdge2)
cp2 <- cpm(gdge2)

write.table(rp2, file.path(p, "atac_peak_counts.tsv"),
  quote = FALSE, sep = "\t")

write.table(cVf6, file.path(p, "FLU_vs_22C_week4.tsv"),
  quote = FALSE, sep = "\t")

write.table(cVf4, file.path(p, "FLU_vs_22C_week2.tsv"),
  quote = FALSE, sep = "\t")

