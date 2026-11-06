library(edgeR)

source("scripts/lib/common.R")


# Differential expression against matched controls and the initial time point.

args <- commandArgs(trailingOnly = TRUE)
count_file <- if (length(args)) args[1] else ref("rna_counts.csv.gz")

gcm <- as.matrix(read.csv(count_file, row.names = 1, check.names = FALSE))
s <- samples("rnaseq")
gcm <- gcm[, s$library]
colnames(gcm) <- s$analysis

p <- outdir("rnaseq/de")

grps <- gsub("\\d$", "", colnames(gcm))

df <- data.frame(row.names = NULL,
  Sample = colnames(gcm),
  Replicate = rep(1:3, 13),
  Group =  grps,
  Time = substr(grps, 1, 1),
  Condition = factor(substr(grps, 3, 3))
)
design <- model.matrix(~0 + Condition + Time, data = df)
design2 <- model.matrix(~0+Group, data = df)

y <- DGEList(gcm, group = grps)
y <- y[filterByExpr(y, groups = grps), , keep = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)

fit <- glmQLFit(y, design)
fit2 <- glmQLFit(y, design2)

Avg_vs_fluc <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0, 0))
Con_vs_fluc <- glmQLFTest(fit, contrast = c(0, 1, -1, 0, 0, 0, 0))
Con_vs_Avg <- glmQLFTest(fit, contrast = c(1, 0, -1, 0, 0, 0, 0))
summary(decideTests(Avg_vs_fluc))
summary(decideTests(Con_vs_fluc))
summary(decideTests(Con_vs_Avg))

write.table(topTags(Avg_vs_fluc, n = Inf, p.value = 0.05),
  "results/rnaseq/de/Avg_vs_fluc.tsv", sep = "\t", quote = FALSE)
write.table(topTags(Con_vs_fluc, n = Inf, p.value = 0.05),
  "results/rnaseq/de/Con_vs_fluc.tsv", sep = "\t", quote = FALSE)
write.table(topTags(Con_vs_Avg, n = Inf, p.value = 0.05),
  "results/rnaseq/de/Con_vs_avg.tsv", sep = "\t", quote = FALSE)

cons <- makeContrasts(levels = design2,
  AvF_Wk3 = Group3WF - Group3WA,
  AvF_Wk4 = Group4WF - Group4WA,
  AvF_Wk5 = Group5WF - Group5WA,
  AvF_Wk6 = Group6WF - Group6WA,
  SvF_Wk3 = Group3WF - Group3WS,
  SvF_Wk4 = Group4WF - Group4WS,
  SvF_Wk5 = Group5WF - Group5WS,
  SvF_Wk6 = Group6WF - Group6WS,
  SvA_Wk3 = Group3WA - Group3WS,
  SvA_Wk4 = Group4WA - Group4WS,
  SvA_Wk5 = Group5WA - Group5WS,
  SvA_Wk6 = Group6WA - Group6WS,
  W0vF3 = Group2WS - Group3WF,
  W0vF4 = Group2WS - Group4WF,
  W0vF5 = Group2WS - Group5WF,
  W0vF6 = Group2WS - Group6WF,
  W0vA3 = Group2WS - Group3WA,
  W0vA4 = Group2WS - Group4WA,
  W0vA5 = Group2WS - Group5WA,
  W0vA6 = Group2WS - Group6WA,
  W0vS3 = Group2WS - Group3WS,
  W0vS4 = Group2WS - Group4WS,
  W0vS5 = Group2WS - Group5WS,
  W0vS6 = Group2WS - Group6WS
)

all_comps <- lapply(1:ncol(cons), function(x) glmQLFTest(fit2, contrast = cons[, x]))
de_mat <- matrix(0, nrow = nrow(y), ncol = ncol(cons),
  dimnames = list(rownames(y), colnames(cons)))
for (i in seq_along(all_comps)) {
  g <- topTags(all_comps[[i]], n = Inf, p.value = 0.05)
  de_mat[rownames(g), i] <- as.data.frame(g)$logFC
  write.table(g, paste0("results/rnaseq/de/", colnames(cons)[i], ".tsv"), sep = "\t", quote = FALSE)
}


# Overall AVG versus FLU comparison.

gcm <- read_matrix(ref("rna_counts_noRRNA.tsv.gz"))
colnames(gcm) <- gsub("JQ_", "", colnames(gcm), fixed = TRUE)
grps <- gsub("\\d$", "", colnames(gcm))


df <- data.frame(row.names = NULL,
  Sample = colnames(gcm),
  Replicate = rep(1:3, 13),
  Group =  grps,
  Week = substr(grps, 1, 1),
  Condition = factor(substr(grps, 3, 3))
)

design <- model.matrix(~0 + Condition + Week, data = df)

y <- DGEList(gcm, group = grps)
y <- y[filterByExpr(y, groups = grps), , keep = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)

gFit <- glmFit(y, model.matrix(~Group, data = df))
res <- residuals(gFit, type = "deviance")
nrgc <- RUVr(y$counts, 1:nrow(y$counts), 2, res)
y <- DGEList(nrgc$normalizedCounts, group = grps)

y <- calcNormFactors(y)
y <- estimateDisp(y)

fit <- glmQLFit(y, design)

Avg_vs_fluc <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0, 0))
summary(decideTests(Avg_vs_fluc))

Avg_vs_fluc_de <- as.data.frame(topTags(Avg_vs_fluc, n = Inf, p.value = 0.05))

write.table(Avg_vs_fluc_de, file.path(p, "AVG_vs_FLU.tsv"),
  quote = FALSE, sep = "\t")

