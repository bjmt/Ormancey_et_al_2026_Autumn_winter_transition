source("scripts/lib/common.R")

p <- outdir("rnaseq/patterns")


# Classify expression patterns using Pearson correlations and model statistics.
# Use --refit to start again from the RNA-seq counts.

args <- commandArgs(trailingOnly = TRUE)

tags_all2 <- read_tsv(ref("rna_pattern_statistics.tsv.gz"))
gene2sym <- structure(tags_all2$Symbol, names = tags_all2$GeneID)

if ("--refit" %in% args) {
  library(edgeR)
  library(RUVSeq)

  gcm <- as.matrix(read.csv(ref("rna_counts.csv.gz"),
    row.names = 1, check.names = FALSE))

  colnames(gcm) <- substr(colnames(gcm), 4, 7)

  grp <- substr(colnames(gcm), 1, 3)

  df <- data.frame(row.names = NULL,
    Sample = colnames(gcm),
    Replicate = rep(1:3, 13),
    Group = grp,
    Time = substr(grp, 1, 1),
    Condition = substr(grp, 3, 3)
  )

  gdge <- DGEList(gcm, group = grp)
  gdge <- gdge[filterByExpr(gdge, group = grp), , keep.lib.sizes = FALSE]
  gdge <- calcNormFactors(gdge)
  gdge <- estimateDisp(gdge)
  gFit <- glmFit(gdge, model.matrix(~Group, data = df))
  res <- residuals(gFit, type = "deviance")
  nrgc <- RUVr(gdge$counts, 1:nrow(gdge$counts), 2, res)
  gdge2 <- DGEList(nrgc$normalizedCounts, group = grp)
  gdge2 <- calcNormFactors(gdge2)

  gcm2 <- gdge2$counts

  cpms <- cpm(gdge2)


  dge <- DGEList(gcm2, group = grp)
  dge <- calcNormFactors(dge)

  nf <- dge$samples


  splines <- poly(1:5, degree = 4)

  patterns <- list(
    linear = splines[, 1],
    transient = -splines[, 2],
    early = -c(splines[, 2][1:3], rep(splines[, 2][3], 2))
  )
  patterns$late <- c(rep(patterns$early[1], 3), patterns$early[2:3])
  patterns2 <- do.call(cbind, lapply(patterns, function(x) rep(x, each = 3)))
  patMod <- model.matrix(~linear+transient+early+late, data = as.data.frame(patterns2))

  gcm_s <- gcm2[, df$Condition == "S" | df$Time == 2]
  gcm_a <- gcm2[, df$Condition == "A" | df$Time == 2]
  gcm_f <- gcm2[, df$Condition == "F" | df$Time == 2]

  dge_f <- DGEList(gcm_f, group = paste0("W", rep(2:6, each = 3)))
  dge_f$samples$norm.factors <- nf$norm.factors[df$Condition == "F" | df$Time == 2]
  dge_f <- estimateDisp(dge_f, patMod)
  fit_f <- glmQLFit(dge_f, patMod, robust = TRUE)
  test_f <- glmQLFTest(fit_f, coef = 2:5)
  tags_f <- as.data.frame(topTags(test_f, Inf))

  dge_s <- DGEList(gcm_s, group = paste0("W", rep(2:6, each = 3)))
  dge_s$samples$norm.factors <- nf$norm.factors[df$Condition == "S" | df$Time == 2]
  dge_s <- estimateDisp(dge_s, patMod)
  fit_s <- glmQLFit(dge_s, patMod, robust = TRUE)
  test_s <- glmQLFTest(fit_s, coef = 2:5)
  tags_s <- as.data.frame(topTags(test_s, Inf))

  dge_a <- DGEList(gcm_a, group = paste0("W", rep(2:6, each = 3)))
  dge_a$samples$norm.factors <- nf$norm.factors[df$Condition == "A" | df$Time == 2]
  dge_a <- estimateDisp(dge_a, patMod)
  fit_a <- glmQLFit(dge_a, patMod, robust = TRUE)
  test_a <- glmQLFTest(fit_a, coef = 2:5)
  tags_a <- as.data.frame(topTags(test_a, Inf))


  get_cor <- function(g, w) {
    g <- cpms[g, df$Condition == w | df$Time == 2]
    cor(g, patterns2)[1, ]^2
  }

  get_cor2 <- function(g, w) {
    g <- cpms[g, df$Condition == w | df$Time == 2]
    cor(g, patterns2)[1, ]
  }

  v_get_cor <- function(g, w, p) vapply(g, function(x) get_cor(x, w)[p], numeric(1))
  tags_all <- data.frame(row.names = NULL,

    GeneID = rownames(tags_f),
    Symbol = gene2sym[rownames(tags_f)],

    Fluctuating_linear_r2 = v_get_cor(rownames(tags_f), "F", "linear"),
    Fluctuating_transient_r2 = v_get_cor(rownames(tags_f), "F", "transient"),
    Fluctuating_early_r2 = v_get_cor(rownames(tags_f), "F", "early"),
    Fluctuating_late_r2 = v_get_cor(rownames(tags_f), "F", "late"),
    Fluctuating_pvalue = tags_f$PValue,
    Fluctuating_fdr = tags_f$FDR,

    Average_linear_r2 = v_get_cor(rownames(tags_f), "A", "linear"),
    Average_transient_r2 = v_get_cor(rownames(tags_f), "A", "transient"),
    Average_early_r2 = v_get_cor(rownames(tags_f), "A", "early"),
    Average_late_r2 = v_get_cor(rownames(tags_f), "A", "late"),
    Average_pvalue = tags_a[rownames(tags_f), ]$PValue,
    Average_fdr = tags_a[rownames(tags_f), ]$FDR,

    Constant_linear_r2 = v_get_cor(rownames(tags_f), "S", "linear"),
    Constant_transient_r2 = v_get_cor(rownames(tags_f), "S", "transient"),
    Constant_early_r2 = v_get_cor(rownames(tags_f), "S", "early"),
    Constant_late_r2 = v_get_cor(rownames(tags_f), "S", "late"),
    Constant_pvalue = tags_s[rownames(tags_f), ]$PValue,
    Constant_fdr = tags_s[rownames(tags_f), ]$FDR

  )

  v_get_cor2 <- function(g, w, p) vapply(g, function(x) get_cor2(x, w)[p], numeric(1))
  tags_all2 <- data.frame(row.names = NULL,

    GeneID = rownames(tags_f),
    Symbol = gene2sym[rownames(tags_f)],

    Fluctuating_linear_r = v_get_cor2(rownames(tags_f), "F", "linear"),
    Fluctuating_transient_r = v_get_cor2(rownames(tags_f), "F", "transient"),
    Fluctuating_early_r = v_get_cor2(rownames(tags_f), "F", "early"),
    Fluctuating_late_r = v_get_cor2(rownames(tags_f), "F", "late"),
    Fluctuating_pvalue = tags_f$PValue,
    Fluctuating_fdr = tags_f$FDR,

    Average_linear_r = v_get_cor2(rownames(tags_f), "A", "linear"),
    Average_transient_r = v_get_cor2(rownames(tags_f), "A", "transient"),
    Average_early_r = v_get_cor2(rownames(tags_f), "A", "early"),
    Average_late_r = v_get_cor2(rownames(tags_f), "A", "late"),
    Average_pvalue = tags_a[rownames(tags_f), ]$PValue,
    Average_fdr = tags_a[rownames(tags_f), ]$FDR,

    Constant_linear_r = v_get_cor2(rownames(tags_f), "S", "linear"),
    Constant_transient_r = v_get_cor2(rownames(tags_f), "S", "transient"),
    Constant_early_r = v_get_cor2(rownames(tags_f), "S", "early"),
    Constant_late_r = v_get_cor2(rownames(tags_f), "S", "late"),
    Constant_pvalue = tags_s[rownames(tags_f), ]$PValue,
    Constant_fdr = tags_s[rownames(tags_f), ]$FDR

  )


} else {
  tags_all <- tags_all2

  r_cols <- grep("_r$", names(tags_all))

  for (i in r_cols) {
    tags_all[[i]] <- tags_all[[i]]^2
  }

  names(tags_all)[r_cols] <- paste0(names(tags_all)[r_cols], "2")
}


# Classify each condition and annotate the direction of change.

for (i in grep("_r2$", names(tags_all))) {
  tags_all[[i]][is.na(tags_all[[i]])] <- 0
}

for (i in grep("_r$", names(tags_all2))) {
  tags_all2[[i]][is.na(tags_all2[[i]])] <- 0
}

th <- .8

tags_class <- with(tags_all, data.frame(row.names = NULL,
  Gene = GeneID,
  Symbol = Symbol,
  Fluctuating = "",
  Average = "",
  Constant = ""
))

for (w in c("Average", "Fluctuating", "Constant")) {
  tags_class[[w]][tags_all[[paste0(w, "_linear_r2")]] >= th] <- paste0(
    tags_class[[w]][tags_all[[paste0(w, "_linear_r2")]] >= th], "C")
  tags_class[[w]][tags_all[[paste0(w, "_transient_r2")]] >= th & tags_class[[w]] != "C"] <- paste0(
    tags_class[[w]][tags_all[[paste0(w, "_transient_r2")]] >= th & tags_class[[w]] != "C"], "T")
  tags_class[[w]][tags_all[[paste0(w, "_early_r2")]] >= th & tags_class[[w]] != "C"] <- paste0(
    tags_class[[w]][tags_all[[paste0(w, "_early_r2")]] >= th & tags_class[[w]] != "C"], "E")
  tags_class[[w]][tags_all[[paste0(w, "_late_r2")]] >= th & tags_class[[w]] != "C"] <- paste0(
    tags_class[[w]][tags_all[[paste0(w, "_late_r2")]] >= th & tags_class[[w]] != "C"], "L")
  tags_class[[w]][tags_class[[w]] == ""] <- "N"

  message(w)
  print(table(tags_class[[w]]))
}

tags_class$Symbol <- tags_all$Symbol


tags_class$Fluctuating_r2 <- mapply(function(x, y) {
  switch(x,
    C = tags_all$Fluctuating_linear_r2[tags_all$GeneID == y],
    T = tags_all$Fluctuating_transient_r2[tags_all$GeneID == y],
    E = tags_all$Fluctuating_early_r2[tags_all$GeneID == y],
    L = tags_all$Fluctuating_late_r2[tags_all$GeneID == y],
    0
  )
}, tags_class$Fluctuating, tags_class$Gene)

tags_class$Average_r2 <- mapply(function(x, y) {
  switch(x,
    C = tags_all$Average_linear_r2[tags_all$GeneID == y],
    T = tags_all$Average_transient_r2[tags_all$GeneID == y],
    E = tags_all$Average_early_r2[tags_all$GeneID == y],
    L = tags_all$Average_late_r2[tags_all$GeneID == y],
    0
  )
}, tags_class$Average, tags_class$Gene)

tags_class$Constant_r2 <- mapply(function(x, y) {
  switch(x,
    C = tags_all$Constant_linear_r2[tags_all$GeneID == y],
    T = tags_all$Constant_transient_r2[tags_all$GeneID == y],
    E = tags_all$Constant_early_r2[tags_all$GeneID == y],
    L = tags_all$Constant_late_r2[tags_all$GeneID == y],
    0
  )
}, tags_class$Constant, tags_class$Gene)

tags_class2 <- tags_class

for (i in 3:5) {
  labels <- c(
    C = "Linear", T = "Transient", E = "Early", L = "Late", N = "Stable"
  )

  tags_class2[[i]] <- labels[tags_class[[i]]]
}

write_tsv(tags_class2, file.path(p, "Pattern_classifications.tsv"))


tags_class3 <- tags_class2
tags_class3$Fluctuating_expr <- NA_character_
tags_class3$Average_expr <- NA_character_
tags_class3$Constant_expr <- NA_character_

for (i in 1:nrow(tags_class3)) {
  tags_class3$Fluctuating_expr[i] <- switch(tags_class3$Fluctuating[i],
    Linear = if (tags_all2$Fluctuating_linear_r[i] > 0) "Increasing" else "Decreasing",
    Early = if (tags_all2$Fluctuating_early_r[i] > 0) "Increasing" else "Decreasing",
    Late = if (tags_all2$Fluctuating_late_r[i] > 0) "Increasing" else "Decreasing",
    Transient = if (tags_all2$Fluctuating_transient_r[i] > 0) "Increasing" else "Decreasing",
    Stable = "Neutral")
  tags_class3$Average_expr[i] <- switch(tags_class3$Average[i],
    Linear = if (tags_all2$Average_linear_r[i] > 0) "Increasing" else "Decreasing",
    Early = if (tags_all2$Average_early_r[i] > 0) "Increasing" else "Decreasing",
    Late = if (tags_all2$Average_late_r[i] > 0) "Increasing" else "Decreasing",
    Transient = if (tags_all2$Average_transient_r[i] > 0) "Increasing" else "Decreasing",
    Stable = "Neutral")
  tags_class3$Constant_expr[i] <- switch(tags_class3$Constant[i],
    Linear = if (tags_all2$Constant_linear_r[i] > 0) "Increasing" else "Decreasing",
    Early = if (tags_all2$Constant_early_r[i] > 0) "Increasing" else "Decreasing",
    Late = if (tags_all2$Constant_late_r[i] > 0) "Increasing" else "Decreasing",
    Transient = if (tags_all2$Constant_transient_r[i] > 0) "Increasing" else "Decreasing",
    Stable = "Neutral")
}

tags_class3$Fluctuating_p <- tags_all2$Fluctuating_pvalue
tags_class3$Average_p <- tags_all2$Average_pvalue
tags_class3$Constant_p <- tags_all2$Constant_pvalue

tags_class3$Fluctuating_r2[tags_class3$Fluctuating_p >= 0.01] <- 0
tags_class3$Fluctuating_expr[tags_class3$Fluctuating_p >= 0.01] <- "Neutral"
tags_class3$Fluctuating[tags_class3$Fluctuating_p >= 0.01] <- "Stable"

tags_class3$Average_r2[tags_class3$Average_p >= 0.01] <- 0
tags_class3$Average_expr[tags_class3$Average_p >= 0.01] <- "Neutral"
tags_class3$Average[tags_class3$Average_p >= 0.01] <- "Stable"

tags_class3$Constant_r2[tags_class3$Constant_p >= 0.01] <- 0
tags_class3$Constant_expr[tags_class3$Constant_p >= 0.01] <- "Neutral"
tags_class3$Constant[tags_class3$Constant_p >= 0.01] <- "Stable"


write_tsv(tags_class3, file.path(p, "Pattern_classifications_v2.tsv"))


# Select cold-responsive genes with stable expression in the constant condition.

classes <- c("Early", "Linear", "Late")
directions <- c("Increasing", "Decreasing")

cls_d <- list()

for (direction in directions) {
  for (cl in classes) {
    keep <- (
      (tags_class3$Fluctuating == cl | tags_class3$Average == cl) &
      (tags_class3$Fluctuating_expr == direction | tags_class3$Average_expr == direction) &
      tags_class3$Constant_expr == "Neutral"
    )

    keep[is.na(keep)] <- FALSE

    label <- paste(
      if (direction == "Increasing") "Inc" else "Dec",
      "cold", tolower(cl), sep = "_"
    )

    cls_d[[label]] <- data.frame(
      Genes = tags_class3$Gene[keep],
      Class = label
    )
  }
}

cls_d <- do.call(rbind, cls_d)

expression_genes <- read_tsv(ref("gene_lengths.tsv"))$Gene
cls_d <- cls_d[cls_d$Genes %in% expression_genes, ]

cls_d <- cls_d[order(cls_d$Class), ]
cls_d <- cls_d[!duplicated(cls_d$Genes), ]

rownames(cls_d) <- NULL

write_tsv(cls_d, file.path(p, "cold_responsive_clusters.tsv"))
