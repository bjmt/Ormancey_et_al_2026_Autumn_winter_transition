source("scripts/lib/common.R")

p <- outdir("gene_ontology/figures")


# Plot GO terms for genes with higher expression in AVG or FLU.

for (condition in c("AVG", "FLU")) {
  terms <- read_tsv(ref(paste0("go/", condition, ".tsv")))
  terms <- terms[terms$Category == "GOTERM_BP_DIRECT", ]
  terms <- terms[order(terms$PValue), ][1:5, ]

  write_tsv(terms, file.path(p, paste0(condition, "_GO.tsv")))

  labels <- vapply(terms$Term, function(term) {
    paste(strwrap(term, width = 28), collapse = "\n")
  }, character(1))

  colour <- if (condition == "AVG") "#E69F00" else "#0070c0"

  pdf(file.path(p, paste0(condition, "_GO.pdf")), width = 6.5, height = 4)

  par(mar = c(4, 13, 3, 1), mgp = c(2.5, .5, 0))

  barplot(rev(-log10(terms$PValue)),
    names.arg = rev(labels),
    horiz = TRUE,
    las = 1,
    col = colour,
    border = NA,
    cex.names = .7,
    xlab = "-log10(P-value)",
    main = paste("Higher in", condition)
  )

  dev.off()
}


# Plot the top 20 enriched biological processes for the accessibility clusters.

for (cluster in c("C1", "C2")) {
  terms <- read.csv(ref(paste0("go/", cluster, ".csv")), check.names = FALSE)
  terms <- terms[terms$source == "GO:BP", ][1:20, ]

  write_tsv(terms, file.path(p, paste0(cluster, "_GO.tsv")))

  labels <- vapply(terms$term_name, function(term) {
    paste(strwrap(term, width = 43), collapse = "\n")
  }, character(1))

  pdf(file.path(p, paste0(cluster, "_GO.pdf")), width = 9, height = 8)

  par(mar = c(4, 19, 3, 1), mgp = c(2.5, .5, 0))

  barplot(rev(terms$negative_log10_of_adjusted_p_value),
    names.arg = rev(labels),
    horiz = TRUE,
    las = 1,
    col = "black",
    border = NA,
    cex.names = .65,
    xlab = "-log10(Adjusted P-value)",
    main = sub("C", "A", cluster)
  )

  dev.off()
}
