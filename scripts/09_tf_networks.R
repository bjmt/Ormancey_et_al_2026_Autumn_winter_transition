library(universalmotif)
library(Biostrings)
library(dplyr)
library(readr)
library(igraph)

source("scripts/lib/common.R")


# Scan motifs in the cluster sequences.

gtab <- rna_expression()
prom2peak <- read_tsv(ref("prom2peak.tsv"))
fmots <- read_meme(ref("motifs.meme"))

c1s <- readDNAStringSet(ref("sequences/C1.fa.gz"))
c2s <- readDNAStringSet(ref("sequences/C2.fa.gz"))

sscan1 <- scan_sequences(fmots, c1s, RC = TRUE, no.overlaps = TRUE)
sscan2 <- scan_sequences(fmots, c2s, RC = TRUE, no.overlaps = TRUE)

p <- outdir("networks")

write_tsv(as.data.frame(sscan1), file.path(p, "A1_motif_hits.tsv"))
write_tsv(as.data.frame(sscan2), file.path(p, "A2_motif_hits.tsv"))


# Assign TFs to motifs and retain positively correlated TF-target pairs.

tom <- read_tsv(ref("tomtom.tsv"), comment.char = "#")

tom2 <- tom[
  tom$`q-value` < 0.05 & tom$Target_ID %in% rownames(gtab),
]

tom3 <- as.data.frame(tom2[, 1:2])

expr <- log2(1 + gtab)

gcor <- function(tf, gene) {
  cor(expr[tf, ], expr[gene, ])
}

prom2peak2 <- prom2peak[prom2peak$name %in% c(sscan1$sequence, sscan2$sequence), ]

motif_peaks <- function(hits) {
  result <- list()

  for (motif in paste0("M", 1:7)) {
    peaks <- unique(hits$sequence[hits$motif == motif])

    result[[motif]] <- data.frame(
      Motif = rep(motif, length(peaks)),
      Peak = peaks
    )
  }

  do.call(rbind, result)
}


tf2g1 <- motif_peaks(sscan1)
tf2g2 <- motif_peaks(sscan2)


tf2g1 <- dplyr::left_join(tf2g1, tom3, by = c("Motif" = "Query_ID"),
  relationship = "many-to-many")

tf2g2 <- dplyr::left_join(tf2g2, tom3, by = c("Motif" = "Query_ID"),
  relationship = "many-to-many")

colnames(tf2g1)[3] <- "TF"
colnames(tf2g2)[3] <- "TF"

tf2g1 <- dplyr::left_join(tf2g1, prom2peak2[, c("name", "Gene")], by = c("Peak" = "name"))
tf2g2 <- dplyr::left_join(tf2g2, prom2peak2[, c("name", "Gene")], by = c("Peak" = "name"))

tf2g1 <- tf2g1[!is.na(tf2g1$Gene), ]
tf2g2 <- tf2g2[!is.na(tf2g2$Gene), ]

tf2g1 <- tf2g1[tf2g1$Gene %in% rownames(gtab), ]
tf2g2 <- tf2g2[tf2g2$Gene %in% rownames(gtab), ]

tf2g1$PCC <- mapply(gcor, tf2g1$TF, tf2g1$Gene)
tf2g2$PCC <- mapply(gcor, tf2g2$TF, tf2g2$Gene)

tf2g1f <- tf2g1[tf2g1$PCC >= .5, ]
tf2g2f <- tf2g2[tf2g2$PCC >= .5, ]



# Plot the two networks and export their membership tables.

colsM <- structure(colorBlindBlack8[-1], names = paste0("M", 1:7))

networks <- list(
  COOLING = tf2g1f,
  CHILLING = tf2g2f
)

for (name in names(networks)) {
  membership <- networks[[name]]

  if (name == "COOLING") {
    g <- graph.data.frame(membership[, c("TF", "Gene", "Motif")])
  } else {
    g <- graph.data.frame(membership[, c("TF", "Gene")])
  }

  g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)

  tf2M <- structure(membership$Motif, names = membership$TF)

  V(g)$color <- colsM[tf2M[names(V(g))]]
  V(g)$color[is.na(V(g)$color)] <- "black"

  pdf(file.path(p, paste0(name, "_network.pdf")), width = 9, height = 7)

  plot(g,
    layout = layout_with_fr(g) * .03,
    rescale = FALSE,
    edge.arrow.size = 0,
    edge.width = .1,
    edge.color = "gray25",
    vertex.frame.color = NA,
    vertex.label = NA,
    vertex.size = log10(degree(g, mode = "all")) * 2
  )

  dev.off()

  edges <- as_data_frame(g)[, 1:2]
  colnames(edges) <- c("TF", "Target")

  keep <- paste0(membership$TF, membership$Gene) %in%
    paste0(edges$TF, edges$Target)

  write_tsv(membership[keep, ], file.path(p, paste0(name, "_network.tsv")))
}

