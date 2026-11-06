#!/usr/bin/env bash
# Motif discovery, enrichment and TF-motif assignments.

set -euo pipefail

cd "$(dirname "$0")/.."


# Decompress the cluster and background sequences.

mkdir -p results/motifs/sequences
mkdir -p results/motifs/streme
mkdir -p results/motifs/sea
mkdir -p results/motifs/centrimo
mkdir -p results/motifs/known_tf

for name in C1 C2 C3 Bkg5; do

  gzip -dc "data/reference/sequences/${name}.fa.gz" \
    > "results/motifs/sequences/${name}.fa"

done


# Discover motifs in each accessibility cluster.

for cluster in C1 C2 C3; do

  streme \
    --oc "results/motifs/streme/${cluster}" \
    --align center \
    --p "results/motifs/sequences/${cluster}.fa" \
    --n results/motifs/sequences/Bkg5.fa

done


# Calculate motif enrichment in each cluster.

for cluster in C1 C2 C3; do

  sea \
    --oc "results/motifs/sea/${cluster}" \
    --align center \
    --m data/reference/discovered_motifs.meme \
    --p "results/motifs/sequences/${cluster}.fa" \
    --n results/motifs/sequences/Bkg5.fa

  centrimo \
    --optimize-score \
    --noseq \
    --neg results/motifs/sequences/Bkg5.fa \
    --oc "results/motifs/centrimo/${cluster}" \
    "results/motifs/sequences/${cluster}.fa" \
    data/reference/discovered_motifs.meme

done


# Enrich known TF motifs in each accessibility cluster.

for cluster in C1 C2 C3; do

  sea \
    --p "results/motifs/sequences/${cluster}.fa" \
    --m data/reference/JASPAR2024_Arabidopsis.meme \
    --n results/motifs/sequences/Bkg5.fa \
    --oc "results/motifs/known_tf/${cluster}" \
    --align center \
    --noseqs \
    --thresh 1 \
    --pvalue

done


# Match the motifs used in the networks to known Arabidopsis TF motifs.

tomtom \
  -oc results/motifs/tomtom \
  data/reference/motifs.meme \
  data/reference/PlantTFDB_Arabidopsis.meme
