#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p data/processed


fetch() {
  local series=$1 remote=$2 dest=$3

  if [[ -s "data/processed/$dest" ]]; then
    return
  fi

  curl --fail --location --retry 3 \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/${series}/suppl/${remote}" \
    -o "data/processed/${dest}.gz"

  gzip -dc "data/processed/${dest}.gz" > "data/processed/${dest}.tmp"
  mv "data/processed/${dest}.tmp" "data/processed/$dest"
  rm "data/processed/${dest}.gz"
}


fetch GSE283016 GSE283016_GeneRPKM.tsv.gz GeneRPKM.tsv
fetch GSE283016 GSE283016_GeneCounts.tsv.gz GeneCounts.tsv
fetch GSE283017 GSE283017_atac_peak_rpkm.tsv.gz atac_peak_rpkm.tsv
