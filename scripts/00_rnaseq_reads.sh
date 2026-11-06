#!/usr/bin/env bash
# Adapter trimming, STAR alignment, StringTie quantification and prepDE counts.

set -euo pipefail

cd "$(dirname "$0")/.."


# Genome files and software settings.

: "${GENOME_FASTA:?Set GENOME_FASTA to the TAIR10 FASTA (uncompressed)}"
: "${ARAPORT_GFF:?Set ARAPORT_GFF to the Araport11 GFF3 (uncompressed)}"
: "${PREPDE:?Set PREPDE to the StringTie prepDE.py3 script}"
threads=${THREADS:-6}
star_index=${STAR_INDEX:-data/genome/STAR-tair10}


mkdir -p "$star_index"
mkdir -p results/rnaseq/trimmed
mkdir -p results/rnaseq/aligned
mkdir -p results/rnaseq/stringtie
mkdir -p results/rnaseq/logs


# Build the STAR index if needed.

if [[ ! -s "$star_index/Genome" ]]; then
  STAR \
    --runMode genomeGenerate \
    --runThreadN "$threads" \
    --genomeDir "$star_index" \
    --genomeFastaFiles "$GENOME_FASTA"
fi


# Trim, align and quantify each paired-end library.

printf '' > results/rnaseq/gtfs.txt

while IFS=$'\t' read -r sample library _analysis _condition _week _rep \
    r1 r2 _submitted1 _submitted2; do

  if [[ "$sample" == sample ]]; then
    continue
  fi

  # Adapter trimming.

  cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -j "$threads" \
    -m 15 \
    -o "results/rnaseq/trimmed/${library}_1.fq.gz" \
    -p "results/rnaseq/trimmed/${library}_2.fq.gz" \
    "data/raw/rnaseq/$r1" "data/raw/rnaseq/$r2" \
    > "results/rnaseq/logs/${library}.cutadapt.txt"


  # Align to TAIR10.

  STAR \
    --genomeDir "$star_index" \
    --runThreadN "$threads" \
    --readFilesCommand "gunzip -c" \
    --readFilesIn "results/rnaseq/trimmed/${library}_1.fq.gz" \
      "results/rnaseq/trimmed/${library}_2.fq.gz" \
    --outFileNamePrefix "results/rnaseq/aligned/${library}." \
    --outSAMstrandField intronMotif \
    --outMultimapperOrder Random \
    --outSAMmultNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 10000 \
    --limitOutSAMoneReadBytes 10000000


  # Quantify the annotated transcripts.

  stringtie \
    "results/rnaseq/aligned/${library}.Aligned.sortedByCoord.out.bam" \
    -G "$ARAPORT_GFF" \
    -e --rf \
    -p "$threads" \
    -A "results/rnaseq/stringtie/${library}.gene_abund.txt" \
    -o "results/rnaseq/stringtie/${library}.gtf"

  printf '%s\t%s\n' "$library" "results/rnaseq/stringtie/${library}.gtf" \
    >> results/rnaseq/gtfs.txt

done < metadata/rnaseq_samples.tsv


# Convert the transcript quantifications to count matrices.

python3 "$PREPDE" \
  -l 151 \
  -i results/rnaseq/gtfs.txt \
  -g results/rnaseq/gene_count_matrix.csv \
  -t results/rnaseq/transcript_count_matrix.csv
