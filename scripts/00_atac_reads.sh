#!/usr/bin/env bash
# Process the 16 ATAC-seq libraries and count reads in peaks.

set -euo pipefail

cd "$(dirname "$0")/.."


# Genome files and software settings.

: "${GENOME_FASTA:?Set GENOME_FASTA to the TAIR10 FASTA (uncompressed)}"
threads=${THREADS:-6}
bt_index=${BOWTIE2_INDEX:-data/genome/bowtie2-tair10}


mkdir -p "$(dirname "$bt_index")"
mkdir -p results/atac/trimmed
mkdir -p results/atac/bams
mkdir -p results/atac/beds
mkdir -p results/atac/bigWigs
mkdir -p results/atac/peaks
mkdir -p results/atac/counts
mkdir -p results/atac/logs


# Build the Bowtie 2 index if needed.

if [[ ! -s "${bt_index}.1.bt2" && ! -s "${bt_index}.1.bt2l" ]]; then
  bowtie2-build "$GENOME_FASTA" "$bt_index"
fi


# Process each paired-end library.

while IFS=$'\t' read -r sample library analysis _condition _week _rep \
    r1 r2 _submitted1 _submitted2; do

  if [[ "$sample" == sample ]]; then
    continue
  fi

  # Adapter trimming.

  fastp \
    -w "$threads" \
    -a CTGTCTCTTATACACATCT \
    -i "data/raw/atac/$r1" \
    -I "data/raw/atac/$r2" \
    -o "results/atac/trimmed/${library}_1.fq.gz" \
    -O "results/atac/trimmed/${library}_2.fq.gz" \
    -h "results/atac/logs/${library}.html" \
    -j "results/atac/logs/${library}.json"


  # Align to TAIR10 and sort by read name.

  bowtie2 \
    --very-sensitive \
    -X 2000 \
    --dovetail \
    -p "$threads" \
    -x "$bt_index" \
    -1 "results/atac/trimmed/${library}_1.fq.gz" \
    -2 "results/atac/trimmed/${library}_2.fq.gz" \
    2> "results/atac/logs/${library}.bowtie2.txt" \
    | samtools sort -n -@ "$threads" -o "results/atac/bams/${library}.name.bam"


  # Correct mate information, sort by coordinate and remove duplicates.

  samtools fixmate -m "results/atac/bams/${library}.name.bam" \
    "results/atac/bams/${library}.fixmate.bam"

  samtools sort -@ "$threads" -o "results/atac/bams/${library}.sort.bam" \
    "results/atac/bams/${library}.fixmate.bam"

  samtools markdup -r "results/atac/bams/${library}.sort.bam" \
    "results/atac/bams/${library}.dedup.bam"

  samtools index "results/atac/bams/${library}.dedup.bam"


  # Retain nuclear alignments and remove blacklist overlaps.

  samtools view -b -F 3852 -q 10 "results/atac/bams/${library}.dedup.bam" 1 2 3 4 5 \
    | bedtools intersect -v -abam stdin -b data/reference/blacklist.bed \
    | samtools sort -n -@ "$threads" -o "results/atac/bams/${library}.nuclear.name.bam"


  # Retain fragments up to 120 bp, then take the 5-prime base of each mate.

  bedtools bamtobed -bedpe -i "results/atac/bams/${library}.nuclear.name.bam" \
    | awk '
        BEGIN {
          OFS = "\t"
        }

        $1 == $4 && $6 - $2 <= 120 && $6 > $2 {
          print $1, $2, $3, $7 "/1", $8, $9
          print $4, $5, $6, $7 "/2", $8, $10
        }
      ' \
    | awk '
        BEGIN {
          OFS = "\t"
        }

        $6 == "+" {
          print $1, $2, $2 + 1, $4, $5, $6
        }

        $6 == "-" {
          print $1, $3 - 1, $3, $4, $5, $6
        }
      ' \
    | bedtools sort -i stdin > "results/atac/beds/${library}.nfr0_120.bed"


  # Expand the one-base read starts to 100 bp for signal tracks.

  bedtools slop -g data/reference/chrom.sizes -l 50 -r 49 \
    -i "results/atac/beds/${library}.nfr0_120.bed" \
    | bedtools sort -i stdin \
    | bedtools genomecov -bg -g data/reference/chrom.sizes -i stdin \
    > "results/atac/bigWigs/${library}.bedGraph"

  bedGraphToBigWig "results/atac/bigWigs/${library}.bedGraph" data/reference/chrom.sizes \
    "results/atac/bigWigs/${library}.bw"


  # Call peaks and count reads in the peak regions.

  macs3 callpeak \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    -p 0.05 \
    -f BED \
    --keep-dup all \
    -g 1.191e8 \
    -t "results/atac/beds/${library}.nfr0_120.bed" \
    -n "$library" \
    --outdir results/atac/peaks

  bedtools coverage -counts -a data/reference/atac_peaks.bed \
    -b "results/atac/beds/${library}.nfr0_120.bed" \
    | awk '
        BEGIN {
          OFS = "\t"
        }

        {
          print $4, $NF
        }
      ' > "results/atac/counts/${analysis}.tsv"

done < metadata/atac_samples.tsv


# Combine the counts in the analysis sample order.

python3 scripts/lib/merge_atac_counts.py
