#!/bin/bash

READS_DIR="../00_RawData"
OUT_DIR="bwt_out"
INDEX="bwt_index/At.TAIR10"

mkdir -p "$OUT_DIR"

for R1 in ${READS_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="${READS_DIR}/${SAMPLE}_2.fastq.gz"

    echo "Processing $SAMPLE ..."

    # Step 1: Align with Bowtie2
    bowtie2 --threads 8 -x "$INDEX" \
        -q \
        -1 "$R1" \
        -2 "$R2" \
    | samtools view -@ 8 -bS - \
    | samtools sort -@ 8 -o "${OUT_DIR}/${SAMPLE}.sorted.bam"

    # Step 2: Mark & remove duplicates
    samtools markdup -r -@ 8 "${OUT_DIR}/${SAMPLE}.sorted.bam" "${OUT_DIR}/${SAMPLE}.dedup.bam"

    # Step 3: Remove mitochondrial & chloroplast reads
    samtools idxstats "${OUT_DIR}/${SAMPLE}.dedup.bam" | cut -f 1 > "${OUT_DIR}/all_chroms.txt"
    # Keep all except organelles (modify names if needed)
    grep -Ev "Mt|Pt|chloroplast|mitochondria|chrM" "${OUT_DIR}/all_chroms.txt" > "${OUT_DIR}/keep_chroms.txt"
    samtools view -@ 8 -b "${OUT_DIR}/${SAMPLE}.dedup.bam" -o "${OUT_DIR}/${SAMPLE}.filtered.bam" -L <(cat "${OUT_DIR}/keep_chroms.txt")

    # Step 4: Index final BAM
    samtools index -@ 8 "${OUT_DIR}/${SAMPLE}.filtered.bam"

    echo "✅ Finished $SAMPLE → ${SAMPLE}.filtered.bam"
done

