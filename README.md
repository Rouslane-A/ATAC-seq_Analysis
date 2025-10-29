## 🧬 ATAC-seq Analysis Pipeline (Arabidopsis thaliana)

### 📘 Overview

This repository provides a reproducible **ATAC-seq analysis workflow** based on the dataset from
**Jégu et al., 2017**, which explored the role of chromatin remodeling in *Arabidopsis thaliana* seedling morphogenesis.

The tutorial walks through each step — from downloading raw reads to peak calling and visualization — using widely adopted bioinformatics tools (Bowtie2, Samtools, MACS2, and UCSC utilities).

---

### 🔬 Workflow Summary

The complete analysis includes the following steps:

1. **Download ATAC-seq FASTQ files** from ENA (`PRJNA351855`)
2. **Quality control** using `FastQC`
3. **Reference genome setup** (Araport or EnsemblPlants)
4. **Genome indexing** with `Bowtie2`
5. **Alignment** and post-processing with `Samtools`
6. **Organellar read removal** (chloroplast and mitochondrial)
7. **Peak calling** with `MACS2`
8. **Signal track generation** (BedGraph → BigWig)
9. **Visualization** using IGV or UCSC Genome Browser

---

### ⚙️ Dependencies

| Tool                                                                 | Version | Description                   |
| -------------------------------------------------------------------- | ------- | ----------------------------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | ≥0.11   | Quality control of raw reads  |
| [Cutadapt](https://cutadapt.readthedocs.io/)                         | ≥3.0    | Adapter trimming              |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/)                | ≥2.4    | Sequence alignment            |
| [Samtools](http://www.htslib.org/)                                   | ≥1.10   | BAM/SAM manipulation          |
| [MACS2](https://github.com/macs3-project/MACS)                       | ≥2.2    | Peak calling                  |
| [Bedtools](https://bedtools.readthedocs.io/)                         | ≥2.29   | Genomic arithmetic            |
| [UCSC Tools](https://hgdownload.soe.ucsc.edu/admin/exe/)             | —       | `bedClip`, `bedGraphToBigWig` |

---

### 🧩 Installation

Clone this repository:

```bash
git clone https://github.com/<your-username>/atacseq-arabidopsis.git
cd atacseq-arabidopsis
```

Create and activate a conda environment (recommended):

```bash
conda create -n atacseq_env fastqc cutadapt bowtie2 samtools macs2 bedtools -c bioconda
conda activate atacseq_env
```

---

### 🚀 Usage

#### 1️⃣ Download ATAC-seq FASTQ files

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR473/002/SRR4733912/SRR4733912_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR473/002/SRR4733912/SRR4733912_2.fastq.gz
```

#### 2️⃣ Quality Check

```bash
mkdir 01_AnalysisStepOne
fastqc -o Quality_ATAC *.gz
```

#### 3️⃣ Genome Preparation

```bash
ln -s /path/to/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa At_Genome
bowtie2-build At_Genome bwt_index/At.TAIR10
```

#### 4️⃣ Alignment

```bash
bowtie2 --threads 8 -x bwt_index/At.TAIR10 \
  -1 SRR4733912_1.fastq.gz -2 SRR4733912_2.fastq.gz \
  -S bwt_out/SRR4733912.sam
```

#### 5️⃣ Convert, Sort, and Index

```bash
samtools view -bS SRR4733912.sam | samtools sort -o SRR4733912.sorted.bam
samtools index SRR4733912.sorted.bam
```

#### 6️⃣ Remove Organellar Reads

```bash
samtools idxstats SRR4733912.sorted.bam | cut -f1 | grep -v Mt | grep -v Pt \
  | xargs samtools view -b SRR4733912.sorted.bam > SRR4733912.sorted.noorg.bam
```

#### 7️⃣ Peak Calling

```bash
macs2 callpeak -t SRR4733912.sorted.noorg.bam -q 0.05 --broad -f BAMPE \
  -n SRR4733912 -B --trackline --outdir peaks/
```

#### 8️⃣ Generate BigWig for Visualization

```bash
bioawk -c fastx '{print $name, length($seq)}' At_Genome > At_chr.sizes
bedtools slop -i SRR4733912_treat_pileup.bdg -g At_chr.sizes -b 0 \
  | bedClip stdin At_chr.sizes SRR4733912_treat_pileup.clipped.bdg
sort -k1,1 -k2,2n SRR4733912_treat_pileup.clipped.bdg > SRR4733912_treat_pileup.sorted.bdg
bedGraphToBigWig SRR4733912_treat_pileup.sorted.bdg At_chr.sizes SRR4733912.bw
```

---

### 📊 Outputs

| File                          | Description                       |
| ----------------------------- | --------------------------------- |
| `SRR4733912.sorted.noorg.bam` | Aligned reads (filtered)          |
| `SRR4733912_peaks.broadPeak`  | Broad peaks from MACS2            |
| `SRR4733912_treat_pileup.bdg` | Signal intensity bedGraph         |
| `SRR4733912.bw`               | BigWig for IGV/UCSC visualization |

---

### 🧠 References

* Jégu T. *et al.* (2017). **The Arabidopsis SWI/SNF Protein BAF60 Mediates Seedling Growth Control by Modulating Chromatin Accessibility**. *Plant Cell*, 29(9): 2214–2231.
* [MACS2 documentation](https://macs3-project.github.io/MACS/)
* [Bowtie2 documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)


## 🚧 Next Steps — Expanding This Project

Here’s how you can extend this repository into a **modern bioinformatics pipeline**:

---

### 🔁 1. Build a **Nextflow Pipeline**

Make this workflow reproducible and scalable.

**Key modules:**

* `download` — fetch raw FASTQs from ENA
* `qc` — FastQC reports
* `align` — Bowtie2 alignment
* `filter` — remove organellar reads
* `peak_calling` — MACS2
* `track_gen` — generate BigWig tracks

📦 Output: reproducible workflow supporting any organism or dataset.
💡 Use Docker/Singularity containers for tool dependencies.

---

### 📊 2. Develop an **Interactive Dashboard**

Use **Streamlit** or **Dash** to visualize:

* QC summaries (FastQC/MultiQC)
* Alignment metrics (SAMtools idxstats)
* Peak distribution plots (chromosome-level)
* Signal intensity tracks (via pyBigWig)
* Volcano/MA plots for differential peaks (future extension)

---


### 🧩 3. Integrate with **Downstream Analysis**

* Differential accessibility (e.g., using **DiffBind**, **csaw**, or **DESeq2**)
* Motif enrichment (via **HOMER** or **MEME Suite**)
* Co-analysis with RNA-seq data to link open chromatin with expression

---

