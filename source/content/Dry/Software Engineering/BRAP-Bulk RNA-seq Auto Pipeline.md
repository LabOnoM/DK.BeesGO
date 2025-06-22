---
title: STAR–Salmon–StringTie–DEXSeq RNA-seq Analysis Pipeline
draft: false
tags:
  - Dry
  - Linux
  - Console
  - RNA
  - High-throughput_Sequencing
---

This entry documents a fully automated RNA-seq processing pipeline that performs read alignment, transcript quantification, splicing analysis, transcriptome assembly, and statistical analysis. The pipeline dynamically adapts to input data quality and structure, ensuring robust and reproducible results.

---
# Short reads mapping
## 📋 Workflow Overview

This pipeline performs:

- ✅ FASTQ quality control and trimming
- ✅ STAR-based genome alignment
- ✅ Salmon transcript quantification (2 rounds)
- ✅ Transcriptome assembly via StringTie
- ✅ Transcriptome refinement via SQANTI3
- ✅ circRNA detection (CIRCexplorer2)
- ✅ Alternative splicing (LeafCutter)
- ✅ Exon-level quantification (DEXSeq)
- ✅ Final analysis and visualization in R
## 🗺️ Workflow Diagram
```mermaid
graph TD

  

subgraph Preprocessing

A1[Start Script and Load Config] --> A2[Download genome.fa and genome.gtf]

A2 --> A3[Filter GTF - cellranger mkgtf]

A3 --> A4[Extract transcriptome - genomeTx.fa]

A4 --> A5[Build Salmon index 1 - Ensembl GTF]

end

  

subgraph FASTQ Processing

B1[Scan FASTQ files] --> B2[Run FastQC]

B2 --> B3["Detect read quality and adapter content"]

B3 --> B4{Is QC pass?}

B4 -- Yes --> B5[Use raw FASTQ]

B4 -- No or auto --> B6[Run trim_galore and re-FastQC]

B1 --> B7["Detect Phred encoding (score)"]

B7 --> B6

B6 --> B8["Use trimmed FASTQ"]

end

  

subgraph Alignment and Quantification Round 1

B5 --> C1[STAR mapping]

B8 --> C1

C1 --> C2[Aligned.sortedByCoord.out.bam]

C1 --> C3[Chimeric.out.junction]

C2 --> C4[samtools index]

C2 --> C5[regtools junctions extract → .junc]

C4 --> C6[featureCounts]

A5 --> D1[Salmon quant 1 - Ensembl GTF]

B5 --> D1

B8 --> D1

D1 --> D2["Detect library type (strandness) from Salmon log"]

D2 --> C6

end

  

subgraph Alternative Splicing Analysis

C5 --> L1[LeafCutter: junction clustering]

end

  

subgraph circRNA Detection

C3 --> M1[CIRCexplorer2 parse]

M1 --> M2[CIRCexplorer2 annotate]

end

  

subgraph Transcriptome Assembly and Refinement

C2 --> E1[StringTie per-sample using STAR BAM]

E1 --> E2[StringTie merge - stdout.gtf]

E2 --> E3[SQANTI3 QC → stdout_final.gtf]

E3 --> E4[Extract transcriptome - genomeTx_stdout.fa]

E3 --> E5[Convert GTF to GFF for DEXSeq]

end

  

subgraph Alignment and Quantification Round 2

E4 --> F1[Build Salmon index 2 - Merged GTF]

F1 --> F2["Salmon quant 2 - Merged GTF (overwrites round 1)"]

B5 --> F2

B8 --> F2

end

  

subgraph DEXSeq Exon-Level Counting

E5 --> Z1[Prepare exon bins - genome_stdoutFinal.gff]

C2 --> Z2[Use STAR BAM for exon counting]

D2 --> Z3[Use strand info from Salmon log]

Z1 --> Z4[Run dexseq_count.py]

Z2 --> Z4

Z3 --> Z4

Z4 --> Z5[Generate exon count matrix]

end

  

subgraph Downstream Analysis

F2 --> Z6[MultiQC input ← final Salmon output]

Z5 --> Z6

Z6 --> Z7[MultiQC and Final Report]

Z7 --> Z8[Run 000.Analysis.R for all downstream analysis]

end
```

---

## 🔁 Step-by-Step Description

### 1. Environment Setup

- Activates `SQANTI3.env` Conda environment
- Parses settings from `ExConfiguration_bashScript.txt`

### 2. Genome Preparation

- Downloads genome FASTA and GTF from Ensembl
- Optionally adds exogenous sequences (e.g., GFP)
- Filters GTF using `cellranger mkgtf`
- Builds:
  - STAR index
  - Salmon transcriptome index (with decoy)
  - Pfam HMM database

### 3. FASTQ Processing

- Detects FASTQ files and normalizes extensions
- Runs **FastQC** on raw reads
- Detects:
  - **Phred encoding** from first 5k reads
  - **Read quality** and **adapter contamination**
- Triggers `trim_galore` if needed
- Reruns FastQC on trimmed reads

### 4. Alignment & Quantification (Round 1)

- Aligns reads using **STAR**
- Outputs:
  - `*.bam` (sorted and indexed)
  - `*.Chimeric.out.junction` (for circRNA)
  - `*.junc` (from `regtools` for LeafCutter)
- Quantifies transcripts using **Salmon** (based on Ensembl GTF)
- Detects **library strandness** from `salmon_quant.log`
- Applies strandness to `featureCounts`

### 5. Splicing and circRNA

- **LeafCutter** analyzes alternative splicing from `.junc` files
- **CIRCexplorer2** annotates circRNAs using STAR chimeric junctions

### 6. Transcriptome Assembly & Refinement

- **StringTie** assembles transcripts per sample using STAR BAMs
- Merges GTFs into a master transcriptome
- **SQANTI3** filters and corrects transcript models
- Generates:
  - `stdout_final.gtf`
  - `genomeTx_stdout.fa`
  - `genome_stdoutFinal.gff` for DEXSeq

### 7. Quantification (Round 2)

- Deletes first Salmon output
- Builds new Salmon index using `stdout_final.gtf`
- Re-runs **Salmon quantification (round 2)** — overwrites round 1
- These outputs are used for MultiQC and downstream analysis

### 8. DEXSeq Exon Counting

- Uses:
  - Refined GFF (`genome_stdoutFinal.gff`)
  - STAR BAMs
  - Strandness from Salmon
- Runs `dexseq_count.py` to generate `*.exon.txt` files

### 9. MultiQC Report

- Runs `multiqc .` to collect:
  - FastQC, Salmon, STAR, featureCounts, trim_galore, etc.
  - Includes only **final** Salmon outputs

### 10. Final R Analysis

- Executes: `Rscript 000.Analysis.R`
- Performs:
  - Normalization
  - PCA, RLE, clustering
  - DEG or splicing analysis
  - HTML and table outputs

---

## 🔍 Dynamic Features

| Property            | Detected From              | Applied In                  |
|---------------------|----------------------------|-----------------------------|
| Phred encoding      | First 5k reads (ASCII)     | `trim_galore`               |
| Read quality        | FastQC                     | Triggers trimming           |
| Adapter presence    | FastQC                     | Triggers trimming           |
| Library strandness  | `salmon_quant.log`         | `featureCounts`, DEXSeq     |

---

## 📦 Output Summary

| Output                        | Description                                |
|-------------------------------|--------------------------------------------|
| `*.bam`, `*.bai`              | STAR-aligned and indexed BAM files         |
| `Salmon_quants/`             | Final transcript quantification (TPM, etc) |
| `*.counts`                   | Gene-level quantification from featureCounts |
| `*.exon.txt`                 | DEXSeq exon bin counts                     |
| `LeafCutter.Output/`         | Splicing analysis results                  |
| `*_CIRCexplorer2.ce`         | circRNA annotations                        |
| `stdout_final.gtf`           | Refined transcriptome annotation           |
| `multiqc_report.html`        | Aggregated QC and alignment metrics        |
| `00.Analysis_Report.html`    | Custom final R analysis report             |

---

## 🚀 Getting Started

1. Prepare the input directory with:
   - `ExConfiguration_bashScript.txt`
   - Raw or trimmed FASTQ files
2. Run:

   ```bash
   bash 00.STAR.sh <path_to_your_project>
   ```

3. Inspect outputs:
   - BAMs, counts, and QC
   - HTML reports
   - Final analysis results from `000.Analysis.R`

---

## 📚 Dependencies

- `STAR`, `Salmon`, `StringTie`, `samtools`, `featureCounts`
- `trim_galore`, `FastQC`, `regtools`, `CIRCexplorer2`
- `SQANTI3`, `LeafCutter`, `DEXSeq`, `MultiQC`
- `R` with custom script: `000.Analysis.R`

---

## 🧠 Notes

- The pipeline is **idempotent**: already-completed steps are skipped unless missing.
- Highly suitable for both **canonical and custom transcriptome exploration**.
- Ideal for both bulk RNA-seq and pseudo-bulk analyses.

---
# Comprehensive analysis in R
