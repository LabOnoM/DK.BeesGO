---
title: "BRAP: Bulk RNA-seq Auto-Pipeline"
draft: false
tags:
  - Dry
  - Linux
  - Console
  - RNA
  - High-throughput_Sequencing
---

This entry documents a fully automated **B**ulk **R**NA-seq **A**uto **P**ipeline (**BRAP**) that performs read alignment, transcript quantification, splicing analysis, transcriptome assembly, and statistical analysis. The pipeline dynamically adapts to input data quality and structure, ensuring robust and reproducible results.

---

## 1. Overview

The diagram below summarizes how the main scripts interact. Clicking **Run** in the Shiny interface writes configuration files and spawns the shell workflow which finishes by generating the final report.

```mermaid
graph TD
  subgraph "00.launcher.sh"
    L["Open Shiny on random port"]
  end
  L --> S["Shiny app"]
  subgraph "Shiny app"
    S --> G[global.R]
    S --> U[ui.R]
    S --> V[server.R]
  end
  subgraph "server.R after Run"
    V --> C["Write config files"]
    C --> B["Call 00.STAR.sh"]
  end
  B --> STAR[00.STAR.sh]
  subgraph "00.STAR.sh"
    STAR --> Rcall["Rscript 000.Analysis.R"]
  end
  Rcall --> R[000.Analysis.R]
  subgraph "R_Source functions"
    R --> F0[Function00_ParameterSetting.R]
    R --> F1[Function01_STAR_FeatureCounts_RUVSeq.R]
    R --> F2[Function02_DESeq2_PreProcessing.R]
    R --> F3[Function03_IsoformSwitchAnalyzer.R]
    R --> F4[Function04_DESeq2_GSEAenrich_Heatmap_PCA_MultiComparsion.R]
    R --> F5[Function05_DEXseq.R]
    R --> F6[Function06_CircRNA.R]
    R --> F7[Function07_Summary.R]
  end
  F7 --> Report[00.Analysis_Report.html]
```

## 2. Short reads mapping

The main bash script performs read mapping and quantification. A detailed workflow for this step is provided below (unchanged from the previous version).
### 2.1. Workflow Overview

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
#### 2.1.1 Workflow Diagram
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

### 2.2. Step-by-Step Description

#### 2.2.1 Environment Setup

- Activates `SQANTI3.env` Conda environment
- Parses settings from `ExConfiguration_bashScript.txt`

#### 2.2.2 Genome Preparation

- Downloads genome FASTA and GTF from Ensembl
- Optionally adds exogenous sequences (e.g., GFP)
- Filters GTF using `cellranger mkgtf`
- Builds:
  - STAR index
  - Salmon transcriptome index (with decoy)
  - Pfam HMM database

#### 2.2.3 FASTQ Processing

- Detects FASTQ files and normalizes extensions
- Runs **FastQC** on raw reads
- Detects:
  - **Phred encoding** from first 5k reads
  - **Read quality** and **adapter contamination**
- Triggers `trim_galore` if needed
- Reruns FastQC on trimmed reads

#### 2.2.4 Alignment & Quantification (Round 1)

- Aligns reads using **STAR**
- Outputs:
  - `*.bam` (sorted and indexed)
  - `*.Chimeric.out.junction` (for circRNA)
  - `*.junc` (from `regtools` for LeafCutter)
- Quantifies transcripts using **Salmon** (based on Ensembl GTF)
- Detects **library strandness** from `salmon_quant.log`
- Applies strandness to `featureCounts`

#### 2.2.5 Splicing and circRNA

- **LeafCutter** analyzes alternative splicing from `.junc` files
- **CIRCexplorer2** annotates circRNAs using STAR chimeric junctions

#### 2.2.6 Transcriptome Assembly & Refinement

- **StringTie** assembles transcripts per sample using STAR BAMs
- Merges GTFs into a master transcriptome
- **SQANTI3** filters and corrects transcript models
- Generates:
  - `stdout_final.gtf`
  - `genomeTx_stdout.fa`
  - `genome_stdoutFinal.gff` for DEXSeq

#### 2.2.7 Quantification (Round 2)

- Deletes first Salmon output
- Builds new Salmon index using `stdout_final.gtf`
- Re-runs **Salmon quantification (round 2)** — overwrites round 1
- These outputs are used for MultiQC and downstream analysis

#### 2.2.8 DEXSeq Exon Counting

- Uses:
  - Refined GFF (`genome_stdoutFinal.gff`)
  - STAR BAMs
  - Strandness from Salmon
- Runs `dexseq_count.py` to generate `*.exon.txt` files

#### 2.2.9 MultiQC Report

- Runs `multiqc .` to collect:
  - FastQC, Salmon, STAR, featureCounts, trim_galore, etc.
  - Includes only **final** Salmon outputs

#### 2.2.10 Final R Analysis

- Executes: `Rscript 000.Analysis.R`
- Performs:
  - Normalization
  - PCA, RLE, clustering
  - DEG or splicing analysis
  - HTML and table outputs

---

### 2.3. Dynamic Features

| Property            | Detected From              | Applied In                  |
|---------------------|----------------------------|-----------------------------|
| Phred encoding      | First 5k reads (ASCII)     | `trim_galore`               |
| Read quality        | FastQC                     | Triggers trimming           |
| Adapter presence    | FastQC                     | Triggers trimming           |
| Library strandness  | `salmon_quant.log`         | `featureCounts`, DEXSeq     |

---

### 2.4. Output Summary

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

### 2.5. Getting Started

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

### 2.6. Dependencies

- `STAR`, `Salmon`, `StringTie`, `samtools`, `featureCounts`
- `trim_galore`, `FastQC`, `regtools`, `CIRCexplorer2`
- `SQANTI3`, `LeafCutter`, `DEXSeq`, `MultiQC`
- `R` with custom script: `000.Analysis.R`

---

### 2.7. Notes

- The pipeline is **idempotent**: already-completed steps are skipped unless missing.
- Highly suitable for both **canonical and custom transcriptome exploration**.
- Ideal for both bulk RNA-seq and pseudo-bulk analyses.

---

## 3.  Comprehensive analysis in R

The R script `000.Analysis.R` orchestrates downstream processing once the bash pipeline is finished. It sequentially calls helper functions located under `R_Source/`.

### 3.1 P0_ParameterSetting
The workflow initializes parameters for all downstream steps.
```mermaid
graph TD
  R["000.Analysis.R"] --> F0
  subgraph Function00_ParameterSetting.R
    F0 --> A["Create result folders"]
    A --> B["Load ExDesign.csv"]
    B --> C["Load ExComparison.csv"]
    C --> D["Connect to Ensembl"]
    D --> E["Save Parameters.Rdata"]
  end
  E --> Next["P1_FeatureCounts_import_RUVSeq"]
```
### 3.2 P1_FeatureCounts_import_RUVSeq
This step imports counts, builds annotations and applies RUVSeq correction.
```mermaid
graph TD
  Prev["P0_ParameterSetting"] --> Start
  subgraph Function01_STAR_FeatureCounts_RUVSeq.R
    Start --> A["Load Parameters.Rdata"]
    A --> B{"Annotation.Rdata?"}
    B -- No --> C["Extract gene list"]
    C --> D["Query Ensembl"]
    D --> E["Fill Entrez via NCBI"]
    E --> F["UniProt/GO/KEGG"]
    F --> G["Save Annotation.Rdata"]
    B -- Yes --> G
    G --> H["Build OrgDb & BSgenome?"]
    H --> I{"CountsFrom"}
    I -- FeatureCounts --> J["Read *.STAR.counts"]
    J --> K["Gene lengths"]
    I -- Tximport --> L["Import Salmon quant"]
    L --> M["Design matrix"]
    M --> N["Create SwitchAnalyzeRlist"]
    N --> K
    K --> O["Filter low expr."]
    O --> P[SeqExpressionSet]
    P --> Q["RLE & PCA before"]
    Q --> R["Between-lane norm"]
    R --> S["RUVs correction"]
    S --> T["RLE & PCA after"]
    T --> U["Save STARFeatureRawData.Rdata"]
  end
  U --> Next["P2_DESeq2_PreProcessing"]
```
### 3.3 P2_DESeq2_PreProcessing
This function normalizes counts and prepares DESeq2 objects.
```mermaid
graph TD
  Prev["P1_FeatureCounts_import_RUVSeq"] --> Start
  subgraph Function02_DESeq2_PreProcessing.R
    Start --> A["Load Annotation.Rdata"]
    A --> B["Load Parameters.Rdata"]
    B --> C["Load STARFeatureRawData.Rdata"]
    C --> D["Create DESeq2 dataset"]
    D --> E["Run DESeq"]
    E --> F["Calculate TPM and log2TPM"]
    F --> G["Export Excel"]
    G --> H["Save DESeq2_Processing.Rdata"]
  end
  H --> Next["P3_IsoformSwitchAnalyzeR"]
```
### 3.4 P3_IsoformSwitchAnalyzeR
This step detects isoform switches.
```mermaid
graph TD
  Prev["P2_DESeq2_PreProcessing"] --> Start
  subgraph Function03_IsoformSwitchAnalyzer.R
    Start --> A["Load DESeq2_Processing"]
    A --> B["Prepare SwitchAnalyzeRlist"]
    B --> C["isoformSwitchTestDRIMSeq"]
    C --> D["Save aSwitchList.Rdata"]
  end
  D --> Next["P4_DESeq2_GSEA_MultiComparsion"]
```
### 3.5 P4_DESeq2_GSEA_MultiComparsion
This function performs differential expression and enrichment analyses.
```mermaid
graph TD
  Prev["P3_IsoformSwitchAnalyzeR"] --> Start
  subgraph Function04_DESeq2_GSEAenrich_Heatmap_PCA_MultiComparsion.R
    Start --> A["Load DESeq2_Processing"]
    A --> B["Run DESeq2"]
    B --> C["clusterProfiler enrichment"]
    C --> D["Heatmaps & PCA"]
    D --> E["Save DESeq2_Results.Rdata"]
  end
  E --> Next["P5_DEXseqSummary"]
```
### 3.6 P5_DEXseqSummary
This function computes exon usage with DEXSeq.
```mermaid
graph TD
  Prev["P4_DESeq2_GSEA_MultiComparsion"] --> Start
  subgraph Function05_DEXseq.R
    Start --> A["Load exon counts"]
    A --> B["Prepare exon bins"]
    B --> C["Run DEXSeq"]
    C --> D["Save DEXseqResOut.Rdata"]
  end
  D --> Next["P6_CircRNA"]
```
### 3.7 P6_CircRNA
This function summarizes circRNA evidence.
```mermaid
graph TD
  Prev["P5_DEXseqSummary"] --> Start
  subgraph Function06_CircRNA.R
    Start --> A["Load circRNA data"]
    A --> B["Identify circRNAs"]
    B --> C["Motif analysis"]
    C --> D["Save circRNA_Results.Rdata"]
  end
  D --> Next["P7_WebSummary"]
```
### 3.8 P7_WebSummary
This step generates the HTML report.
```mermaid
graph TD
  Prev["P6_CircRNA"] --> Start
  subgraph Function07_Summary.R
    Start --> A["Gather all results"]
    A --> B["Render Analysis_Report.rmd"]
    B --> C["Save 00.Analysis_Report.html"]
  end
  C --> Output["00.Analysis_Report.html"]
```

## 4. Structure of Analysis_Report.html
The HTML report generated by `P7_WebSummary` contains:
- **0 Description of Workflow**
- **1 Sample Sequencing Statistics** with reference info, QC, Salmon statistics and RUVSeq correction
- **2 Analysis** covering differential expression, enrichment, isoform switches and circular RNA
- **3 Support** and **4 Citation** sections
- **5 Supplementary** with source code and running log

---

# 5. Example results from this workflow

1. https://d3dcaz4rv8jgb4.cloudfront.net/ from [Zheng Y, Wang Z, Weng Y, Sitosari H, He Y, Zhang X, Shiotsu N, Fukuhara Y, Ikegame M, Okamura H. **Gingipain regulates isoform switches of PD-L1 in macrophages infected with Porphyromonas gingivalis**. _Scientific reports_. **2025** Mar 26;15(1):10462.](https://www.nature.com/articles/s41598-025-94954-7)
2. https://dndy5us1uro9a.cloudfront.net/BulkRNAseq/00.Analysis_Report.html from [Weng Y, Wang Z, Sitosari H, Ono M, Okamura H, Oohashi T. **_O_‐GlcNAcylation regulates osteoblast differentiation through the morphological changes in mitochondria, cytoskeleton, and endoplasmic reticulum**. _BioFactors_. **2025** Jan;51(1):e2131.](https://iubmb.onlinelibrary.wiley.com/doi/abs/10.1002/biof.2131)