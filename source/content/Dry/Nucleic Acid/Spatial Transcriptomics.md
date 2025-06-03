---
title: Spatial Transcriptomics
tags:
  - RNA-seq
  - Dry
  - Transcriptomics
  - Spatial
  - Biotech
---

Spatial transcriptomics is a cutting-edge methodology that enables researchers to measure and map gene expression within intact tissue sections, preserving the spatial context of gene activity. Unlike traditional RNA sequencing techniques that require tissue dissociation—thereby losing information about the original cellular environment—spatial transcriptomics maintains the architecture of the tissue, allowing for a more comprehensive understanding of cellular function and interaction within their native microenvironments.

---

### 🧬 Core Principles

At its core, spatial transcriptomics combines high-throughput gene expression profiling with spatial information, allowing for the localization of mRNA transcripts within tissue sections. This spatially resolved approach provides insights into how gene expression varies across different regions of a tissue, which is crucial for understanding complex biological processes such as development, disease progression, and tissue organization.

---

### 🔬 Technological Approaches

Spatial transcriptomics encompasses a variety of techniques, broadly categorized into:

- **Imaging-Based Methods**: These techniques are successors of [[in situ hybridization]] *in situ* hybridization, such as single-molecule fluorescence in situ hybridization (smFISH), multiplexed error-robust fluorescence in situ hybridization (MERFISH), and sequential fluorescence *in situ* hybridization (seqFISH), utilize fluorescent probes to visualize and quantify RNA molecules directly within tissue sections. They offer high spatial resolution, often at the single-cell or subcellular level, but are typically limited to analyzing a predefined set of genes.
- **Sequencing-Based Methods**: Approaches like 10x Genomics' Visium and BGI's Stereo-seq involve capturing RNA from tissue sections using spatially barcoded arrays or beads, followed by next-generation sequencing. These methods enable unbiased, whole-transcriptome profiling while retaining spatial information, although they may have lower spatial resolution compared to imaging-based techniques.

---

### 🧪 Applications

Spatial transcriptomics has broad applications across various fields of biomedical research, including:

- **Cancer Research**: Understanding tumor heterogeneity, microenvironment interactions, and identifying spatially distinct biomarkers.
- **Neuroscience**: Mapping gene expression across different brain regions to study neural development and function.
- **Developmental Biology**: Investigating gene expression patterns during embryogenesis and tissue differentiation.
- **Immunology**: Exploring the spatial organization of immune cells within tissues to understand immune responses.

---

## Summary

|Feature|Imaging-Based|Sequencing-Based|
|---|---|---|
|**Technique**|Microscopy with fluorescent probes|Spatial barcoding and RNA sequencing|
|**Resolution**|High (single-cell/subcellular)|Variable (multi-cell to near single-cell)|
|**Gene Detection**|Targeted (hundreds to thousands)|Unbiased (whole transcriptome)|
|**Throughput**|Moderate|High|
|**Data Complexity**|High (image processing required)|High (sequencing data analysis)|
|**Best Use Cases**|Detailed spatial mapping of known genes|Discovery of novel expression patterns
In summary, spatial transcriptomics represents a significant advancement in molecular biology, providing a powerful tool to study gene expression within the spatial context of intact tissues. By integrating spatial information with transcriptomic data, researchers can gain deeper insights into the complex interplay between cellular function and tissue architecture.