# DK.BeesGO Repository Summary

This document presents a synthesized overview of the DK.BeesGO project, compiling structured information directly from the internal LLM-Wiki (`.wiki/`).

## 1. Project Mission & Overview
DK.BeesGO is a digital knowledge repository designed to bridge wet-lab empirical experiments with dry-lab computational and bioinformatics workflows [[DK.BeesGO]]. By establishing bidirectional links, DK.BeesGO maps connections between biological material entities, experimental protocols, analytical data, and downstream pipelines.

The deployment of the knowledge base is automated via GitHub Actions, rendering an interactive Obsidian notebook using Quartz [[readme]].

## 2. Structural Paradigm: Parallel Wet/Dry Design
A core conceptual feature of DK.BeesGO is its parallel directory structure [[overview]]. Topics are organized under cellular or molecular categories (such as `Nucleic Acid`, `Protein`, `Lipid`), which are intentionally mirrored under both:
- **`Wet/`**: Housing experimental protocols, reagents, and material metadata.
- **`Dry/`**: Housing analytical scripts, computational pipelines, and software dependencies.

This structure allows researchers to easily trace biological samples from initial bench assays to final bioinformatics analyses.

## 3. Curated Components & Interactive Tools
DK.BeesGO integrates several resources for bioscience education and research:
- **ShinyApps Choice**: High-quality Shiny applications curated from the bioinformatics community (e.g., `iDEP`, `ShinyGO`, `Heatmapper`) to lower barriers to entry for non-programming users [[DK.BeesGO]].
- **Journal Club**: Multi-turn analysis of recent literature, covering topics like inflammaging, PTSD brain fog, and genomic evolution [[overview]].
- **RNA-seq Tool Kits**: Dedicated bioinformatics pipelines, such as `BRAP`, `SCRIB`, and `GREP1` [[overview]].

## 4. History & Governance
The project baseline was established on **2026-06-08**, marking the full onboarding of the repository with:
- **Git LFS** configuration for binary tracking [[timeline]].
- **re_gent** version control layer initialization for AI agent auditability [[timeline]].
- **AGENTS.md** policies regulating workspace hygiene and wiki-first resolution.

---
*Sources Cited:*
- `[[readme]]` (.wiki/sources/readme.md)
- `[[DK.BeesGO]]` (.wiki/entities/DK.BeesGO.md)
- `[[overview]]` (.wiki/overview.md)
- `[[timeline]]` (.wiki/timeline.md)
