# DK.BeesGO - Mapping Science: From Bench to Bioinformatics

[DK.BeesGO](https://www.bs-gou.com/DK.BeesGO/) is a digital knowledge repository designed to integrate both wet-lab and dry-lab experiments in a unified, interactive platform. By leveraging the bidirectional linking capabilities of Quartz, DK.BeesGO enables a dynamic Graph View that visually maps the interconnections between biological materials, experimental procedures, analytical data, and bioinformatics workflows.

## Executive Summary
This repository houses the Quartz source code, markdown content, and templates for hosting the DK.BeesGO notebook. Biological materials (such as Nucleic Acids, Proteins, and Lipids) are mapped parallelly under both `Wet/` and `Dry/` folders to trace their transition from bench protocols to dry-lab computational analyses.

- **Status**: Onboarded & Active.
- **Key Methods**: Wet-lab molecular protocols (in situ hybridization, reagents), Dry-lab bioinformatics (Spatial Transcriptomics, Image Processing, RNA-seq toolkits).

## Hypothesis Evolution
- **H1 (Mapping Bench to Bioinformatics)**: Biological materials represent a central entity whose physical properties (Wet) and digital representations (Dry) can be bidirectionally mapped to enable holistic experimental audit trails.

## Verified Timeline
- **2026-06-08** `[verified]`: Project cloned and initialized. Verified Git and Git LFS availability. Initialized AROS onboarding, `re_gent` VCS layer, and LLM-Wiki indexing.

## Repository Structure Map
- `.wiki/` (symlinked as `Wiki`): The internal knowledge base (LLM-Wiki) for agent grounding.
- `00.RawData/`: Primary data repository (contains `INDEX.csv`).
- `source/`: Quartz source code and site content.
  - `source/content/`: Main markdown pages for the notebook.
    - `source/content/Wet/`: Wet-lab experimental protocols.
    - `source/content/Dry/`: Dry-lab computational workflows and Jupyter notebooks.
    - `source/content/BSGOU-System/`: Internal discussion files, RNA-seq toolkits, and journal club logs.
    - `source/content/ShinyApps in BSGOU-choice/`: Curated outstanding ShinyApps.
- `AGENTS.md`: Agent governance policies and testing specifications.

## Local Development & Rendering
To build and serve the Quartz notebook locally:
1. Navigate to the `source/` folder.
2. Run `npm install` to install dependencies.
3. Run `npx quartz build --serve`.

For details on the deployment of Obsidian notebooks to GitHub Pages, refer to the original template tutorial: https://dev.to/defenderofbasic/host-your-obsidian-notebook-on-github-pages-for-free-8l1.
