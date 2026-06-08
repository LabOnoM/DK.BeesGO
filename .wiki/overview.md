# Repository Overview

This document provides a high-level conceptual overview of the DK.BeesGO repository structure and scientific organization.

## Parallel Structure Design

To maintain scalability and map biological molecules directly from bench to computer, topics are organized in parallel directories:
- `/Wet/`: Reagents, experimental protocols, biological materials.
- `/Dry/`: Analytical scripts, bioinformatics pipelines, software dependencies.

This parallel structure facilitates tracing how the same biological materials transition across domains.

## Curated Components

- **Journal Club**: Multi-turn analysis of recent papers (e.g. inflammaging, PTSD brain fog, genomic mutations).
- **RNA-seq Tool Kits**: Diagnostic and analytical pipelines including `BRAP`, `SCRIB`, `GREP1`.
- **ShinyApps Choice**: Interactive web environments curated for DeSci and education.

## Related Entities

- [[DK.BeesGO]]
- [[SCHEMA]]
