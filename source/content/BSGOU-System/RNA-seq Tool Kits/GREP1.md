---
title: "GREP1: GEO Retrieval & Extraction Pipeline in One shiny-app"
draft: false
tags:
  - Dry
  - Linux
  - Console
  - R
  - Shiny-App
  - Database
---

GREP1 provides a Shiny interface for retrieving and preparing GEO/SRA sequencing data. The application guides the user through three sequential tasks:

1. **Retrieve SRR information** from user-supplied GSE accession numbers.
2. **Download** `.sra` or TenX BAM files.
3. **Decompress** the archives to FASTQ files for downstream analysis.

## Folder overview
- `global.R` – shared setup loaded by both UI and server components.
- `ui.R` – builds the Shiny layout with modules for each step.
- `server.R` – coordinates the modules and maintains application state.
- `01.ShinyModules/` – server and UI code for the downloader and decompressor.
- `03.R_Source/` – standalone R scripts (`RetrieveGSEinfo.R`, `DownloadSRA.R`,
  `DecompressSRA.R`, `ReArrangeFiles.R`).
- `04.bash_Source/` – shell helpers called by the R scripts.
- `00.launcher.sh` – convenience script to set up the environment and start the
  Shiny server.

## Step 1: Get SRR ID list
The first step fetches GSM and SRR information from GEO for each input GSE accession. A child R process executes `RetrieveGSEinfo.R`, which collects run metadata using RSelenium and the ENA API.

```mermaid
flowchart TD
    subgraph "Downloader_server.R"
        A[User clicks GET ID list]
        A --> B["observeEvent: Step 1"]
        B --> C[Save GSE_IDls & ns to Downloader_server_para.RData]
        C --> D[Kill port 4778 and old process]
        D --> E[Launch Rscript RetrieveGSEinfo.R --WD <WD>]
        E --> Q[Monitor progress / errors]
        Q --> R[Load Downloader_server_rout.RData]
        R --> S[values$out_tb <- out_tb]
        S --> T[para$Step1_done <- 1]
        T --> U["Render DataTable (Step1_ui)"]
    end
    subgraph "RetrieveGSEinfo.R operations"
        F[Load packages & parse WD]
        F --> G[Load GSE_IDls from RData]
        G --> H{Loop over GSE IDs}
        H --> I[Fetch GSM list from GEO]
        I --> J{Large GSM list?}
        J -->|">1000"| K[Parallel foreach]
        J -->|"<=1000"| L[Sequential loop]
        K --> M[Collect SRR IDs]
        L --> M
        M --> N[Scrape run metadata via RSelenium]
        N --> O[Determine DataType via ENA API]
        O --> P["Write 00.GSE_SRR_List.csv\nand Downloader_server_rout.RData"]
    end
    D --> F
    P --> Q
```

## Step 2: Download .sra or TenX BAM files
`DownloadSRA.R` reads the SRR list, fetches TenX BAM files if requested, and runs a shell script to prefetch SRA archives in parallel.

```mermaid
flowchart TD
    subgraph "Downloader_server.R"
        A["User clicks Start Download"]
        B["observeEvent in Downloader_server.R"]
        C["Collect DataType selections\nand rows to download"]
        D["Save Downloader_server_para2.RData"]
        E["Launch Rscript DownloadSRA.R --WD <WD>"]
        A --> B
        B --> C
        C --> D
        D --> E
    end

    subgraph "DownloadSRA.R operations"
        F["Load parameters & previous GSE_SRR lists"]
        G["Skip already downloaded SRR IDs"]
        H["Update DataType selections in 00.GSE_SRR_List.csv"]
        I{"TenX BAM samples?"}
        J["Fetch BAM links via RSelenium\nParallel download & bamtofastq"]
        K["Skip"]
        L["Write interim GSE_SRR_List.csv"]
        M["Invoke 01.GEO_SRA_Download.sh"]
        N["Check read types with vdb-dump"]
        O["Write final GSE_SRR_List.csv"]
        P["Monitor progress / console"]
        Q["Process finishes"]
        R["Load GSE_SRR_List.csv"]
        S["para$Step2_done <- 1"]
        T["Render DataTable (Step2_ui)"]
        E --> F
        F --> G
        G --> H
        H --> I
        I -- yes --> J
        I -- no --> K
        J --> L
        K --> L
        L --> M
    end

    subgraph "01.GEO_SRA_Download.sh"
        M1["Load GSE_SRR_List.csv to get SRR IDs"]
        M2["Init counters and progress files"]
        M3["task(sra_id): prefetch with retries\nvalidate using vdb-validate"]
        M4["GNU parallel -j <core> task ::: SRR IDs"]
        M5["Update .completed_jobs.count via flock"]
        M1 --> M2
        M2 --> M3
        M3 --> M4
        M4 --> M5
    end

    M --> M1
    M5 --> N
    N --> O
    O --> P
    P --> Q
    Q --> R
    R --> S
    S --> T
```

## Step 3: Decompress SRA files
After downloading, `DecompressSRA.R` rearranges the archives and invokes `02.fasterq_dump_gzip.sh` to produce compressed FASTQ files. TenX BAM samples are renamed and tracked in history logs.

```mermaid
flowchart TD
    subgraph "Decompressor_server.R"
        A["User clicks Start Decompress"]
        A --> B
        B["observeEvent(input$Decompress)"]
        B --> C
        C["Collect DataType edits\nfrom values$tout03"]
        C --> D
        D["Write updated GSE_SRR_List.csv"]
        D --> E
        E["Save DecompressSRA_server_para.RData"]
        E --> F
        F["Kill old process if running"]
        F --> G
        G["Launch Rscript DecompressSRA.R --WD <WD>"]
        G --> H
        H["Monitor progress / console"]
        H --> I
        I["Process finishes"]
        I --> J
        J["para$Deco_done <- 1"]
        J --> K
        K["Render tables in Decompressor_ui"]
    end

    G --> L
    subgraph "DecompressSRA.R operations"
        L["Load packages & parse WD"]
        L --> M
        M["Load DecompressSRA_server_para.RData"]
        M --> N
        N["ReArrangeFiles(WD)"]
        N --> N1A["Read GSE_SRR_List.csv"]
        subgraph "ReArrangeFiles.R"
            N1A
            N1A --> N1B
            N1B["Query SRA layout via NCBI"]
            N1B --> N1C
            N1C["Generate prefix & file names"]
            N1C --> N1D
            N1D["Write FileNameMap.csv"]
            N1D --> N1E
            N1E["Write AlignerInput.txt"]
        end
        N --> O
        O["Handle TenX_bam samples\nrename files in parallel"]
        O --> P
        P["Write *_HistoryOrigin_Log.txt"]
        P --> Q
        Q["Create HistoryOrigin_Log.txt\nfor other samples"]
        Q --> R
        R["Run 02.fasterq_dump_gzip.sh"]
        R --> R1A["Parse FileNameMap.csv arrays"]
        subgraph "02.fasterq_dump_gzip.sh"
            R1A
            R1A --> R1B
            R1B["For each SRR_ID"]
            R1B --> R1C
            R1C["Prefetch if missing\nand run fasterq-dump"]
            R1C --> R1D
            R1D["Compress FASTQ with pigz"]
            R1D --> R1E
            R1E["Move outputs to final names"]
        end
        R --> S
        S["Rename scRNA-seq FASTQ files\nusing vdb-dump"]
        S --> T
        T["Print progress Done"]
    end

    T --> I
```

## Launching the app
Run `00.launcher.sh` from this directory to set up the Conda environment and start the Shiny server. The application opens in your default web browser.

## Demo video
A short demonstration of the GREP1 workflow can be seen on YouTube:

<a href="https://www.youtube.com/watch?v=YdIe83-7Yr8" target="_blank" rel="noopener noreferrer">
  <img src="https://img.youtube.com/vi/YdIe83-7Yr8/hqdefault.jpg" alt="YouTube Thumbnail">
</a>
---

If you found this helpful, feel free to comment, share, and follow for more. Your support encourages us to keep creating quality content.