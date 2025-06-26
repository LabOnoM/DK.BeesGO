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
# Workflow overview

## Step 1: GET SRR ID list 

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


---
## Step 2: Download .sra or TenX scRNA-seq BAM files

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



---

# Running Demo in Video

A running demo is shown in a YouTube Video below:

<a href="https://www.youtube.com/watch?v=YdIe83-7Yr8" target="_blank" rel="noopener noreferrer">
  <img src="https://img.youtube.com/vi/YdIe83-7Yr8/hqdefault.jpg" alt="YouTube Thumbnail">
</a>

---

If you found this helpful, feel free to comment, share, and follow for more. Your support encourages us to keep creating quality content.