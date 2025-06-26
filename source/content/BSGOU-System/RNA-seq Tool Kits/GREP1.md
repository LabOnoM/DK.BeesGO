---
title: GEO Retrieval & Extraction Pipeline in One shiny-app (GREP1)
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
    A[User clicks GET ID list / GETSRRIDls] --> B{Working directory set?}
    B -- No --> Z[Show notification: folder required]
    B -- Yes --> C["Reset Step1 state (if any)"]
    C --> D
    D[Save GSE IDs / params to Downloader_server_para.RData] --> E
    E[Kill port 4778 and old process] --> F
    F[Launch RetrieveGSEinfo.R via processx / Rscript] --> G

    subgraph RetrieveGSEinfo.R operations
        G[Load packages & parse WD] --> H[Load GSE_IDls from RData]
        H --> I{Loop over GSE IDs}
        I --> J["Fetch GSM list from GEO (lines 77–117)"]
        J --> K{Large GSM list?}
        K -- ">1000" --> L["Parallel foreach (lines 138–191)"]
        K -- "≤1000" --> M["Sequential loop (lines 192–259)"]
        L --> N[Collect SRR IDs]
        M --> N
        N --> O["Use RSelenium to scrape run metadata (lines 262–343)"]
        O --> P[Determine DataType using ENA API]
        P --> Q["Write 00.GSE_SRR_List.csv & Downloader_server_rout.RData (lines 348–395)"]
    end

    Q --> R[Monitor progress / errors]
    R --> S[Load Downloader_server_rout.RData]
    S --> T[Set values$out_tb and para$Step1_done]
    T --> U[Render Step1 UI with DataTable]

```


---

# Running Demo in Video

A running demo is shown in a YouTube Video below:

[![Real-Time Console Monitor in R ShinyApp? + GEO FASTQ Downloader | GREP1 | Devlog #1](https://img.youtube.com/vi/YdIe83-7Yr8/hqdefault.jpg)](https://youtu.be/YdIe83-7Yr8?si=K1hZyBiYviPoPzWx)

---

If you found this helpful, feel free to comment, share, and follow for more. Your support encourages us to keep creating quality content.