---
title: CARLIN Lineage Analysis on Stereo-seq
draft: false
tags:
  - Dry
  - Python
  - Linux
  - Console
  - Spatial
  - RNA
  - High-throughput_Sequencing
---
Following with the [[Extract Tissue CIDs from Stereo-seq]], this workflow enables the extraction and analysis of CARLIN barcode amplicons from **unmapped Stereo-seq reads**, which are filtered and cleaned by SAW but not aligned to the reference genome (as CARLIN sequences are synthetic).

## 🧭 **1. Overview of the Problem**

- SAW (`saw count`) processes Stereo-seq data using the STAR aligner.
- CARLIN barcodes are **non-genomic**, so STAR cannot align them.
- CARLIN reads are therefore retained in the `--unmapped-fastq` output.
- However, SAW outputs only **R2** for unmapped reads; **R1 (UMI)** must be recovered manually.

## 📥 **2. Generate Unmapped Reads via SAW**

```bash
saw count \
  --unmapped-fastq \
  --clean-reads-fastq \
  ... (other required arguments)
```

✅ This produces a lot of `*_1.pure_unmapped_reads.fq` files. Just same as the output from [--clean-reads-fastq](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/analysis/pipelines/saw-commands#saw-count) option (see [[Extract Tissue CIDs from Stereo-seq#2. The [--clean-reads-fastq](https //stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/analysis/pipelines/saw-commands saw-count) option in `saw` software]]), although the fq was named as `*_1*`, it actually only contains the biological sequences from R2. 

To reduce the computing cost, we only do the analysis on the reads located within the tissue. Therefore, we need to filter and combine the `--unmapped-fastq` resulted files into a single `R2.tissue.unmapped.fastq` by using the `tissue_xy_coords.txt` from [[Extract Tissue CIDs from Stereo-seq#1. What's inside the `*.tissue.gef`?]] with `filter_cleanR1_by_XY_parallel.sh` below:
```bash
#!/bin/bash

# ===== CONFIG =====
XY_LIST="barcodeToPos_tissue.txt"
OUTFILE="R2.tissue.unmapped.fastq"
TMPDIR="tmp_cleanR1"
mkdir -p "$TMPDIR"

echo "📄 Preparing whitelist..."
# Extract and sort unique XY coordinates (Cx, Cy)
cut -f2,3 "$XY_LIST" | sort | uniq > whitelist_xy.txt

# Export for GNU parallel
export TMPDIR
export WL="$(realpath whitelist_xy.txt)"  # Use full path to avoid parallel context issues

# ===== FUNCTION: PROCESS EACH FILE =====
process_cleanR1() {
    fq="$1"
    base=$(basename "$fq" _1.pure_unmapped_reads.fq)
    tmpout="$TMPDIR/${base}.tissue.fq"

    gawk -v WL="$WL" -v OUT="$tmpout" '
        BEGIN {
            while ((getline line < WL) > 0) {
                gsub(/\r/, "", line)  # Remove Windows-style line endings
                valid_xy[line] = 1
            }
        }
        {
            if (substr($0,1,1) == "@") {
                id = $0
                getline seq
                getline plus
                getline qual
                if (match(id, /\|Cx:i:([0-9]+)\|Cy:i:([0-9]+)/, m)) {
                    xy = m[1] "\t" m[2]
                    if (xy in valid_xy) {
                        print id >> OUT
                        print seq >> OUT
                        print plus >> OUT
                        print qual >> OUT
                    }
                }
            }
        }
    ' "$fq"
}
export -f process_cleanR1

# ===== PARALLEL PROCESSING =====
echo "🚀 Filtering reads based on XY coordinates..."
ls *_1.pure_unmapped_reads.fq | parallel --bar -j "$(nproc)" process_cleanR1 {}

# ===== MERGE RESULTS =====
echo "🧬 Merging into $OUTFILE..."
if ls $TMPDIR/*.tissue.fq 1> /dev/null 2>&1; then
    cat $TMPDIR/*.tissue.fq > "$OUTFILE"
    echo "✅ Done. Output written to $OUTFILE"
else
    echo "⚠️  No matching reads found. Output file was not created."
fi
```

## 🔁 **3. Recover Corresponding R1 Reads**

Since SAW doesn’t output R1 unmapped reads, we:

- Extract the **read IDs** from `R2.tissue.unmapped.fastq`
- Match and recover corresponding R1 reads from the **raw `.fq.gz` files**

Python script
```python
import os
import re
import subprocess
import sqlite3
from tqdm import tqdm
import multiprocessing

r2_fastq = "/mnt/md0/22_Pam/MyCustomBank/reads/R2.tissue.unmapped.fastq"
r1_files = sorted([f for f in os.listdir(".") if f.endswith("_1.fq.gz")])
db_path = "r1_reads.db"
output_file = "R1.tissue.unmapped.fastq"

# === Step 1: Load ordered R2 IDs ===
r2_order = []
r2_set = set()
print("📥 Extracting ordered R2 read IDs...")
with open(r2_fastq, "r") as f:
    while True:
        header = f.readline()
        _ = f.readline(); _ = f.readline(); _ = f.readline()
        if not header:
            break
        match = re.match(r"@([^| ]+)", header)
        if match:
            rid = match.group(1)
            r2_order.append(rid)
            r2_set.add(rid)

print(f"✅ Loaded {len(r2_order)} R2 read IDs.")

# === Step 2: Worker function just returns matched reads ===
def collect_r1_reads(r1_file):
    result = []
    cmd = ["zcat", r1_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    while True:
        header = proc.stdout.readline()
        seq = proc.stdout.readline()
        plus = proc.stdout.readline()
        qual = proc.stdout.readline()
        if not qual:
            break
        match = re.match(r"@([^/| ]+)", header)
        if match:
            key = match.group(1)
            if key in r2_set:
                result.append((key, seq, plus, qual))
    return result

# === Step 3: Collect from workers, write to SQLite in main ===
print(f"🚀 Collecting matched R1 reads in parallel (no DB writes yet)...")
matched_records = []
with multiprocessing.Pool(processes=min(8, multiprocessing.cpu_count())) as pool:
    for records in tqdm(pool.imap_unordered(collect_r1_reads, r1_files), total=len(r1_files)):
        matched_records.extend(records)

# === Step 4: Write all matched reads to SQLite (sequentially) ===
print(f"💾 Writing {len(matched_records)} matched R1 reads to {db_path}")
if os.path.exists(db_path):
    os.remove(db_path)

conn = sqlite3.connect(db_path)
cur = conn.cursor()
cur.execute("CREATE TABLE reads (id TEXT PRIMARY KEY, seq TEXT, plus TEXT, qual TEXT)")
cur.executemany("INSERT OR IGNORE INTO reads VALUES (?, ?, ?, ?)", matched_records)
conn.commit()
conn.close()

# === Step 5: Write R1 output in R2 order ===
print(f"📝 Writing final R1 FASTQ to {output_file} in R2 order...")
with sqlite3.connect(db_path) as conn, open(output_file, "w") as out:
    cur = conn.cursor()
    written = 0
    for rid in tqdm(r2_order, desc="Writing ordered R1"):
        cur.execute("SELECT seq, plus, qual FROM reads WHERE id = ?", (rid,))
        row = cur.fetchone()
        if row:
            seq, plus, qual = row
            out.write(f"@{rid}\n{seq}{plus}{qual}")
            written += 1

print(f"✅ Done! Total ordered R1 reads written: {written}")

```
Output:
```bash
(base) ubuntu4@ubuntu4:/mnt/md0/22_Pam/Stereoseq/A02598A4/00.Rawdata/reads$ python3 ./R1_Recover.py
📥 Extracting ordered R2 read IDs...
✅ Loaded 32121315 R2 read IDs.
🚀 Collecting matched R1 reads in parallel (no DB writes yet)...
100%|████████████████████████████████████████████████████████████████████| 16/16 [09:42<00:00, 36.41s/it]
💾 Writing 32121315 matched R1 reads to r1_reads.db
📝 Writing final R1 FASTQ to R1.tissue.unmapped.fastq in R2 order...
Writing ordered R1: 100%|████████████████████████████████| 32121315/32121315 [04:37<00:00, 115636.00it/s]
✅ Done! Total ordered R1 reads written: 32121315

```

**Final output:**  
`R1.tissue.unmapped.fastq` with headers exactly matching `R2.tissue.unmapped.fastq`

## ⚙️ **4. Define a Custom CARLIN Configuration**

**`CustomCfg.json`:**
```json
{
  "type": "Bulk",
  "UMI": {
    "length": 26,
    "location": "L"
  },
  "read_perspective": {
    "ShouldComplement": "N",
    "ShouldReverse": "Y"
  },

  "trim": {
    "Primer5": "ignore",
    "Primer3": "ignore",
    "SecondarySequence": "ignore"
  }
}
```

Properties of reads from Stereo-seq:
- ❗**CARLIN barcodes were captured passively by the polyA-based reverse transcription in the Stereo-seq platform** — not through a targeted PCR with CARLIN-specific primers. So:
	- There are **no defined 5′ or 3′ primers** in the reads.
	- The CARLIN region may appear **anywhere within the R2 sequence**, and
	- Its flanking context is likely **non-deterministic**.
- CID+UMI is 26 bp at the 5′ end of R1
- R2 reads are in reverse orientation relative to the CARLIN amplicon reference

Please check [https://gitlab.com/hormozlab/carlin](https://gitlab.com/hormozlab/carlin) for more details about the **`CustomCfg.json`** options.

## 🚀 **5. Run the CARLIN Analysis Pipeline**
In MatLab:
```matlab
analyze_CARLIN({'./reads/R1.tissue.unmapped.fastq'; './reads/R2.tissue.unmapped.fastq'}, './CustomCfg.json', './Output');
cd("./Output/");
load('./Summary.mat');
load('./Bank.mat');
p_clonal = bank.compute_clonal_pvalue(summary);
writematrix(p_clonal,'P_Con.txt');
p_frequency = bank.compute_frequency_pvalue(summary);
writematrix(p_frequency,'P_Freq.txt');
plot_allele_frequency_CDF(summary, 'A02598A4');
plot_indel_freq_vs_length(summary);
plot_site_decomposition(summary, true, 'A02598A4', '# of Transcripts');
plot_stargate.create(summary);
```

## Reference:
1. [SAW User Manual V8.1](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1)
2. [From CID to ATGC: Decoding Stereo-seq Barcodes](https://www.bs-gou.com/2025/06/08/how-to-encode-barcode-in-stereo-seq.html)
3. [[Extract Tissue CIDs from Stereo-seq]]

## 🧾 TL;DR — Extract CARLIN Amplicons from Unmapped Stereo-seq Reads
This workflow enables extraction and analysis of synthetic CARLIN barcodes from unmapped reads produced by SAW (`saw count`) in Stereo-seq. Because CARLIN sequences are not aligned to the genome, they reside in `--unmapped-fastq` outputs. The workflow filters these unmapped R2 reads to those located within the tissue, then recovers corresponding R1 reads (which contain UMI/CID information). The paired R1/R2 fastqs are then processed using the CARLIN Matlab pipeline with a custom configuration.

### 📄 Required Files
| 📁 File                       | 🧾 Description                          |
|------------------------------|----------------------------------------|
| `*_1.pure_unmapped_reads.fq` | R2-only unmapped reads from SAW        |
| `barcodeToPos_tissue.txt`    | List of XY coordinates within tissue   |
| `*_1.fq.gz`                  | Raw R1 fastq files (compressed)        |
| `CustomCfg.json`             | Custom CARLIN config for analysis      |

### 🛠 Required Scripts
| ⚙️ Script                   | 🔧 Purpose                                     |
|----------------------------|-----------------------------------------------|
| `filter_cleanR1_by_XY_parallel.sh` | Filters R2 reads based on tissue XYs    |
| `R1_Recover.py`            | Recovers matching R1 reads from raw fastqs    |
| `analyze_CARLIN` (Matlab)  | Runs CARLIN barcode detection and stats       |

### 🧭 Workflow Summary
1. Run `saw count` with `--unmapped-fastq` to extract unmapped R2 reads.
2. Filter these R2 reads to tissue-relevant reads using `barcodeToPos_tissue.txt`.
3. Merge filtered reads into `R2.tissue.unmapped.fastq`.
4. Use `R1_Recover.py` to extract matching R1 reads from raw gzipped fastqs.
5. Save result as `R1.tissue.unmapped.fastq` with correct read order.
6. Analyze the paired R1/R2 fastqs using the CARLIN Matlab pipeline with `CustomCfg.json`.

### 📂 Output Files
| 📤 Output                    | 📌 Content                                         |
|-----------------------------|---------------------------------------------------|
| `R2.tissue.unmapped.fastq`  | Filtered R2 reads from within the tissue          |
| `R1.tissue.unmapped.fastq`  | Recovered and ordered R1 reads with UMI/CID info |
| `Summary.mat`, `Bank.mat`   | Output from CARLIN barcode calling                |
| `P_Con.txt`, `P_Freq.txt`   | Clonal and frequency p-values from CARLIN         |
| `*.png` or figures          | Diagnostic plots from CARLIN analysis             |

### 🧬 Workflow Diagram
~~~mermaid
flowchart TD
  subgraph Input_Files
    A1["📁 *_1.pure_unmapped_reads.fq"]
    A2["📁 barcodeToPos_tissue.txt"]
    A3["📁 *_1.fq.gz"]
    A4["📁 CustomCfg.json"]
  end

  subgraph Processing_Steps
    B1["⚙️ filter_cleanR1_by_XY_parallel.sh"]
    B2["📄 R2.tissue.unmapped.fastq"]
    B3["⚙️ R1_Recover.py"]
    B4["📄 R1.tissue.unmapped.fastq"]
  end

  subgraph Outputs
    C1["⚙️ analyze_CARLIN"]
    C2["📤 Summary.mat"]
    C3["📤 Bank.mat"]
    C4["📤 P_Con.txt / P_Freq.txt"]
    C5["📤 *.png plots"]
  end

  A1 --> B1
  A2 --> B1
  B1 --> B2
  B2 --> B3
  A3 --> B3
  B3 --> B4
  B4 --> C1
  B2 --> C1
  A4 --> C1
  C1 --> C2
  C1 --> C3
  C1 --> C4
  C1 --> C5
~~~

---

If you found this helpful, feel free to comment, share, and follow for more. Your support encourages us to keep creating quality content.