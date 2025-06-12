---
title: CARLIN Lineage Analysis from Stereo-seq Data
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

### 🧭 **1. Overview of the Problem**

- SAW (`saw count`) processes Stereo-seq data using the STAR aligner.
- CARLIN barcodes are **non-genomic**, so STAR cannot align them.
- CARLIN reads are therefore retained in the `--unmapped-fastq` output.
- However, SAW outputs only **R2** for unmapped reads; **R1 (UMI)** must be recovered manually.

### 📥 **2. Generate Unmapped Reads via SAW**

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
import multiprocessing
from tqdm import tqdm

r2_fastq = "/mnt/md0/22_Pam/MyCustomBank/reads/R2.tissue.unmapped.fastq"
output_file = "R1.tissue.unmapped.fastq"

# === Step 1: Build read_key → full_header mapping from R2 ===
def load_r2_read_map():
    read_map = {}
    with open(r2_fastq, "r") as f:
        while True:
            header = f.readline()
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                break
            match = re.match(r"@([^| ]+)", header)
            if match:
                read_map[match.group(1)] = header.strip()
    return read_map

read_map = load_r2_read_map()
print(f"✅ Loaded {len(read_map)} read IDs from R2")

# === Step 2: Worker function to process a single file ===
def process_r1_file(file_path_and_map):
    file_path, read_map = file_path_and_map
    matched = 0
    buffer = []
    filename = os.path.basename(file_path)

    cmd = ["zcat", file_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=8192)

    pbar = tqdm(desc=filename, unit="read", position=0, leave=False)
    while True:
        header = proc.stdout.readline()
        seq = proc.stdout.readline()
        plus = proc.stdout.readline()
        qual = proc.stdout.readline()
        if not qual:
            break
        pbar.update(1)

        match = re.match(r"@([^/| ]+)", header)
        if match:
            key = match.group(1)
            if key in read_map:
                new_header = read_map[key] + "\n"
                buffer.append(new_header)
                buffer.append(seq)
                buffer.append(plus)
                buffer.append(qual)
                matched += 1
    pbar.close()
    return matched, buffer

# === Step 3: Writer function to collect and write results ===
def writer(queue, total_files):
    with open(output_file, "w") as out:
        completed = 0
        pbar = tqdm(total=total_files, desc="📝 Writing output", position=0)
        while completed < total_files:
            result = queue.get()
            if result is None:
                completed += 1
            else:
                out.writelines(result)
                pbar.update(0)
        pbar.close()

# === Step 4: Multiprocessing orchestration ===
if __name__ == "__main__":
    manager = multiprocessing.Manager()
    read_map_shared = manager.dict(read_map)
    queue = multiprocessing.Queue()
    r1_files = sorted([f for f in os.listdir(".") if f.endswith("_1.fq.gz")])
    total_files = len(r1_files)

    writer_proc = multiprocessing.Process(target=writer, args=(queue, total_files))
    writer_proc.start()

    def wrapped(file):
        matched, buffer = process_r1_file((file, read_map))
        queue.put(buffer)
        queue.put(None)
        return matched

    print("🚀 Processing R1 files in parallel with processes...\n")
    with multiprocessing.Pool(processes=min(multiprocessing.cpu_count(), 16)) as pool:
        results = list(pool.map(wrapped, r1_files))

    writer_proc.join()
    print(f"\n✅ Finished. Total R1 reads matched and written: {sum(results)}")

```

**Final output:**  
`R1.tissue.unmapped.fastq` with headers exactly matching `R2.tissue.unmapped.fastq`

## ⚙️ **4. Define a Custom CARLIN Configuration**

**`CustomCfg.json`:**
```json
{
  "type": "Bulk",
  "UMI.length": 26,
  "UMI.location": "L",
  "read_perspective": {
    "ShouldComplement": "N",
    "ShouldReverse": "Y"
  },
  "trim": {
    "Primer5": "exact",
    "Primer3": "malformed",
    "SecondarySequence": "ignore"
  }
}

```

Based on the observations:
- CID+UMI is 26 bp at the 5′ end of R1
- R2 reads are in reverse orientation relative to the CARLIN amplicon reference

Please check [https://gitlab.com/hormozlab/carlin](https://gitlab.com/hormozlab/carlin) for more details about the **`CustomCfg.json`** options.

🚀 **5. Run the CARLIN Analysis Pipeline**
In MatLab:
```matlab
analyze_CARLIN({'./reads/R1.tissue.unmapped.fastq', './reads/R2.tissue.unmapped.fastq'}, CFG_TYPE='./CustomCfg.json, './Output');
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