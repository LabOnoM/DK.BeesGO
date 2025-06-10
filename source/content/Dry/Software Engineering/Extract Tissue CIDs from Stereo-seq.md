---
title: Extract Tissue CIDs from Stereo-seq
tags:
  - Dry
  - Python
  - Linux
  - Console
  - Spatial
  - RNA
  - High-throughput_Sequencing
---
## 1. What's inside the `*.tissue.gef`?

Let's first explore the resulted `<SN>.tissue.gef` in Python:
```python
import h5py
import numpy as np

def explore_gef(h5file):
	def visitor(name, obj):
		if isinstance(obj, h5py.Dataset):
		print(f"{name}: shape={obj.shape}, dtype={obj.dtype}")
	with h5py.File(h5file, 'r') as f:
		f.visititems(visitor)
		
gef_file = "./<SN>.tissue.gef"

with h5py.File('./<SN>.tissue.gef', 'r') as f:
	def visit(name):
		print(name)
	f.visit(visit)
```
Output:
```bash
contour
contour/tissueContour
geneExp
geneExp/bin1
geneExp/bin1/exon
geneExp/bin1/expression
geneExp/bin1/gene
geneExp/bin10
geneExp/bin10/exon
geneExp/bin10/expression
geneExp/bin10/gene
geneExp/bin100
geneExp/bin100/exon
geneExp/bin100/expression
geneExp/bin100/gene
geneExp/bin150
geneExp/bin150/exon
geneExp/bin150/expression
geneExp/bin150/gene
geneExp/bin20
geneExp/bin20/exon
geneExp/bin20/expression
geneExp/bin20/gene
geneExp/bin200
geneExp/bin200/exon
geneExp/bin200/expression
geneExp/bin200/gene
geneExp/bin5
geneExp/bin5/exon
geneExp/bin5/expression
geneExp/bin5/gene
geneExp/bin50
geneExp/bin50/exon
geneExp/bin50/expression
geneExp/bin50/gene
stat
stat/gene
wholeExp
wholeExp/bin1
wholeExp/bin10
wholeExp/bin100
wholeExp/bin150
wholeExp/bin20
wholeExp/bin200
wholeExp/bin5
wholeExp/bin50
wholeExpExon
wholeExpExon/bin1
wholeExpExon/bin10
wholeExpExon/bin100
wholeExpExon/bin150
wholeExpExon/bin20
wholeExpExon/bin200
wholeExpExon/bin5
wholeExpExon/bin50
```

Let's check where `<SN>.tissue.gef` file stores the XY coordinates continuous in Python:
```python
explore_gef(gef_file)
```
Output:
```bash
contour/tissueContour: shape=(), dtype=object
geneExp/bin1/exon: shape=(116772117,), dtype=uint8
geneExp/bin1/expression: shape=(116772117,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', 'u1')]
geneExp/bin1/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin10/exon: shape=(93972549,), dtype=uint16
geneExp/bin10/expression: shape=(93972549,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
geneExp/bin10/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin100/exon: shape=(32675685,), dtype=uint16
geneExp/bin100/expression: shape=(32675685,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
geneExp/bin100/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin150/exon: shape=(22060836,), dtype=uint16
geneExp/bin150/expression: shape=(22060836,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
geneExp/bin150/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin20/exon: shape=(78271410,), dtype=uint16
geneExp/bin20/expression: shape=(78271410,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
geneExp/bin20/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin200/exon: shape=(15879595,), dtype=uint16
geneExp/bin200/expression: shape=(15879595,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
geneExp/bin200/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin5/exon: shape=(103779528,), dtype=uint8
geneExp/bin5/expression: shape=(103779528,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', 'u1')]
geneExp/bin5/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
geneExp/bin50/exon: shape=(52649679,), dtype=uint16
geneExp/bin50/expression: shape=(52649679,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
geneExp/bin50/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
stat/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('MIDcount', '<u4'), ('E10', '<f4')]
wholeExp/bin1: shape=(26460, 26460), dtype=[('MIDcount', 'u1'), ('genecount', '<u2')]
wholeExp/bin10: shape=(2646, 2646), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
wholeExp/bin100: shape=(265, 265), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
wholeExp/bin150: shape=(177, 177), dtype=[('MIDcount', '<u4'), ('genecount', '<u2')]
wholeExp/bin20: shape=(1323, 1323), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
wholeExp/bin200: shape=(133, 133), dtype=[('MIDcount', '<u4'), ('genecount', '<u2')]
wholeExp/bin5: shape=(5292, 5292), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
wholeExp/bin50: shape=(530, 530), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
wholeExpExon/bin1: shape=(26460, 26460), dtype=uint8
wholeExpExon/bin10: shape=(2646, 2646), dtype=uint16
wholeExpExon/bin100: shape=(265, 265), dtype=uint16
wholeExpExon/bin150: shape=(177, 177), dtype=uint32
wholeExpExon/bin20: shape=(1323, 1323), dtype=uint16
wholeExpExon/bin200: shape=(133, 133), dtype=uint32
wholeExpExon/bin5: shape=(5292, 5292), dtype=uint16
wholeExpExon/bin50: shape=(530, 530), dtype=uint16
```

Now, we know that the `<SN>.tissue.gef` file stores the XY coordinates in `geneExp/bin1/expression` on the smallest level. Therefore, we will save them out into `tissue_xy_coords.txt` for later use by using Python:
```python
with h5py.File(gef_file, 'r') as f:
	expression = f['geneExp/bin1/expression']
	x_coords = expression['x']
	y_coords = expression['y']

xy_set = set(zip(x_coords, y_coords))

# Save to tab-separated list
with open('tissue_xy_coords.txt', 'w') as out:
	for x, y in sorted(xy_set):
		out.write(f"{x}\t{y}\n")
```

Check the output `tissue_xy_coords.txt`:
```bash
ubuntu4@ubuntu4:~$ head tissue_xy_coords.txt
5713	26441
5787	26358
5788	26318
5788	26319
5788	26321
5788	26329
5789	26316
5789	26327
5790	26305
5790	26310
```

## 2. The  [--clean-reads-fastq](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/analysis/pipelines/saw-commands#saw-count) option in `saw` software

Actually, with the [--clean-reads-fastq](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/analysis/pipelines/saw-commands#saw-count) option in saw software, we can get the read ID along with the XY coordinates and the biological sequence originally from Read2, as shown below:
```bash
ubuntu@ubuntu:~$ head <ID>_<LaneNo.>_<SampleNo.>_1.clean_reads.fq
@E150018299L1C004R03402098317|Cx:i:19892|Cy:i:17727 E8BBE63524CF BEAB3
GCCAATAGTATAGTGTGGTGTGCTTTTACGTGATGGCGAGTGGGCAGCGGGCGGTGGGCTGTACACAGCCGTCTGTCCTTTGAATCTCAATCTGCCTGCG
+
9A>CCCCACCCC>CACB?CCC>ACCCCCA>C@<AB8<ABB?>B;6BC7?/8=B>ACA=>CC@BAB<BC=CCCACABB<CCCC1CCAC8C/A;CCCCC@CA
@E150018299L1C004R03402098337|Cx:i:13026|Cy:i:14469 1C20554DC029 A5B9F
GCGTGGCACTCAGGCGGGGCCCTGGGAGCGCTGCGGGCACGGGGTGGCCGGCAGGACGCGGGCTGGATGGCTCTGGCCGCGCCAGGAGGAGGCCGACCTG
+
1B1AC>ABC<CCCBB8<C0>>1C;A95C:>9>@869<2C?@C;C1C?A=ACA@CC6BC=CA>6CAA:BCC?BCCCA<1@CCA:@CA1<A)CC93-CCCCC
@E150018299L1C004R03402098484|Cx:i:17752|Cy:i:16736 8E7513DF377D FD27F
GCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAA
```

Therefore, we can use the `tissue_xy_coords.txt` to filter the `*_1.clean_reads.fq` files directly into a merged `R2.tissue.fq.gz` file by using the `filter_cleanR1_by_XY_parallel.sh` below:

```bash
#!/bin/bash

# ===== CONFIG =====
XY_LIST="tissue_xy_coords.txt"
OUTFILE="R2.tissue.fq.gz"
TMPDIR="tmp_cleanR1"
mkdir -p "$TMPDIR"

echo "📄 Preparing whitelist..."
# Extract and sort unique XY coordinates (Cx, Cy)
cut -f1,2 "$XY_LIST" | sort | uniq > whitelist_xy.txt

# Export for GNU parallel
export TMPDIR
export WL="$(realpath whitelist_xy.txt)"  # Use full path to avoid parallel context issues

# ===== FUNCTION: PROCESS EACH FILE =====
process_cleanR1() {
    fq="$1"
    base=$(basename "$fq" .clean_reads.fq)
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
ls *_1.clean_reads.fq | parallel --bar -j "$(nproc)" process_cleanR1 {}

# ===== MERGE RESULTS =====
echo "🧬 Merging into $OUTFILE..."
if ls $TMPDIR/*.tissue.fq 1> /dev/null 2>&1; then
    cat $TMPDIR/*.tissue.fq | gzip > "$OUTFILE"
    echo "✅ Done. Output written to $OUTFILE"
else
    echo "⚠️  No matching reads found. Output file was not created."
fi

```
Output:
```bash
ubuntu4@ubuntu4:~$ zcat R2.tissue.fq.gz|head -n 12
@E150018299L1C001R00100000906/2
GTCTTAGGAAGACAATGTGAAGTTGGAGGTCGAGAGACCTCTGATAAGGTCGCCATGCCTCTCAGTACGTCAGCAGAAGCAGTGGTATCAACGCAGAGTA
+
DDDDDCCDDDDDDDDDCDDDDDDDDDCDCDDDCDCDCDDDDDDCDDCCCCDDDDCDDDDDCDDCDDDDDDDCCDCDCCDDCDCDDDCDDCCDDDCDCDDB
@E150018299L1C001R00100001311/2
CCTGATCAGTCTTAGGAAGACAAAGAAAAGTGGTGCCTAAGGCCTCACCTGATAAGGTCGCCATGCCTCTCAGTACGTCAGCAGAAGCAGTGGTATCAAC
+
BCDDCDCCCDCDDDCBAACCCDC@CCCCBDCDCCDCCDCCCBCCDBDBCDCCCCCDCCCDBBCCDCDDCCB>CDCBCCDCBCCC=<CCCCCBBACCC?BC
@E150018299L1C001R00100002321/2
ACAAAGCACCGAGGAAAAAAGTAATGAGTCTGATAAGGTCGCCATGCCTCTCAGTACGTCAGCAGAAGCAGTGGTATCAACGCAGAGTACATTCCGTTGA
+
CDDDCCDCCCDDDDDCCDCDDCCCDCCDCCDDCDCCDDDDDCCDDDCCDDDCDDDCDDDDCDDCDCCDDCCDDDDDDDCDDDDCDCDCCCDDCDCCDCDD
```

To here, we got what we want.

## 3. If you need un-clean reads

According to [SAW User Manual V8.1](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1): 

![Workflow-of-CID-mapping](https://raw.githubusercontent.com/LabOnoM/DK.BeesGO/master/source/content/00.Images/Workflow-of-CID-mapping.png)

### 3.1. [CID mapping](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/algorithms/gene-expression-algorithms#cid-mapping)
CID mapping requires FASTQs and a chip mask file, recording position information for sequencing reads.

Check the amount of Ns in Coordinate IDs:
- If there is 1 N base in CID, the N base will be replaced by `A/T/C/G`.
- CIDs without N bases will be directly poured into the match pool.
- CIDs with more than one N base will be discarded.

In the absense of an N base, if a CID does not match to any positions, each base of the CID will replaced by the other three types iterately until a successful match. After the mapping, only the unique CID match will be retained for subsequent steps.

### 3.2. [RNA filtering](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/algorithms/gene-expression-algorithms#rna-filtering)

Before genome alignment, it is necessary to confirm that the reads entering the next step are cDNA sequences. Ideally, this part of the data only contains cDNA fragments, but there may be cases where the fragmented cDNA is too short or some non-cDNA fragments are detected.

Reads will be discarded if any of the following conditions are triggered:

- with a length of less than 30 after cutting out adapter sequences,
- mapped to DNB sequences,
- with a length of less than 30 after cutting out the poly-A sequence.

The above three are collectively known as "non-relevant short reads".

### 3.3. [MID filtering](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/algorithms/gene-expression-algorithms#mid-filtering)

MID sequences will be filtered out if they match any of the following:

- having more than one base with quality <= Q10,
- having N bases >= 1.

In summary, the [--clean-reads-fastq](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1/analysis/pipelines/saw-commands) option lets the SAW software output the Clean Reads (before RNA alignment) in FASTQ format, which have undergone **CID mapping**, **RNA filtering**, and **MID filtering**. 

Therefore, if you also need the reads before **CID mapping**, **RNA filtering**, and **MID filtering** (unclean reads), you can get if by following the below steps: 

### 3.4. Extract Raw Reads by Barcodes

As we mentioned in in blog of [From CID to ATGC: Decoding Stereo-seq Barcodes](https://www.bs-gou.com/2025/06/08/how-to-encode-barcode-in-stereo-seq.html), we can use the [STOmics/ST_BarcodeMap](https://github.com/wong-ziyi/ST_BarcodeMap) to convert the CID numbers into the actual ATGC sequences, as shown below:
```bash
./ST_BarcodeMap-0.0.1 --in A02598A4.barcodeToPos.h5 --out barcodeToPos.txt --action 3
```
Output:
```bash
ubuntu@ubuntu:~$ head barcodeToPos.txt
TGCTCTATCGCAACCCATGCTCCAG	26457	26459
GGATTACCACGTGTCCATTTACCCC	26456	26459
AATTGAAGCCGACACTGTATGGGGG	26453	26459
GCATATCGACCTCCTTCGATCCCTG	26447	26459
GTCGATGTGCTTGCAATCATGAGCT	26445	26459
AAGGTTTGCTTGAGCACGGGGCCGA	26444	26459
GAGCCTATGAACTTTCGCTCCGAAA	26443	26459
AAACTTAACAGCCGCCTGTCCTTTC	26442	26459
GCTCTCGGCGTCATCCATTTACTAC	26441	26459
GTCGTTTGTCCACGATAGCAACATT	26437	26459
```



```bash
awk 'NR==FNR {xy[$1"\t"$2]=1; next} ($2"\t"$3) in xy' tissue_xy_coords.txt ./A02598A4/00.Rawdata/mask/A02598A4.barcodeToPos.txt > A02598A4.barcodeToPos_tissue.txt

```

```bash
./ST_BarcodeMap-0.0.1 --in A02598A4.barcodeToPos.h5 --out barcodes.txt --action 3

awk 'NR==FNR {xy[$1"\t"$2]=1; next} ($2"\t"$3) in xy' tissue_xy_coords.txt ./A02598A4/00.Rawdata/mask/A02598A4.barcodeToPos.txt > A02598A4.barcodeToPos_tissue.txt
cut -f1 A02598A4.barcodeToPos_tissue.txt > barcodes_in_tissue.txt
```

## Reference:
1. [SAW User Manual V8.1](https://stereotoolss-organization.gitbook.io/saw-user-manual-v8.1)
2. [From CID to ATGC: Decoding Stereo-seq Barcodes](https://www.bs-gou.com/2025/06/08/how-to-encode-barcode-in-stereo-seq.html)