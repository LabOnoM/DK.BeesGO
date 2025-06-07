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
```Python
import h5py
import numpy as np

def explore_gef(h5file):
	def visitor(name, obj):
		if isinstance(obj, h5py.Dataset):
		print(f"{name}: shape={obj.shape}, dtype={obj.dtype}")
	with h5py.File(h5file, 'r') as f:
		f.visititems(visitor)
		
gef_file = "./A02598A4.tissue.gef"

with h5py.File('./A02598A4.tissue.gef', 'r') as f:
	def visit(name):
		print(name)
	f.visit(visit)
```
>contour
>contour/tissueContour
>geneExp
>geneExp/bin1
>geneExp/bin1/exon
>geneExp/bin1/expression
>geneExp/bin1/gene
>geneExp/bin10
>geneExp/bin10/exon
>geneExp/bin10/expression
>geneExp/bin10/gene
>geneExp/bin100
>geneExp/bin100/exon
>geneExp/bin100/expression
>geneExp/bin100/gene
>geneExp/bin150
>geneExp/bin150/exon
>geneExp/bin150/expression
>geneExp/bin150/gene
>geneExp/bin20
>geneExp/bin20/exon
>geneExp/bin20/expression
>geneExp/bin20/gene
>geneExp/bin200
>geneExp/bin200/exon
>geneExp/bin200/expression
>geneExp/bin200/gene
>geneExp/bin5
>geneExp/bin5/exon
>geneExp/bin5/expression
>geneExp/bin5/gene
>geneExp/bin50
>geneExp/bin50/exon
>geneExp/bin50/expression
>geneExp/bin50/gene
>stat
>stat/gene
>wholeExp
>wholeExp/bin1
>wholeExp/bin10
>wholeExp/bin100
>wholeExp/bin150
>wholeExp/bin20
>wholeExp/bin200
>wholeExp/bin5
>wholeExp/bin50
>wholeExpExon
>wholeExpExon/bin1
>wholeExpExon/bin10
>wholeExpExon/bin100
>wholeExpExon/bin150
>wholeExpExon/bin20
>wholeExpExon/bin200
>wholeExpExon/bin5
>wholeExpExon/bin50

```Python
explore_gef(gef_file)
```
>contour/tissueContour: shape=(), dtype=object
>geneExp/bin1/exon: shape=(116772117,), dtype=uint8
>geneExp/bin1/expression: shape=(116772117,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', 'u1')]
>geneExp/bin1/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin10/exon: shape=(93972549,), dtype=uint16
>geneExp/bin10/expression: shape=(93972549,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
>geneExp/bin10/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin100/exon: shape=(32675685,), dtype=uint16
>geneExp/bin100/expression: shape=(32675685,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
>geneExp/bin100/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin150/exon: shape=(22060836,), dtype=uint16
>geneExp/bin150/expression: shape=(22060836,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
>geneExp/bin150/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin20/exon: shape=(78271410,), dtype=uint16
>geneExp/bin20/expression: shape=(78271410,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
>geneExp/bin20/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin200/exon: shape=(15879595,), dtype=uint16
>geneExp/bin200/expression: shape=(15879595,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
>geneExp/bin200/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin5/exon: shape=(103779528,), dtype=uint8
>geneExp/bin5/expression: shape=(103779528,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', 'u1')]
>geneExp/bin5/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>geneExp/bin50/exon: shape=(52649679,), dtype=uint16
>geneExp/bin50/expression: shape=(52649679,), dtype=[('x', '<i4'), ('y', '<i4'), ('count', '<u2')]
vgeneExp/bin50/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('offset', '<u4'), ('count', '<u4')]
>stat/gene: shape=(57584,), dtype=[('geneID', 'S64'), ('geneName', 'S64'), ('MIDcount', '<u4'), ('E10', '<f4')]
>wholeExp/bin1: shape=(26460, 26460), dtype=[('MIDcount', 'u1'), ('genecount', '<u2')]
>wholeExp/bin10: shape=(2646, 2646), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
>wholeExp/bin100: shape=(265, 265), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
>wholeExp/bin150: shape=(177, 177), dtype=[('MIDcount', '<u4'), ('genecount', '<u2')]
>wholeExp/bin20: shape=(1323, 1323), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
>wholeExp/bin200: shape=(133, 133), dtype=[('MIDcount', '<u4'), ('genecount', '<u2')]
>wholeExp/bin5: shape=(5292, 5292), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
>wholeExp/bin50: shape=(530, 530), dtype=[('MIDcount', '<u2'), ('genecount', '<u2')]
>wholeExpExon/bin1: shape=(26460, 26460), dtype=uint8
>wholeExpExon/bin10: shape=(2646, 2646), dtype=uint16
>wholeExpExon/bin100: shape=(265, 265), dtype=uint16
>wholeExpExon/bin150: shape=(177, 177), dtype=uint32
>wholeExpExon/bin20: shape=(1323, 1323), dtype=uint16
>wholeExpExon/bin200: shape=(133, 133), dtype=uint32
>wholeExpExon/bin5: shape=(5292, 5292), dtype=uint16
>wholeExpExon/bin50: shape=(530, 530), dtype=uint16

```Python
with h5py.File(gef_file, 'r') as f:
	expression = f['geneExp/bin1/expression']
	x_coords = expression['x']
	y_coords = expression['y']
```

```Python
xy_set = set(zip(x_coords, y_coords))

# Save to tab-separated list
with open('tissue_xy_coords.txt', 'w') as out:
	for x, y in sorted(xy_set):
		out.write(f"{x}\t{y}\n")
```

```bash
cut -f1 A02598A4.barcodeToPos_tissue.txt > barcodes_in_tissue.txt
```
