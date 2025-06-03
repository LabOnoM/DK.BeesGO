---
title: FASTQ files Cheat Sheet - Linux OS
tags:
  - High-throughput_Sequencing
  - FileFormat
  - Linux
  - bash
  - Console
---
Get read sequence by searching the read ID:

```bash
zgrep -A1 'E150021049L1C006R010900881884' E150021049_L01_57_1.fq.gz
```

Check few lines of fastq.gz:

```bash
gunzip -c my_fastq.fastq.gz | head -n 10
```