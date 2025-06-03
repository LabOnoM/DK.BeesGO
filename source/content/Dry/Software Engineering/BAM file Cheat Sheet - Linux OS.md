---
title: BAM file Cheat Sheet - Linux OS
tags:
  - High-throughput_Sequencing
  - FileFormat
  - Linux
  - Console
---
Check tag in bam file:

```bash
samtools view -f 4 A02598A4.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam | head
```

- `-f 4` selects reads **flagged as unmapped**
- If this returns reads → your BAM **includes unmapped reads**
- If nothing returns → unmapped reads were **filtered out before or during BAM generation**