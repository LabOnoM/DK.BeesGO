---
title: "BRAP: Bulk RNA-seq Auto-Pipeline"
draft: true
tags:
  - Evolution
  - Genetics
  - RNA
  - Splicing
---
# Reference:
[https://pubmed.ncbi.nlm.nih.gov/39966351/](https://pubmed.ncbi.nlm.nih.gov/39966351/)
# Summary
This study investigates the physiological effects of a human-specific amino acid substitution (I197V) in the NOVA1 protein, which is essential for neural development and function. The researchers generated humanized mice carrying the I197V variant and compared them to wild-type mice to analyze the molecular and behavioral consequences. While the substitution had minimal impact on NOVA1's RNA binding capacity, it led to specific effects on alternative splicing and vocalization patterns in the humanized mice. The findings suggest that this human-specific NOVA1 variant may have been part of an ancient evolutionary selective sweep in a common ancestral population of Homo sapiens, potentially contributing to the development of spoken language through differential RNA regulation during brain development.

# Key Points:
1. NOVA1 is a neuronal RNA-binding protein essential for survival and development in mice and humans
2. A single amino acid change (I197V) in NOVA1 is unique to modern humans
3. Humanized mice carrying the I197V variant showed specific effects on alternative splicing and vocalization patterns compared to wild-type mice
4. The I197V substitution likely underwent a strong evolutionary selective sweep in the common ancestral population of modern humans
5. The findings suggest a potential role for the human-specific NOVA1 variant in the evolution of human language and communication

# Logic Flow

```mermaid
flowchart TB

%% SECTION: EVOLUTIONARY CONTEXT
subgraph "🔬 Evolutionary Context"
A1["Modern humans possess unique vocal abilities"]
A2["NOVA1 has a human-specific I197V variant"]
A3["Unknown functional impact of I197V"]
A1 --> A2 --> A3
end

%% SECTION: RESEARCH GOAL
subgraph "🎯 Research Question"
B1["Does NOVA1 I197V affect brain function or vocal behavior?"]
B2["Could it play a role in evolution of human language?"]
A3 --> B1 --> B2
end

%% SECTION: METHODS
subgraph "🧪 Experimental Design"
C1["CRISPR-engineered humanized Nova1^hu/hu mice"]
C2["Compared with Nova1^wt/wt controls"]
C3["Multi-layer analysis:\n• RNA-seq\n• CLIP\n• Behavior"]
B2 --> C1 --> C2 --> C3
end

%% SECTION: MOLECULAR FINDINGS (with branching)
subgraph "🧬 Molecular Observations"
C3 --> D1["NOVA1 I197V has minimal effect on RNA binding"]
D1 --> D2a["UCAU binding motifs conserved (CLIP)"]
D1 --> D2b["Genome-wide CLIP peaks largely unchanged"]
D1 --> D2c["In vitro Kd values show no binding difference"]
D1 --> D2d["Suggests no loss of RNA-binding function"]
D1 --> D3["But... some peaks show small differences → Explore splicing"]
end

%% SECTION: EXPRESSION & LOCALIZATION
subgraph "🧠 Gene Expression"
E1["NOVA1 highly expressed in midbrain (esp. PAG & Amb)"]
E2["Protein levels & spatial expression unchanged"]
D2d --> E1 --> E2
end

%% SECTION: ALTERNATIVE SPLICING (with branching)
subgraph "🧩 Alternative Splicing"
D3 --> F1["720 alternative splicing events altered in Nova1^hu/hu"]
F1 --> F2a["41% of spliced transcripts have NOVA1 CLIP peaks"]
F1 --> F2b["Splicing enrichment in:\n• Synaptic genes\n• Chromatin regulators"]
F1 --> F2c["27 behavior-related genes altered"]
F2c --> F3["4 vocalization-related genes differentially spliced:\nAuts2, Nrxn2, Myh14, Srpx2"]
end

%% SECTION: BEHAVIORAL OBSERVATIONS - PUPS
subgraph "🐭 Pup Vocalization (P7)"
G1["USVs recorded during isolation"]
G1 --> G2a["Total call number similar"]
G1 --> G2b["Changes in pitch jump syllable ratios"]
G1 --> G2c["Higher frequency and altered complexity"]
F3 --> G1
end

%% SECTION: BEHAVIORAL OBSERVATIONS - ADULTS
subgraph "🐭 Adult Vocalization"
H1["Courtship USVs recorded in adult males"]
H1 --> H2a["Total call number unchanged"]
H1 --> H2b["Altered pitch in long 's' syllables"]
H1 --> H2c["Increased variance in high-frequency pitch jumps"]
G2c --> H1
end

%% SECTION: EVOLUTIONARY ANALYSIS (with branches)
subgraph "📈 Evolutionary Significance"
A2 --> I1["I197V nearly fixed in all modern humans"]
I1 --> I2a["Tajima’s D = -2.48 → strong purifying selection"]
I1 --> I2b["CLUES2 & ARGweaver: selection coefficient S ≈ 19"]
I1 --> I2c["Supports ancient selective sweep (pre-Out-of-Africa)"]
end

%% SECTION: CONCLUSIONS
subgraph "🏁 Conclusions"
J1["NOVA1 I197V does not impair core RNA-binding"]
J2["But subtly alters splicing of brain-expressed genes"]
J3["Especially impacts vocalization-related splicing"]
J4["Leads to altered vocal behaviors in both pups and adults"]
J5["Suggests potential evolutionary contribution to human language development"]
H2c --> J1 --> J2 --> J3 --> J4 --> J5
I2c --> J5
end
```

