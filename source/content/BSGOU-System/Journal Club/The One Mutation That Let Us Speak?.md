---
title: "BRAP: Bulk RNA-seq Auto-Pipeline"
draft: false
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

%% SECTION 1: EVOLUTIONARY BACKGROUND
subgraph "🔬 Evolutionary Background"
A1["Modern humans have unique vocal abilities"]
A2["NOVA1 has human-specific I197V variant"]
A3["Functional impact of I197V unknown"]
A1 --> A2 --> A3
end

  

%% SECTION 2: RESEARCH OBJECTIVE
subgraph "🎯 Research Objective"
B1["Hypothesis: I197V affects neural RNA regulation and vocalization"]
B2["Goal: Understand evolutionary role of NOVA1 I197V in speech-related traits"]
A3 --> B1 --> B2
end

  

%% SECTION 3: METHODS
subgraph "🧪 Experimental Design"
C1["CRISPR knock-in → $$Nova1^{hu/hu}$$ mice (I197V)"]
C2["Compare with $$Nova1^{wt/wt}$$ mice"]
C3["Multi-modal analysis:<br>• Gene expression<br>• RNA-binding (CLIP)<br>• Splicing<br>• Behavior"]
B2 --> C1 --> C2 --> C3
end

  

%% SECTION 4: RNA-BINDING RESULTS (MOLECULAR)
subgraph "🧬 NOVA1 RNA Binding Characteristics"
D1["I197V shows minimal effect on RNA binding"]
D1 --> D2a["✅ UCAU motifs enriched in both genotypes (CLIP)"]
D1 --> D2b["📍 Peak locations largely unchanged across brain"]
D1 --> D2c["🧪 Gel-shift: Binding affinity (Kd) unchanged"]
D1 --> D2d["→ Suggests I197V preserves RNA-binding function"]
end

C3 --> D1

%% SECTION 5: GENE EXPRESSION & LOCALIZATION
subgraph "🧠 Expression & Localization"
E1["NOVA1 highly expressed in midbrain (PAG, Amb)"]
E2["Expression levels & spatial patterns unchanged"]
D2b --> E1 --> E2
end

  

%% SECTION 6: ALTERNATIVE SPLICING
subgraph "🧩 Alternative Splicing Regulation"
F1["720 differential splicing events in $$Nova1^{hu/hu}$$"]
F1 --> F2a["📎 41% have NOVA1 binding peaks"]
F1 --> F2b["🔧 Enrichment in synapse, morphogenesis, chromatin processes"]
F1 --> F2c["🧠 27 behavior-related genes altered"]
F2c --> F3["📣 4 vocalization-related genes altered:<br>Auts2, Nrxn2, Myh14, Srpx2"]
end

D2d --> F1

%% SECTION 7: BEHAVIOR - PUP VOCALIZATION
subgraph "🐭 Pup Vocalization (P7)"
G1["USVs recorded at postnatal day 7"]
G1 --> G2a["🔢 Total call number unchanged"]
G1 --> G2b["🎶 Syllable composition shifted<br>→ fewer pitch jumps"]
G1 --> G2c["📈 Higher proportion of high-frequency calls"]
end

F3 --> G1

%% SECTION 8: BEHAVIOR - ADULT VOCALIZATION
subgraph "🐭 Adult Vocalization"
H1["Courtship USVs recorded from adult males"]
H1 --> H2a["⏳ Call count similar"]
H1 --> H2b["📉 Long 's' syllables have lower frequency"]
H1 --> H2c["🎛️ More variance in pitch jumps"]
end

G2c --> H1

%% SECTION 9: EVOLUTIONARY GENETICS
subgraph "📈 Evolutionary Analysis"
I1["I197V nearly fixed in all modern humans"]
I2["Tajima’s D = –2.48 → strong purifying selection"]
I3["CLUES2 + ARGweaver: selection coefficient S ≈ 19"]
I4["**Conclusion**: Ancient selective sweep before human population split"]
end

B2 --> I1 --> I2 --> I3 --> I4

%% SECTION 10: CONCLUSIONS
subgraph "🏁 **Conclusion**"
J1["✔️ NOVA1 I197V retains RNA-binding function"]
J2["🧩 Alters specific splicing programs"]
J3["🎤 Modifies vocal behavior in pups & adults"]
J4["🧬 Supports possible evolutionary role in development of vocal communication"]
end

H2c --> J1 --> J2 --> J3 --> J4

I4 --> J4
```


