---
bibliography: references.bib
csl: nature-spaced.csl
link-citations: true
---

%% Title Evolution and Cancer Mutations in TP53

%% Subtitle PEAD7902 — Assignment 1

%% Author Min Li (3036744858)

%% Date October 2026

# Summary

TP53 is the most frequently mutated tumour suppressor. We aligned six vertebrate p53 orthologues, quantified residue conservation, and mapped 4,225 somatic TP53 mutations from 3,839 tumours onto the alignment and the p53 core-domain structure. The DNA-binding domain was the most conserved region (mean score 0.88 versus 0.72 elsewhere). Five codons (R273, R248, R175, G245 and R282) carried 19.8% of all mutations; all five were invariant across the six species and located at the DNA interface (R248, R273) or in the buried core (R175, G245, R282). Conservation and mutation density were positively correlated (Spearman ρ = 0.51). Cancer mutations therefore concentrate on evolutionarily invariant residues, whereas conservation alone does not determine hotspot rank.

# Key points

- The DBD is the conserved core: 57% of its residues are invariant versus 14% elsewhere.
- R273, R248, R175, G245 and R282 carry 19.8% of all TP53 mutations.
- All five are invariant across vertebrates.
- Two mechanisms recur: DNA-contact loss (R248, R273) and core destabilisation (R175, G245, R282).
- Conservation predicts functional importance, not mutation frequency.

# Results

## Conservation across six vertebrate orthologues

Six p53 orthologues (human, cattle, mouse, chicken, frog, zebrafish; UniProt [@uniprot2025]) were aligned with Clustal Omega [@sievers2011; @madeira2022] (Figure 1). Identity to human falls from ~80% in mammals to ~50% in chicken, frog and zebrafish, yet the DNA-binding domain (DBD, 102–292) [@cho1994; @joerger2008; @kastenhuber2017] stays 68–91% identical (Figure 2, Table 1); mean conservation is 0.88 in the DBD versus 0.72 elsewhere.

![](../../assets/r/fig1_alignment.png){width=15cm}

%% Image Caption **Figure 1.** Multiple sequence alignment of six p53 orthologues. Cells are shaded by column conservation score (darker = more conserved) and warm-grey cells are alignment gaps. The DBD is visibly darker than the terminal regions. Six vertebrate orthologues were aligned with Clustal Omega.

![](../../assets/r/fig2_cross_species.png){width=15cm}

%% Image Caption **Figure 2.** Conservation across species and domains. (a) Clustal Omega guide tree with branch lengths. (b) Identity to human p53 over the full length and within the DBD. (c) Conservation score per domain (boxes: median and quartiles; dots: individual residues).

%% Table Caption **Table 1. Sequences used and identity to human p53.**

| Species | UniProt | Length (aa) | Identity, full length | Identity, DBD (102–292) |
|---|---|---|---|---|
| *Homo sapiens* (human) | P04637 | 393 | 100.0% | 100.0% |
| *Bos taurus* (cattle) | P67939 | 386 | 81.6% | 91.1% |
| *Mus musculus* (mouse) | P02340 | 390 | 78.0% | 88.5% |
| *Gallus gallus* (chicken) | P10360 | 367 | 54.9% | 71.7% |
| *Xenopus laevis* (frog) | P07193 | 363 | 52.9% | 68.1% |
| *Danio rerio* (zebrafish) | P79734 | 373 | 50.4% | 71.2% |

## Five mutation hotspots

Among 4,225 somatic TP53 mutations from 3,839 tumours (32 TCGA Pan-Cancer Atlas cohorts [@hoadley2018], cBioPortal [@cerami2012; @debruijn2023]; 65% missense), the five most frequently mutated missense codons (R273, R248, R175, G245 and R282) carry 19.8% of all and 30.7% of missense mutations (Figure 3, Table 2). Each mutates to few alternatives: R175H (90%), R282W (96%), R248Q/W (54/40%), R273C/H (47/38%), G245S (46%).

![](../../assets/r/fig3_hotspots.png){width=15cm}

%% Image Caption **Figure 3.** The five hotspots and their substitution spectra. (a) Matrix of hotspot codons (columns) against mutant residues (rows, grouped by residue class); cell shading and the number give the count of each substitution, and the column labels give the codon total and its percentage of all TP53 mutations. (b) Consequence classes of all 4,225 TP53 mutations.

%% Table Caption **Table 2. Five major TP53 hotspots (TCGA Pan-Cancer Atlas, cBioPortal).**

| Codon | Main substitutions (% of codon) | n | % of all TP53 mutations | Conservation | Structural role | Structural evidence (1TUP) |
|---|---|---|---|---|---|---|
| R273 | C 47, H 38, L 10 | 267 | 6.3% | 1.00 | DNA contact | 2.5 Å to DNA |
| R248 | Q 54, W 40, L 4 | 224 | 5.3% | 1.00 | DNA contact | 2.3 Å to DNA |
| R175 | H 90, G 5, C 4 | 167 | 4.0% | 1.00 | structural (buried) | 2% accessible surface |
| G245 | S 46, D 21, V 19, C 11 | 91 | 2.2% | 1.00 | structural (buried) | 0% accessible surface |
| R282 | W 96, G 3 | 89 | 2.1% | 1.00 | structural (buried) | 2% accessible surface |

## Hotspots coincide with invariant residues

All five are invariant across the six vertebrates and lie in the DBD; five random codons would all be invariant with *P* = 0.005. Mutation density follows conservation: 13.8 missense events per DBD residue versus 0.4 outside (Figures 4, 5). In 1TUP [@cho1994], R248 and R273 contact the DNA backbone (2.3–2.5 Å), whereas R175, G245 and R282 are buried (0–2% accessible surface) in the core (Figure 6).

![](../../assets/r/fig4_landscape.png){width=15cm}

%% Image Caption **Figure 4.** Conservation and mutation density along p53. Upper panel: per-residue conservation with the DBD shaded. Lower panel: missense mutations per codon (log scale); hotspots in red. Domains are drawn below.

![](../../assets/r/fig5_conservation_vs_mutation.png){width=15cm}

%% Image Caption **Figure 5.** Conservation predicts where mutations accumulate. (a) Distribution of conservation scores; all five hotspots fall in the invariant class (137 residues). (b) Conservation versus missense density (Spearman ρ = 0.51, *P* = 3×10⁻²⁷); each point is one codon, grey for residues outside the DBD, teal for DBD residues and red for the five hotspot codons. (c) Mutations per residue in conserved versus variable halves and inside versus outside the DBD.

![](../../assets/r/fig6_structure.png){width=15cm}

%% Image Caption **Figure 6.** Structural context of the hotspots. (a) DNA-interface view of the p53 core domain (PDB 1TUP) with the DNA duplex; sticks mark the DNA-contact hotspots R248 and R273 (red). (b) Rotated view with the DNA removed, showing the buried structural hotspots R175, G245 and R282 (blue); Zn²⁺ is the purple sphere. (c) Relative solvent accessibility of each hotspot, which separates the buried core residues from the DNA-contacting ones.

# Discussion

## What conservation reveals about disease

Conservation flags residues that cannot change without cost, and tumours select them: purifying selection across vertebrates and positive selection in tumours converge on one interface and one folding core. Yet 137 residues are invariant and only a few are hotspots, so conservation alone does not determine which invariant residues become hotspots.

## How large databases enable discovery of functional residues

Large databases turn observations into statistical evidence. cBioPortal, used here, aggregates thousands of tumours, letting recurrence at one codon be separated from passenger background; COSMIC adds curated somatic calls [@sondka2024], and the IARC TP53 database [@bouaoun2016] adds germline variants and functional assays, which can provide independent cross-checks. In this dataset, mutations cluster at a few codons, these codons are invariant across species, and the affected residues contact DNA or stabilise the core. Cross-referencing frequency, conservation, structure and assay data therefore identifies functionally important residues and separates drivers from passengers. Limitations remain: cohort composition and sequencing depth bias tumour representation, and indels are under-called.

# Methods

Six reviewed UniProtKB entries (P04637, P02340, P67939, P10360, P07193, P79734) were downloaded through the UniProt REST API [@uniprot2025] and aligned with Clustal Omega 1.2.4 through the EMBL-EBI REST service [@sievers2011; @madeira2022], producing a 409-column alignment.

Conservation was quantified for each alignment column as the Shannon entropy [@shannon1948], with score = 1 − H/log₂20; Jensen–Shannon divergence [@lin1991] and identity to human were also computed, and columns were mapped onto human residue numbering.

All TP53 (Entrez 7157) mutations from the 32 TCGA Pan-Cancer Atlas studies were retrieved through the cBioPortal REST API [@cerami2012; @debruijn2023], yielding 4,225 unique events in 3,839 samples; hotspots were ranked by missense count.

Enrichment among invariant codons was tested with a hypergeometric test, hotspot and non-hotspot conservation scores were compared with a Mann–Whitney U test, and the relation between conservation and mutation density was assessed with Spearman correlation.

The core domain with the most DNA contacts was selected from PDB entry 1TUP [@cho1994]; Shrake–Rupley solvent accessibility [@shrake1973] was computed in Biopython, with residue maxima from Tien et al. [@tien2013].

Domain boundaries follow a UniProt-anchored segmentation: DBD 102–292 and TET 325–356 are taken from UniProt [@uniprot2025]; TAD1 (1–40) and TAD2 (41–60) follow the conventional schematic, PRD (61–101) spans the proline-rich region and the adjacent segment up to the DBD start, and linker (293–324) and CTD (357–393) fill the remaining regions. The 1TUP core domain [@cho1994] is the structural counterpart of the DBD. Analyses used Python 3.9.6 with Biopython 1.85 and NumPy 2.0.2, and R 4.6.1 with ggplot2 4.0.3 and patchwork 1.3.2; Figure 6a,b were rendered with open-source PyMOL 3.2.0.2. The cBioPortal API was accessed on 25 October 2026.

# Data Availability

All data analysed here are publicly available: protein sequences from UniProtKB (P04637, P02340, P67939, P10360, P07193, P79734); somatic mutations from the cBioPortal REST API (TCGA Pan-Cancer Atlas, 32 studies; accessed 25 October 2026); and the structure from the PDB (1TUP). Derived data (alignment, per-residue conservation, codon mutation counts, hotspot and structure metrics, and the source data for Figures 5 and 6) are available in the GitHub repository https://github.com/WhyLIM/HKU-PEAD7902. The cBioPortal raw JSON cache (about 3.4 MB) is not included; run scripts/04_fetch_mutations.py to re-download it.

# Code Availability

Analysis code is available in the GitHub repository https://github.com/WhyLIM/HKU-PEAD7902. Full reproduction steps and a file guide are in the repository README.

# Declarations

Generative artificial intelligence tools were used in the process, including but not limited to experimental design, code writing and computational methods. All steps were reviewed and verified manually by the author.

%% References Title References
