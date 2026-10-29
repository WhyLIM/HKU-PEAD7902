# TP53 conservation and cancer hotspots

Coursework for **HKU PEAD7902 (Assignment 1)**. The analysis asks whether the codons that mutate most often in human cancer are the residues that evolution refuses to change.

Six vertebrate p53 orthologues (human, cattle, mouse, chicken, frog, zebrafish) were aligned with Clustal Omega. Per-column conservation (Shannon entropy) was mapped onto human p53, and 4,225 somatic TP53 mutations from 3,839 tumours (32 TCGA Pan-Cancer Atlas studies, retrieved through cBioPortal) were counted per codon. Five codons — R273, R248, R175, G245 and R282 — carry 19.8% of all mutations and all five are invariant across the six species. The structure of the p53 core domain (PDB 1TUP) shows that R248 and R273 contact DNA, whereas R175, G245 and R282 are buried in the folding core.

## Repository layout

```
scripts/          analysis, figure and report-export code (01-10)
data/
  sequences/      six UniProt FASTA files
  mutations/      merged, de-duplicated mutation table (TSV)
  structure/      1TUP.pdb
results/
  alignment/      Clustal Omega alignment, guide tree, job id
  tables/         conservation, codon counts, hotspots, structure metrics
  *.json          summary statistics
assets/r/         the six report figures (PNG)
references/       RIS / EndNote XML / Markdown reference list
TP53_Assignment_Report.md      English report
TP53_Assignment_Report_ZH.md   Chinese report
requirements.txt
```

## Data sources

| Data | Source | Access |
|---|---|---|
| Protein sequences | UniProtKB (P04637, P02340, P67939, P10360, P07193, P79734) | REST API |
| Somatic mutations | cBioPortal, TCGA Pan-Cancer Atlas (32 studies) | REST API, accessed 2026-10-25 |
| Structure | PDB 1TUP | RCSB PDB |

The cBioPortal raw JSON cache (~3.4 MB) is **not** included; `scripts/04_fetch_mutations.py` downloads it again (the API is rate-limited, so the script retries and caches).

## Requirements

- Python 3.9.6 with `requirements.txt` (Biopython, pandas, NumPy, SciPy, requests)
- R 4.6.1 with ggplot2, patchwork, ape, dplyr, tidyr, readr, scales, ggrepel, jsonlite, ragg and png
- Optional: open-source PyMOL 3.2.0.2, only to re-render Figure 6a,b; the pre-rendered PNGs in `results/pymol/` are included

## Reproducing the analysis

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

.venv/bin/python scripts/01_fetch_sequences.py   # UniProt FASTA
.venv/bin/python scripts/02_run_msa.py           # Clustal Omega via EBI
.venv/bin/python scripts/03_conservation.py      # per-column conservation
.venv/bin/python scripts/04_fetch_mutations.py   # cBioPortal (cached)
.venv/bin/python scripts/05_hotspots.py          # hotspots + statistics
.venv/bin/python scripts/06_structure_analysis.py  # 1TUP distances and SASA
PYOPENGL_PLATFORM=osmesa .venv-pymol/bin/python scripts/07_render_fig6_pymol.py
Rscript scripts/08_make_figures_r.R              # six figures in assets/r/
```

Figure 6a,b are rendered before the R figure script because panel c of Figure 6 reads `results/pymol/fig6_view_*.png`. The pre-rendered PNGs are included, so the PyMOL step is optional for a quick rebuild.

Report utilities (not needed for the figures):

```bash
.venv/bin/python scripts/09_make_references.py      # references/references.ris / EndNote XML
.venv/bin/python scripts/10_build_word_sources.py   # Markdown sources for md-to-word
```

## Outputs

| Figure | File | Content |
|---|---|---|
| 1 | `assets/r/fig1_alignment.png` | six-species alignment coloured by conservation |
| 2 | `assets/r/fig2_cross_species.png` | guide tree, identity, per-domain conservation |
| 3 | `assets/r/fig3_hotspots.png` | hotspot substitution matrix and consequence classes |
| 4 | `assets/r/fig4_landscape.png` | conservation and mutation density along p53 |
| 5 | `assets/r/fig5_conservation_vs_mutation.png` | conservation versus mutation density |
| 6 | `assets/r/fig6_structure.png` | hotspot environment in 1TUP and solvent accessibility |

## Not included

- `data/mutations/raw/` (cBioPortal JSON cache; re-download with script 04)
- `results/pymol/*.pse` (PyMOL sessions)
- virtual environments and package caches

The Word exports are included under `word/`; the figure PDFs are in `assets/r/pdf/`.

## Notes

- The mutation table is de-duplicated by sample and protein change; hotspots are ranked by missense count.
- Conservation is `1 - H / log2(20)`; invariant columns score 1.0.
- Domain boundaries use UniProt DBD 102–292 and TET 325–356; the remaining boundaries follow the conventional schematic.
- The reports use 15 DOI-verified references (see `references/`).