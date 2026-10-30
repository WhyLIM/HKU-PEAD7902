#!/usr/bin/env python3
"""Step 3 - per-residue evolutionary conservation from the Clustal Omega alignment.

For every alignment column the script computes
  * Shannon entropy  H = -sum p_i log2 p_i   (0 = perfectly conserved, 4.32 = uniform)
  * conservation score  = 1 - H / log2(20)
  * Jensen-Shannon divergence from the alignment background composition
  * fraction of the other five species identical to the human residue
and maps each column onto the human p53 residue numbering (1-393).

Outputs:
  results/tables/conservation.csv
  results/tables/pairwise_identity.csv
  results/tables/alignment_columns.csv
"""
import json
import math
import pathlib

import pandas as pd
from Bio import AlignIO

ROOT = pathlib.Path(__file__).resolve().parents[1]
ALN = ROOT / "results" / "alignment" / "tp53_clustalo.aln"
TBL = ROOT / "results" / "tables"
TBL.mkdir(parents=True, exist_ok=True)

HUMAN = "sp|P04637|P53_HUMAN"
AA = "ACDEFGHIKLMNPQRSTVWY"

DOMAINS = [
    (row.domain, int(row.start), int(row.end))
    for row in pd.read_csv(ROOT / "data" / "domain_map.csv").itertuples()
]

DBD_START, DBD_END = next((start, end) for name, start, end in DOMAINS if name == "DBD")


def domain_of(pos):
    for name, start, end in DOMAINS:
        if start <= pos <= end:
            return name
    return "?"


def main():
    aln = AlignIO.read(ALN, "clustal")
    ids = [r.id for r in aln]
    seqs = {r.id: str(r.seq) for r in aln}
    human = seqs[HUMAN]
    ncol = aln.get_alignment_length()

    # background amino-acid composition of the alignment
    counts = {a: 0 for a in AA}
    for s in seqs.values():
        for ch in s:
            if ch in counts:
                counts[ch] += 1
    total = sum(counts.values())
    bg = {a: counts[a] / total for a in AA}

    columns = []
    human_pos = 0
    for col in range(ncol):
        chars = [seqs[i][col] for i in ids]
        residues = [c for c in chars if c in AA]
        if human[col] in AA:
            human_pos += 1
            pos = human_pos
        else:
            pos = None

        n = len(residues)
        gap_frac = 1 - n / len(ids)
        if n:
            freq = {a: residues.count(a) / n for a in AA}
            entropy = -sum(p * math.log2(p) for p in freq.values() if p > 0)
            # Jensen-Shannon divergence (base 2) between column and background
            m = {a: 0.5 * (freq[a] + bg[a]) for a in AA}
            jsd = 0.5 * sum(freq[a] * math.log2(freq[a] / m[a]) for a in AA if freq[a] > 0)
            jsd += 0.5 * sum(bg[a] * math.log2(bg[a] / m[a]) for a in AA if bg[a] > 0)
            ident = ([c for i, c in zip(ids, chars) if i != HUMAN and c in AA].count(human[col])
                     / max(1, sum(1 for i, c in zip(ids, chars) if i != HUMAN and c in AA)))
        else:
            entropy, jsd, ident = float("nan"), float("nan"), float("nan")

        columns.append(dict(
            alignment_column=col + 1, human_pos=pos, aa=human[col],
            n_residues=n, gap_fraction=round(gap_frac, 3),
            entropy=entropy, jsd=jsd,
            conservation=(1 - entropy / math.log2(20)) if n else float("nan"),
            identity_to_human=ident,
        ))

    df = pd.DataFrame(columns)
    dfc = df[df.human_pos.notna()].copy()
    dfc["human_pos"] = dfc["human_pos"].astype(int)
    dfc["domain"] = dfc["human_pos"].map(domain_of)
    dfc.to_csv(TBL / "conservation.csv", index=False)
    df.to_csv(TBL / "alignment_columns.csv", index=False)

    # pairwise identity over the whole alignment and over the DBD only
    rows = []
    dbd_cols = [c for c in range(ncol) if df.iloc[c].human_pos and DBD_START <= df.iloc[c].human_pos <= DBD_END]
    for other in ids:
        if other == HUMAN:
            continue
        pairs_all = [(human[c], seqs[other][c]) for c in range(ncol)
                     if human[c] in AA and seqs[other][c] in AA]
        pairs_dbd = [(human[c], seqs[other][c]) for c in dbd_cols
                     if human[c] in AA and seqs[other][c] in AA]
        rows.append(dict(
            species=other.split("|")[2],
            identity_full=round(100 * sum(a == b for a, b in pairs_all) / len(pairs_all), 1),
            identity_dbd=round(100 * sum(a == b for a, b in pairs_dbd) / len(pairs_dbd), 1),
            aligned_length=len(pairs_all),
        ))
    pd.DataFrame(rows).to_csv(TBL / "pairwise_identity.csv", index=False)

    summary = dict(
        alignment_length=int(ncol),
        human_length=int(dfc.shape[0]),
        mean_conservation=float(dfc.conservation.mean()),
        mean_conservation_dbd=float(dfc[dfc.domain == "DBD"].conservation.mean()),
        mean_conservation_outside_dbd=float(dfc[dfc.domain != "DBD"].conservation.mean()),
        conserved_positions_ge_0_8=int((dfc.conservation >= 0.8).sum()),
    )
    (ROOT / "results" / "conservation_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(dfc.groupby("domain").conservation.agg(["mean", "count"]).round(3))


if __name__ == "__main__":
    main()
