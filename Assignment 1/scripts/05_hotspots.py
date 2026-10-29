#!/usr/bin/env python3
"""Step 5 - identify TP53 missense mutation hotspots and relate them to conservation.

Inputs
  data/mutations/tp53_mutations_pancan.tsv   (step 4)
  results/tables/conservation.csv            (step 3)

Outputs
  results/tables/hotspots.csv                top-5 hotspot codons + conservation
  results/tables/hotspot_substitutions.csv   substitution spectrum at each hotspot
  results/tables/codon_mutation_counts.csv   per-codon mutation counts
  results/hotspot_summary.json               headline numbers used in the report
"""
import json
import math
import pathlib
import re

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
TBL = ROOT / "results" / "tables"

PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z*])$")


def hypergeom_sf(k, N, K, n):
    """P(X >= k) for X ~ Hypergeometric(N, K, n)."""
    total = math.comb(N, n)
    return sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(k, min(K, n) + 1)) / total


def main():
    mut = pd.read_csv(ROOT / "data" / "mutations" / "tp53_mutations_pancan.tsv", sep="\t")
    mut["proteinChange"] = mut["proteinChange"].fillna("")
    parsed = mut["proteinChange"].str.extract(PATTERN)
    parsed.columns = ["wt", "pos", "mut"]
    mut = pd.concat([mut, parsed], axis=1)
    n_total = len(mut)
    n_samples = mut["sampleId"].nunique()
    n_studies = mut["studyId"].nunique()
    n_parsed = int(mut["pos"].notna().sum())

    mut["pos"] = pd.to_numeric(mut["pos"], errors="coerce")
    coding = mut[mut["pos"].notna()].copy()
    coding["pos"] = coding["pos"].astype(int)
    missense = coding[(coding["mut"].notna()) & (coding["mut"] != "*") & (coding["mut"] != coding["wt"])]

    counts = (missense.groupby(["pos", "wt", "mut"]).size().rename("n").reset_index()
              .sort_values("n", ascending=False))
    per_codon = (missense.groupby(["pos", "wt"]).size().rename("n_missense").reset_index()
                 .sort_values("n_missense", ascending=False))

    conserv = pd.read_csv(TBL / "conservation.csv")[["human_pos", "conservation", "domain"]]
    per_codon = per_codon.merge(conserv, left_on="pos", right_on="human_pos", how="left")

    # ---- five hotspots: most frequently mutated codons (missense) -------------
    top5 = per_codon.head(5).copy()
    top5["rank"] = range(1, 6)
    top5["pct_of_all_mutations"] = 100 * top5["n_missense"] / n_total
    top5["pct_of_missense"] = 100 * top5["n_missense"] / len(missense)
    # percentile rank with ties counted as fully conserved (max statistic)
    top5["conservation_pctile"] = [
        float((conserv["conservation"] <= c).mean() * 100) for c in top5["conservation"]
    ]

    subs = []
    for _, row in top5.iterrows():
        s = counts[counts["pos"] == row["pos"]].copy()
        s["pct_of_codon"] = 100 * s["n"] / row["n_missense"]
        subs.append(s)
    subs = pd.concat(subs)
    subs["hotspot"] = subs["pos"].map(lambda p: f"{top5.set_index('pos').loc[p, 'wt']}{p}")
    subs.to_csv(TBL / "hotspot_substitutions.csv", index=False)

    top5.to_csv(TBL / "hotspots.csv", index=False)
    per_codon.to_csv(TBL / "codon_mutation_counts.csv", index=False)

    # ---- conservation vs mutation density ------------------------------------
    # background set: every codon with an annotated wild-type residue and a
    # conservation score; mutation counts are filled with zeros where unobserved.
    background = conserv.copy()
    background = background.merge(per_codon[["pos", "n_missense"]], left_on="human_pos",
                                  right_on="pos", how="left")
    background["n_missense"] = background["n_missense"].fillna(0).astype(int)
    background["n_all"] = background["human_pos"].map(
        coding.groupby("pos").size()).fillna(0).astype(int)

    rho_missense = stats.spearmanr(background["conservation"], background["n_missense"])
    rho_all = stats.spearmanr(background["conservation"], background["n_all"])

    # mutation density (per residue) in conserved vs variable halves
    median_cons = background["conservation"].median()
    conserved_half = background[background["conservation"] >= median_cons]
    variable_half = background[background["conservation"] < median_cons]

    # hotspots vs the rest (Mann-Whitney; n = 5, reported with the caveat)
    hotspot_scores = top5["conservation"].tolist()
    rest_scores = background[~background["human_pos"].isin(top5["pos"])]["conservation"].tolist()
    mwu = stats.mannwhitneyu(hotspot_scores, rest_scores, alternative="greater")

    # enrichment: are the five hotspots among the most conserved 20% of codons?
    cutoff = background["conservation"].quantile(0.80)
    k = int((top5["conservation"] >= cutoff).sum())
    p_enrich = hypergeom_sf(k, len(background), int((background["conservation"] >= cutoff).sum()), 5)

    summary = dict(
        n_mutations=n_total,
        n_samples=int(n_samples),
        n_studies=int(n_studies),
        n_coding_parseable=n_parsed,
        n_missense=int(len(missense)),
        pct_missense=round(100 * len(missense) / n_total, 1),
        n_nonsense=int((coding["mut"] == "*").sum()),
        pct_nonsense=round(100 * (coding["mut"] == "*").sum() / n_total, 1),
        top5_codons=[f"{r.wt}{r.pos}" for r in top5.itertuples()],
        top5_counts=[int(v) for v in top5["n_missense"]],
        top5_pct_of_missense=round(float(top5["pct_of_missense"].sum()), 1),
        top5_pct_of_all=round(float(top5["pct_of_all_mutations"].sum()), 1),
        top5_mean_conservation=round(float(top5["conservation"].mean()), 3),
        all_mean_conservation=round(float(background["conservation"].mean()), 3),
        conserved20pct_cutoff=round(float(cutoff), 3),
        hotspots_in_top20pct=k,
        hypergeom_p=p_enrich,
        mwu_p=float(mwu.pvalue),
        spearman_conservation_vs_missense=dict(rho=float(rho_missense.correlation),
                                               p=float(rho_missense.pvalue)),
        spearman_conservation_vs_all=dict(rho=float(rho_all.correlation), p=float(rho_all.pvalue)),
        mutations_per_residue_conserved_half=round(float(conserved_half["n_all"].sum() / len(conserved_half)), 2),
        mutations_per_residue_variable_half=round(float(variable_half["n_all"].sum() / len(variable_half)), 2),
        mutations_per_residue_dbd=round(float(background[background.domain == "DBD"]["n_all"].sum()
                                               / (background.domain == "DBD").sum()), 2),
        mutations_per_residue_outside_dbd=round(float(background[background.domain != "DBD"]["n_all"].sum()
                                                      / (background.domain != "DBD").sum()), 2),
        n_invariant_residues=int((background["conservation"] >= 0.999).sum()),
        dbd_share_of_all_mutations=round(100 * float(background[background.domain == "DBD"]["n_all"].sum())
                                         / float(background["n_all"].sum()), 1),
        dbd_share_of_residues=round(100 * float((background.domain == "DBD").sum()) / len(background), 1),
    )
    (ROOT / "results" / "hotspot_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("\n", top5.to_string(index=False))


if __name__ == "__main__":
    main()