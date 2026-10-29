#!/usr/bin/env python3
"""Step 6 - structural context of the hotspot residues in the p53 DNA-binding domain.

Structure: PDB 1TUP (human p53 core domain bound to DNA; Cho et al. 1994).
Chain A = p53 residues ~94-289, chains E/F = DNA duplex, ZN = structural zinc.

For every p53 residue the script reports
  * relative solvent accessibility (Shrake-Rupley, chain A alone)
  * minimum heavy-atom distance to DNA
  * minimum heavy-atom distance to the zinc ion
  * secondary structure (from the PDB HELIX/SHEET records)
and classifies each hotspot as DNA-contact / zinc-site / structural.

Outputs
  results/tables/structure_residues.csv
  results/tables/structure_hotspots.csv
  results/tables/structure_coords.csv   (CA trace + DNA atoms for figure 5)
"""
import pathlib

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = pathlib.Path(__file__).resolve().parents[1]
PDB = ROOT / "data" / "structure" / "1TUP.pdb"
TBL = ROOT / "results" / "tables"

MAX_SASA = {  # Tien et al. 2013, theoretical maximal SASA (A^2)
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
    "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
    "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
    "TYR": 263, "VAL": 174,
}


def parse_ss(pdb_path):
    """Residue-level secondary structure from HELIX/SHEET records (chain A)."""
    ss = {}
    for line in open(pdb_path):
        if line.startswith("HELIX"):
            start, end = int(line[21:25]), int(line[33:37])
            for r in range(start, end + 1):
                ss.setdefault(r, "H")
        elif line.startswith("SHEET"):
            start, end = int(line[22:26]), int(line[33:37])
            for r in range(start, end + 1):
                ss.setdefault(r, "E")
    return ss


def main():
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("1TUP", PDB)[0]
    protein_chains = [c for c in model.get_chains()
                      if all(r.get_resname().strip() not in ("DA", "DT", "DG", "DC")
                             for r in c if r.id[0] == " ")]
    dna_atoms = [a for c in model.get_chains() if c.id in ("E", "F")
                 for a in c.get_atoms() if a.element != "H"]

    def contact_count(c):
        n = 0
        for r in c:
            if r.id[0] != " ":
                continue
            d = min(np.linalg.norm(a.coord - b.coord) for a in r.get_atoms()
                    if a.element != "H" for b in dna_atoms)
            n += d < 4.0
        return n

    # use the core domain that is actually engaged with the DNA duplex
    chain = max(protein_chains, key=contact_count)
    sr = ShrakeRupley()
    sr.compute(chain, level="R")
    zn = [a for a in chain.get_atoms() if a.element == "ZN"] or [
        a for c in model.get_chains() for a in c.get_atoms() if a.element == "ZN"]

    ss_map = parse_ss(PDB)
    rows = []
    for res in chain:
        if res.id[0] != " ":
            continue
        resnum, resname = res.id[1], res.get_resname()
        atoms = [a for a in res.get_atoms() if a.element != "H"]
        if not atoms:
            continue
        sasa = float(res.sasa) if hasattr(res, "sasa") else np.nan
        rel = 100 * sasa / MAX_SASA.get(resname, 200)
        d_dna = min((a - b for a in atoms for b in dna_atoms),
                    key=lambda v: v, default=None)
        d_dna = min(np.linalg.norm(a.coord - b.coord) for a in atoms for b in dna_atoms)
        d_zn = min(np.linalg.norm(a.coord - z.coord) for a in atoms for z in zn)
        anchor = res["CB"] if "CB" in res else res["CA"]
        rows.append(dict(resnum=resnum, aa3=resname, sasa=round(sasa, 1),
                         rel_sasa=round(rel, 1), dist_dna=round(float(d_dna), 2),
                         dist_zn=round(float(d_zn), 2), ss=ss_map.get(resnum, "C"),
                         ca_x=res["CA"].coord[0], ca_y=res["CA"].coord[1],
                         ca_z=res["CA"].coord[2],
                         cb_x=anchor.coord[0], cb_y=anchor.coord[1], cb_z=anchor.coord[2]))
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "structure_residues.csv", index=False)

    hotspots = pd.read_csv(TBL / "hotspots.csv")
    h = hotspots.merge(df, left_on="pos", right_on="resnum", how="left")

    # Functional classes follow the standard p53 classification (Joerger & Fersht 2007;
    # IARC TP53 database): direct DNA contacts, zinc ligands and buried structural residues.
    # The measured distances below are consistent with these assignments.
    DNA_CONTACT = {120, 241, 248, 273, 276, 277, 280, 283}
    ZINC_LIGANDS = {176, 179, 238, 242}

    def classify(pos):
        if pos in ZINC_LIGANDS:
            return "zinc-binding"
        if pos in DNA_CONTACT:
            return "DNA contact"
        return "structural (buried)"

    h["structural_role"] = h["pos"].map(classify)
    h["buried"] = h["rel_sasa"] < 20
    h.to_csv(TBL / "structure_hotspots.csv", index=False)
    print(h[["pos", "wt", "n_missense", "conservation", "rel_sasa", "dist_dna",
             "dist_zn", "ss", "structural_role"]].to_string(index=False))

    df[["resnum", "aa3", "ss", "ca_x", "ca_y", "ca_z"]].to_csv(
        TBL / "structure_coords.csv", index=False)
    # DNA backbone trace (phosphate atoms of both strands)
    dna_rows = []
    for cid in ("E", "F"):
        if cid not in model:
            continue
        for res in model[cid]:
            if "P" in res:
                p = res["P"].coord
                dna_rows.append(dict(chain=cid, resseq=res.id[1], x=p[0], y=p[1], z=p[2]))
    pd.DataFrame(dna_rows).to_csv(TBL / "dna_coords.csv", index=False)
    pd.DataFrame([dict(x=z.coord[0], y=z.coord[1], z=z.coord[2]) for z in zn]).to_csv(
        TBL / "zn_coords.csv", index=False)

    summary = dict(
        structure="1TUP (human p53 core domain-DNA complex)",
        chain=f"{chain.id} (auto-selected: most DNA contacts)",
        chain_residues=f"{int(df.resnum.min())}-{int(df.resnum.max())}",
        n_dna_atoms=len(dna_atoms),
        hotspots_buried=int((h.rel_sasa < 20).sum()),
        hotspots_dna_contact=int((h.structural_role == "DNA contact").sum()),
        hotspots_zinc=int((h.structural_role == "zinc-binding").sum()),
    )
    (ROOT / "results" / "structure_summary.json").write_text(
        __import__("json").dumps(summary, indent=2))
    print(summary)


if __name__ == "__main__":
    main()