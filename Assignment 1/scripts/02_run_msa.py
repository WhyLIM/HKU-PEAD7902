#!/usr/bin/env python3
"""Step 2 - run a multiple sequence alignment with Clustal Omega (EBI REST service).

The six sequences are submitted as one FASTA file; the service returns a
Clustal-format alignment (with column numbering) and the guide tree.

Run:  python scripts/02_run_msa.py
"""
import pathlib
import time

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
SEQ = ROOT / "data" / "sequences"
OUT = ROOT / "results" / "alignment"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"


def main():
    fasta = "".join((SEQ / f"{a}.fasta").read_text() for a in
                    ["P04637", "P02340", "P67939", "P10360", "P07193", "P79734"])
    (OUT / "all_tp53.fasta").write_text(fasta)

    r = requests.post(f"{BASE}/run", data={
        "email": "student@connect.hku.hk",
        "stype": "protein",
        "sequence": fasta,
    }, timeout=120)
    r.raise_for_status()
    job = r.text.strip()
    print("job", job)

    for _ in range(60):
        status = requests.get(f"{BASE}/status/{job}", timeout=60).text.strip()
        print(" ", status)
        if status == "FINISHED":
            break
        if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(status)
        time.sleep(5)
    else:
        raise TimeoutError("Clustal Omega job did not finish")

    aln = requests.get(f"{BASE}/result/{job}/aln-clustal_num", timeout=60).text
    (OUT / "tp53_clustalo.aln").write_text(aln)
    tree = requests.get(f"{BASE}/result/{job}/phylotree", timeout=60).text
    (OUT / "tp53_clustalo.nwk").write_text(tree)
    (OUT / "job_id.txt").write_text(job + "\n")
    print("saved alignment and guide tree")


if __name__ == "__main__":
    main()
