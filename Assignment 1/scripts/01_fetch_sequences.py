#!/usr/bin/env python3
"""Step 1 - download TP53 ortholog protein sequences from UniProtKB (reviewed entries).

Run:  python scripts/01_fetch_sequences.py
"""
import pathlib
import time

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
SEQ = ROOT / "data" / "sequences"
SEQ.mkdir(parents=True, exist_ok=True)

ACCESSIONS = {
    "P04637": "Homo sapiens (human)",
    "P02340": "Mus musculus (mouse)",
    "P67939": "Bos taurus (cattle)",
    "P10360": "Gallus gallus (chicken)",
    "P07193": "Xenopus laevis (African clawed frog)",
    "P79734": "Danio rerio (zebrafish)",
}

for acc, name in ACCESSIONS.items():
    out = SEQ / f"{acc}.fasta"
    if out.exists() and out.stat().st_size > 50:
        print(f"cached  {acc}  {name}")
        continue
    for attempt in range(5):
        try:
            r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60)
            if r.ok and r.text.startswith(">"):
                out.write_text(r.text)
                print(f"fetched {acc}  {name}")
                break
        except Exception as exc:
            print(f"  retry {attempt + 1}: {exc}")
        time.sleep(3)
    else:
        raise RuntimeError(f"could not download {acc}")
