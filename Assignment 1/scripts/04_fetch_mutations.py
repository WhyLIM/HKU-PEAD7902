#!/usr/bin/env python3
"""Fetch all somatic TP53 mutations from TCGA Pan-Cancer Atlas studies via the cBioPortal REST API.

Outputs
-------
data/mutations/raw/<studyId>.json : raw API responses
data/mutations/tp53_mutations_pancan.tsv : merged, de-duplicated mutation table
"""
import concurrent.futures as cf
import json
import pathlib
import sys
import time

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "mutations" / "raw"
RAW.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data" / "mutations" / "tp53_mutations_pancan.tsv"
API = "https://www.cbioportal.org/api"
TP53_ENTREZ = 7157


def get(url, params=None, tries=5):
    for attempt in range(tries):
        try:
            r = requests.get(url, params=params, timeout=120)
            if r.status_code == 200:
                return r.json()
            print(f"  HTTP {r.status_code}: {url}", file=sys.stderr)
        except Exception as exc:  # network flakiness is common on this endpoint
            print(f"  retry {attempt+1}/{tries}: {exc}", file=sys.stderr)
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}")


def main():
    studies = get(f"{API}/studies", {"projection": "SUMMARY", "pageSize": 2000})
    study_ids = sorted(
        s["studyId"] for s in studies if s["studyId"].endswith("_tcga_pan_can_atlas_2018")
    )
    print(f"{len(study_ids)} TCGA Pan-Cancer Atlas studies")

    def fetch_study(study):
        cache = RAW / f"{study}.json"
        if cache.exists() and cache.stat().st_size > 2:
            try:
                return study, json.loads(cache.read_text())
            except (OSError, json.JSONDecodeError):
                pass
        try:
            data = get(
                f"{API}/molecular-profiles/{study}_mutations/mutations",
                {
                    "sampleListId": f"{study}_sequenced",
                    "entrezGeneId": TP53_ENTREZ,
                    "projection": "SUMMARY",
                },
                tries=8,
            )
        except Exception as exc:
            print(f"FAILED {study}: {exc}", file=sys.stderr)
            return study, None
        cache.write_text(json.dumps(data))
        return study, data

    rows = []
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        for study, data in pool.map(fetch_study, study_ids):
            if data is None:
                print(f"{study:45s} not retrieved", flush=True)
                continue
            print(f"{study:45s} {len(data):5d} TP53 mutations", flush=True)
            rows.extend(data)

    seen, merged = set(), []
    for m in rows:
        key = (m["sampleId"], m.get("proteinChange"), m.get("startPosition"), m.get("variantAllele"))
        if key in seen:
            continue
        seen.add(key)
        merged.append(m)

    cols = [
        "sampleId", "patientId", "studyId", "proteinChange", "mutationType",
        "variantType", "chr", "startPosition", "referenceAllele", "variantAllele",
        "tumorAltCount", "tumorRefCount", "hgvspShort",
    ]
    with OUT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m in merged:
            fh.write("\t".join(str(m.get(c, "")) for c in cols) + "\n")
    print(f"\nWrote {len(merged)} unique TP53 mutations to {OUT}")


if __name__ == "__main__":
    main()