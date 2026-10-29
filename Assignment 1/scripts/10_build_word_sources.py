#!/usr/bin/env python3
"""Build the md-to-word sources directly from the two Markdown reports so that
the Word files contain exactly the report content (no extra abstract, title or
table notes)."""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]

KEYS = {
    1: "uniprot2025", 2: "sievers2011", 3: "madeira2022", 4: "cho1994",
    5: "joerger2008", 6: "kastenhuber2017", 7: "hoadley2018", 8: "cerami2012",
    9: "debruijn2023", 10: "sondka2024", 11: "bouaoun2016", 12: "shannon1948",
    13: "lin1991", 14: "shrake1973", 15: "tien2013",
}

CITE = re.compile(r"\[(\d+(?:[–\-,]\d+)?(?:,\d+)*)\]")


def keys_for(expr: str) -> list[str]:
    nums = []
    for part in expr.split(","):
        if "–" in part or "-" in part:
            a, b = re.split(r"[–-]", part)
            nums.extend(range(int(a), int(b) + 1))
        else:
            nums.append(int(part))
    return [KEYS[n] for n in nums]


def convert_citations(text: str) -> str:
    def repl(m):
        keys = keys_for(m.group(1))
        return "[@" + "; @".join(keys) + "]"
    return CITE.sub(repl, text)


def build(src: pathlib.Path, out: pathlib.Path, lang: str, width: str):
    lines = src.read_text().splitlines()
    title = next(l[2:].strip() for l in lines if l.startswith("# "))
    if lang == "en":
        front = ('%% Title ' + title + "\n\n"
                 "%% Subtitle PEAD7902 — Assignment 1\n\n"
                 "%% Author Min Li (3036744858)\n\n"
                 "%% Date October 2026")
        csl = "nature-spaced.csl"
    else:
        front = ('%% 文档标题 ' + title + "\n\n"
                 "%% 副标题 PEAD7902 — 作业一\n\n"
                 "%% 作者 Min Li（3036744858）\n\n"
                 "%% 日期 2026 年 10 月")
        csl = "gb7714-spaced.csl"

    start = next(i for i, l in enumerate(lines)
                 if l.startswith("## Summary") or l.startswith("## 摘要"))
    end = next(i for i, l in enumerate(lines) if l.startswith("## 参考文献") or l.startswith("## References"))

    body = []
    i = start
    while i < end:
        line = lines[i]
        s = line.strip()
        if not s:
            body.append("")
            i += 1
            continue
        if s.startswith("### "):
            body.append("## " + s[4:])
            i += 1
            continue
        if s.startswith("## "):
            body.append("# " + s[3:])
            i += 1
            continue
        m = re.match(r"!\[[^]]*\]\(([^)]+)\)", s)
        if m:
            path = m.group(1).replace("assets/", "../../assets/")
            body.append(f"![]({path}){{width={width}}}")
            i += 1
            # caption immediately after the image (skip blank lines)
            j = i
            while j < end and not lines[j].strip():
                j += 1
            if j < end and lines[j].strip().startswith(("**Figure", "**图")):
                body.append("")
                body.append("%% Image Caption " + lines[j].strip())
                i = j + 1
            continue
        if s.startswith(("**Table", "**表")):
            j = i + 1
            while j < end and not lines[j].strip():
                j += 1
            if j < end and lines[j].strip().startswith("|"):
                body.append("%% Table Caption " + s)
                i += 1
                continue
        body.append(convert_citations(s).replace("`", ""))
        i += 1

    text = ("---\n"
            "bibliography: references.bib\n"
            f"csl: {csl}\n"
            "link-citations: true\n"
            "---\n\n"
            + front + "\n\n"
            + "\n".join(body).rstrip()
            + "\n\n%% References Title " + ("References" if lang == "en" else "参考文献") + "\n")
    out.write_text(text)
    print("wrote", out, len(text), "chars")


if __name__ == "__main__":
    build(ROOT / "TP53_Assignment_Report.md", ROOT / "word/en/manuscript.md", "en", "15cm")
    build(ROOT / "TP53_Assignment_Report_ZH.md", ROOT / "word/zh/manuscript.md", "zh", "14cm")