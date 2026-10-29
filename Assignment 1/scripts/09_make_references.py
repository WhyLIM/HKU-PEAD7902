#!/usr/bin/env python3
"""Build EndNote-importable reference files (RIS and EndNote XML) from the
references used in the report, after verifying them against Crossref.

Each record below was checked against api.crossref.org by DOI on 2026-10-26;
title, journal, year, volume, issue, pages and author order follow the
publisher metadata. Two errors in the first draft were corrected:
  * Hoadley et al. 2018 had a wrong DOI (10.1016/j.cell.2018.03.046 belongs to
    an unrelated Cell paper); the correct DOI is 10.1016/j.cell.2018.03.022.
  * Joerger & Fersht 2008 was cited with a non-existent title; the published
    title is "Structural biology of the tumor suppressor p53".

Long consortium author lists are truncated to the first ten authors, which is
the usual convention for citation lists; the RIS/XML records keep those authors
and EndNote styles add "et al." automatically.

Outputs
  references/references.ris
  references/references_endnote.xml
  references/references.md
"""
import pathlib
import xml.sax.saxutils as sx

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "references"
OUT.mkdir(exist_ok=True)

REFS = [
    dict(
        key="uniprot2025",
        authors=["UniProt Consortium"],
        title="UniProt: the Universal Protein Knowledgebase in 2025",
        journal="Nucleic Acids Research",
        abbr="Nucleic Acids Res",
        year=2025, volume="53", issue="D1", pages="D609-D617",
        doi="10.1093/nar/gkae1010",
        note="Corporate author; Crossref lists the individual curators.",
    ),
    dict(
        key="sievers2011",
        authors=["Sievers, Fabian", "Wilm, Andreas", "Dineen, David", "Gibson, Toby J",
                 "Karplus, Kevin", "Li, Weizhong", "Lopez, Rodrigo", "McWilliam, Hamish",
                 "Remmert, Michael", "Soding, Johannes", "Thompson, Julie D",
                 "Higgins, Desmond G"],
        title="Fast, scalable generation of high-quality protein multiple sequence "
              "alignments using Clustal Omega",
        journal="Molecular Systems Biology", abbr="Mol Syst Biol",
        year=2011, volume="7", issue="1", pages="539",
        doi="10.1038/msb.2011.75",
        note="Article number 539 (journal has no page range for this item).",
    ),
    dict(
        key="madeira2022",
        authors=["Madeira, Fabio", "Pearce, Matt", "Tivey, Adrian R N", "Basutkar, Prasad",
                 "Lee, Joon", "Edbali, Ossama", "Madhusoodanan, Nandana",
                 "Kolesnikov, Anton", "Lopez, Rodrigo"],
        title="Search and sequence analysis tools services from EMBL-EBI in 2022",
        journal="Nucleic Acids Research", abbr="Nucleic Acids Res",
        year=2022, volume="50", issue="W1", pages="W276-W279",
        doi="10.1093/nar/gkac240",
    ),
    dict(
        key="cerami2012",
        authors=["Cerami, Ethan", "Gao, Jianjiong", "Dogrusoz, Ugur", "Gross, Benjamin E",
                 "Sumer, Selcuk Onur", "Aksoy, Bulent Arman", "Jacobsen, Anders",
                 "Byrne, Caitlin J", "Heuer, Michael L", "Larsson, Erik", "Antipin, Yevgeniy",
                 "Reva, Boris", "Goldberg, Arthur P", "Sander, Chris", "Schultz, Nikolaus"],
        title="The cBio cancer genomics portal: an open platform for exploring "
              "multidimensional cancer genomics data",
        journal="Cancer Discovery", abbr="Cancer Discov",
        year=2012, volume="2", issue="5", pages="401-404",
        doi="10.1158/2159-8290.CD-12-0095",
    ),
    dict(
        key="debruijn2023",
        authors=["de Bruijn, Ino", "Kundra, Ritika", "Mastrogiacomo, Brooke",
                 "Tran, Thinh Ngoc", "Sikina, Luke", "Mazor, Tali", "Li, Xiang",
                 "Ochoa, Angelica", "Zhao, Gaofei", "Lai, Bryan"],
        title="Analysis and visualization of longitudinal genomic and clinical data from "
              "the AACR Project GENIE Biopharma Collaborative in cBioPortal",
        journal="Cancer Research", abbr="Cancer Res",
        year=2023, volume="83", issue="23", pages="3861-3867",
        doi="10.1158/0008-5472.CAN-23-0816",
        note="73 authors in total; first ten listed.",
    ),
    dict(
        key="hoadley2018",
        authors=["Hoadley, Katherine A", "Yau, Christina", "Hinoue, Toshinori",
                 "Wolf, Denise M", "Lazar, Alexander J", "Drill, Esther", "Shen, Ronglai",
                 "Taylor, Alison M", "Cherniack, Andrew D", "Thorsson, Vesteinn"],
        title="Cell-of-origin patterns dominate the molecular classification of 10,000 "
              "tumors from 33 types of cancer",
        journal="Cell", abbr="Cell",
        year=2018, volume="173", issue="2", pages="291-304.e6",
        doi="10.1016/j.cell.2018.03.022",
        note="Large consortium paper; first ten authors listed.",
    ),
    dict(
        key="bouaoun2016",
        authors=["Bouaoun, Liacine", "Sonkin, Dmitriy", "Ardin, Maude", "Hollstein, Monica",
                 "Byrnes, Graham", "Zavadil, Jiri", "Olivier, Magali"],
        title="TP53 variations in human cancers: new lessons from the IARC TP53 database "
              "and genomics data",
        journal="Human Mutation", abbr="Hum Mutat",
        year=2016, volume="37", issue="9", pages="865-876",
        doi="10.1002/humu.23035",
    ),
    dict(
        key="cho1994",
        authors=["Cho, Yunje", "Gorina, Svetlana", "Jeffrey, Philip D", "Pavletich, Nikola P"],
        title="Crystal structure of a p53 tumor suppressor-DNA complex: understanding "
              "tumorigenic mutations",
        journal="Science", abbr="Science",
        year=1994, volume="265", issue="5170", pages="346-355",
        doi="10.1126/science.8023157",
    ),
    dict(
        key="joerger2008",
        authors=["Joerger, Andreas C", "Fersht, Alan R"],
        title="Structural biology of the tumor suppressor p53",
        journal="Annual Review of Biochemistry", abbr="Annu Rev Biochem",
        year=2008, volume="77", issue="1", pages="557-582",
        doi="10.1146/annurev.biochem.77.060806.091238",
    ),
    dict(
        key="kastenhuber2017",
        authors=["Kastenhuber, Edward R", "Lowe, Scott W"],
        title="Putting p53 in context",
        journal="Cell", abbr="Cell",
        year=2017, volume="170", issue="6", pages="1062-1078",
        doi="10.1016/j.cell.2017.08.028",
    ),
    dict(
        key="sondka2024",
        authors=["Sondka, Zbyslaw", "Dhir, Nidhi Bindal", "Carvalho-Silva, Denise", "Jupe, Steven"],
        title="COSMIC: a curated database of somatic variants and clinical data for cancer",
        journal="Nucleic Acids Research", abbr="Nucleic Acids Res",
        year=2024, volume="52", issue="D1", pages="D1210-D1217",
        doi="10.1093/nar/gkad986",
    ),
    dict(
        key="shrake1973",
        authors=["Shrake, A", "Rupley, J A"],
        title="Environment and exposure to solvent of protein atoms. Lysozyme and insulin",
        journal="Journal of Molecular Biology", abbr="J Mol Biol",
        year=1973, volume="79", issue="2", pages="351-371",
        doi="10.1016/0022-2836(73)90011-9",
    ),
    dict(
        key="tien2013",
        authors=["Tien, Matthew Z", "Meyer, Austin G", "Sydykova, Dariya K",
                 "Spielman, Stephanie J", "Wilke, Claus O"],
        title="Maximum allowed solvent accessibilities of residues in proteins",
        journal="PLoS ONE", abbr="PLoS One",
        year=2013, volume="8", issue="11", pages="e80635",
        doi="10.1371/journal.pone.0080635",
    ),
    dict(
        key="shannon1948",
        authors=["Shannon, Claude E"],
        title="A mathematical theory of communication",
        journal="Bell System Technical Journal", abbr="Bell Syst Tech J",
        year=1948, volume="27", issue="3", pages="379-423",
        doi="10.1002/j.1538-7305.1948.tb01338.x",
    ),
    dict(
        key="lin1991",
        authors=["Lin, J"],
        title="Divergence measures based on the Shannon entropy",
        journal="IEEE Transactions on Information Theory", abbr="IEEE Trans Inf Theory",
        year=1991, volume="37", issue="1", pages="145-151",
        doi="10.1109/18.61115",
    ),
]


def split_pages(p):
    return p.split("-", 1) if "-" in p else (p, "")


def ris():
    lines = []
    for r in REFS:
        sp, ep = split_pages(r["pages"])
        lines.append("TY  - JOUR")
        lines.append(f"TI  - {r['title']}")
        for a in r["authors"]:
            lines.append(f"AU  - {a}")
        lines.append(f"T2  - {r['journal']}")
        lines.append(f"JA  - {r['abbr']}")
        lines.append(f"VL  - {r['volume']}")
        lines.append(f"IS  - {r['issue']}")
        lines.append(f"SP  - {sp}")
        if ep:
            lines.append(f"EP  - {ep}")
        lines.append(f"PY  - {r['year']}")
        lines.append(f"DO  - {r['doi']}")
        lines.append(f"UR  - https://doi.org/{r['doi']}")
        if r.get("note"):
            lines.append(f"N1  - {r['note']}")
        lines.append("ER  -")
        lines.append("")
    return "\n".join(lines)


def endnote_xml():
    recs = []
    for i, r in enumerate(REFS, 1):
        authors = "".join(f"        <author>{sx.escape(a)}</author>\n" for a in r["authors"])
        recs.append(f"""  <record>
    <rec-number>{i}</rec-number>
    <ref-type name="Journal Article">17</ref-type>
    <contributors>
      <authors>
{authors}      </authors>
    </contributors>
    <titles>
      <title>{sx.escape(r['title'])}</title>
      <secondary-title>{sx.escape(r['journal'])}</secondary-title>
    </titles>
    <periodical>
      <full-title>{sx.escape(r['journal'])}</full-title>
      <abbr-1>{sx.escape(r['abbr'])}</abbr-1>
    </periodical>
    <pages>{sx.escape(r['pages'])}</pages>
    <volume>{r['volume']}</volume>
    <number>{r['issue']}</number>
    <dates>
      <year>{r['year']}</year>
    </dates>
    <electronic-resource-num>{r['doi']}</electronic-resource-num>
    <urls>
      <related-urls>
        <url>https://doi.org/{r['doi']}</url>
      </related-urls>
    </urls>
    <notes>{sx.escape(r.get('note', ''))}</notes>
  </record>""")
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<xml>\n<records>\n' + "\n".join(recs) + '\n</records>\n</xml>\n')


def initials(name):
    """'Sievers, Fabian' -> 'Sievers F'; 'Tivey, Adrian R N' -> 'Tivey ARN'."""
    if "," not in name:
        return name
    family, given = [p.strip() for p in name.split(",", 1)]
    letters = "".join(w[0] for w in given.replace("-", " ").split() if w)
    return f"{family} {letters}"


def markdown():
    out = ["# References (verified against Crossref by DOI, 2026-10-28)", ""]
    for i, r in enumerate(REFS, 1):
        authors = [initials(a) for a in r["authors"]]
        if len(authors) > 6:
            au = ", ".join(authors[:6]) + ", et al."
        else:
            au = ", ".join(authors) + "."
        pages = r["pages"].replace("-", "\u2013")
        out.append(f"{i}. {au} {r['title']}. *{r['abbr']}*. {r['year']};"
                   f"{r['volume']}({r['issue']}):{pages}. doi:{r['doi']}")
        out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    (OUT / "references.ris").write_text(ris(), encoding="utf-8")
    (OUT / "references_endnote.xml").write_text(endnote_xml(), encoding="utf-8")
    (OUT / "references.md").write_text(markdown(), encoding="utf-8")
    print(f"wrote {len(REFS)} records to {OUT}")