#!/usr/bin/env python3
"""
Convert ClinVar TSV/VCF into a compact CSV usable by the variant-effect
prediction workflow.

Output columns:
    chrom,pos,ref,alt,gene_symbol,clinical_significance,hgvs_p

- Keeps only Pathogenic and Benign (collapses Likely_* into those).
- Works with variant_summary.txt.gz (tab-delimited) or ClinVar VCF.
"""

import argparse, re
from pathlib import Path
import pandas as pd
import numpy as np

# Map ClinVar terms to simplified labels
KEEP_MAP = {
    "pathogenic": "Pathogenic",
    "likely_pathogenic": "Pathogenic",
    "benign": "Benign",
    "likely_benign": "Benign"
}

def normalize_clinsig(s, collapse_likely=True):
    """Normalize ClinVar clinical significance string."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return None
    parts = re.split(r"[|/,;]", str(s))
    for p in parts:
        key = p.strip().lower().replace(" ", "_")
        if key in KEEP_MAP:
            if collapse_likely:
                return KEEP_MAP[key]
            if key in ["pathogenic","benign"]:
                return key.capitalize()
    return None

def tsv_to_csv(path: Path, out_csv: Path) -> Path:
    """Convert ClinVar TSV (e.g., variant_summary.txt.gz) to clean CSV."""
    compression = "infer"  # pandas auto-detects .gz
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, compression=compression)
    cols = {c.lower(): c for c in df.columns}

    def pick(*cands):
        for c in cands:
            if c.lower() in cols:
                return cols[c.lower()]
        return None

    col_chrom = pick("chromosome","chrom","chr")
    col_pos   = pick("start","pos","position")
    col_ref   = pick("referenceallele","ref")
    col_alt   = pick("alternateallele","alt","altallele")
    col_gene  = pick("genesymbol","gene_symbol","gene","symbol")
    col_sig   = pick("clinicalsignificance","clinsig","clinical_significance")
    col_hgvs_p= pick("hgvsp","hgvs_p","protein_change","protein_hgvs")

    if not all([col_chrom,col_pos,col_ref,col_alt,col_gene,col_sig]):
        raise ValueError("Missing required columns in TSV.")

    df["clinical_significance"] = df[col_sig].apply(normalize_clinsig)
    df = df[df["clinical_significance"].isin(["Pathogenic","Benign"])].copy()

    out = pd.DataFrame({
        "chrom": df[col_chrom].astype(str),
        "pos": pd.to_numeric(df[col_pos], errors="coerce").astype("Int64"),
        "ref": df[col_ref].astype(str),
        "alt": df[col_alt].astype(str),
        "gene_symbol": df[col_gene].astype(str),
        "clinical_significance": df["clinical_significance"],
        "hgvs_p": df[col_hgvs_p].astype(str) if col_hgvs_p else ""
    }).dropna(subset=["pos"])

    out.to_csv(out_csv, index=False)
    return out_csv

def vcf_to_csv(path: Path, out_csv: Path) -> Path:
    """Convert ClinVar VCF to clean CSV (requires cyvcf2)."""
    from cyvcf2 import VCF
    rows = []
    for rec in VCF(str(path)):
        info = rec.INFO or {}
        label = normalize_clinsig(info.get("CLNSIG") or info.get("CLNSIGCONF") or "")
        if label not in ["Pathogenic","Benign"]:
            continue
        gene = ""
        gi = info.get("GENEINFO")
        if gi:
            gene = str(gi).split("|")[0].split(":")[0]
        hgvs_p = ""
        for key in ["HGVSP","HGVSp","HGVS_P","HGVS"]:
            if key in info and info[key]:
                hgvs_p = str(info[key]); break
        chrom = str(rec.CHROM); pos = int(rec.POS); ref = rec.REF
        for alt in rec.ALT or []:
            rows.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "gene_symbol": gene, "clinical_significance": label, "hgvs_p": hgvs_p
            })
    out = pd.DataFrame.from_records(rows)
    if not out.empty:
        out.to_csv(out_csv, index=False)
    return out_csv

def main(inp: Path, out: Path):
    sfx = "".join(inp.suffixes).lower()
    if any(x in sfx for x in [".tsv",".txt"]):
        tsv_to_csv(inp, out)
        print("Wrote CSV:", out)
    elif any(x in sfx for x in [".vcf",".vcf.gz",".gz"]):
        vcf_to_csv(inp, out)
        print("Wrote CSV:", out)
    else:
        raise ValueError("Unknown input type. Use TSV/TXT or VCF/VCF.GZ")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Convert ClinVar TSV/VCF into CSV for ML pipeline")
    ap.add_argument("--input", type=Path, required=True, help="Path to ClinVar TSV or VCF")
    ap.add_argument("--output", type=Path, required=True, help="Output CSV path")
    args = ap.parse_args()
    main(args.input, args.output)

