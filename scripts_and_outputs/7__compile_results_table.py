#!/usr/bin/env python3
"""
Compile per-site selection results from CODEML, FUBAR, and MEME into a summary TSV.

Usage:
  python3 compile_results_table.py \
      --alignment results_ATB_new/ATB_diphtheriae_subset.fa \
      --outdir    results_ATB_new/diphtheriae \
      --group     diphtheriae \
      --output    selection_summary_diphtheriae.tsv

Expects CODEML output at <outdir>/M8/codeml.out and <outdir>/M8/rst,
FUBAR JSON at <alignment>.FUBAR.json, MEME JSON at <alignment>.MEME.json.
"""

import argparse
import json
import re
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

SCRIPT = Path(__file__).name
BEB_THRESHOLD   = 0.80
FUBAR_THRESHOLD = 0.80

COLS = [
    "group", "site", "ref_seq", "ref_residue", "alt_residues",
    "codeml_beb_pp", "codeml_posterior_mean_ω", "codeml_posterior_mean_ω_sd",
    "fubar_pp_pos", "fubar_bayes_factor", "fubar_dS", "fubar_dN", "fubar_omega",
    "meme_pval", "meme_branches",
]


def parse_beb_pp(codeml_out):
    sites = {}
    text = Path(codeml_out).read_text()
    start = text.find("Bayes Empirical Bayes")
    if start == -1:
        return sites
    for line in text[start:].split("\n"):
        m = re.match(r"\s*(\d+)\s+\S\s+([\d.]+)(\**)", line)
        if m:
            sites[int(m.group(1))] = (float(m.group(2)), m.group(3))
    return sites


def parse_beb_rst(rst_path):
    sites = {}
    text = Path(rst_path).read_text()
    start = text.find("Bayes Empirical Bayes (BEB)")
    if start == -1:
        return sites
    for line in text[start:].split("\n"):
        m = re.match(r"\s*(\d+)\s+\S\s+[\d.\s()+\-]+\s+([\d.]+)\s+\+-\s+([\d.]+)", line)
        if m:
            sites[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return sites


def parse_fubar_json(json_path):
    d = json.loads(Path(json_path).read_text())
    headers = [h[0] for h in d["MLE"]["headers"]]
    pp_i    = headers.index("Prob[alpha<beta]")
    bf_i    = headers.index("BayesFactor[alpha<beta]")
    alpha_i = headers.index("alpha")
    beta_i  = headers.index("beta")
    sites = {}
    for i, row in enumerate(d["MLE"]["content"]["0"], 1):
        alpha = row[alpha_i]
        beta  = row[beta_i]
        sites[i] = {
            "pp":    row[pp_i],
            "bf":    row[bf_i],
            "dS":    alpha,
            "dN":    beta,
            "omega": beta / alpha if alpha > 0 else float("inf"),
        }
    return sites


def parse_meme_json(json_path):
    d = json.loads(Path(json_path).read_text())
    headers = [h[0] for h in d["MLE"]["headers"]]
    pval_i = headers.index("p-value")
    nb_i   = headers.index("# branches under selection")
    sites = {}
    for i, row in enumerate(d["MLE"]["content"]["0"], 1):
        sites[i] = {"pval": row[pval_i], "branches": int(row[nb_i])}
    return sites


def codon_to_aa(codon):
    codon = codon.replace("-", "N")
    try:
        aa = str(Seq(codon).translate(table=11))
        return aa if aa not in ("*", "X") else None
    except Exception:
        return None


def get_alt_aas(aln_dict, site, ref_aa):
    nt_start = (site - 1) * 3
    aas = set()
    for seq in aln_dict.values():
        codon = seq[nt_start:nt_start + 3]
        if len(codon) < 3:
            continue
        aa = codon_to_aa(codon)
        if aa and aa != ref_aa:
            aas.add(aa)
    return ",".join(sorted(aas)) if aas else "NA"


def compile_results(aln, outdir, group_label):
    beb_pp  = parse_beb_pp(outdir / "M8" / "codeml.out")
    beb_rst = parse_beb_rst(outdir / "M8" / "rst")
    fubar   = parse_fubar_json(outdir / f"{aln.name}.FUBAR.json")
    meme_f  = outdir / f"{aln.name}.MEME.json"
    meme    = parse_meme_json(meme_f) if meme_f.exists() else {}

    aln_dict = {r.id: str(r.seq) for r in SeqIO.parse(aln, "fasta")}
    ref_id   = next(iter(aln_dict))
    ref_seq  = aln_dict[ref_id]

    sig = (set(k for k, v in beb_pp.items() if v[0] >= BEB_THRESHOLD) |
           set(k for k, v in fubar.items()  if v["pp"] >= FUBAR_THRESHOLD))

    rows = []
    for site in sorted(sig):
        nt_start  = (site - 1) * 3
        ref_codon = ref_seq[nt_start:nt_start + 3] if nt_start + 3 <= len(ref_seq) else "---"
        ref_aa    = codon_to_aa(ref_codon) or "-"
        alts      = get_alt_aas(aln_dict, site, ref_aa)

        bp, bstar = beb_pp.get(site, (float("nan"), ""))
        pw, psd   = beb_rst.get(site, (float("nan"), float("nan")))
        fv        = fubar.get(site, {})
        mv        = meme.get(site, {})

        def fmt(v, n=4):
            return f"{v:.{n}f}" if v == v else "NA"

        rows.append({
            "group":                      group_label,
            "site":                       site,
            "ref_seq":                    ref_id,
            "ref_residue":                ref_aa,
            "alt_residues":               alts,
            "codeml_beb_pp":              fmt(bp) + bstar,
            "codeml_posterior_mean_ω":    fmt(pw),
            "codeml_posterior_mean_ω_sd": fmt(psd),
            "fubar_pp_pos":               fmt(fv.get("pp", float("nan"))),
            "fubar_bayes_factor":         fmt(fv.get("bf", float("nan"))),
            "fubar_dS":                   fmt(fv.get("dS", float("nan"))),
            "fubar_dN":                   fmt(fv.get("dN", float("nan"))),
            "fubar_omega":                fmt(fv.get("omega", float("nan"))),
            "meme_pval":                  fmt(mv.get("pval", float("nan"))),
            "meme_branches":              mv.get("branches", "NA"),
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description="Compile CODEML/FUBAR/MEME results into a summary TSV")
    parser.add_argument("--alignment", required=True, help="Codon alignment FASTA used for analyses")
    parser.add_argument("--outdir",    required=True, help="Directory containing CODEML model subdirs")
    parser.add_argument("--group",     default="group", help="Group label for output table")
    parser.add_argument("--output",    default="selection_summary.tsv", help="Output TSV path")
    args = parser.parse_args()

    t0  = time.time()
    aln = Path(args.alignment)
    outdir = Path(args.outdir)

    print(f"[{SCRIPT}] Compiling results for group '{args.group}' ...", flush=True)
    rows = compile_results(aln, outdir, args.group)

    out_path = Path(args.output)
    with open(out_path, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in COLS) + "\n")

    print(f"[{SCRIPT}] {len(rows)} significant sites → {out_path} ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
