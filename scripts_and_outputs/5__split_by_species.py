#!/usr/bin/env python3
"""
Split a combined codon alignment into per-species FASTAs.
Uses the ATB file list TSV for species assignments.
Unknown-species samples are assigned to the species of their nearest
labelled neighbour in a FastTree GTR phylogeny of the alignment.

Usage:
  python3 5__split_by_species.py \
      --alignment results_ATB_new/ATB_new_codon_aln_derep_clean.fa \
      --filelist  ATB_file_list.all.latest.tsv.gz \
      --tree      results_ATB_new/ATB_new.tree \
      --outdir    results_ATB_new \
      --species   diphtheriae ulcerans
"""

import argparse
import csv
import gzip
import time
from pathlib import Path

from Bio import Phylo, SeqIO

SCRIPT = Path(__file__).name

def step(msg):
    step._t = time.time()
    print(f"[{SCRIPT}] {msg}", flush=True)

def done(msg=""):
    elapsed = time.time() - step._t
    suffix = f" — {msg}" if msg else ""
    print(f"[{SCRIPT}]   done ({elapsed:.1f}s){suffix}", flush=True)

step._t = time.time()
_script_start = time.time()


def load_species(filelist):
    species = {}
    opener = gzip.open if str(filelist).endswith(".gz") else open
    with opener(filelist, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species[row["sample"]] = row["sylph_species"]
    return species


def acc_from_id(seq_id):
    return seq_id.split(".contig")[0]


def assign_unknowns_from_tree(tree, species_map):
    """
    For tips whose species is 'unknown' or absent, find the nearest tip
    (by sum of branch lengths) that has a known species and use that.
    Returns an updated copy of species_map.
    """
    assigned = dict(species_map)
    tips = tree.get_terminals()

    known_tips   = [t for t in tips
                    if assigned.get(acc_from_id(t.name), "unknown") != "unknown"
                    and acc_from_id(t.name) in assigned]
    unknown_tips = [t for t in tips if t not in set(known_tips)]

    print(f"[{SCRIPT}]   Known: {len(known_tips)}, Unknown (to assign): {len(unknown_tips)}")

    for unk in unknown_tips:
        best_sp, best_dist = None, float("inf")
        for kn in known_tips:
            try:
                dist = tree.distance(unk, kn)
            except Exception:
                continue
            if dist < best_dist:
                best_dist = dist
                best_sp = assigned[acc_from_id(kn.name)]
        acc = acc_from_id(unk.name)
        assigned[acc] = best_sp
        sp_short = best_sp.split()[-1] if best_sp else "?"
        print(f"[{SCRIPT}]     {unk.name} → {sp_short} (dist={best_dist:.6f})")

    return assigned


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, help="Codon alignment FASTA")
    parser.add_argument("--filelist",  required=True, help="ATB file list TSV(.gz)")
    parser.add_argument("--tree",      required=True, help="Newick tree (rooted)")
    parser.add_argument("--outdir",    required=True, help="Output directory")
    parser.add_argument("--prefix",    default="ATB", help="Output file prefix")
    parser.add_argument("--species",   nargs="+",
                        default=["diphtheriae", "ulcerans"],
                        help="Species epithet(s) to output (e.g. diphtheriae ulcerans)")
    args = parser.parse_args()

    print(f"[{SCRIPT}] Starting", flush=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_species = load_species(args.filelist)
    species_map = {acc: (sp if sp else "unknown") for acc, sp in raw_species.items()}

    step("Assigning unknown-species tips from tree neighbours ...")
    tree = Phylo.read(args.tree, "newick")
    species_map = assign_unknowns_from_tree(tree, species_map)
    done()

    buckets = {}
    for rec in SeqIO.parse(args.alignment, "fasta"):
        acc = acc_from_id(rec.id)
        sp = species_map.get(acc, "unknown")
        sp_short = sp.split()[-1] if sp and sp != "unknown" else "unknown"
        buckets.setdefault(sp_short, []).append(rec)

    for sp_short, recs in sorted(buckets.items()):
        if args.species and sp_short not in args.species:
            print(f"[{SCRIPT}]   Skipping {sp_short} ({len(recs)} seqs) — not in --species list")
            continue
        out = outdir / f"{args.prefix}_{sp_short}_subset.fa"
        SeqIO.write(recs, out, "fasta")
        print(f"[{SCRIPT}]   {sp_short}: {len(recs)} sequences → {out}")

    print(f"[{SCRIPT}] Finished in {time.time() - _script_start:.1f}s")


if __name__ == "__main__":
    main()
