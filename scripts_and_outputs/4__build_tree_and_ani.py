#!/usr/bin/env python3
"""
Build a FastTree GTR phylogeny from a codon alignment, run skani all-vs-all ANI
on the corresponding genomes, and plot the tree alongside the ANI heatmap.

Usage:
  python3 4__build_tree_and_ani.py \
      --alignment results_ATB_new/ATB_new_codon_aln_derep_clean.fa \
      --genomes   ATB_genomes \
      --filelist  ATB_file_list.all.latest.tsv.gz \
      --outdir    results_ATB_new \
      --prefix    ATB_new \
      --outgroup  SAMN07109093
"""

import argparse
import csv
import gzip
import subprocess
import time
from pathlib import Path

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

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SP_COLORS = {
    "Corynebacterium diphtheriae":        "#2166ac",
    "Corynebacterium ulcerans":           "#d6604d",
    "Corynebacterium silvaticum":         "#4dac26",
    "Corynebacterium pseudotuberculosis": "#984ea3",
    "Corynebacterium rouxii":             "#ff7f00",
}
DEFAULT_COLOR = "#888888"


def load_species(filelist):
    species = {}
    opener = gzip.open if str(filelist).endswith(".gz") else open
    with opener(filelist, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species[row["sample"]] = row["sylph_species"]
    return species


def build_fasttree(aln, outdir, prefix, threads=1):
    treefile = outdir / f"{prefix}.tree"
    step("Building FastTree GTR phylogeny ...")
    env = {**__import__("os").environ, "OMP_NUM_THREADS": str(threads)}
    with open(treefile, "w") as fh:
        subprocess.run(["FastTree", "-nt", "-gtr", str(aln)],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True, env=env)
    done(str(treefile))
    return treefile


def run_skani(aln, genomes_dir, outdir, threads=1):
    accs = set()
    for rec in SeqIO.parse(aln, "fasta"):
        accs.add(rec.id.split(".contig")[0])

    genome_files = {}
    for fa in Path(genomes_dir).rglob("*.fa"):
        acc = fa.stem.split(".fa")[0]
        if acc in accs:
            genome_files[acc] = str(fa)
    print(f"[{SCRIPT}] {len(genome_files)}/{len(accs)} genome files found for skani")

    list_file = outdir / "skani_query_list.txt"
    with open(list_file, "w") as f:
        for path in genome_files.values():
            f.write(path + "\n")

    matrix_file = outdir / "skani_matrix.txt"
    step("Running skani all-vs-all ANI ...")
    subprocess.run(["skani", "triangle", "-l", str(list_file),
                    "-o", str(matrix_file), "-t", str(threads)],
                   check=True, capture_output=True)
    done(str(matrix_file))
    return matrix_file


def load_skani_matrix(matrix_file):
    with open(matrix_file) as f:
        lines = f.read().splitlines()
    n = int(lines[0])
    labels, rows = [], []
    for line in lines[1:]:
        parts = line.split("\t")
        labels.append(Path(parts[0]).stem)
        rows.append([float(x) if x != "N/A" else np.nan for x in parts[1:]])
    mat = np.full((n, n), np.nan)
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            mat[i][j] = v
            mat[j][i] = v
    np.fill_diagonal(mat, 100.0)
    return pd.DataFrame(mat, index=labels, columns=labels)


def plot(treefile, df_ani, species, outdir, prefix, outgroup):
    tree = Phylo.read(treefile, "newick")
    if outgroup:
        og = next((c for c in tree.get_terminals() if c.name.startswith(outgroup)), None)
        if og is None:
            raise ValueError(f"Outgroup '{outgroup}' not found in tree")
        tree.root_with_outgroup(og)
    else:
        tree.root_at_midpoint()
    tree.ladderize()

    tip_order = [t.name for t in tree.get_terminals()]
    tip_accs  = [t.split(".contig")[0] for t in tip_order]
    df_ordered = df_ani.reindex(index=tip_accs, columns=tip_accs)

    def tip_label(clade):
        if not clade.name:
            return ""
        acc = clade.name.split(".contig")[0]
        sp = species.get(acc, "unknown")
        return f"{clade.name} ({sp.replace('Corynebacterium ', 'C. ')})"

    fig = plt.figure(figsize=(18, 14))
    ax_tree = fig.add_axes([0.01, 0.05, 0.30, 0.88])
    ax_heat = fig.add_axes([0.45, 0.05, 0.46, 0.88])
    ax_cbar = fig.add_axes([0.93, 0.05, 0.02, 0.88])

    Phylo.draw(tree, axes=ax_tree, do_show=False, label_func=tip_label)
    for text in ax_tree.texts:
        t = text.get_text().strip()
        if not t:
            continue
        acc = t.split(" (")[0].split(".contig")[0]
        sp = species.get(acc, "unknown")
        text.set_color(SP_COLORS.get(sp, DEFAULT_COLOR))
        text.set_fontsize(5)
    ax_tree.set_title("FastTree GTR", fontsize=9)
    ax_tree.set_xlabel("branch length", fontsize=7)
    ax_tree.tick_params(labelsize=6)

    ani_vals = df_ordered.values
    vmin = np.nanmin(ani_vals)
    im = ax_heat.imshow(ani_vals, aspect="auto", cmap="RdYlGn",
                        vmin=vmin, vmax=100, origin="upper")
    labels = [f"{a} ({species.get(a,'?').replace('Corynebacterium ','C. ')})"
              for a in tip_accs]
    ax_heat.set_xticks(range(len(tip_order)))
    ax_heat.set_xticklabels(labels, rotation=90, fontsize=4)
    ax_heat.set_yticks([])
    ax_heat.set_yticklabels([])
    ax_heat.set_title("skani ANI (%)", fontsize=9)

    plt.colorbar(im, cax=ax_cbar)
    ax_cbar.tick_params(labelsize=7)
    ax_cbar.set_ylabel("ANI (%)", fontsize=7)

    patches = [mpatches.Patch(color=c, label=sp.replace("Corynebacterium ", "C. "))
               for sp, c in SP_COLORS.items()]
    patches.append(mpatches.Patch(color=DEFAULT_COLOR, label="unknown"))
    ax_tree.legend(handles=patches, fontsize=6, loc="lower left")

    out = outdir / f"{prefix}_tree_vs_ANI.pdf"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    done(str(out))

    # Also save rerooted tree
    Phylo.write(tree, outdir / f"{prefix}.tree", "newick")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--genomes",   required=True)
    parser.add_argument("--filelist",  required=True, help="ATB file list TSV(.gz)")
    parser.add_argument("--outdir",    required=True)
    parser.add_argument("--prefix",    default="ATB")
    parser.add_argument("--outgroup",  default="", help="Accession prefix to root on")
    parser.add_argument("--threads",   type=int, default=1, help="Threads for FastTree and skani")
    args = parser.parse_args()

    print(f"[{SCRIPT}] Starting", flush=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species = load_species(args.filelist)
    treefile = build_fasttree(args.alignment, outdir, args.prefix, args.threads)
    matrix_file = run_skani(args.alignment, args.genomes, outdir, args.threads)
    df_ani = load_skani_matrix(matrix_file)
    step("Plotting tree vs ANI heatmap ...")
    plot(treefile, df_ani, species, outdir, args.prefix, args.outgroup)
    print(f"[{SCRIPT}] Finished in {time.time() - _script_start:.1f}s")


if __name__ == "__main__":
    main()
