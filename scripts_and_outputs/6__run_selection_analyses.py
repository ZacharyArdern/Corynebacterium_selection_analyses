#!/usr/bin/env python3
"""
run_selection_analyses.py
Steps (from aligned, dereplicated codon alignment):
  1. Infer phylogeny with IQ-TREE2
  2. Run CODEML site models M0, M1a, M2a, M7, M8
  3. Run FUBAR (HyPhy)
  4. Run MEME (HyPhy)
  5. Parse and compile results into summary TSV

Usage:
  python3 run_selection_analyses.py --alignment seqs_codon_aln_derep_clean.fa --outdir results/ [--prefix prefix]

Requirements:
  - iqtree2 in PATH
  - codeml (PAML) in PATH
  - HyPhy installed in conda env 'hyphy' (conda run -n hyphy hyphy)
"""

import argparse
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HYPHY_CMD = ["hyphy"]
CODEML = "codeml"
IQTREE = "iqtree2"
CODON_TABLE = 11
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

# ── Step 1: IQ-TREE ──────────────────────────────────────────────────────────

def run_iqtree(aln, outdir, prefix, threads=1):
    treefile = outdir / f"{prefix}.treefile"
    if treefile.exists():
        return treefile
    subprocess.run([
        IQTREE, "-s", str(aln), "-T", str(threads), "-m", "TEST",
        "--prefix", str(outdir / prefix), "-B", "1000", "--quiet"
    ], check=True)
    return treefile


# ── Step 2: CODEML ───────────────────────────────────────────────────────────

MODELS = {"M0": 0, "M1a": 1, "M2a": 2, "M7": 7, "M8": 8}

def _run_one_codeml(model, nssite, aln, treefile, outdir):
    mdir = outdir / model
    mdir.mkdir(exist_ok=True)
    out = mdir / "codeml.out"
    if out.exists() and out.stat().st_size > 0:
        return
    (mdir / "codeml.ctl").write_text(f"""\
seqfile  = {aln.resolve()}
treefile = {treefile.resolve()}
outfile  = {out.resolve()}
noisy    = 0
verbose  = 0
runmode  = 0
seqtype  = 1
CodonFreq = 2
model    = 0
NSsites  = {nssite}
icode    = 0
fix_kappa = 0
kappa    = 2
fix_omega = 0
omega    = 0.5
getSE    = 0
RateAncestor = 1
cleandata = 0
""")
    subprocess.run([CODEML], cwd=str(mdir),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_codeml(aln, treefile, outdir, threads=4):
    with ThreadPoolExecutor(max_workers=min(len(MODELS), threads)) as ex:
        futures = {ex.submit(_run_one_codeml, m, n, aln, treefile, outdir): m
                   for m, n in MODELS.items()}
        for f in futures:
            f.result()  # re-raise any exceptions


# ── Step 3: FUBAR ────────────────────────────────────────────────────────────

def run_fubar(aln, treefile, outdir, threads=1):
    json_out = outdir / f"{aln.name}.FUBAR.json"
    if json_out.exists():
        return json_out
    log = outdir / "fubar.log"
    with open(log, "w") as f:
        subprocess.run(
            HYPHY_CMD + ["fubar", "--alignment", str(aln), "--tree", str(treefile),
                         "--output", str(json_out), "--cpu", str(threads)],
            stdout=f, stderr=f, check=True
        )
    return json_out


# ── Step 4: MEME ─────────────────────────────────────────────────────────────

def run_meme(aln, treefile, outdir, threads=1):
    json_out = outdir / f"{aln.name}.MEME.json"
    log = outdir / "meme.log"
    if json_out.exists():
        return json_out
    with open(log, "w") as f:
        subprocess.run(
            HYPHY_CMD + ["meme", "--alignment", str(aln), "--tree", str(treefile),
                         "--output", str(json_out), "--cpu", str(threads)],
            stdout=f, stderr=f, check=True
        )
    return json_out




# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Run selection analyses on a codon alignment")
    parser.add_argument("--alignment", required=True, help="Codon alignment FASTA (dereplicated, clean headers)")
    parser.add_argument("--outdir",    required=True, help="Output directory")
    parser.add_argument("--prefix",    default="analysis", help="Output file prefix for IQ-TREE and CODEML outputs")
    parser.add_argument("--group",     default="group", help="Group label for output table")
    parser.add_argument("--output",    default="selection_summary.tsv", help="Summary TSV output path")
    parser.add_argument("--threads",   type=int, default=4, help="Threads for IQ-TREE, CODEML (model parallelism), and HyPhy")
    args = parser.parse_args()

    aln    = Path(args.alignment)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[{SCRIPT}] Starting (group: {args.group}, threads: {args.threads})", flush=True)

    step("[1] IQ-TREE")
    treefile = run_iqtree(aln, outdir, args.prefix, args.threads)
    done()

    step("[2] CODEML (5 models in parallel)")
    run_codeml(aln, treefile, outdir, args.threads)
    done()

    step("[3] FUBAR")
    run_fubar(aln, treefile, outdir, args.threads)
    done()

    step("[4] MEME")
    run_meme(aln, treefile, outdir, args.threads)
    done()

    step("[5] Compiling results")
    compile_script = Path(__file__).parent / "7__compile_results_table.py"
    subprocess.run([
        sys.executable, str(compile_script),
        "--alignment", str(aln),
        "--outdir",    str(outdir),
        "--group",     args.group,
        "--output",    args.output,
    ], check=True)
    done()

    print(f"[{SCRIPT}] Finished in {time.time() - _script_start:.1f}s")


if __name__ == "__main__":
    main()
