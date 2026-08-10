#!/usr/bin/env bash
set -euo pipefail

SCRIPTS="scripts_and_outputs"
RESULTS="results_ATB_new"
PREFIX="ATB_new"
OUTGROUP="SAMN07109093"
THREADS=8          # default; override with --threads N on the command line

# Parse --threads argument
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ── Dependency checks ────────────────────────────────────────────────────────
echo "=== Checking dependencies ==="
MISSING=0

check() {
    local name=$1; shift
    if "$@" &>/dev/null 2>&1; then
        printf "  %-22s ✓\n" "$name"
    else
        printf "  %-22s ✗  (not found)\n" "$name"
        MISSING=$((MISSING + 1))
    fi
}

check python3          python3 --version
check diamond          diamond --version
check orfipy           orfipy --version
check mafft            mafft --version
check pal2nal          pal2nal.pl --help
check FastTree         FastTree -help
check skani            skani --version
check iqtree2          iqtree2 --version
check codeml           bash -c "command -v codeml"
check hyphy            hyphy --version

if [ "$MISSING" -gt 0 ]; then
    echo ""
    echo "ERROR: $MISSING required program(s) not found. Aborting." >&2
    exit 1
fi
echo ""

# ── Helper ───────────────────────────────────────────────────────────────────
_STEP_T=0
step() {
    _STEP_T=$SECONDS
    echo ""
    echo "=== [$(date '+%H:%M:%S')] $1 ==="
}
done_step() {
    echo "=== done in $((SECONDS - _STEP_T))s ==="
}

echo "=== ATB tox selection pipeline — $(date) ==="

# ── Step 1 ───────────────────────────────────────────────────────────────────
step "1/6  Download ATB genomes"
python "$SCRIPTS/1__get_seqs_from_ATB.py"
done_step

# ── Step 2 ───────────────────────────────────────────────────────────────────
step "2/6  orfipy + DIAMOND blastp tox search"
python "$SCRIPTS/2__search_tox_gene.py" --threads "$THREADS"
done_step

# ── Step 3 ───────────────────────────────────────────────────────────────────
step "3/6  Align, filter, dereplicate"
python "$SCRIPTS/3__align_process_seqs.py" \
    --coords ATB_tox_coords.tsv \
    --genomes ATB_genomes \
    --outdir "$RESULTS" \
    --prefix "$PREFIX" \
    --threads "$THREADS"
done_step

# ── Step 4 ───────────────────────────────────────────────────────────────────
step "4/6  FastTree phylogeny + skani ANI"
python "$SCRIPTS/4__build_tree_and_ani.py" \
    --alignment "$RESULTS/${PREFIX}_codon_aln_derep_clean.fa" \
    --genomes ATB_genomes \
    --filelist ATB_file_list.all.latest.tsv.gz \
    --outdir "$RESULTS" \
    --prefix "$PREFIX" \
    --outgroup "$OUTGROUP" \
    --threads "$THREADS"
done_step

# ── Step 5 ───────────────────────────────────────────────────────────────────
step "5/6  Split alignment by species"
python "$SCRIPTS/5__split_by_species.py" \
    --alignment "$RESULTS/${PREFIX}_codon_aln_derep_clean.fa" \
    --filelist ATB_file_list.all.latest.tsv.gz \
    --tree "$RESULTS/${PREFIX}.tree" \
    --outdir "$RESULTS" \
    --prefix "$PREFIX" \
    --species diphtheriae ulcerans
done_step

# ── Step 6+7 (parallel) ──────────────────────────────────────────────────────
# Script 6 runs IQ-TREE, CODEML, FUBAR, MEME, then calls script 7
# (7__compile_results_table.py) internally to produce the summary TSV.
step "6-7/7  Selection analyses + results table (diphtheriae + ulcerans in parallel)"
python "$SCRIPTS/6__run_selection_analyses.py" \
    --alignment "$RESULTS/${PREFIX}_diphtheriae_subset.fa" \
    --outdir "$RESULTS/diphtheriae" \
    --prefix "${PREFIX}_diphtheriae" \
    --group diphtheriae \
    --output "$SCRIPTS/selection_summary_diphtheriae.tsv" \
    --threads "$((THREADS / 2))" &

python "$SCRIPTS/6__run_selection_analyses.py" \
    --alignment "$RESULTS/${PREFIX}_ulcerans_subset.fa" \
    --outdir "$RESULTS/ulcerans" \
    --prefix "${PREFIX}_ulcerans" \
    --group ulcerans \
    --output "$SCRIPTS/selection_summary_ulcerans.tsv" \
    --threads "$((THREADS / 2))" &

wait
done_step

echo ""
echo "=== Pipeline complete — $(date) ==="
echo ""

echo "(No large temporary files to clean up.)"
