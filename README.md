# ATB Toxin Selection Analysis — Scripts and Outputs

## Quick start

```bash
git clone https://github.com/ZacharyArdern/Corynebacterium_selection_analyses
cd Corynebacterium_selection_analyses
conda create -n selection_analysis -c bioconda -c conda-forge \
    python diamond orfipy mafft pal2nal fasttree skani iqtree paml hyphy \
    biopython matplotlib pandas numpy
conda activate selection_analysis
bash run_pipeline.sh --threads 8
```

---

## Input files

| File | Description |
|---|---|
| `tox.fas` | Nucleotide sequences of tox alleles; tox_2 is extracted as the DIAMOND search reference. Downloaded automatically from https://gitlab.pasteur.fr/BEBP/diphtoscan/-/raw/main/diphtoscan/data/tox/sequences/tox.fas |

---

## Required programs

| Program | Purpose | Version used |
|---|---|---|
| DIAMOND | blastp search of 6FT-translated ORFs against tox_2 protein reference | 2.1.25 |
| orfipy | 6-frame translation of genome assemblies for ORF extraction | 0.0.4 |
| MAFFT | Protein multiple sequence alignment | 7.526 |
| pal2nal | Back-translation of protein alignment to codon alignment | v14 |
| FastTree | Rapid ML phylogeny for tree-based species assignment and visualisation | 2.1.11 |
| skani | Whole-genome ANI estimation | 0.3.2 |
| IQ-TREE2 | Maximum likelihood phylogeny inference for selection analyses | 2.3.6 |
| PAML (codeml) | Site model selection analyses (M0, M1a, M2a, M7, M8) | 4.10.10 |
| HyPhy | FUBAR and MEME selection analyses | 2.5.101 |

### Installation

All programs can be installed via conda. Create a dedicated environment and install all tools in one step:

```bash
conda create -n selection_analysis -c bioconda -c conda-forge \
    python diamond orfipy mafft pal2nal fasttree skani iqtree paml hyphy \
    biopython matplotlib pandas numpy
conda activate selection_analysis
```

---

## Running scripts

The full pipeline can be run end-to-end with:

```bash
bash scripts_and_outputs/run_pipeline.sh
```

Or run individual steps manually:

```bash
# 1. Download ATB genomes (Corynebacterium)
python scripts_and_outputs/1__get_seqs_from_ATB.py

# 2. Search genomes for tox gene using orfipy + DIAMOND blastp
python scripts_and_outputs/2__search_tox_gene.py

# 3. Extract, align, filter (≥500 aa), and dereplicate sequences
python scripts_and_outputs/3__align_process_seqs.py \
    --coords ATB_tox_coords.tsv --genomes ATB_genomes \
    --outdir results_ATB_new --prefix ATB_new

# 4. Build FastTree phylogeny and skani ANI matrix; plot
python scripts_and_outputs/4__build_tree_and_ani.py \
    --alignment results_ATB_new/ATB_new_codon_aln_derep_clean.fa \
    --genomes ATB_genomes --filelist ATB_file_list.all.latest.tsv.gz \
    --outdir results_ATB_new --prefix ATB_new --outgroup SAMN07109093

# 5. Split alignment by species (assigns unknowns from tree context)
python scripts_and_outputs/5__split_by_species.py \
    --alignment results_ATB_new/ATB_new_codon_aln_derep_clean.fa \
    --filelist ATB_file_list.all.latest.tsv.gz \
    --tree results_ATB_new/ATB_new.tree \
    --outdir results_ATB_new --species diphtheriae ulcerans

# 6. Run selection analyses (run separately for each group)
python scripts_and_outputs/6__run_selection_analyses.py \
    --alignment results_ATB_new/ATB_diphtheriae_subset.fa \
    --outdir results_ATB_new/diphtheriae --prefix ATB_diphtheriae \
    --group diphtheriae --output scripts_and_outputs/selection_summary_diphtheriae.tsv

python scripts_and_outputs/6__run_selection_analyses.py \
    --alignment results_ATB_new/ATB_ulcerans_subset.fa \
    --outdir results_ATB_new/ulcerans --prefix ATB_ulcerans \
    --group ulcerans --output scripts_and_outputs/selection_summary_ulcerans.tsv

# (Optional) Re-run table compilation only, without re-running analyses
python scripts_and_outputs/7__compile_results_table.py \
    --alignment results_ATB_new/ATB_diphtheriae_subset.fa \
    --outdir results_ATB_new/diphtheriae \
    --group diphtheriae --output scripts_and_outputs/selection_summary_diphtheriae.tsv
```

---

## Pipeline overview

This analysis was run on 7th–8th August 2026. Note that ATB results may change as the
database is updated.

On this date, ~9,258 Corynebacterium genus genomes were obtained from ATB. Sequences
were extracted by 6-frame translation of all genomes with orfipy followed by DIAMOND
blastp against the tox_2 protein reference (1,761 hits). Sequences shorter than 500 aa
were excluded. After alignment and dereplication, 65 unique sequences remained.

A FastTree GTR phylogeny and skani whole-genome ANI matrix were used to verify
species assignments. Samples with unknown species in the ATB file list were assigned
to the nearest labelled neighbour in the tox gene tree. The 65 sequences were split
into 31 *C. diphtheriae* and 29 *C. ulcerans* for separate selection analyses.

Selection analyses were run on each species subset: IQ-TREE2 phylogeny inference
followed by CODEML site models (M0, M1a, M2a, M7, M8), FUBAR, and MEME.

---

## Scripts

### 1. `1__get_seqs_from_ATB.py`
Downloads the ATB genome file list from OSF, filters to Corynebacterium genus entries,
and downloads and unpacks all genome tar.xz batch archives into `ATB_genomes/`.

Output: `ATB_genomes/` directory, `ATB_file_list.all.latest.tsv.gz`

---

### 2. `2__search_tox_gene.py`
Downloads `tox.fas` from the diphtoscan GitLab repository, extracts the tox_2 protein
sequence, builds a DIAMOND database, then 6-frame translates all ATB genomes in batches
using orfipy (≥400 aa ORFs, bacterial codon table) and runs DIAMOND blastp to identify
tox gene hits. ORF coordinates are converted back to genome nucleotide positions.

Output: `ATB_tox_coords.tsv` (hit coordinates: contig, sseqid, pident, length, gstart,
gend, sstart, send, slen, qcovhsp, evalue)

---

### 3. `3__align_process_seqs.py`
Extracts tox nucleotide sequences from genome files using tox hit coordinates, then
prepares them for selection analysis.

Steps:
1. Extract sequences from genome FASTA files using tox hit coordinates (reverse complement for minus-strand hits)
2. Translate sequences (genetic code 11); remove any with internal stop codons
3. Remove sequences shorter than 500 aa (1,500 nt)
4. Align proteins with MAFFT
5. Back-translate protein alignment to codon alignment with pal2nal
6. Dereplicate: remove identical sequences and sequences fully embedded in longer ones

Input: `ATB_tox_coords.tsv` + `ATB_genomes/` (1,761 hits)
Output: codon alignment of unique sequences (65 sequences after filtering and dereplication)

Usage:
```
python3 3__align_process_seqs.py --coords ATB_tox_coords.tsv --genomes ATB_genomes \
    --outdir results_ATB_new --prefix ATB_new
```

Can also accept a pre-extracted nucleotide FASTA directly via `--input <fasta>` instead of `--coords`/`--genomes`, skipping the extraction step.

---

### 4. `4__build_tree_and_ani.py`
Builds a FastTree GTR phylogeny from the codon alignment, runs skani all-vs-all ANI
on the corresponding genome assemblies, and plots the tree alongside the ANI heatmap
for visual concordance checking.

Input: codon alignment, genome directory, ATB file list
Output: `{prefix}.tree` (Newick, rerooted), `{prefix}_tree_vs_ANI.pdf`

Usage:
```
python3 4__build_tree_and_ani.py \
    --alignment results_ATB_new/ATB_new_codon_aln_derep_clean.fa \
    --genomes ATB_genomes --filelist ATB_file_list.all.latest.tsv.gz \
    --outdir results_ATB_new --prefix ATB_new --outgroup SAMN07109093  # replace with an appropriate outgroup accession
```

---

### 5. `5__split_by_species.py`
Splits the combined codon alignment into per-species FASTAs. Uses the ATB file list
for species assignments. Samples with unknown species are assigned to the species of
their nearest labelled neighbour in the tox gene phylogeny (by sum of branch lengths).

Input: codon alignment + ATB file list + FastTree newick tree
Output: one FASTA per requested species (e.g. `ATB_diphtheriae_subset.fa`, `ATB_ulcerans_subset.fa`)

Usage:
```
python3 5__split_by_species.py \
    --alignment results_ATB_new/ATB_new_codon_aln_derep_clean.fa \
    --filelist ATB_file_list.all.latest.tsv.gz \
    --tree results_ATB_new/ATB_new.tree \
    --outdir results_ATB_new --species diphtheriae ulcerans
```

---

### 6. `6__run_selection_analyses.py`
Runs the selection analysis pipeline on a codon alignment. Run separately for each
species group. All output goes to `--outdir`.

Steps:
1. Infer ML phylogeny with IQ-TREE2 (model selection, UFBoot)
2. Run CODEML site models M0, M1a, M2a, M7, M8 in parallel
3. Run FUBAR (HyPhy) — Bayesian per-site dN/dS estimation
4. Run MEME (HyPhy) — episodic selection test
5. Call `compile_results_table.py` to produce the summary TSV

Usage:
```
python3 6__run_selection_analyses.py --alignment <codon_aln.fa> --outdir <outdir> \
    --prefix <prefix> --group <label> --output <summary.tsv>
```

---

### `7__compile_results_table.py`
Parses CODEML (M8 BEB), FUBAR, and MEME output files and compiles significant sites
into a summary TSV. Called automatically by script 6, but can also be run standalone
to regenerate the table without re-running analyses.

Significance thresholds: CODEML BEB pp ≥ 0.80, FUBAR pp_pos ≥ 0.80, MEME p ≤ 0.05

Usage:
```
python3 compile_results_table.py --alignment <codon_aln.fa> --outdir <outdir> \
    --group <label> --output <summary.tsv>
```

---

## Outputs

| File | Description |
|---|---|
| `scripts_and_outputs/selection_summary_diphtheriae.tsv` | Significant sites from 31 *C. diphtheriae* sequences |
| `scripts_and_outputs/selection_summary_ulcerans.tsv` | Significant sites from 29 *C. ulcerans* sequences |

### TSV columns
| Column | Description |
|---|---|
| group | Analysis group (diphtheriae / ulcerans) |
| site | Amino acid position (1-based) |
| ref_seq | Reference sequence ID |
| ref_residue | Amino acid at this site in the reference sequence |
| alt_residues | Alternative amino acids observed at this site (NA if invariant) |
| codeml_beb_pp | CODEML BEB posterior probability of positive selection (M8) |
| codeml_posterior_mean_ω | CODEML posterior mean ω |
| codeml_posterior_mean_ω_sd | SD of posterior mean ω |
| fubar_pp_pos | FUBAR posterior probability that dN > dS |
| fubar_bayes_factor | FUBAR Bayes Factor for positive selection |
| fubar_dS | FUBAR posterior mean dS |
| fubar_dN | FUBAR posterior mean dN |
| fubar_omega | FUBAR dN/dS |
| meme_pval | MEME p-value for episodic selection |
| meme_branches | Number of branches under selection (MEME) |
