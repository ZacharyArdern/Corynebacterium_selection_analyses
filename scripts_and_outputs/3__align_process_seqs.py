#!/usr/bin/env python3
"""
align_process_seqs.py
Steps:
  0. (optional) Extract tox sequences from genome files using blastx coords TSV
  1. Translate nucleotide sequences to protein (codon table 11)
  2. Align proteins with MAFFT
  3. Filter: remove sequences with internal stop codons or < 500 aa (1500 nt)
  4. Back-translate protein alignment to codon alignment with pal2nal
  5. Dereplicate: remove identical nt sequences and sequences embedded in longer ones

Usage:
  # From pre-extracted FASTA:
  python3 align_process_seqs.py --input seqs.fa --outdir results/ [--prefix prefix]

  # From blastx coords + genome directory:
  python3 align_process_seqs.py --coords ATB_tox_coords.tsv --genomes ATB_genomes/ \
      --outdir results/ [--prefix prefix]
"""

import argparse
import csv
import gzip
import subprocess
import sys
import time
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

CODON_TABLE = 11
PAL2NAL = "pal2nal.pl"
MAFFT = "mafft"
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


def extract_from_coords(coords_tsv, genomes_dir):
    """Extract nucleotide sequences from genome files using blastx hit coordinates."""
    genomes_dir = Path(genomes_dir)

    # Index all genome files by accession
    print("Indexing genome files ...")
    genome_files = {}
    for fa in list(genomes_dir.rglob("*.fa")) + list(genomes_dir.rglob("*.fa.gz")):
        acc = fa.stem.split(".fa")[0]
        genome_files[acc] = fa
    print(f"  Found {len(genome_files)} genome files")

    # Parse coords and extract
    records = []
    missing = []
    with open(coords_tsv) as f:
        for line in f:
            fields = line.strip().split("\t")
            contig_id = fields[0]  # e.g. SAMD00732958.contig00020
            pident, length = float(fields[2]), int(fields[3])
            qstart, qend = int(fields[4]), int(fields[5])
            acc = contig_id.split(".contig")[0]

            fa_path = genome_files.get(acc)
            if fa_path is None:
                missing.append(acc)
                continue

            opener = gzip.open if str(fa_path).endswith(".gz") else open
            with opener(fa_path, "rt") as fh:
                genome = {r.id: r for r in SeqIO.parse(fh, "fasta")}

            # Find contig — genome headers may have extra info after space
            contig = None
            for rid, rec in genome.items():
                if rid == contig_id or rid.startswith(contig_id + " "):
                    contig = rec
                    break
            if contig is None:
                # Try matching just the contig part
                contig_suffix = contig_id.split(".", 1)[1] if "." in contig_id else contig_id
                for rid, rec in genome.items():
                    if contig_suffix in rid:
                        contig = rec
                        break
            if contig is None:
                missing.append(contig_id)
                continue

            # Extract region, reverse complement if on minus strand
            start, end = min(qstart, qend) - 1, max(qstart, qend)
            seq = contig.seq[start:end]
            if qstart > qend:
                seq = seq.reverse_complement()

            rec = SeqRecord(seq, id=contig_id, description="")
            records.append(rec)

    print(f"  Extracted {len(records)} sequences ({len(missing)} not found)")
    return records


def translate_seqs(records, table=CODON_TABLE):
    proteins = []
    for rec in records:
        aa_str = str(rec.seq.translate(table=table)).rstrip("*")
        new = rec[:]
        new.seq = Seq(aa_str)
        proteins.append(new)
    return proteins


def has_internal_stop(rec, table=CODON_TABLE):
    aa = str(rec.seq.translate(table=table))
    return "*" in aa[:-1]


def run_mafft(input_fa, output_fa, threads=1):
    with open(output_fa, "w") as f:
        subprocess.run([MAFFT, "--auto", "--quiet", "--thread", str(threads), str(input_fa)],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)


def run_pal2nal(prot_aln, nuc_fa, output_fa):
    result = subprocess.run(
        [PAL2NAL, str(prot_aln), str(nuc_fa), "-codontable", str(CODON_TABLE), "-output", "fasta"],
        capture_output=True, text=True, check=True
    )
    with open(output_fa, "w") as f:
        f.write(result.stdout)


def dereplicate(records):
    seen = {}
    for rec in records:
        s = str(rec.seq).replace("-", "")
        if s not in seen:
            seen[s] = rec
    unique = list(seen.values())

    ungapped = [str(r.seq).replace("-", "") for r in unique]
    unique_sorted = sorted(zip(ungapped, unique), key=lambda x: len(x[0]), reverse=True)
    ungapped = [s for s, _ in unique_sorted]
    unique   = [r for _, r in unique_sorted]

    keep = []
    for i, (s, rec) in enumerate(zip(ungapped, unique)):
        embedded = any(len(ungapped[j]) > len(s) and s in ungapped[j]
                       for j in range(len(ungapped)) if j != i)
        if not embedded:
            keep.append(rec)
    return keep


def main():
    parser = argparse.ArgumentParser(description="Align and process nucleotide sequences for selection analysis")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input",  help="Pre-extracted nucleotide FASTA")
    group.add_argument("--coords", help="blastx coords TSV from 2_search_tox_gene.py")
    parser.add_argument("--genomes", help="Genome directory (required with --coords)")
    parser.add_argument("--outdir",  required=True, help="Output directory")
    parser.add_argument("--prefix",  default="seqs", help="Output file prefix")
    parser.add_argument("--threads", type=int, default=1, help="Threads for MAFFT")
    args = parser.parse_args()

    if args.coords and not args.genomes:
        parser.error("--genomes is required when using --coords")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    px = args.prefix

    print(f"[{SCRIPT}] Starting", flush=True)

    # Step 0: Extract sequences if coords provided
    if args.coords:
        step("Step 0: Extracting sequences from genome files ...")
        records = extract_from_coords(args.coords, args.genomes)
        done(f"{len(records)} sequences extracted")
    else:
        records = list(SeqIO.parse(args.input, "fasta"))
        print(f"[{SCRIPT}] Input: {len(records)} sequences")

    # Step 1 & 2: Translate + align proteins with MAFFT
    step("Step 1-2: Translating and aligning proteins with MAFFT ...")
    proteins = translate_seqs(records)
    prot_fa = outdir / f"{px}_proteins.fa"
    SeqIO.write(proteins, prot_fa, "fasta")
    prot_aln = outdir / f"{px}_proteins_aln.fa"
    run_mafft(prot_fa, prot_aln, args.threads)
    done(f"{sum(1 for _ in SeqIO.parse(prot_aln, 'fasta'))} sequences aligned")

    # Step 3: Filter
    clean_nuc = [r for r in records if not has_internal_stop(r)]
    n_stops = len(records) - len(clean_nuc)
    clean_nuc = [r for r in clean_nuc if len(r.seq) >= 1500]
    n_short = len(records) - n_stops - len(clean_nuc)
    print(f"[{SCRIPT}] Step 3: Filtered {n_stops} internal stops, {n_short} short (<500 aa) → {len(clean_nuc)} remaining")

    clean_prot = translate_seqs(clean_nuc)
    clean_prot_fa = outdir / f"{px}_proteins_clean.fa"
    SeqIO.write(clean_prot, clean_prot_fa, "fasta")
    clean_prot_aln = outdir / f"{px}_proteins_clean_aln.fa"
    step("Step 3 (cont): Re-aligning filtered proteins ...")
    run_mafft(clean_prot_fa, clean_prot_aln, args.threads)
    done()

    clean_nuc_fa = outdir / f"{px}_nuc_clean.fa"
    SeqIO.write(clean_nuc, clean_nuc_fa, "fasta")

    # Step 4: Back-translate to codon alignment with pal2nal
    step("Step 4: Back-translating to codon alignment with pal2nal ...")
    codon_aln = outdir / f"{px}_codon_aln.fa"
    run_pal2nal(clean_prot_aln, clean_nuc_fa, codon_aln)
    done(f"{sum(1 for _ in SeqIO.parse(codon_aln, 'fasta'))} sequences")

    # Step 5: Dereplicate
    step("Step 5: Dereplicating ...")
    codon_records = list(SeqIO.parse(codon_aln, "fasta"))
    derep = dereplicate(codon_records)
    done(f"{len(codon_records)} → {len(derep)} unique non-embedded sequences")

    derep_fa = outdir / f"{px}_codon_aln_derep.fa"
    SeqIO.write(derep, derep_fa, "fasta")

    clean_fa = outdir / f"{px}_codon_aln_derep_clean.fa"
    with open(clean_fa, "w") as f:
        for rec in derep:
            f.write(f">{rec.id}\n{rec.seq}\n")

    print(f"[{SCRIPT}] Final alignment: {clean_fa}")
    print(f"[{SCRIPT}] Finished in {time.time() - _script_start:.1f}s")


if __name__ == "__main__":
    main()
