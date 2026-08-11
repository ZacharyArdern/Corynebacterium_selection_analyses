#!/usr/bin/env python3
"""
Extract tox_2 protein from tox.fas, 6FT-translate all ATB genomes with orfipy
in parallel batches, then run DIAMOND blastp to find tox hits. Converts hit
coordinates back to genome nucleotide positions in blastx-equivalent format
for script 3.

tox.fas source:
https://gitlab.pasteur.fr/BEBP/diphtoscan/-/raw/main/diphtoscan/data/tox/sequences/tox.fas
"""
import gzip
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import SeqIO

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

TOX_FAS_URL   = "https://gitlab.pasteur.fr/BEBP/diphtoscan/-/raw/main/diphtoscan/data/tox/sequences/tox.fas"
TOX_FAS       = "tox.fas"
TOX2_PROT_OUT = "tox_2_protein.fa"
DMND_DB       = "tox_2.dmnd"
GENOMES_DIR   = "ATB_genomes"
ORF_PEP       = "tmp_orfipy_pep.fa"
ORF_BED       = "tmp_orfipy.bed"
BLAST_RAW     = "tmp_blastp_raw.tsv"
BLAST_OUT     = "ATB_tox_coords.tsv"
BATCH_SIZE     = 300
MIN_ORF_NT     = 1200  # 400 aa minimum
ORFIPY_TABLE   = 9     # orfipy index 9 = NCBI transl_table=11 (Bacterial)


def download_tox_fas():
    if Path(TOX_FAS).exists():
        print(f"[{SCRIPT}] {TOX_FAS} already exists, skipping download")
        return
    step(f"Downloading {TOX_FAS} ...")
    urllib.request.urlretrieve(TOX_FAS_URL, TOX_FAS)
    done()


def extract_tox2_protein():
    for rec in SeqIO.parse(TOX_FAS, "fasta"):
        if rec.id == "tox_2":
            prot_str = str(rec.seq.translate(table=11)).rstrip("*")
            with open(TOX2_PROT_OUT, "w") as f:
                f.write(f">tox_2\n{prot_str}\n")
            print(f"[{SCRIPT}] Extracted tox_2 protein ({len(prot_str)} aa) -> {TOX2_PROT_OUT}")
            return
    sys.exit("ERROR: tox_2 not found in tox.fas")


def build_diamond_db():
    if Path(DMND_DB).exists():
        print(f"[{SCRIPT}] DIAMOND db already exists: {DMND_DB}")
        return
    step("Building DIAMOND db ...")
    subprocess.run(
        ["diamond", "makedb", "--in", TOX2_PROT_OUT, "--db", DMND_DB, "--quiet"],
        check=True, capture_output=True
    )
    done(DMND_DB)


ORFIPY_OUTDIR = Path("orfipy_batches")

def _process_batch(batch_id, genome_fastas):
    """Concatenate a batch of genomes to a temp file, run orfipy on it, clean up."""
    ORFIPY_OUTDIR.mkdir(exist_ok=True)
    batch_fa  = ORFIPY_OUTDIR / f"tmp_batch_{batch_id}.fa"
    pep_name  = f"tmp_batch_{batch_id}_pep.fa"
    bed_name  = f"tmp_batch_{batch_id}.bed"

    # Write concatenated batch to disk (avoids /dev/stdin issues on HPC)
    with open(batch_fa, "wb") as out:
        for fa in genome_fastas:
            opener = gzip.open if str(fa).endswith(".gz") else open
            with opener(fa, "rb") as fh:
                out.write(fh.read())

    subprocess.run(
        ["orfipy", str(batch_fa),
         "--pep", pep_name,
         "--bed", bed_name,
         "--table", str(ORFIPY_TABLE),
         "--between-stops",
         "--min", str(MIN_ORF_NT),
         "--strand", "b",
         "--procs", "1",
         "--chunk-size", "500"],
        check=True, stderr=subprocess.DEVNULL, cwd=str(ORFIPY_OUTDIR)
    )
    batch_fa.unlink()
    return ORFIPY_OUTDIR / pep_name, ORFIPY_OUTDIR / bed_name


def run_orfipy_batched():
    if Path(ORF_PEP).exists() and Path(ORF_BED).exists():
        n = sum(1 for line in open(ORF_PEP) if line.startswith(">"))
        print(f"[{SCRIPT}] ORF files already exist ({n:,} ORFs), skipping orfipy")
        return

    genome_fastas = sorted(
        list(Path(GENOMES_DIR).rglob("*.fa")) +
        list(Path(GENOMES_DIR).rglob("*.fa.gz"))
    )
    batches = [genome_fastas[i:i + BATCH_SIZE]
               for i in range(0, len(genome_fastas), BATCH_SIZE)]
    step(f"6FT-translating {len(genome_fastas)} genomes — "
         f"{len(batches)} batches of {BATCH_SIZE}, {ORFIPY_WORKERS} parallel, "
         f"min {MIN_ORF_NT // 3} aa ...")

    # Skip batches whose output files already exist (resumable)
    pending = [(i, batch) for i, batch in enumerate(batches)
               if not (ORFIPY_OUTDIR / f"tmp_batch_{i}_pep.fa").exists()]
    if len(pending) < len(batches):
        print(f"[{SCRIPT}]   {len(batches)-len(pending)} batches already done, "
              f"running {len(pending)} remaining ...", flush=True)

    pep_parts, bed_parts = [], []
    completed = len(batches) - len(pending)
    with ThreadPoolExecutor(max_workers=ORFIPY_WORKERS) as ex:
        futures = {ex.submit(_process_batch, i, batch): i
                   for i, batch in pending}
        for fut in as_completed(futures):
            pep, bed = fut.result()
            pep_parts.append(pep)
            bed_parts.append(bed)
            completed += 1
            if completed % ORFIPY_WORKERS == 0 or completed == len(batches):
                print(f"[{SCRIPT}]   {completed}/{len(batches)} batches done ...", flush=True)

    # Include any already-completed batches
    for i in range(len(batches)):
        p = ORFIPY_OUTDIR / f"tmp_batch_{i}_pep.fa"
        b = ORFIPY_OUTDIR / f"tmp_batch_{i}.bed"
        if p not in pep_parts and p.exists():
            pep_parts.append(p)
            bed_parts.append(b)

    # Merge in order so output is deterministic
    pep_parts.sort(key=lambda p: int(p.stem.split("_")[2]))
    bed_parts.sort(key=lambda p: int(p.stem.split("_")[2]))

    with open(ORF_PEP, "wb") as fout:
        for p in pep_parts:
            fout.write(p.read_bytes()); p.unlink()
    with open(ORF_BED, "wb") as fout:
        for b in bed_parts:
            fout.write(b.read_bytes()); b.unlink()

    n = sum(1 for line in open(ORF_PEP) if line.startswith(">"))
    done(f"{n:,} ORFs -> {ORF_PEP}")


def run_blastp(threads=8):
    step("Running DIAMOND blastp (ORFs vs tox_2) ...")
    subprocess.run(
        ["diamond", "blastp",
         "--db", DMND_DB,
         "--query", ORF_PEP,
         "--out", BLAST_RAW,
         "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
         "qstart", "qend", "sstart", "send", "slen", "qcovhsp", "evalue",
         "--max-target-seqs", "1",
         "--threads", str(threads),
         "-b", "2", "-c", "1", "--masking", "0", "--quiet"],
        check=True, capture_output=True
    )
    n = sum(1 for _ in open(BLAST_RAW))
    done(f"{n} raw hits")


def convert_coords():
    """Convert blastp ORF hits + BED coords to blastx-equivalent genome nt coords."""
    step("Converting ORF coordinates to genome nt positions ...")
    bed = {}
    with open(ORF_BED) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 6:
                continue
            contig, orf_start, orf_end, raw_name, _, strand = parts[:6]
            orf_name = raw_name.split(";")[0].removeprefix("ID=")
            bed[orf_name] = (contig, int(orf_start), int(orf_end), strand)

    written = 0
    with open(BLAST_RAW) as fin, open(BLAST_OUT, "w") as fout:
        for line in fin:
            parts = line.rstrip().split("\t")
            orf_name, sseqid, pident, length, qa, qb, sstart, send, slen, qcovhsp, evalue = parts
            qa, qb = int(qa), int(qb)
            if orf_name not in bed:
                continue
            contig, orf_start, orf_end, strand = bed[orf_name]
            if strand == "+":
                gstart = orf_start + (qa - 1) * 3 + 1
                gend   = orf_start + qb * 3
            else:
                gstart = orf_end - (qa - 1) * 3
                gend   = orf_end - qb * 3 + 1
            fout.write("\t".join([contig, sseqid, pident, length,
                                   str(gstart), str(gend),
                                   sstart, send, slen, qcovhsp, evalue]) + "\n")
            written += 1
    done(f"{written} hits -> {BLAST_OUT}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=8,
                        help="Total threads available (orfipy workers = threads//4, DIAMOND = threads)")
    args = parser.parse_args()

    global ORFIPY_WORKERS
    ORFIPY_WORKERS = max(1, args.threads // 4)

    print(f"[{SCRIPT}] Starting (threads={args.threads}, orfipy_workers={ORFIPY_WORKERS})", flush=True)
    download_tox_fas()
    extract_tox2_protein()
    build_diamond_db()
    run_orfipy_batched()
    run_blastp(args.threads)
    convert_coords()
    for tmp in (ORF_PEP, ORF_BED, BLAST_RAW):
        Path(tmp).unlink(missing_ok=True)
    print(f"[{SCRIPT}] Finished in {time.time() - _script_start:.1f}s")


if __name__ == "__main__":
    main()
