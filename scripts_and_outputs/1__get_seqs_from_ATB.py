#!/usr/bin/env python3
"""
Download the latest ATB file list from OSF, filter to Corynebacterium genomes,
download the tar.xz archives, and unpack them.
"""
import csv
import gzip
import os
import sys
import tarfile
import time
import urllib.parse
import urllib.request
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

FILE_LIST_URL = "https://osf.io/download/3xs6h/"
FILE_LIST_OUT = "ATB_file_list.all.latest.tsv.gz"
OUTDIR        = "ATB_genomes"


def download_file_list():
    if Path(FILE_LIST_OUT).exists():
        print(f"[{SCRIPT}] {FILE_LIST_OUT} already exists, skipping download")
        return
    step("Fetching ATB file list ...")
    final_url = urllib.request.urlopen(FILE_LIST_URL).url
    urllib.request.urlretrieve(final_url, FILE_LIST_OUT)
    done(f"saved to {FILE_LIST_OUT}")


def get_coryne_batches():
    batches = {}
    with gzip.open(FILE_LIST_OUT, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "Corynebacterium" not in row.get("sylph_species", ""):
                continue
            url = row.get("tar_xz_url", "").strip()
            fname = row.get("filename_in_tar_xz", "").strip()
            if url and fname:
                batches.setdefault(url, []).append(Path(OUTDIR) / fname)
    print(f"[{SCRIPT}] Found {len(batches)} unique tar.xz batches for Corynebacterium")
    return sorted(batches.items())


def get_filename_from_url(url):
    parsed = urllib.parse.urlparse(url)
    params = urllib.parse.parse_qs(parsed.query)
    disposition = params.get("response-content-disposition", [""])[0]
    for part in disposition.split(";"):
        part = part.strip()
        if part.startswith("filename="):
            return part.split("=", 1)[1].strip('"')
    return parsed.path.split("/")[-1]


def download_and_unpack(url, outdir, wanted):
    """Download tar.xz and extract only the files listed in `wanted` (set of tar member names)."""
    final_url = urllib.request.urlopen(url).url  # follow OSF redirect
    filename = get_filename_from_url(final_url)
    dest = Path(outdir) / filename
    if not dest.exists():
        urllib.request.urlretrieve(final_url, dest)
    with tarfile.open(dest, "r:xz") as tar:
        for member in tar.getmembers():
            if member.name in wanted:
                tar.extract(member, path=outdir, filter="data")
    dest.unlink()  # remove archive after unpacking


def main():
    print(f"[{SCRIPT}] Starting", flush=True)
    Path(OUTDIR).mkdir(exist_ok=True)

    download_file_list()
    batches = get_coryne_batches()

    n_skipped = 0
    for i, (url, expected_paths) in enumerate(batches, 1):
        if all(p.exists() for p in expected_paths):
            n_skipped += 1
            continue
        step(f"[{i}/{len(batches)}] Downloading {Path(url).name} ...")
        # Tar member names match filename_in_tar_xz (e.g. "batch.41/SAMN00000.fa")
        wanted = {str(p.relative_to(OUTDIR)) for p in expected_paths}
        try:
            download_and_unpack(url, OUTDIR, wanted)
            done()
        except Exception as e:
            print(f"ERROR: {e}", file=sys.stderr)

    if n_skipped:
        print(f"[{SCRIPT}] Skipped {n_skipped}/{len(batches)} batches (already unpacked)")
    total = time.time() - _script_start
    print(f"[{SCRIPT}] Finished in {total:.1f}s")


if __name__ == "__main__":
    main()
