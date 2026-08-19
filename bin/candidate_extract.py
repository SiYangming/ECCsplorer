#!/usr/bin/env python3
"""
candidate_extract.py — Extract eccDNA candidate regions by intersecting
split-read regions with peak regions from all/DiscordantRead mappings.

Extracted from ECCsplorer/lib/eccMapper.py:eccMapper.extract_candidate_regions.

Input:
    --sr       : Split-read regions BED (from haarz)
    --peak_all : Peak regions BED (all reads)
    --peak_dr  : Peak regions BED (discordant reads)
Output:
    --output   : Candidate regions BED (high-conf preferred, fallback to low-conf)

Parameters:
    --max_len  : Maximum eccDNA length (default: 35000)
    --min_len  : Minimum eccDNA length (default: 100)
    --merge    : Merge close regions distance (default: 1000)
"""

import argparse
import sys
import os
import subprocess
import shlex
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract eccDNA candidate regions via bedtools intersect"
    )
    parser.add_argument("--sr", required=True, help="Split-read regions BED file")
    parser.add_argument("--peak_all", required=True, help="Peak regions (all reads) BED file")
    parser.add_argument("--peak_dr", required=True, help="Peak regions (discordant reads) BED file")
    parser.add_argument("--output", required=True, help="Output candidate regions BED file")
    parser.add_argument("--hiconf_out", required=True, help="High-confidence (3/3) candidate regions BED")
    parser.add_argument("--lowconf_out", required=True, help="Low-confidence (2/3) candidate regions BED")
    parser.add_argument("--max_len", type=int, default=35000, help="Maximum eccDNA length")
    parser.add_argument("--min_len", type=int, default=100, help="Minimum eccDNA length")
    parser.add_argument("--merge", type=int, default=1000, help="Merge distance for regions")
    parser.add_argument("--bedtools", default="bedtools", help="Path to bedtools executable")
    return parser.parse_args()


def run_bedtools(cmd, error_msg="bedtools command failed"):
    """Run a bedtools command (may contain shell pipes) and return True if it succeeded."""
    proc = subprocess.run(
        cmd, shell=True, capture_output=True, text=True
    )
    if proc.returncode != 0 and proc.stderr.strip():
        print(f"WARNING: {error_msg}: {proc.stderr.strip()}", file=sys.stderr)
    return proc.returncode == 0


def check_file(filepath, label):
    """Ensure a file exists and is non-empty."""
    if not os.path.isfile(filepath):
        print(f"ERROR: {label} file '{filepath}' does not exist.", file=sys.stderr)
        sys.exit(1)


def filter_by_length(filepath, prefix, max_len, min_len):
    """Filter BED file by region length and return path to filtered file."""
    filtered = f"{prefix}_lenfilt.bed"
    cmd = (
        f"awk '($3-$2)<={max_len} && ($3-$2)>={min_len}' {filepath} > {filtered}"
    )
    subprocess.run(cmd, shell=True, capture_output=True)
    return filtered if os.path.isfile(filtered) else filepath


def main():
    args = parse_args()

    check_file(args.sr, "SR regions")
    check_file(args.peak_all, "Peak all")
    check_file(args.peak_dr, "Peak DR")

    prefix = os.path.splitext(os.path.basename(args.output))[0]
    tmpdir = tempfile.mkdtemp(prefix="eccsplorer_candidates_")

    # Sanitize BED files (ensure start < end)
    for fpath, label in [
        (args.sr, "SR"),
        (args.peak_all, "peak_all"),
        (args.peak_dr, "peak_dr"),
    ]:
        sanitized = os.path.join(tmpdir, f"{label}_clean.bed")
        subprocess.run(
            f"awk '$2<$3' {fpath} > {sanitized}",
            shell=True, capture_output=True
        )
        if os.path.isfile(sanitized) and os.path.getsize(sanitized) > 0:
            setattr(args, f"_{label}_clean", sanitized)
        else:
            setattr(args, f"_{label}_clean", fpath)

    # ---- High-confidence: 3/3 (SR ∩ peak_all ∩ peak_dr) ----
    hiconf_file = os.path.join(tmpdir, f"{prefix}_hiconf.bed")
    cmd_hiconf = (
        f"bedtools intersect -u -a {args._SR_clean} -b {args._peak_all_clean} | "
        f"bedtools intersect -u -a stdin -b {args._peak_dr_clean} > {hiconf_file}"
    )
    run_bedtools(cmd_hiconf, "High-confidence intersect failed")

    # ---- Low-confidence: 2/3 combos ----
    # SR ∩ all, not DR
    sr_all_file = os.path.join(tmpdir, f"{prefix}_lowconf_SR-all.bed")
    cmd_sr_all = (
        f"bedtools intersect -u -a {args._SR_clean} -b {args._peak_all_clean} | "
        f"bedtools intersect -v -a stdin -b {args._peak_dr_clean} > {sr_all_file}"
    )
    run_bedtools(cmd_sr_all, "SR-all intersect failed")

    # SR ∩ DR, not all
    sr_dr_file = os.path.join(tmpdir, f"{prefix}_lowconf_SR-DR.bed")
    cmd_sr_dr = (
        f"bedtools intersect -u -a {args._SR_clean} -b {args._peak_dr_clean} | "
        f"bedtools intersect -v -a stdin -b {args._peak_all_clean} > {sr_dr_file}"
    )
    run_bedtools(cmd_sr_dr, "SR-DR intersect failed")

    # all ∩ DR, not SR
    dr_all_file = os.path.join(tmpdir, f"{prefix}_lowconf_DR-all.bed")
    cmd_dr_all = (
        f"bedtools intersect -u -a {args._peak_all_clean} -b {args._peak_dr_clean} | "
        f"bedtools intersect -v -a stdin -b {args._SR_clean} > {dr_all_file}"
    )
    run_bedtools(cmd_dr_all, "DR-all intersect failed")

    # Concatenate low-conf regions
    lowconf_file = os.path.join(tmpdir, f"{prefix}_lowconf.bed")
    subprocess.run(
        f"cat {sr_all_file} {sr_dr_file} {dr_all_file} | sort -k1,1 -k2,2n > {lowconf_file}",
        shell=True, capture_output=True
    )

    # Decide which to use: high-conf if non-empty, else low-conf
    hiconf_ok = os.path.isfile(hiconf_file) and os.path.getsize(hiconf_file) > 0
    lowconf_ok = os.path.isfile(lowconf_file) and os.path.getsize(lowconf_file) > 0

    # Filter by length and write BOTH outputs (replicates eccMapper extract_candidate_regions
    # which emits {pre}_hiconf-ECC-REGIONS.bed and {pre}_lowconf-ECC-regions.bed)
    if hiconf_ok:
        final = filter_by_length(hiconf_file, prefix, args.max_len, args.min_len)
        subprocess.run(f"cp {final} {args.hiconf_out}", shell=True, capture_output=True)
        print("Using high-confidence (3/3) candidate regions.", file=sys.stderr)
    else:
        open(args.hiconf_out, 'w').close()

    if lowconf_ok:
        final = filter_by_length(lowconf_file, prefix, args.max_len, args.min_len)
        subprocess.run(f"cp {final} {args.lowconf_out}", shell=True, capture_output=True)
    else:
        open(args.lowconf_out, 'w').close()

    if hiconf_ok:
        candidate_source = hiconf_file
    elif lowconf_ok:
        candidate_source = lowconf_file
    else:
        print("WARNING: No candidate regions found.", file=sys.stderr)
        open(args.output, 'w').close()
        sys.exit(0)

    # Copy to output (unified candidates, hiconf preferred)
    subprocess.run(f"cp {candidate_source} {args.output}", shell=True, capture_output=True)

    # Cleanup temp directory
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
