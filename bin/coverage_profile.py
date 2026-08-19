#!/usr/bin/env python3
"""
coverage_profile.py — Extract per-base coverage for each eccDNA candidate region.

Extracted from ECCsplorer/lib/eccMapper.py:analyze_candidate_region.

For each candidate region in the candidates BED file:
    1. Extract sequence with bedtools getfasta (--ref)
    2. Calculate per-base coverage with bedtools coverage -d against all alignments
    3. Merge coverage columns from multiple BAMs into a single summary TSV

Input:
    --candidates  : BED file of candidate regions
    --alignments  : Comma-separated list of alignment BAM/BED files
    --ref         : Reference genome FASTA file
Output:
    --outdir      : Output directory for per-candidate coverage TSV files
                    (each file: <candidate_id>_coverage.tsv)

Parameters:
    --names       : Comma-separated list of sample names for alignment columns
    --bedtools    : Path to bedtools executable
"""

import argparse
import sys
import os
import subprocess
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract per-base coverage profiles for eccDNA candidates"
    )
    parser.add_argument(
        "--candidates", required=True,
        help="BED file of candidate regions"
    )
    parser.add_argument(
        "--alignments", required=True,
        help="Comma-separated list of alignment BAM/BED files"
    )
    parser.add_argument(
        "--ref", required=True,
        help="Reference genome FASTA file"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for per-candidate coverage TSV files"
    )
    parser.add_argument(
        "--names", default=None,
        help="Comma-separated list of sample names (default: use filenames)"
    )
    parser.add_argument(
        "--bedtools", default="bedtools",
        help="Path to bedtools executable"
    )
    return parser.parse_args()


def run_coverage_d(region_str, alignment, bedtools="bedtools", max_retries=5):
    """Run bedtools coverage -d for a region against an alignment file."""
    import time
    for attempt in range(max_retries):
        try:
            result = subprocess.check_output(
                f"echo '{region_str}' | {bedtools} coverage -d -a stdin -b {alignment}",
                shell=True, universal_newlines=True, stderr=subprocess.DEVNULL
            )
            return result
        except subprocess.CalledProcessError:
            if attempt < max_retries - 1:
                time.sleep(5)
            else:
                print(
                    f"WARNING: bedtools coverage failed for region '{region_str}' "
                    f"against '{alignment}' after {max_retries} attempts.",
                    file=sys.stderr
                )
                return "chr\t0\t0\t0\t0\n"
    return "chr\t0\t0\t0\t0\n"


def main():
    args = parse_args()

    # Parse alignment files and names
    alignment_files = [f.strip() for f in args.alignments.split(',') if f.strip()]
    if args.names:
        names = [n.strip() for n in args.names.split(',') if n.strip()]
        if len(names) != len(alignment_files):
            print(
                f"ERROR: {len(names)} names but {len(alignment_files)} alignment files.",
                file=sys.stderr
            )
            sys.exit(1)
    else:
        names = [
            os.path.splitext(os.path.basename(f))[0]
            for f in alignment_files
        ]

    # Read candidates
    with open(args.candidates) as f:
        candidates = [line.strip().split('\t') for line in f if line.strip()]

    if not candidates:
        print("WARNING: No candidates in BED file.", file=sys.stderr)
        sys.exit(0)

    os.makedirs(args.outdir, exist_ok=True)

    # Process each candidate
    for cand in candidates:
        chrom, start, end = cand[0], cand[1], cand[2]
        candidate_id = f"{chrom}:{start}-{end}"
        region_str = f"{chrom}\t{start}\t{end}"

        # Create candidate subdirectory
        cand_dir = os.path.join(args.outdir, candidate_id)
        os.makedirs(cand_dir, exist_ok=True)

        # Extract sequence with bedtools getfasta
        fasta_file = os.path.join(cand_dir, f"{candidate_id}.fasta")
        if not os.path.isfile(fasta_file) or os.path.getsize(fasta_file) == 0:
            with open(fasta_file, 'w') as out_fasta:
                try:
                    seq = subprocess.check_output(
                        f"{args.bedtools} getfasta -fi {args.ref} -bed /dev/stdin",
                        shell=True, universal_newlines=True,
                        input=f"{chrom}\t{start}\t{end}"
                    )
                    out_fasta.write(seq)
                except subprocess.CalledProcessError as e:
                    print(
                        f"WARNING: bedtools getfasta failed for {candidate_id}: {e}",
                        file=sys.stderr
                    )

        # Calculate per-base coverage for each alignment
        summary_path = os.path.join(cand_dir, f"{candidate_id}_coverage.tsv")
        if os.path.isfile(summary_path) and os.path.getsize(summary_path) > 0:
            print(f"SKIP: {candidate_id} coverage exists.", file=sys.stderr)
            continue

        all_columns = []
        for i, aln_file in enumerate(alignment_files):
            cov_text = run_coverage_d(region_str, aln_file, args.bedtools)
            # Parse coverage data
            lines = cov_text.strip().split('\n')
            parsed = [line.split('\t') for line in lines if line.strip()]
            if not parsed:
                parsed = [['chr', '0', '0', '0', '0']]

            cov_arr = np.array(parsed)

            if i == 0:
                # First alignment: keep chr, start, end, pos columns
                all_columns.append(cov_arr[:, :4])
            # Append coverage column (no per-column header; header built once below)
            all_columns.append(cov_arr[:, 4:5])

        # Build summary: chr, start, end, pos, [name1_val, name2_val, ...]
        header_row = np.array([['chr', 'start', 'end', 'pos'] + list(names[:len(all_columns)-1])])
        body = all_columns[0]
        for col in all_columns[1:]:
            body = np.hstack([body, col])
        final = np.vstack([header_row, body]) if len(body) > 1 else header_row

        # Save
        np.savetxt(summary_path, final, delimiter='\t', fmt='%s')
        print(f"OK: {candidate_id} ({len(final)-1} positions)", file=sys.stderr)


if __name__ == "__main__":
    main()
