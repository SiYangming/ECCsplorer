#!/usr/bin/env python3
"""
peak_detect.py — Detect peaks in windowed coverage data using scipy.signal.find_peaks.

Extracted from ECCsplorer/lib/eccMapper.py:get_rough_coverage (peak calculation logic).

Input:
    --coverage  : TSV file with columns [chr, start, end, coverage_value]
                  (output from bedtools coverage -mean over genome windows)
Output:
    --output    : BED file with peak region boundaries [chr, peak_start, peak_end]
                  (using scipy peak_widths at rel_height=0.9)

Parameters:
    --threshold : Minimum peak height (default: None = auto)
    --distance  : Minimum distance between peaks in windows (default: 20)
"""

import argparse
import sys
import os
import numpy as np
from scipy.signal import find_peaks, peak_widths


def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect peaks in windowed coverage data"
    )
    parser.add_argument(
        "--coverage", required=True,
        help="Windowed coverage TSV (chr, start, end, coverage_value)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output BED file with peak regions"
    )
    parser.add_argument(
        "--threshold", type=float, default=None,
        help="Minimum peak height (default: None = auto from data)"
    )
    parser.add_argument(
        "--distance", type=int, default=20,
        help="Minimum distance between peaks in windows (default: 20)"
    )
    parser.add_argument(
        "--prominence", type=float, default=1.0,
        help="Minimum peak prominence (default: 1.0)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Read coverage data
    try:
        cov_raw = np.loadtxt(args.coverage, dtype=str, delimiter='\t', ndmin=2)
    except Exception as e:
        print(f"ERROR: Failed to read coverage file '{args.coverage}': {e}", file=sys.stderr)
        sys.exit(1)

    if cov_raw.shape[0] == 0:
        print("WARNING: Empty coverage file, no peaks to detect.", file=sys.stderr)
        # Write empty output
        open(args.output, 'w').close()
        sys.exit(0)

    # Coverage values are in column 3 (0-indexed)
    coverage_values = cov_raw[:, 3].astype(float)

    # Find peaks
    peaks, properties = find_peaks(
        coverage_values,
        threshold=args.threshold,
        distance=args.distance,
        prominence=args.prominence
    )

    if len(peaks) == 0:
        print("WARNING: No peaks detected.", file=sys.stderr)
        open(args.output, 'w').close()
        sys.exit(0)

    # Calculate peak widths (region boundaries)
    widths_results = peak_widths(coverage_values, peaks, rel_height=0.9)

    # Build regions: left boundary from width[2], right boundary from width[3]
    regions = []
    for i in range(len(peaks)):
        left_idx = int(widths_results[2][i])
        right_idx = int(widths_results[3][i])
        if left_idx < 0:
            left_idx = 0
        if right_idx >= cov_raw.shape[0]:
            right_idx = cov_raw.shape[0] - 1

        chrom = cov_raw[left_idx, 0]
        start = cov_raw[left_idx, 1]
        end = cov_raw[right_idx, 2]
        regions.append([chrom, start, end])

    # Save peak regions
    np.savetxt(args.output, np.array(regions), delimiter='\t', fmt='%s')


if __name__ == "__main__":
    main()
