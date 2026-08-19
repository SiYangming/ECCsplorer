#!/usr/bin/env python3
"""
clu_prepare.py - Prepare reads for RepeatExplorer2 (seqclust) clustering.

Extracted from ECCsplorer/lib/eccPrepare.py:eccSampleConcat.prepare_for_clustering.

Input:
    --r1          : R1 reads (fasta/fastq, optionally gzipped), dataset A (treatment)
    --r2          : R2 reads (fasta/fastq, optionally gzipped), dataset A
    --c1          : (optional) R1 reads, dataset B (control)
    --c2          : (optional) R2 reads, dataset B
    --pre_a       : readID prefix dataset A (default: TR)
    --pre_b       : readID prefix dataset B (default: CO)
Output:
    --out         : REPEATEXPLORER_READY.fa (concatenated interlaced reads)
Options:
    --read-count  : number of paired-end reads to subsample (default: max available)
    --seed        : random seed (default: 42)
"""

import argparse
import math
import os
import random
import sys

import numpy as np
from Bio import SeqIO


def load_records(path, fmt_hint=None):
    """Yield Biopython records from fasta/fastq (gzip supported)."""
    if path.endswith('.gz'):
        import gzip
        handle = gzip.open(path, 'rt')
    else:
        handle = open(path, 'r')
    # detect format from extension
    if fmt_hint is None:
        lower = path.lower()
        if lower.endswith('.fq') or lower.endswith('.fastq') or lower.endswith('.fq.gz') or lower.endswith('.fastq.gz'):
            fmt = 'fastq'
        else:
            fmt = 'fasta'
    else:
        fmt = fmt_hint
    return SeqIO.parse(handle, fmt)


def get_best_read_length(in_file):
    """Optimal read length minimizing cumulative base loss (integer search)."""
    lengths = [len(rec.seq) for rec in load_records(in_file)]
    if not lengths:
        raise ValueError(f"No reads found in {in_file}")
    arr = np.array(lengths, dtype=float)

    def cumulative_bases(set_read_length):
        lost = 0.0
        for l in arr:
            if l < set_read_length:
                lost += l
            elif l > set_read_length:
                lost += (l - set_read_length)
        return lost

    # integer argmin over plausible read lengths
    lo, hi = 1, int(arr.max()) + 2
    best = min(range(lo, hi), key=cumulative_bases)
    return int(best)


def get_max_read_count(r1, r2):
    """Number of read pairs with both mates >= best_read_length."""
    return sum(1 for a, b in zip(load_records(r1), load_records(r2)))


def prexing_reads(r1, r2, prefix, best_len, read_count_set, seed):
    """Subsample, truncate, prefix read IDs, write interlaced FASTA."""
    random.seed(seed)
    np.random.seed(seed)
    rng = np.random.default_rng(seed)

    total = sum(1 for _ in load_records(r1))
    n_use = min(read_count_set, total)
    pick = set(rng.choice(total, size=n_use, replace=False).tolist())

    tmp = f"REPEATEXPLORER_{prefix}.fa.tmp"
    i = 0
    with open(tmp, 'w') as out:
        for rec1, rec2 in zip(load_records(r1), load_records(r2)):
            if len(rec1.seq) >= best_len and len(rec2.seq) >= best_len:
                if i in pick:
                    rec1.seq = rec1.seq[0:best_len]
                    rec1.id = f"{prefix}_{rec1.id}_#0/1"
                    rec1.description = ''
                    rec2.seq = rec2.seq[0:best_len]
                    rec2.id = f"{prefix}_{rec2.id}_#0/2"
                    rec2.description = ''
                    SeqIO.write(rec1, out, 'fasta')
                    SeqIO.write(rec2, out, 'fasta')
                i += 1
    return tmp


def main():
    parser = argparse.ArgumentParser(description='Prepare reads for RepeatExplorer2 clustering')
    parser.add_argument('--r1', required=True, help='R1 reads (dataset A / treatment)')
    parser.add_argument('--r2', required=True, help='R2 reads (dataset A / treatment)')
    parser.add_argument('--c1', default=None, help='R1 reads (dataset B / control, optional)')
    parser.add_argument('--c2', default=None, help='R2 reads (dataset B / control, optional)')
    parser.add_argument('--pre_a', default='TR', help='ReadID prefix dataset A (default: TR)')
    parser.add_argument('--pre_b', default='CO', help='ReadID prefix dataset B (default: CO)')
    parser.add_argument('--out', required=True, help='Output REPEATEXPLORER_READY.fa')
    parser.add_argument('--read-count', type=str, default=None,
                        help='Subsample read count (int or "auto"=0.1x genome coverage; default: max available)')
    parser.add_argument('--genome-size', type=int, default=None, help='Genome size in bp (for --read-count auto)')
    parser.add_argument('--seed', type=int, default=12, help='Random seed (default: 12, matches original PIPELINE_SEED)')
    args = parser.parse_args()

    read_sets = [[args.r1, args.r2, args.pre_a]]
    if args.c1 and args.c2:
        read_sets.append([args.c1, args.c2, args.pre_b])

    # 1) optimal read length (minimum across all read files)
    all_lens = []
    for rs in read_sets:
        all_lens.append(get_best_read_length(rs[0]))
        all_lens.append(get_best_read_length(rs[1]))
    best_len = min(all_lens)
    print(f'Optimal usable read length: {best_len}bp', file=sys.stderr)

    # 2) max available read pairs (minimum across sets)
    max_counts = [get_max_read_count(rs[0], rs[1]) for rs in read_sets]
    max_use = min(max_counts)
    if max_use == 0:
        sys.exit('Read count is 0. Cannot prepare reads for clustering.')
    if args.read_count == 'auto':
        if not args.genome_size:
            sys.exit('--read-count auto requires --genome-size')
        n_use = int(math.floor(args.genome_size * 0.1 / best_len))  # REPEX_FOLD_COV=0.1
    elif args.read_count:
        n_use = int(args.read_count)
    else:
        n_use = max_use
    n_use = min(n_use, max_use)
    print(f'Using {n_use} paired-end reads @ {best_len}bp', file=sys.stderr)

    # 3) prexing (interlaced) + concatenate
    temps = [prexing_reads(rs[0], rs[1], rs[2], best_len, n_use, args.seed) for rs in read_sets]
    with open(args.out, 'w') as out:
        for t in temps:
            with open(t, 'r') as f:
                for line in f:
                    out.write(line)
            os.remove(t)
    print(f'Wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
