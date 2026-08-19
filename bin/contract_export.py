#!/usr/bin/env python3
"""
contract_export.py — Rewrite slim outputs into the original ECCsplorer
eccpipe_results/ directory tree and file naming (output contract compatibility).

Input args (optional, may be omitted):
  --hiconf, --lowconf, --candidates, --normalized, --sequences,
  --blast_m6, --clu_candidates, --clu_table, --clu_list,
  --map_html, --clu_html, --comp_html, --prefix
Output: ./eccpipe_results/{mapping_results,clustering_results,comparative_results,...}
"""
import argparse
import os
import shutil


def ensure(path):
    os.makedirs(path, exist_ok=True)
    return path


def copy_if(src, dst):
    if src and os.path.isfile(src) and os.path.getsize(src) > 0:
        shutil.copy(src, dst)


def main():
    parser = argparse.ArgumentParser(description='Export slim outputs to eccpipe_results tree')
    parser.add_argument('--hiconf', default=None)
    parser.add_argument('--lowconf', default=None)
    parser.add_argument('--candidates', default=None)
    parser.add_argument('--normalized', default=None)
    parser.add_argument('--sequences', default=None)
    parser.add_argument('--blast_m6', default=None)
    parser.add_argument('--clu_candidates', default=None)
    parser.add_argument('--clu_table', default=None)
    parser.add_argument('--clu_list', default=None)
    parser.add_argument('--map_html', default=None)
    parser.add_argument('--clu_html', default=None)
    parser.add_argument('--comp_html', default=None)
    parser.add_argument('--prefix', default='TR')
    args = parser.parse_args()

    root = 'eccpipe_results'
    mapping = ensure(os.path.join(root, 'mapping_results'))
    clustering = ensure(os.path.join(root, 'clustering_results'))
    comparative = ensure(os.path.join(root, 'comparative_results'))
    cand_dir = ensure(os.path.join(mapping, 'candidates'))

    # mapping_results
    copy_if(args.hiconf, os.path.join(mapping, f"{args.prefix}_hiconf-ECC-REGIONS.bed"))
    copy_if(args.lowconf, os.path.join(mapping, f"{args.prefix}_lowconf-ECC-regions.bed"))
    copy_if(args.candidates, os.path.join(mapping, f"{args.prefix}_candidates.bed"))
    if args.normalized:
        copy_if(args.normalized, os.path.join(mapping, f"{args.prefix}_summary_region-coverages_normalized.csv"))
    if args.sequences:
        copy_if(args.sequences, os.path.join(mapping, f"{args.prefix}_ECC-SEQUENCES.fasta"))
    if args.blast_m6:
        copy_if(args.blast_m6, os.path.join(cand_dir, f"{args.prefix}_blast.m6"))
    if args.map_html:
        copy_if(args.map_html, os.path.join(mapping, 'eccMap_summary.html'))

    # clustering_results
    copy_if(args.clu_candidates, os.path.join(clustering, 'COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv'))
    copy_if(args.clu_table, os.path.join(clustering, 'comparative_cluster_table.csv'))
    copy_if(args.clu_list, os.path.join(clustering, 'comp_cl_tab_eccCandidates_list.csv'))
    if args.clu_html:
        copy_if(args.clu_html, os.path.join(clustering, 'eccCL_summary.html'))

    # comparative_results
    if args.comp_html:
        copy_if(args.comp_html, os.path.join(comparative, 'eccComp_summary.html'))

    print(f'contract_export: wrote eccpipe_results tree under {os.path.abspath(root)}')


if __name__ == '__main__':
    main()
