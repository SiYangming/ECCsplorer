#!/usr/bin/env python3
"""
html_report.py — Generate HTML summary report for ECCsplorer eccDNA candidates.

Extracted from ECCsplorer/lib/eccHTML_templates.py (HTML templates)
and ECCsplorer/lib/eccMapper.py:eccMapper.html_summarize.

Input:
    --candidates  : TSV file with candidate data
                    Columns: chr, start, end, TR.all, TR.SR, TR.DR, CO.all, CO.SR, CO.DR,
                             enrich.all, id, sequence (optional)
    --blast       : Directory containing per-candidate BLAST result files (*_blast.m6)
    --outdir      : Output directory (also scanned for per-candidate plot files)
Output:
    --output      : HTML report file

Parameters:
    --prefix      : Prefix for data set labels (default: "sample")
    --mapped_reads: Comma-separated mapped read counts
    --mapped_bases: Comma-separated mapped base counts
"""

import argparse
import sys
import os
import re


# HTML template fragments (from eccHTML_templates.py)
HTML_STYLE = """<style>
body { font-family: Arial, Helvetica, sans-serif; margin: 20px; }
h2 { color: #333; border-bottom: 2px solid #333; padding-bottom: 5px; }
h3 { color: #555; }
table.dataframe { border-collapse: collapse; width: 95%; margin: 10px 0; }
table.dataframe td, table.dataframe th { border: 1px solid #ddd; padding: 6px; text-align: left; }
table.dataframe tr.firstline { background-color: #f2f2f2; font-weight: bold; }
table.dataframe tr:nth-child(even) { background-color: #f9f9f9; }
.seqcell { font-family: monospace; font-size: 10px; max-width: 300px; }
.status-ok { color: green; font-weight: bold; }
.status-warn { color: orange; font-weight: bold; }
.status-err { color: red; font-weight: bold; }
</style>"""

HTML_SUMMARY_HEADER = """<h2>ECCsplorer — eccMap Summary</h2>
<p>
<table border="0" class="dataframe">
<tr class="firstline">
    <th>Data Set</th>
    <th>Mapped Reads</th>
    <th>Mapped Bases</th>
</tr>
<tr><td>{pre1}</td><td>{reads1}</td><td>{bases1}</td></tr>
<tr><td>{pre2}</td><td>{reads2}</td><td>{bases2}</td></tr>
</table>
<br>
</p>"""

HTML_CANDIDATE_HEADER = """<h3>ECCsplorer — Candidate Overview</h3>
<p>
<table border="0" class="dataframe">
<tbody>
<tr class="firstline">
    <th>Candidate</th>
    <th>Mapping Plots</th>
    <th>Enrichment Score</th>
    <th>Length</th>
    <th>Position</th>
    <th>BLAST Best Hit</th>
    <th>Estimated Sequence</th>
</tr>"""

HTML_CANDIDATE_ROW = """<tr>
    <td>{cand}</td>
    <td><a href="{plot}"><img src="{plot}" alt="{cand}" height="300px"></a></td>
    <td>{enr}</td>
    <td>{length}</td>
    <td>{pos}</td>
    <td>{blast}</td>
    <td class="seqcell"><div style="height:300px; overflow:auto;">{seq}</div></td>
</tr>"""

HTML_CANDIDATE_ROW_NO_PLOT = """<tr>
    <td>{cand}</td>
    <td><span class="status-warn">No plot</span></td>
    <td>{enr}</td>
    <td>{length}</td>
    <td>{pos}</td>
    <td>{blast}</td>
    <td class="seqcell"><div style="height:300px; overflow:auto;">{seq}</div></td>
</tr>"""

HTML_FOOTER = """</tbody></table></p><br><br>"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate HTML summary report for ECCsplorer candidates"
    )
    parser.add_argument("--candidates", required=True,
                        help="TSV file with candidate region data")
    parser.add_argument("--blast", default=None,
                        help="Directory containing per-candidate BLAST result files")
    parser.add_argument("--outdir", default=None,
                        help="Directory with per-candidate plot files (default: same as candidates dir)")
    parser.add_argument("--output", required=True,
                        help="Output HTML file path")
    parser.add_argument("--prefix", default="sample",
                        help="Prefix for data set labels")
    parser.add_argument("--mapped_reads", default=None,
                        help="Comma-separated mapped read counts")
    parser.add_argument("--mapped_bases", default=None,
                        help="Comma-separated mapped base counts")
    parser.add_argument("--manhattan_plot", default=None,
                        help="Path to Manhattan plot image (relative to output)")
    return parser.parse_args()


def parse_blast_best_hit(blast_file):
    """Get best BLAST hit (highest bitscore) from BLAST m6 file."""
    if not blast_file or not os.path.isfile(blast_file):
        return ""
    best_score = -1
    best_hit = ""
    try:
        with open(blast_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                try:
                    score = float(parts[11])  # bitscore
                    if score > best_score:
                        best_score = score
                        best_hit = parts[1]  # sseqid
                except ValueError:
                    continue
    except Exception:
        pass
    if best_hit:
        return best_hit
    return ""


def parse_sequence(fasta_file):
    """Read sequence from fasta file."""
    if not fasta_file or not os.path.isfile(fasta_file):
        return ""
    try:
        with open(fasta_file) as f:
            seq = ""
            for line in f:
                if not line.startswith(">"):
                    seq += line.strip()
            return seq[:1000]  # Truncate for display
    except Exception:
        return ""


def main():
    args = parse_args()

    # Determine directories
    outdir = args.outdir or os.path.dirname(args.candidates)
    blast_dir = args.blast or outdir

    # Parse mapped stats
    mapped_reads = [0, 0]
    mapped_bases = [0, 0]
    if args.mapped_reads:
        parts = args.mapped_reads.split(',')
        for i, p in enumerate(parts[:2]):
            try:
                mapped_reads[i] = int(p.strip())
            except ValueError:
                pass
    if args.mapped_bases:
        parts = args.mapped_bases.split(',')
        for i, p in enumerate(parts[:2]):
            try:
                mapped_bases[i] = int(p.strip())
            except ValueError:
                pass

    # Prefixes
    pre1 = f"{args.prefix}_TR"
    pre2 = f"{args.prefix}_CO"

    # Read candidate data
    try:
        with open(args.candidates) as f:
            header = f.readline().strip().split('\t')
            candidates = []
            for line in f:
                parts = line.strip().split('\t')
                if parts:
                    candidates.append(dict(zip(header, parts)))
    except Exception as e:
        print(f"ERROR: Failed to read candidates file: {e}", file=sys.stderr)
        sys.exit(1)

    # Generate HTML
    html_parts = ["<!DOCTYPE html>\n<html>\n<head>"]
    html_parts.append(f"<title>ECCsplorer Report — {args.prefix}</title>")
    html_parts.append(HTML_STYLE)
    html_parts.append("</head>\n<body>")

    # Summary section
    manhattan_link = args.manhattan_plot or ""
    if manhattan_link and os.path.isfile(os.path.join(os.path.dirname(args.output), manhattan_link)):
        html_parts.append(
            f'<p><a href="{manhattan_link}">'
            f'<img src="{manhattan_link}" alt="Manhattan plot" height="450px"></a></p>'
        )
    elif manhattan_link:
        html_parts.append(
            f'<p><a href="{manhattan_link}">'
            f'<img src="{manhattan_link}" alt="Manhattan plot" height="450px"></a></p>'
        )

    html_parts.append(
        HTML_SUMMARY_HEADER.format(
            pre1=pre1, reads1=f"{mapped_reads[0]:,}", bases1=f"{mapped_bases[0]:,}",
            pre2=pre2, reads2=f"{mapped_reads[1]:,}", bases2=f"{mapped_bases[1]:,}"
        )
    )

    # Candidate table
    html_parts.append(HTML_CANDIDATE_HEADER)

    for cand in candidates:
        cand_id = cand.get('id', f"{cand.get('chr','?')}:{cand.get('start','?')}-{cand.get('end','?')}")
        enr = cand.get('enrich.all', 'N/A')
        length = cand.get('length', str(int(cand.get('end', 0)) - int(cand.get('start', 0))))
        pos = f"{cand.get('chr','?')}:{cand.get('start','?')}-{cand.get('end','?')}"
        seq = cand.get('sequence', '')[:1000]

        # Blast hit
        blast_path = os.path.join(blast_dir, f"{cand_id}_blast.m6")
        blast_hit = parse_blast_best_hit(blast_path)
        if blast_hit:
            blast_rel = os.path.relpath(blast_path, os.path.dirname(args.output))
            blast_display = f'<a href="{blast_rel}">{blast_hit}</a>'
        else:
            blast_display = "N/A"

        # Plot
        plot_extensions = ['.png', '.jpeg', '.jpg', '.bmp', '.tiff', '.pdf']
        plot_path = None
        for ext in plot_extensions:
            cand_dir = os.path.join(outdir, cand_id)
            test_path = os.path.join(cand_dir, f"{cand_id}_multiplot{ext}")
            if os.path.isfile(test_path):
                plot_path = test_path
                break
        if plot_path:
            plot_rel = os.path.relpath(plot_path, os.path.dirname(args.output))
            html_parts.append(
                HTML_CANDIDATE_ROW.format(
                    cand=cand_id, plot=plot_rel,
                    enr=enr, length=length, pos=pos,
                    blast=blast_display, seq=seq
                )
            )
        else:
            html_parts.append(
                HTML_CANDIDATE_ROW_NO_PLOT.format(
                    cand=cand_id,
                    enr=enr, length=length, pos=pos,
                    blast=blast_display, seq=seq
                )
            )

    html_parts.append(HTML_FOOTER)
    html_parts.append("</body>\n</html>")

    # Write output
    with open(args.output, 'w') as f:
        f.write('\n'.join(html_parts))

    print(f"HTML report generated: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
