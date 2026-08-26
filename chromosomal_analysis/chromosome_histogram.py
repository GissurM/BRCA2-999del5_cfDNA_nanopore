import os
import pysam
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path
import re

import pandas as pd

import argparse

parser = argparse.ArgumentParser(description="Extract chromosome thread counts from BAM files and plot histograms.")
parser.add_argument('--input_dir', default=None, help='Directory containing BAM files')
parser.add_argument('--out_dir', default='bam_chrom_histograms', help='Directory to save output histograms')
parser.add_argument(
    '--counts_dir',
    default=r'/mnt/d/BRCA2-misc_files/hg38_1-24',
    help='Directory containing chromosome_counts_table_*.csv files',
)
args = parser.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

combined_counts = Counter()

# Define standard chromosomes
standard_chroms_ordered = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM', 'chrMT']
standard_chroms = set(standard_chroms_ordered)

if args.input_dir:
    for file in os.listdir(args.input_dir):
        if file.endswith('.bam'):
            bam_path = os.path.join(args.input_dir, file)
            sample_name = os.path.splitext(file)[0]
            chrom_counts = Counter()
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                for read in bam.fetch(until_eof=True):
                    if not read.is_unmapped:
                        chrom = bam.get_reference_name(read.reference_id)
                        if chrom in standard_chroms:
                            chrom_counts[chrom] += 1
                            combined_counts[chrom] += 1
            # Plot for this BAM
            plt.figure(figsize=(10, 6))
            chroms = [c for c in standard_chroms_ordered if c in chrom_counts]
            counts = [chrom_counts[c] / 1e5 for c in chroms]  # Scale to hundred thousands
            plt.bar(chroms, counts)
            plt.xlabel('Chromosome')
            plt.ylabel('Number of threads (×10⁵)')
            plt.title(f'Threads per Chromosome: {sample_name}')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(args.out_dir, f'{sample_name}_chrom_hist.png'))
            plt.close()

# Combined histogram
if combined_counts:
    plt.figure(figsize=(10, 6))
    chroms = [c for c in standard_chroms_ordered if c in combined_counts]
    counts = [combined_counts[c] / 1e5 for c in chroms]  # Scale to hundred thousands
    plt.bar(chroms, counts)
    plt.xlabel('Chromosome')
    plt.ylabel('Number of threads (×10⁵)')
    plt.title('Threads per Chromosome: All samples Combined')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'combined_chrom_hist.png'))
    plt.close()
if args.input_dir:
    print(f"Histograms saved to {args.out_dir}")


def _read_percentage_tables(counts_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read chromosome count tables and split by BRCA2/ctrl groups."""
    files = sorted(counts_dir.glob('chromosome_counts_table_*.csv'))
    if not files:
        raise FileNotFoundError(f'No chromosome_counts_table_*.csv files found in {counts_dir}')

    rows: list[dict] = []
    for file_path in files:
        sample = file_path.stem.replace('chromosome_counts_table_', '')
        sample_lower = sample.lower()
        if sample_lower.startswith('brca2-'):
            group = 'BRCA2'
        elif sample_lower.startswith('ctrl-'):
            group = 'ctrl'
        else:
            continue

        df = pd.read_csv(file_path, sep=';', dtype=str)
        if not {'Chromosome', 'Percentage'}.issubset(df.columns):
            continue

        tmp = df[['Chromosome', 'Percentage']].copy()
        tmp['chrom'] = tmp['Chromosome'].astype(str).str.strip().str.upper()
        tmp['chrom_num'] = tmp['chrom'].str.extract(r'^CHR(\d+)$', expand=False)
        tmp = tmp[tmp['chrom_num'].notna()].copy()
        tmp['chrom_num'] = tmp['chrom_num'].astype(int)
        tmp = tmp[(tmp['chrom_num'] >= 1) & (tmp['chrom_num'] <= 22)]

        # Input uses decimal comma (e.g., 9,12)
        tmp['pct'] = pd.to_numeric(
            tmp['Percentage'].astype(str).str.replace(',', '.', regex=False),
            errors='coerce',
        )
        tmp = tmp[tmp['pct'].notna()].copy()

        for _, row in tmp.iterrows():
            rows.append(
                {
                    'sample': sample,
                    'group': group,
                    'chrom_num': int(row['chrom_num']),
                    'pct': float(row['pct']),
                }
            )

    if not rows:
        raise RuntimeError('No valid chromosome percentage rows were parsed from count tables.')

    all_df = pd.DataFrame(rows)
    brca2_df = all_df[all_df['group'] == 'BRCA2'].copy()
    ctrl_df = all_df[all_df['group'] == 'ctrl'].copy()
    if brca2_df.empty or ctrl_df.empty:
        raise RuntimeError('Both BRCA2 and ctrl groups are required in chromosome_count tables.')
    return brca2_df, ctrl_df


def _plot_group_chromosome_percentages(counts_dir: Path, out_dir: Path) -> None:
    """Plot BRCA2 vs ctrl chromosome percentage profiles from count-table CSVs."""
    brca2_df, ctrl_df = _read_percentage_tables(counts_dir)

    brca2_mean = brca2_df.groupby('chrom_num', as_index=True)['pct'].mean()
    ctrl_mean = ctrl_df.groupby('chrom_num', as_index=True)['pct'].mean()

    chroms = [c for c in range(1, 23) if c in brca2_mean.index and c in ctrl_mean.index]
    if not chroms:
        raise RuntimeError('No overlapping autosomes (chr1-22) found for BRCA2 and ctrl groups.')

    x = list(range(1, len(chroms) + 1))
    brca2_y = [float(brca2_mean.loc[c]) for c in chroms]
    ctrl_y = [float(ctrl_mean.loc[c]) for c in chroms]

    # Match the example style: grid background, overlapping lines, and visible markers.
    plt.style.use('ggplot')
    plt.figure(figsize=(12.8, 6.0))

    # BRCA2 first (red squares), drawn slightly underneath.
    plt.plot(
        x,
        brca2_y,
        color='red',
        marker='s',
        linewidth=2.6,
        markersize=6.8,
        alpha=0.98,
        label=f'BRCA2 (n={brca2_df["sample"].nunique()})',
        zorder=2,
    )

    # ctrl second (blue dots), drawn on top for overlap readability.
    plt.plot(
        x,
        ctrl_y,
        color='blue',
        marker='o',
        linewidth=2.0,
        markersize=6.2,
        alpha=0.72,
        label=f'ctrl (n={ctrl_df["sample"].nunique()})',
        zorder=3,
    )

    plt.xticks(x, [str(c) for c in chroms])
    plt.xlabel('Chromosome')
    plt.ylabel('% of cfDNA fragments')
    plt.title('Chromosomal cfDNA fragment burden by group')
    plt.grid(True, alpha=0.4)
    plt.legend(loc='upper right', frameon=True)
    plt.tight_layout()

    out_path = out_dir / 'chromosome_percentage_profile_ctrl_vs_BRCA2.png'
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Comparative profile saved to {out_path}')


try:
    _plot_group_chromosome_percentages(Path(args.counts_dir), Path(args.out_dir))
except Exception as exc:
    print(f'Could not generate comparative chromosome percentage profile: {exc}')
