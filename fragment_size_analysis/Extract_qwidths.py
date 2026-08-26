#!/usr/bin/env python3
"""
Universal cfDNA detailed size extraction from BAM files.

This script is designed to be reused across datasets and feeds downstream
PCA/KS-MWU scripts by writing one distribution CSV and one count CSV per
fragment class.

Outputs per fragment class:
- <fragment>_size_distributions.csv
- <fragment>_size_counts.csv

Both files include:
- bam_file
- group
- total_fragments_all
- fragment size bin columns (e.g. mononucleosomal_125bp)
"""

import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

# ============================================================================
# USER SETTINGS (edit these defaults for your environment)
# ============================================================================
DEFAULT_BAM_DIRECTORY = "/mnt/d/coronary_cfDNA"
DEFAULT_OUTPUT_DIRECTORY = "/mnt/d/coronary_cfDNA/cfDNA_detailed_size_analysis"

# If provided, use metadata sample->group mapping instead of regex inference.
DEFAULT_SAMPLE_METADATA_CSV = None
DEFAULT_METADATA_SAMPLE_COL = "sample"
DEFAULT_METADATA_GROUP_COL = "group"

# Group inference regex rules used when metadata mapping is absent.
GROUP_RULES = [
    {"label": "stemi", "pattern": r"^cor_.*barcode[_-]?(0[1-5])\b"},
    {"label": "recovery", "pattern": r"^rec_.*barcode[_-]?(0[6-9]|10)\b"},
    {"label": "control", "pattern": r"^con_.*barcode[_-]?(1[1-5])\b"},
]

# Exclude BAM names containing any of these substrings.
EXCLUDED_BAM_NAME_SUBSTRINGS = ["combined"]

BIN_SIZE = 5
SIZE_RANGES = {
    "mononucleosomal": (25, 260),
    "dinucleosomal": (260, 480),
    "trinucleosomal": (480, 650),
    "hmw": (650, 100000),
}
# ============================================================================


def normalize_runtime_path(path_str: str) -> str:
    """Normalize path style between Windows and WSL/POSIX paths."""
    path_str = str(path_str)

    if os.name == "nt" and path_str.startswith("/mnt/"):
        parts = path_str.split("/")
        if len(parts) >= 4 and len(parts[2]) == 1:
            drive = parts[2].upper() + ":"
            tail = parts[3:]
            return os.path.join(drive + os.sep, *tail)

    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", path_str):
        drive = path_str[0].lower()
        rest = path_str[2:].replace("\\", "/").lstrip("/")
        return f"/mnt/{drive}/{rest}"

    return path_str


def parse_args():
    parser = argparse.ArgumentParser(description="Universal cfDNA size-bin extraction from BAM files")
    parser.add_argument("--bam-dir", default=DEFAULT_BAM_DIRECTORY, help="BAM directory (recursive)")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIRECTORY, help="Output directory")
    parser.add_argument("--bin-size", type=int, default=BIN_SIZE, help="Bin size in bp")
    parser.add_argument(
        "--max-reads-per-bam",
        type=int,
        default=1000000,
        help="Maximum reads to process per BAM for speed",
    )
    parser.add_argument(
        "--sample-metadata-csv",
        default=DEFAULT_SAMPLE_METADATA_CSV,
        help="Optional CSV with sample->group mapping",
    )
    parser.add_argument("--metadata-sample-col", default=DEFAULT_METADATA_SAMPLE_COL)
    parser.add_argument("--metadata-group-col", default=DEFAULT_METADATA_GROUP_COL)
    return parser.parse_args()


def ensure_output_directory(out_dir: str):
    os.makedirs(out_dir, exist_ok=True)


def get_bam_files(directory):
    """Find BAM files recursively and drop likely aggregate files."""
    bam_files = glob.glob(os.path.join(directory, "**", "*.bam"), recursive=True)
    bam_files = sorted(set(bam_files))

    filtered = []
    for file_path in bam_files:
        name = os.path.basename(file_path).lower()
        if any(token.lower() in name for token in EXCLUDED_BAM_NAME_SUBSTRINGS):
            continue
        filtered.append(file_path)

    return filtered


def load_group_map(metadata_csv, sample_col, group_col):
    if metadata_csv is None:
        return None
    if not os.path.exists(metadata_csv):
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_csv}")

    meta_df = pd.read_csv(metadata_csv)
    if sample_col not in meta_df.columns or group_col not in meta_df.columns:
        raise ValueError(f"Metadata CSV must contain '{sample_col}' and '{group_col}' columns")

    mapping = {}
    for _, row in meta_df[[sample_col, group_col]].dropna().iterrows():
        key = str(row[sample_col]).strip()
        mapping[key] = str(row[group_col]).strip().lower()
    return mapping


def infer_group_from_name(sample_name):
    normalized = os.path.splitext(os.path.basename(str(sample_name)))[0].lower()
    for rule in GROUP_RULES:
        if re.search(rule["pattern"], normalized):
            return rule["label"]
    return "unknown"


def assign_group(sample_name, group_map=None):
    sample_key = str(sample_name).strip()
    if group_map is not None and sample_key in group_map:
        return group_map[sample_key], "metadata"
    return infer_group_from_name(sample_key), "inferred"


def extract_fragment_sizes(bam_file, max_reads):
    """Extract query_length values from BAM alignments for cfDNA fragment sizes."""
    print(f"Processing BAM: {os.path.basename(bam_file)}")

    fragment_sizes = []
    total_reads = 0

    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch(until_eof=True):
            total_reads += 1
            fragment_size = read.query_length
            if fragment_size and 25 <= fragment_size <= 1500:
                fragment_sizes.append(fragment_size)

            if total_reads % 100000 == 0:
                print(f"  Reads processed: {total_reads}; valid fragments: {len(fragment_sizes)}")

            if total_reads >= max_reads:
                break

    if fragment_sizes:
        print(
            f"  Final: reads={total_reads}, valid={len(fragment_sizes)}, "
            f"min={min(fragment_sizes)}, median={np.median(fragment_sizes):.1f}, max={max(fragment_sizes)}"
        )
    else:
        print(f"  Final: reads={total_reads}, no valid fragments in 25-1500bp")

    return fragment_sizes


def create_size_bins(fragment_sizes, min_size, max_size, bin_size=5):
    """Create bin centers and counts for requested size range."""
    bin_edges = np.arange(min_size, max_size + bin_size, bin_size)
    bin_centers = bin_edges[:-1] + bin_size // 2

    filtered = [size for size in fragment_sizes if min_size <= size <= max_size]
    if not filtered:
        bin_counts = np.zeros(len(bin_centers))
    else:
        bin_counts, _ = np.histogram(filtered, bins=bin_edges)

    return bin_centers, bin_counts


def process_all_bams(bam_files, bin_size, max_reads, group_map=None):
    percentage_data = {frag_type: [] for frag_type in SIZE_RANGES.keys()}
    count_data = {frag_type: [] for frag_type in SIZE_RANGES.keys()}
    sample_names = []
    sample_groups = []
    sample_group_source = []
    sample_total_fragments = []

    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).replace(".bam", "")
        group, group_source = assign_group(sample_name, group_map=group_map)

        sample_names.append(sample_name)
        sample_groups.append(group)
        sample_group_source.append(group_source)

        print("\n" + "=" * 70)
        print(f"Sample: {sample_name} | group={group} ({group_source})")
        print("=" * 70)

        fragment_sizes = extract_fragment_sizes(bam_file, max_reads=max_reads)

        if not fragment_sizes:
            sample_total_fragments.append(0)
            for frag_type in SIZE_RANGES.keys():
                min_size, max_size = SIZE_RANGES[frag_type]
                bin_centers, _ = create_size_bins([], min_size, max_size, bin_size)
                zeros = np.zeros(len(bin_centers))
                count_data[frag_type].append(zeros)
                percentage_data[frag_type].append(zeros)
            continue

        total_valid = len(fragment_sizes)
        sample_total_fragments.append(total_valid)

        for frag_type, (min_size, max_size) in SIZE_RANGES.items():
            bin_centers, bin_counts = create_size_bins(fragment_sizes, min_size, max_size, bin_size)
            count_data[frag_type].append(bin_counts)

            percentages = (bin_counts / total_valid) * 100 if total_valid > 0 else np.zeros(len(bin_counts))
            percentage_data[frag_type].append(percentages)

            in_range = int(np.sum(bin_counts))
            print(f"  {frag_type}: {in_range} fragments ({100 * in_range / total_valid:.1f}% of sample), bins={len(bin_centers)}")

    return {
        "percentage_data": percentage_data,
        "count_data": count_data,
        "sample_names": sample_names,
        "sample_groups": sample_groups,
        "sample_group_source": sample_group_source,
        "sample_total_fragments": sample_total_fragments,
    }


def save_size_distribution_data(results, output_dir, bin_size):
    ensure_output_directory(output_dir)

    percentage_data = results["percentage_data"]
    count_data = results["count_data"]
    sample_names = results["sample_names"]
    sample_groups = results["sample_groups"]
    sample_group_source = results["sample_group_source"]
    sample_total_fragments = results["sample_total_fragments"]

    for frag_type in SIZE_RANGES.keys():
        min_size, max_size = SIZE_RANGES[frag_type]
        bin_edges = np.arange(min_size, max_size + bin_size, bin_size)
        bin_centers = bin_edges[:-1] + bin_size // 2

        column_names = [f"{frag_type}_{int(center)}bp" for center in bin_centers]

        pct_matrix = np.array(percentage_data[frag_type])
        cnt_matrix = np.array(count_data[frag_type])

        df_pct = pd.DataFrame(pct_matrix, columns=column_names)
        df_pct.insert(0, "bam_file", sample_names)
        df_pct.insert(1, "group", sample_groups)
        df_pct.insert(2, "group_source", sample_group_source)
        df_pct.insert(3, "total_fragments_all", sample_total_fragments)

        df_cnt = pd.DataFrame(cnt_matrix, columns=column_names)
        df_cnt.insert(0, "bam_file", sample_names)
        df_cnt.insert(1, "group", sample_groups)
        df_cnt.insert(2, "group_source", sample_group_source)
        df_cnt.insert(3, "total_fragments_all", sample_total_fragments)
        df_cnt["total_fragments_in_range"] = np.sum(cnt_matrix, axis=1)

        pct_file = os.path.join(output_dir, f"{frag_type}_size_distributions.csv")
        cnt_file = os.path.join(output_dir, f"{frag_type}_size_counts.csv")
        df_pct.to_csv(pct_file, index=False)
        df_cnt.to_csv(cnt_file, index=False)

        print(f"Saved {frag_type} percentage data: {pct_file}")
        print(f"Saved {frag_type} count data: {cnt_file}")


def create_summary_plot(results, output_dir, bin_size):
    """Create simple mean-per-group summary curves for each fragment class."""
    percentage_data = results["percentage_data"]
    sample_groups = results["sample_groups"]

    unique_groups = sorted(set(sample_groups))
    color_map = {
        "stemi": "#d62728",
        "recovery": "#1f77b4",
        "control": "#2ca02c",
        "unknown": "#7f7f7f",
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle("cfDNA Fragment Size Distributions by Fragment Type", fontsize=16, fontweight="bold")

    for idx, (frag_type, (min_size, max_size)) in enumerate(SIZE_RANGES.items()):
        row = idx // 2
        col = idx % 2
        ax = axes[row, col]

        bin_edges = np.arange(min_size, max_size + bin_size, bin_size)
        bin_centers = bin_edges[:-1] + bin_size // 2
        data_matrix = np.array(percentage_data[frag_type])

        for group in unique_groups:
            indices = [i for i, g in enumerate(sample_groups) if g == group]
            if not indices:
                continue
            mean_curve = np.mean(data_matrix[indices], axis=0)
            ax.plot(bin_centers, mean_curve, label=f"{group} (n={len(indices)})", linewidth=2, color=color_map.get(group, None))

        ax.set_xlabel("Fragment Size (bp)")
        ax.set_ylabel("Percentage of sample (%)")
        ax.set_title(frag_type)
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.tight_layout()
    out_path = os.path.join(output_dir, "fragment_size_distributions_summary.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Summary plot saved: {out_path}")


def main():
    args = parse_args()

    bam_dir = normalize_runtime_path(args.bam_dir)
    out_dir = normalize_runtime_path(args.output_dir)
    metadata_csv = normalize_runtime_path(args.sample_metadata_csv) if args.sample_metadata_csv else None

    print("=== Universal cfDNA Detailed Size Extraction ===")
    print(f"BAM directory: {bam_dir}")
    print(f"Output directory: {out_dir}")
    print(f"Bin size: {args.bin_size} bp")
    print(f"Metadata CSV: {metadata_csv}")

    group_map = load_group_map(
        metadata_csv,
        sample_col=args.metadata_sample_col,
        group_col=args.metadata_group_col,
    )

    bam_files = get_bam_files(bam_dir)
    if not bam_files:
        raise FileNotFoundError(f"No BAM files found under {bam_dir}")

    print(f"Found {len(bam_files)} BAM files")
    for file_path in bam_files[:20]:
        print(f"  - {os.path.basename(file_path)}")
    if len(bam_files) > 20:
        print(f"  ... and {len(bam_files) - 20} more")

    results = process_all_bams(
        bam_files=bam_files,
        bin_size=args.bin_size,
        max_reads=args.max_reads_per_bam,
        group_map=group_map,
    )

    save_size_distribution_data(results, out_dir, bin_size=args.bin_size)
    create_summary_plot(results, out_dir, bin_size=args.bin_size)

    print("\n=== Extraction Complete ===")
    print(f"Samples processed: {len(results['sample_names'])}")
    print(f"Output directory: {out_dir}")


if __name__ == "__main__":
    main()
