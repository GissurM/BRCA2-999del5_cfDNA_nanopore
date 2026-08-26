#!/usr/bin/env python3
"""
Generic TelSeq Results Parser
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


def parse_range_token(token: str) -> List[str]:
    token = token.strip()
    if not token:
        return []
    if "-" not in token:
        return [token.zfill(2)]
    left, right = token.split("-", 1)
    left_i = int(left)
    right_i = int(right)
    if right_i < left_i:
        left_i, right_i = right_i, left_i
    return [str(i).zfill(2) for i in range(left_i, right_i + 1)]


def parse_group_map(raw: str) -> Dict[str, str]:
    """Parse barcode-to-group mapping.

    Format example:
      BRCA2:01-06,13-18;Control:07-12,19-24
    """
    mapping: Dict[str, str] = {}
    if not raw.strip():
        return mapping

    for block in raw.split(";"):
        block = block.strip()
        if not block:
            continue
        if ":" not in block:
            raise ValueError(f"Invalid group-map block '{block}'. Use Group:codes format.")
        group, ranges = block.split(":", 1)
        group = group.strip()
        if not group:
            raise ValueError(f"Invalid group name in block '{block}'.")

        codes: List[str] = []
        for tok in ranges.split(","):
            codes.extend(parse_range_token(tok))
        for code in codes:
            mapping[code] = group

    return mapping


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Parse TelSeq output into analysis-ready CSV")
    parser.add_argument(
        "--input-file",
        default="./data/all_samples_telomere_results.txt",
        help="TelSeq output text file path",
    )
    parser.add_argument(
        "--output-csv",
        default="./results/telseq_correct_summary.csv",
        help="Output summary CSV path",
    )
    parser.add_argument(
        "--group-map",
        default="BRCA2:01-06,13-18;Control:07-12,19-24",
        help="Barcode grouping map: Group:codes;Group:codes with codes like 01-06,13,15",
    )
    parser.add_argument(
        "--barcode-regex",
        default=r"barcode(\d+)",
        help="Regex used to extract barcode from ReadGroup text",
    )
    parser.add_argument(
        "--skip-readgroup-prefix",
        default="dd7255dd3f1acd4f81c77bc70fabc9f91ae28526_dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
        help="Optional ReadGroup prefix to skip (set empty string to disable)",
    )
    parser.add_argument(
        "--keep-unknown",
        action="store_true",
        help="Keep samples not present in --group-map as Group=Unknown",
    )
    return parser


def parse_telseq_correctly(args: argparse.Namespace) -> pd.DataFrame:
    """Parse TelSeq results with configurable mappings and paths."""

    input_file = Path(args.input_file)
    output_file = Path(args.output_csv)
    group_map = parse_group_map(args.group_map)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    results = []

    with open(input_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # Skip header and empty lines
            if line.startswith("ReadGroup") or line.strip() == "":
                continue

            # Look for lines with actual sequence data and barcode
            if "barcode" in line and not (
                args.skip_readgroup_prefix and line.strip().startswith(f"{args.skip_readgroup_prefix}\t")
            ):

                # Split by tabs
                parts = line.strip().split("\t")

                # Skip if not enough parts or if it contains UNKNOWN
                if len(parts) < 7 or "UNKNOWN" in line:
                    continue

                try:
                    # Extract barcode from ReadGroup (first column)
                    barcode_match = re.search(args.barcode_regex, parts[0], flags=re.IGNORECASE)
                    if not barcode_match:
                        continue

                    barcode = barcode_match.group(1).zfill(2)

                    # Extract the numeric values
                    # Format: ReadGroup, Library, Sample, Total, Mapped, Duplicates, LENGTH_ESTIMATE
                    total_reads = int(parts[3])
                    mapped_reads = int(parts[4])
                    duplicates = int(parts[5])
                    length_estimate = float(parts[6])

                    # Only process samples with actual data
                    if total_reads > 0 and mapped_reads > 0:

                        # Determine group based on barcode
                        group = group_map.get(barcode, "Unknown")
                        if group == "Unknown" and not args.keep_unknown:
                            continue

                        results.append({
                            "Barcode": barcode,
                            "Group": group,
                            "Total_Reads": total_reads,
                            "Mapped_Reads": mapped_reads,
                            "Duplicate_Reads": duplicates,
                            "Telomere_Length_Estimate": length_estimate,
                            "Mapping_Rate": round(mapped_reads / total_reads * 100, 2),
                            "Duplicate_Rate": round(duplicates / total_reads * 100, 2) if total_reads > 0 else 0,
                        })

                except (ValueError, IndexError) as e:
                    print(f"Error parsing line: {e}")
                    continue

    # Convert to DataFrame
    df = pd.DataFrame(results)

    if len(df) == 0:
        print("No valid data found!")
        return df

    # Sort by barcode number
    df["Barcode_Num"] = df["Barcode"].astype(int)
    df = df.sort_values("Barcode_Num").drop("Barcode_Num", axis=1).reset_index(drop=True)

    # Save to CSV
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)

    print(f"Correct summary saved to: {output_file}")
    print(f"Successfully extracted data for {len(df)} samples with valid telomere data")

    # Print comprehensive analysis
    print("\n" + "="*80)
    print("TELOMERE LENGTH ANALYSIS - CORRECTED RESULTS")
    print("="*80)

    print(f"\nSample Distribution:")
    group_counts = df["Group"].value_counts()
    for group, count in group_counts.items():
        print(f"  {group}: {count} samples")

    print(f"\nTelomere Length Statistics:")
    print(f"  Overall Mean: {df['Telomere_Length_Estimate'].mean():.3f}")
    print(f"  Overall Median: {df['Telomere_Length_Estimate'].median():.3f}")
    print(f"  Standard Deviation: {df['Telomere_Length_Estimate'].std():.3f}")
    print(f"  Range: {df['Telomere_Length_Estimate'].min():.3f} - {df['Telomere_Length_Estimate'].max():.3f}")

    print(f"\nGroup Comparison:")
    for group in sorted(df["Group"].unique()):
        group_data = df[df["Group"] == group]
        lengths = group_data["Telomere_Length_Estimate"]
        print(f"  {group:>7}: Mean={lengths.mean():7.3f}, Median={lengths.median():7.3f}, "
              f"SD={lengths.std():6.3f}, n={len(lengths)}")

    print(f"\nSequencing Quality Metrics:")
    print(f"  Average mapping rate: {df['Mapping_Rate'].mean():.1f}%")
    print(f"  Average duplicate rate: {df['Duplicate_Rate'].mean():.1f}%")
    print(f"  Total reads range: {df['Total_Reads'].min():,} - {df['Total_Reads'].max():,}")
    print(f"  Mapped reads range: {df['Mapped_Reads'].min():,} - {df['Mapped_Reads'].max():,}")

    print(f"\nDetailed Sample Results:")
    print(f"{'Group':>7} {'Barcode':>7} {'Tel_Length':>11} {'Total_Reads':>12} {'Mapped_Reads':>13} {'Mapping%':>9}")
    print("-" * 70)
    for _, row in df.iterrows():
        print(f"{row['Group']:>7} {row['Barcode']:>7} {row['Telomere_Length_Estimate']:>11.3f} "
              f"{row['Total_Reads']:>12,} {row['Mapped_Reads']:>13,} {row['Mapping_Rate']:>8.1f}%")

    # Statistical analysis
    groups = sorted(df["Group"].unique().tolist())
    if len(groups) == 2:
        g1, g2 = groups
        g1_lengths = df[df["Group"] == g1]["Telomere_Length_Estimate"]
        g2_lengths = df[df["Group"] == g2]["Telomere_Length_Estimate"]
        if len(g1_lengths) > 0 and len(g2_lengths) > 0:
            print(f"\nStatistical Comparison:")
            print(f"  {g1} mean: {g1_lengths.mean():.3f} ± {g1_lengths.std():.3f}")
            print(f"  {g2} mean: {g2_lengths.mean():.3f} ± {g2_lengths.std():.3f}")
            print(f"  Difference ({g1} - {g2}): {g1_lengths.mean() - g2_lengths.mean():+.3f}")

            pooled_std = ((g1_lengths.std() ** 2 + g2_lengths.std() ** 2) / 2) ** 0.5
            effect_size = (g1_lengths.mean() - g2_lengths.mean()) / pooled_std if pooled_std > 0 else 0.0
            print(f"  Effect size (Cohen's d): {effect_size:.3f}")

            if abs(effect_size) < 0.2:
                interpretation = "negligible"
            elif abs(effect_size) < 0.5:
                interpretation = "small"
            elif abs(effect_size) < 0.8:
                interpretation = "medium"
            else:
                interpretation = "large"
            print(f"  Effect size interpretation: {interpretation}")

    print(f"\n" + "="*80)
    print("INTERPRETATION AND RECOMMENDATIONS")
    print("="*80)

    # Provide interpretation
    overall_mean = df["Telomere_Length_Estimate"].mean()
    if overall_mean > 10:
        tel_status = "relatively long"
    elif overall_mean > 5:
        tel_status = "moderate length"
    else:
        tel_status = "relatively short"

    print(f"\nTelomere Length Interpretation:")
    print(f"  - Overall telomere lengths are {tel_status} (mean = {overall_mean:.2f})")
    print(f"  - Range indicates substantial variation between samples")

    mapping_rate = df["Mapping_Rate"].mean()
    if mapping_rate > 90:
        qual_status = "excellent"
    elif mapping_rate > 80:
        qual_status = "good"
    elif mapping_rate > 70:
        qual_status = "acceptable"
    else:
        qual_status = "concerning"
    
    print(f"\nData Quality Assessment:")
    print(f"  - Mapping rates are {qual_status} (average {mapping_rate:.1f}%)")
    if mapping_rate < 80:
        print(f"  - Consider investigating low mapping rates")

    return df


if __name__ == "__main__":
    cli_args = build_arg_parser().parse_args()
    df = parse_telseq_correctly(cli_args)
