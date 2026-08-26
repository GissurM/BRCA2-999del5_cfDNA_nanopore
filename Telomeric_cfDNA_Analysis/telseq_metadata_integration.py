#!/usr/bin/env python3
"""
Telomere Analysis with Metadata Integration
Combines telomere length data with metadata covariates in a reusable pipeline.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge TelSeq summary with metadata and report covariate diagnostics")
    parser.add_argument(
        "--telomere-csv",
        default="./results/telseq_correct_summary.csv",
        help="Path to telomere summary CSV",
    )
    parser.add_argument(
        "--metadata-csv",
        default="./data/metadata_age_gender.csv",
        help="Path to metadata CSV",
    )
    parser.add_argument(
        "--output-csv",
        default="./results/telomere_with_metadata.csv",
        help="Output merged CSV path",
    )
    parser.add_argument(
        "--barcode-column-telomere",
        default="Barcode",
        help="Barcode column name in telomere CSV",
    )
    parser.add_argument(
        "--barcode-column-metadata",
        default="Filename",
        help="Source column in metadata used to extract barcode if --metadata-barcode-regex is set",
    )
    parser.add_argument(
        "--metadata-barcode-regex",
        default=r"Barcode(\d+)",
        help="Regex to extract barcode from metadata source column",
    )
    parser.add_argument("--group-column", default="Group", help="Group column in telomere CSV")
    parser.add_argument("--length-column", default="Telomere_Length_Estimate", help="Telomere length column")
    parser.add_argument("--age-column", default="age", help="Age column")
    parser.add_argument("--gender-column", default="gender", help="Gender column")
    parser.add_argument("--bmi-column", default="BMI", help="BMI column")
    parser.add_argument(
        "--write-visualization-script",
        action="store_true",
        help="Also write a standalone visualization script",
    )
    parser.add_argument(
        "--visualization-script-path",
        default="./comprehensive_telomere_visualization.py",
        help="Path for generated visualization script",
    )
    return parser


def load_and_merge_data(args: argparse.Namespace) -> pd.DataFrame:
    """Load telomere data and metadata, then merge them."""

    print("="*80)
    print("TELOMERE DATA INTEGRATION WITH METADATA")
    print("="*80)

    telomere_path = Path(args.telomere_csv)
    metadata_path = Path(args.metadata_csv)
    if not telomere_path.exists():
        raise FileNotFoundError(f"Telomere CSV not found: {telomere_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_path}")

    # Load telomere data
    telomere_df = pd.read_csv(telomere_path)

    # Load metadata
    metadata_df = pd.read_csv(metadata_path)

    print(f"Telomere data: {len(telomere_df)} samples")
    print(f"Metadata: {len(metadata_df)} samples")

    if args.barcode_column_telomere not in telomere_df.columns:
        raise ValueError(f"Missing telomere barcode column: {args.barcode_column_telomere}")
    if args.barcode_column_metadata not in metadata_df.columns:
        raise ValueError(f"Missing metadata barcode source column: {args.barcode_column_metadata}")

    # Extract barcode numbers from metadata filename
    metadata_df["Barcode"] = metadata_df[args.barcode_column_metadata].astype(str).str.extract(args.metadata_barcode_regex)[0]
    if metadata_df["Barcode"].isna().all():
        raise ValueError("No barcodes extracted from metadata. Check --metadata-barcode-regex and source column.")

    # Ensure barcode is string format for matching
    telomere_df["Barcode"] = telomere_df[args.barcode_column_telomere].astype(str).str.extract(r"(\d+)")[0].str.zfill(2)
    metadata_df["Barcode"] = metadata_df["Barcode"].astype(str).str.extract(r"(\d+)")[0].str.zfill(2)

    # Merge datasets
    merged_df = pd.merge(telomere_df, metadata_df, on="Barcode", how="inner")

    print(f"Successfully merged: {len(merged_df)} samples")

    required = [args.group_column, args.length_column, args.age_column, args.gender_column, args.bmi_column]
    missing = [c for c in required if c not in merged_df.columns]
    if missing:
        raise ValueError(f"Merged dataframe missing required columns: {missing}")

    # Display merged data structure
    print(f"\nMerged dataset columns: {list(merged_df.columns)}")
    print(f"\nFirst few rows:")
    print(merged_df[["Barcode", args.group_column, args.length_column, args.age_column, args.gender_column, args.bmi_column]].head())

    return merged_df


def analyze_metadata_distribution(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Analyze the distribution of metadata variables"""

    print(f"\n" + "="*80)
    print("METADATA DISTRIBUTION ANALYSIS")
    print("="*80)

    print(f"\nSample distribution by group:")
    print(df[args.group_column].value_counts())

    print(f"\nGender distribution:")
    gender_by_group = pd.crosstab(df[args.group_column], df[args.gender_column])
    print(gender_by_group)

    print(f"\nAge statistics by group:")
    age_stats = df.groupby(args.group_column)[args.age_column].agg(['count', 'mean', 'std', 'min', 'max']).round(2)
    print(age_stats)

    print(f"\nBMI statistics by group:")
    bmi_stats = df.groupby(args.group_column)[args.bmi_column].agg(['count', 'mean', 'std', 'min', 'max']).round(2)
    print(bmi_stats)

    print(f"\nAge distribution across all samples:")
    print(f"  Range: {df[args.age_column].min()} - {df[args.age_column].max()} years")
    print(f"  Mean: {df[args.age_column].mean():.1f} ± {df[args.age_column].std():.1f} years")

    print(f"\nBMI distribution across all samples:")
    print(f"  Range: {df[args.bmi_column].min()} - {df[args.bmi_column].max()}")
    print(f"  Mean: {df[args.bmi_column].mean():.1f} ± {df[args.bmi_column].std():.1f}")

    # Check for group differences in covariates
    print(f"\n" + "-"*50)
    print("GROUP COMPARISONS FOR COVARIATES")
    print("-"*50)

    unique_groups = sorted(df[args.group_column].dropna().unique().tolist())
    if len(unique_groups) == 2:
        g1, g2 = unique_groups
        g1_age = df[df[args.group_column] == g1][args.age_column]
        g2_age = df[df[args.group_column] == g2][args.age_column]
        age_ttest = stats.ttest_ind(g1_age, g2_age)
        print(f"\nAge comparison:")
        print(f"  {g1}: {g1_age.mean():.1f} ± {g1_age.std():.1f} years")
        print(f"  {g2}: {g2_age.mean():.1f} ± {g2_age.std():.1f} years")
        print(f"  t-test p-value: {age_ttest.pvalue:.4f}")

        g1_bmi = df[df[args.group_column] == g1][args.bmi_column]
        g2_bmi = df[df[args.group_column] == g2][args.bmi_column]
        bmi_ttest = stats.ttest_ind(g1_bmi, g2_bmi)
        print(f"\nBMI comparison:")
        print(f"  {g1}: {g1_bmi.mean():.1f} ± {g1_bmi.std():.1f}")
        print(f"  {g2}: {g2_bmi.mean():.1f} ± {g2_bmi.std():.1f}")
        print(f"  t-test p-value: {bmi_ttest.pvalue:.4f}")

    # Gender comparison (Chi-square)
    chi2, chi2_p = stats.chi2_contingency(gender_by_group)[:2]
    print(f"\nGender distribution comparison:")
    print(f"  Chi-square p-value: {chi2_p:.4f}")

    return df


def correlation_analysis(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Analyze correlations between telomere length and covariates"""

    print(f"\n" + "="*80)
    print("CORRELATION ANALYSIS: TELOMERE LENGTH vs COVARIATES")
    print("="*80)

    # Age correlation
    age_corr, age_p = pearsonr(df[args.age_column], df[args.length_column])
    age_spearman, age_sp = spearmanr(df[args.age_column], df[args.length_column])

    print(f"\nAge vs Telomere Length:")
    print(f"  Pearson r = {age_corr:.3f}, p = {age_p:.4f}")
    print(f"  Spearman ρ = {age_spearman:.3f}, p = {age_sp:.4f}")

    # BMI correlation
    bmi_corr, bmi_p = pearsonr(df[args.bmi_column], df[args.length_column])
    bmi_spearman, bmi_sp = spearmanr(df[args.bmi_column], df[args.length_column])

    print(f"\nBMI vs Telomere Length:")
    print(f"  Pearson r = {bmi_corr:.3f}, p = {bmi_p:.4f}")
    print(f"  Spearman ρ = {bmi_spearman:.3f}, p = {bmi_sp:.4f}")

    # Gender comparison
    gender_vals = sorted(df[args.gender_column].dropna().astype(str).unique().tolist())
    if len(gender_vals) == 2:
        g1, g2 = gender_vals
        g1_tel = df[df[args.gender_column] == g1][args.length_column]
        g2_tel = df[df[args.gender_column] == g2][args.length_column]
        gender_ttest = stats.ttest_ind(g1_tel, g2_tel)

        print(f"\nGender vs Telomere Length:")
        print(f"  {g1} (n={len(g1_tel)}): {g1_tel.mean():.3f} ± {g1_tel.std():.3f}")
        print(f"  {g2} (n={len(g2_tel)}): {g2_tel.mean():.3f} ± {g2_tel.std():.3f}")
        print(f"  t-test p-value: {gender_ttest.pvalue:.4f}")

    print(f"\nGender vs Telomere Length:")
    if len(gender_vals) != 2:
        print("  Skipped: gender comparison requires exactly 2 categories")

    # Group-specific correlations
    print(f"\n" + "-"*50)
    print("GROUP-SPECIFIC CORRELATIONS")
    print("-"*50)

    for group in sorted(df[args.group_column].dropna().unique().tolist()):
        group_data = df[df[args.group_column] == group]

        print(f"\n{group} Group (n={len(group_data)}):")

        # Age correlation within group
        if len(group_data) > 3:  # Need minimum samples for correlation
            g_age_corr, g_age_p = pearsonr(group_data[args.age_column], group_data[args.length_column])
            g_bmi_corr, g_bmi_p = pearsonr(group_data[args.bmi_column], group_data[args.length_column])

            print(f"  Age correlation: r = {g_age_corr:.3f}, p = {g_age_p:.4f}")
            print(f"  BMI correlation: r = {g_bmi_corr:.3f}, p = {g_bmi_p:.4f}")

    return df


def save_merged_dataset(df: pd.DataFrame, args: argparse.Namespace) -> Path:
    """Save the merged dataset for further analysis"""

    output_file = Path(args.output_csv)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)

    print(f"\n" + "="*80)
    print("DATASET SAVED FOR VISUALIZATION")
    print("="*80)
    print(f"Merged dataset saved to: {output_file}")
    print(f"Ready for visualization with {len(df)} samples")

    # Show final dataset summary
    print(f"\nFinal dataset summary:")
    groups = sorted(df[args.group_column].dropna().unique().tolist())
    if len(groups) >= 2:
        print(f"{'Variable':<20} {groups[0]:<15} {groups[1]:<15} {'Total':<10}")
    else:
        print(f"{'Variable':<20} {'Group':<15} {'N/A':<15} {'Total':<10}")
    print("-" * 65)

    for var in [args.age_column, args.bmi_column]:
        g1_vals = df[df[args.group_column] == groups[0]][var] if len(groups) >= 1 else pd.Series(dtype=float)
        g2_vals = df[df[args.group_column] == groups[1]][var] if len(groups) >= 2 else pd.Series(dtype=float)
        total_mean = df[var].mean()
        g1_str = f"{g1_vals.mean():.1f}±{g1_vals.std():.1f}" if len(g1_vals) else "NA"
        g2_str = f"{g2_vals.mean():.1f}±{g2_vals.std():.1f}" if len(g2_vals) else "NA"
        print(f"{var:<20} {g1_str:<15} {g2_str:<15} {total_mean:.1f}")

    # Telomere length
    g1_tel = df[df[args.group_column] == groups[0]][args.length_column] if len(groups) >= 1 else pd.Series(dtype=float)
    g2_tel = df[df[args.group_column] == groups[1]][args.length_column] if len(groups) >= 2 else pd.Series(dtype=float)
    total_tel = df[args.length_column]
    g1_tel_str = f"{g1_tel.mean():.2f}±{g1_tel.std():.2f}" if len(g1_tel) else "NA"
    g2_tel_str = f"{g2_tel.mean():.2f}±{g2_tel.std():.2f}" if len(g2_tel) else "NA"
    print(f"{'Telomere Length':<20} {g1_tel_str:<15} {g2_tel_str:<15} {total_tel.mean():.2f}")

    print(f"\nGender distribution:")
    print(df.groupby([args.group_column, args.gender_column]).size().unstack(fill_value=0))

    return output_file

def create_visualization_script(args: argparse.Namespace) -> Path:
    """Create a reusable standalone visualization script."""

    viz_script = f'''#!/usr/bin/env python3
"""
Comprehensive Telomere Visualization with Metadata (Generic)
"""

from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LinearRegression


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Create comprehensive telomere visualization from merged CSV")
    p.add_argument("--input-csv", default="{args.output_csv}")
    p.add_argument("--output-png", default="./results/comprehensive_telomere_analysis.png")
    p.add_argument("--group-column", default="{args.group_column}")
    p.add_argument("--length-column", default="{args.length_column}")
    p.add_argument("--age-column", default="{args.age_column}")
    p.add_argument("--gender-column", default="{args.gender_column}")
    p.add_argument("--bmi-column", default="{args.bmi_column}")
    return p


def main() -> None:
    args = build_parser().parse_args()
    df = pd.read_csv(args.input_csv)

    plt.style.use("default")
    sns.set_palette("husl")
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["savefig.dpi"] = 300

    groups = sorted(df[args.group_column].dropna().unique().tolist())
    palette = dict(zip(groups, sns.color_palette("tab10", n_colors=max(1, len(groups)))))

    fig = plt.figure(figsize=(20, 16))

    ax1 = plt.subplot(3, 4, 1)
    sns.boxplot(data=df, x=args.group_column, y=args.length_column, ax=ax1)
    sns.stripplot(data=df, x=args.group_column, y=args.length_column, color="black", alpha=0.6, ax=ax1)
    ax1.set_title("Telomere Length by Group")

    ax2 = plt.subplot(3, 4, 2)
    for g in groups:
        sub = df[df[args.group_column] == g]
        ax2.scatter(sub[args.age_column], sub[args.length_column], alpha=0.7, s=60, label=g, color=palette[g])
        if len(sub) > 2:
            z = np.polyfit(sub[args.age_column], sub[args.length_column], 1)
            p = np.poly1d(z)
            corr, pval = stats.pearsonr(sub[args.age_column], sub[args.length_column])
            xs = np.sort(sub[args.age_column].to_numpy())
            ax2.plot(xs, p(xs), linestyle="--", alpha=0.8, color=palette[g], label=f"{{g}}: r={{corr:.2f}}, p={{pval:.3f}}")
    ax2.set_xlabel("Age")
    ax2.set_ylabel("Telomere length")
    ax2.set_title("Age vs Telomere Length")
    ax2.legend(fontsize=8)

    ax3 = plt.subplot(3, 4, 3)
    for g in groups:
        sub = df[df[args.group_column] == g]
        ax3.scatter(sub[args.bmi_column], sub[args.length_column], alpha=0.7, s=60, label=g, color=palette[g])
        if len(sub) > 2:
            z = np.polyfit(sub[args.bmi_column], sub[args.length_column], 1)
            p = np.poly1d(z)
            corr, pval = stats.pearsonr(sub[args.bmi_column], sub[args.length_column])
            xs = np.sort(sub[args.bmi_column].to_numpy())
            ax3.plot(xs, p(xs), linestyle="--", alpha=0.8, color=palette[g], label=f"{{g}}: r={{corr:.2f}}, p={{pval:.3f}}")
    ax3.set_xlabel("BMI")
    ax3.set_ylabel("Telomere length")
    ax3.set_title("BMI vs Telomere Length")
    ax3.legend(fontsize=8)

    ax4 = plt.subplot(3, 4, 4)
    age_reg = LinearRegression()
    age_reg.fit(df[[args.age_column]], df[args.length_column])
    df["age_adjusted_telomere"] = df[args.length_column] - age_reg.predict(df[[args.age_column]])
    sns.boxplot(data=df, x=args.gender_column, y="age_adjusted_telomere", hue=args.group_column, ax=ax4)
    ax4.set_title("Age-adjusted length by gender")

    ax5 = plt.subplot(3, 4, 5)
    sns.histplot(data=df, x=args.age_column, hue=args.group_column, multiple="stack", bins=8, ax=ax5)
    ax5.set_title("Age Distribution")

    ax6 = plt.subplot(3, 4, 6)
    sns.histplot(data=df, x=args.bmi_column, hue=args.group_column, multiple="stack", bins=8, ax=ax6)
    ax6.set_title("BMI Distribution")

    ax7 = plt.subplot(3, 4, 7)
    order = df.sort_values([args.group_column, args.length_column]).reset_index(drop=True)
    colors = [palette[g] for g in order[args.group_column]]
    ax7.bar(np.arange(len(order)), order[args.length_column], color=colors, alpha=0.75)
    ax7.set_title("Individual Sample Results")
    ax7.set_xlabel("Sample")
    ax7.set_ylabel("Telomere length")
    ax7.tick_params(axis="x", labelbottom=False)

    ax8 = plt.subplot(3, 4, 8)
    tmp = df.copy()
    tmp["group_numeric"] = pd.factorize(tmp[args.group_column])[0]
    tmp["gender_numeric"] = pd.factorize(tmp[args.gender_column])[0]
    corr_cols = [args.length_column, args.age_column, args.bmi_column, "gender_numeric", "group_numeric"]
    sns.heatmap(tmp[corr_cols].corr(), annot=True, cmap="coolwarm", center=0, ax=ax8)
    ax8.set_title("Correlation Matrix")

    ax9 = plt.subplot(3, 4, 9)
    sns.boxplot(data=df, x=args.group_column, y="age_adjusted_telomere", ax=ax9)
    sns.stripplot(data=df, x=args.group_column, y="age_adjusted_telomere", color="black", alpha=0.6, ax=ax9)
    ax9.set_title("Age-Adjusted Telomere Length")

    ax10 = plt.subplot(3, 4, 10)
    for g in groups:
        sub = df[df[args.group_column] == g]
        ax10.scatter(sub[args.age_column], sub[args.length_column], s=sub[args.bmi_column] * 3, alpha=0.6, label=g, color=palette[g])
    ax10.set_title("Age vs Telomere (size=BMI)")
    ax10.set_xlabel("Age")
    ax10.set_ylabel("Telomere length")
    ax10.legend(fontsize=8)

    ax11 = plt.subplot(3, 4, 11)
    summary = df.groupby(args.group_column).agg(
        age_mean=(args.age_column, "mean"),
        bmi_mean=(args.bmi_column, "mean"),
        tel_mean=(args.length_column, "mean"),
    ).reset_index()
    x = np.arange(len(summary))
    w = 0.25
    ax11.bar(x - w, summary["age_mean"] / 10.0, w, label="Age/10", alpha=0.7)
    ax11.bar(x, summary["bmi_mean"] / 10.0, w, label="BMI/10", alpha=0.7)
    ax11.bar(x + w, summary["tel_mean"], w, label="Tel length", alpha=0.7)
    ax11.set_xticks(x)
    ax11.set_xticklabels(summary[args.group_column].tolist())
    ax11.set_title("Mean Characteristics by Group")
    ax11.legend(fontsize=8)

    ax12 = plt.subplot(3, 4, 12)
    ax12.axis("off")
    overall_corr, overall_p = stats.pearsonr(df[args.age_column], df[args.length_column])
    txt_lines = [
        "Statistical Summary:",
        f"Samples: {{len(df)}}",
        f"Groups: {{', '.join(groups)}}",
        f"Overall age-length r={{overall_corr:.3f}}, p={{overall_p:.4f}}",
        "",
    ]
    for g in groups:
        sub = df[df[args.group_column] == g]
        txt_lines.append(
            f"{{g}}: mean={{sub[args.length_column].mean():.3f}} ± {{sub[args.length_column].std():.3f}} (n={{len(sub)}})"
        )
    ax12.text(0.02, 0.95, "\n".join(txt_lines), va="top", fontsize=10, family="monospace")

    plt.tight_layout()
    plt.savefig(args.output_png, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Saved: {{args.output_png}}")


if __name__ == "__main__":
    main()
'''

    script_path = Path(args.visualization_script_path)
    script_path.parent.mkdir(parents=True, exist_ok=True)
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(viz_script)

    print(f"\nComprehensive visualization script created: {script_path}")
    print(f"Run with: python {script_path}")

    return script_path


def main() -> pd.DataFrame:
    """Main analysis function"""

    args = build_arg_parser().parse_args()

    # Load and merge data
    merged_df = load_and_merge_data(args)

    # Analyze metadata distribution
    merged_df = analyze_metadata_distribution(merged_df, args)

    # Correlation analysis
    merged_df = correlation_analysis(merged_df, args)

    # Save merged dataset
    output_file = save_merged_dataset(merged_df, args)

    viz_script: Optional[Path] = None
    if args.write_visualization_script:
        viz_script = create_visualization_script(args)

    print(f"\n" + "="*80)
    print("ANALYSIS COMPLETE - READY FOR VISUALIZATION")
    print("="*80)
    print(f"✅ Merged dataset: {output_file}")
    if viz_script is not None:
        print(f"✅ Visualization script: {viz_script}")
    print(f"✅ {len(merged_df)} samples with complete data")

    return merged_df


if __name__ == "__main__":
    df = main()
