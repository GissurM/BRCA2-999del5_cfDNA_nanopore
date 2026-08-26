#!/usr/bin/env python3
"""
Comprehensive Telomere Visualization with Metadata
Creates publication-quality figures from a merged telomere+metadata CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create comprehensive telomere visualizations from merged CSV")
    parser.add_argument(
        "--input-csv",
        default="./results/telomere_with_metadata.csv",
        help="Merged telomere+metadata CSV",
    )
    parser.add_argument(
        "--output-png",
        default="./results/comprehensive_telomere_analysis.png",
        help="Output figure path",
    )
    parser.add_argument(
        "--output-summary-csv",
        default="./results/comprehensive_telomere_summary.csv",
        help="Optional summary statistics CSV path",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive matplotlib window (recommended for servers/CI)",
    )
    parser.add_argument("--group-column", default="Group", help="Group column name")
    parser.add_argument("--length-column", default="Telomere_Length_Estimate", help="Telomere length column")
    parser.add_argument("--age-column", default="age", help="Age column")
    parser.add_argument("--gender-column", default="gender", help="Gender column")
    parser.add_argument("--bmi-column", default="BMI", help="BMI column")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    input_csv = Path(args.input_csv)
    output_png = Path(args.output_png)
    output_summary_csv = Path(args.output_summary_csv) if args.output_summary_csv else None

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    output_png.parent.mkdir(parents=True, exist_ok=True)
    if output_summary_csv is not None:
        output_summary_csv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)
    required = [args.group_column, args.length_column, args.age_column, args.gender_column, args.bmi_column]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")

    if len(df) < 3:
        raise RuntimeError("Need at least 3 samples for this visualization.")

    plt.style.use('default')
    sns.set_palette("husl")
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

    groups = sorted(df[args.group_column].dropna().unique().tolist())
    if len(groups) < 1:
        raise RuntimeError("No groups found in the selected group column.")
    palette = dict(zip(groups, sns.color_palette("tab10", n_colors=max(1, len(groups)))))

    fig = plt.figure(figsize=(20, 16))

    # Plot 1: Basic group comparison
    ax1 = plt.subplot(3, 4, 1)
    sns.boxplot(data=df, x=args.group_column, y=args.length_column, ax=ax1)
    sns.stripplot(data=df, x=args.group_column, y=args.length_column, color='black', alpha=0.6, ax=ax1)
    ax1.set_title('Telomere Length by Group')
    ax1.set_ylabel('Telomere Length')

    # Plot 2: Age vs Telomere Length with group-wise trends
    ax2 = plt.subplot(3, 4, 2)
    for group in groups:
        group_data = df[df[args.group_column] == group]
        ax2.scatter(group_data[args.age_column], group_data[args.length_column], c=[palette[group]], alpha=0.7, s=60, label=group)
        if len(group_data) > 2:
            z = np.polyfit(group_data[args.age_column], group_data[args.length_column], 1)
            p = np.poly1d(z)
            corr, p_val = stats.pearsonr(group_data[args.age_column], group_data[args.length_column])
            xs = np.sort(group_data[args.age_column].to_numpy())
            ax2.plot(xs, p(xs), color=palette[group], linestyle='--', alpha=0.8, label=f'{group}: r={corr:.3f}, p={p_val:.3f}')
    ax2.set_xlabel('Age')
    ax2.set_ylabel('Telomere Length')
    ax2.set_title('Age vs Telomere Length')
    ax2.legend(fontsize=8, loc='best')

    # Plot 3: BMI vs Telomere Length with group-wise trends
    ax3 = plt.subplot(3, 4, 3)
    for group in groups:
        group_data = df[df[args.group_column] == group]
        ax3.scatter(group_data[args.bmi_column], group_data[args.length_column], c=[palette[group]], alpha=0.7, s=60, label=group)
        if len(group_data) > 2:
            z_bmi = np.polyfit(group_data[args.bmi_column], group_data[args.length_column], 1)
            p_bmi = np.poly1d(z_bmi)
            corr_bmi, p_val_bmi = stats.pearsonr(group_data[args.bmi_column], group_data[args.length_column])
            xs = np.sort(group_data[args.bmi_column].to_numpy())
            ax3.plot(xs, p_bmi(xs), color=palette[group], linestyle='--', alpha=0.8, label=f'{group}: r={corr_bmi:.3f}, p={p_val_bmi:.3f}')
    ax3.set_xlabel('BMI')
    ax3.set_ylabel('Telomere Length')
    ax3.set_title('BMI vs Telomere Length')
    ax3.legend(fontsize=8, loc='best')

    # Plot 4: Gender comparison (age-adjusted)
    ax4 = plt.subplot(3, 4, 4)
    age_reg = LinearRegression()
    age_reg.fit(df[[args.age_column]], df[args.length_column])
    df['age_adjusted_telomere'] = df[args.length_column] - age_reg.predict(df[[args.age_column]])
    sns.boxplot(data=df, x=args.gender_column, y='age_adjusted_telomere', hue=args.group_column, ax=ax4)
    ax4.set_title('Age-adjusted Telomere by Gender')
    ax4.set_xlabel('Gender')
    ax4.set_ylabel('Age-adjusted Telomere')

    # Plot 5: Age distribution
    ax5 = plt.subplot(3, 4, 5)
    sns.histplot(data=df, x=args.age_column, hue=args.group_column, multiple='stack', bins=8, ax=ax5)
    ax5.set_title('Age Distribution')

    # Plot 6: BMI distribution
    ax6 = plt.subplot(3, 4, 6)
    sns.histplot(data=df, x=args.bmi_column, hue=args.group_column, multiple='stack', bins=8, ax=ax6)
    ax6.set_title('BMI Distribution')

    # Plot 7: Individual samples
    ax7 = plt.subplot(3, 4, 7)
    ordered = df.sort_values([args.group_column, args.length_column]).reset_index(drop=True)
    colors_ind = [palette[g] for g in ordered[args.group_column]]
    ax7.bar(range(len(ordered)), ordered[args.length_column], color=colors_ind, alpha=0.7)
    ax7.set_title('Individual Sample Results')
    ax7.set_xlabel('Sample')
    ax7.set_ylabel('Telomere Length')
    ax7.tick_params(axis='x', labelbottom=False)

    # Plot 8: Correlation matrix
    ax8 = plt.subplot(3, 4, 8)
    df_corr = df.copy()
    df_corr['gender_numeric'] = pd.factorize(df_corr[args.gender_column])[0]
    df_corr['group_numeric'] = pd.factorize(df_corr[args.group_column])[0]
    corr_vars = [args.length_column, args.age_column, args.bmi_column, 'gender_numeric', 'group_numeric']
    sns.heatmap(df_corr[corr_vars].corr(), annot=True, cmap='coolwarm', center=0, ax=ax8)
    ax8.set_title('Correlation Matrix')

    # Plot 9: Age-adjusted comparison by group
    ax9 = plt.subplot(3, 4, 9)
    sns.boxplot(data=df, x=args.group_column, y='age_adjusted_telomere', ax=ax9)
    sns.stripplot(data=df, x=args.group_column, y='age_adjusted_telomere', color='black', alpha=0.6, ax=ax9)
    ax9.set_title('Age-Adjusted Telomere Length')

    # Plot 10: Detailed scatter with size=BMI
    ax10 = plt.subplot(3, 4, 10)
    for group in groups:
        group_data = df[df[args.group_column] == group]
        ax10.scatter(group_data[args.age_column], group_data[args.length_column], c=[palette[group]], s=group_data[args.bmi_column] * 3, alpha=0.6, label=group)
    ax10.set_xlabel('Age')
    ax10.set_ylabel('Telomere Length')
    ax10.set_title('Age vs Telomere (size=BMI)')
    ax10.legend(fontsize=8)

    # Plot 11: Mean characteristics
    ax11 = plt.subplot(3, 4, 11)
    sample_chars = df.groupby(args.group_column).agg(
        age_mean=(args.age_column, 'mean'),
        bmi_mean=(args.bmi_column, 'mean'),
        tel_mean=(args.length_column, 'mean'),
    ).reset_index()
    x = np.arange(len(sample_chars))
    width = 0.25
    ax11.bar(x - width, sample_chars['age_mean'] / 10, width, label='Age/10', alpha=0.7)
    ax11.bar(x, sample_chars['bmi_mean'] / 10, width, label='BMI/10', alpha=0.7)
    ax11.bar(x + width, sample_chars['tel_mean'], width, label='Telomere Length', alpha=0.7)
    ax11.set_xlabel('Group')
    ax11.set_ylabel('Value')
    ax11.set_title('Mean Characteristics by Group')
    ax11.set_xticks(x)
    ax11.set_xticklabels(sample_chars[args.group_column])
    ax11.legend(fontsize=8)

    # Plot 12: Statistical summary
    ax12 = plt.subplot(3, 4, 12)
    ax12.axis('off')
    overall_corr, overall_p = stats.pearsonr(df[args.age_column], df[args.length_column]) if len(df) > 2 else (np.nan, np.nan)
    summary_lines = [
        'Statistical Summary:',
        '',
        f'Samples: {len(df)}',
        f'Groups: {", ".join(groups)}',
        f'Overall age-length r={overall_corr:.3f}, p={overall_p:.4f}',
        '',
    ]
    for group in groups:
        sub = df[df[args.group_column] == group][args.length_column]
        summary_lines.append(f'{group}: {sub.mean():.2f} ± {sub.std():.2f} (n={len(sub)})')
    ax12.text(0.1, 0.9, "\n".join(summary_lines), transform=ax12.transAxes, fontsize=10, verticalalignment='top', fontfamily='monospace')

    if output_summary_csv is not None:
        summary_rows = []
        for group in groups:
            sub = df[df[args.group_column] == group]
            summary_rows.append(
                {
                    "group": group,
                    "n": int(len(sub)),
                    "telomere_mean": float(sub[args.length_column].mean()),
                    "telomere_sd": float(sub[args.length_column].std()),
                    "age_mean": float(sub[args.age_column].mean()),
                    "bmi_mean": float(sub[args.bmi_column].mean()),
                }
            )
        pd.DataFrame(summary_rows).to_csv(output_summary_csv, index=False)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    if args.no_show:
        plt.close()
    else:
        plt.show()

    print('Comprehensive visualization completed!')
    print(f'Saved as: {output_png}')
    if output_summary_csv is not None:
        print(f'Summary CSV: {output_summary_csv}')


if __name__ == "__main__":
    main()
