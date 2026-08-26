#!/usr/bin/env python3
"""
PCA and boxplot analysis for MethAtlas deconvolution output.

Expected input format (correct CSV):
- Rows: tissues
- Columns: sample IDs (e.g., con_barcode11, cor_barcode01, rec_barcode06)

Example Windows path provided by user:
C:\\Users\\gissu\\meth_atlas\\methatlas_combined_450k_BRCA2_deconv_output.csv

The script accepts either Windows or WSL paths and normalizes automatically.
"""

import argparse
import os
import re
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Circle
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


# ============================================================================
# USER SETTINGS (edit these for new projects)
# ============================================================================
DEFAULT_INPUT_CSV = (
    r"/your/meth/atlas/output_csv"
)
DEFAULT_OUTPUT_DIR = r"/your/output/directory"

# Optional metadata file (CSV) with sample->group mapping.
# Leave as None to infer groups from sample names using GROUP_RULES.
DEFAULT_SAMPLE_METADATA_CSV = None

# Metadata column names used when --sample-metadata-csv is provided.
DEFAULT_METADATA_SAMPLE_COL = "sample"
DEFAULT_METADATA_GROUP_COL = "group"

# If False, rows assigned as "unknown" are dropped before analysis.
KEEP_UNKNOWN_GROUPS = False
DEFAULT_N_PERMUTATIONS = 100000

GROUP_COLOR_MAP = {
    "BRCA2": "#d62728",
    "control": "#1f77b4",
    "unknown": "#7f7f7f",
}
GROUP_ORDER = ["BRCA2", "control", "unknown"]

# Display labels for groups in plots/tables.
GROUP_DISPLAY_NAMES = {
    "BRCA2": "BRCA2",
}

# Sample-name regex rules used when metadata mapping is not provided.
# EDIT THESE RULES for your naming convention.
GROUP_RULES = [
    {"label": "BRCA2", "pattern": r"^brca2-*barcode*"},
    {"label": "control", "pattern": r"^ctrl-*barcode*"},
]

SELECTED_TISSUES = [
    "Neutrophils",
    "Erythrocyte progenitors",
    "Hepatocytes",
    "Adipocytes",
    "Monocytes",
    "Left atrium",
    "Vascular endothelial cells",
]
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


def ensure_output_directory(path: str):
    """Create output directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)


def infer_group_from_sample_name(sample_name: str) -> str:
    """Infer biological group from sample naming convention using GROUP_RULES."""
    normalized = os.path.splitext(os.path.basename(str(sample_name)))[0].lower()

    for rule in GROUP_RULES:
        if re.search(rule["pattern"], normalized):
            return str(rule["label"]).strip().lower()

    return "unknown"


def format_group_label(group: str) -> str:
    """Human-readable group labels."""
    group = str(group).lower()
    if group in GROUP_DISPLAY_NAMES:
        return GROUP_DISPLAY_NAMES[group]
    return group.capitalize()


def get_ordered_groups(groups):
    """Return groups in preferred display order, then any extras alphabetically."""
    observed = set(groups)
    ordered = [group for group in GROUP_ORDER if group in observed]
    remaining = sorted([group for group in observed if group not in GROUP_ORDER])
    return ordered + remaining


def build_group_color_map(groups):
    """Use configured colors and auto-assign colors for unseen groups."""
    color_map = dict(GROUP_COLOR_MAP)
    extras = [g for g in get_ordered_groups(groups) if g not in color_map]
    if extras:
        palette = plt.cm.tab20(np.linspace(0, 1, max(3, len(extras))))
        for idx, group in enumerate(extras):
            rgb = palette[idx % len(palette)][:3]
            color_map[group] = "#%02x%02x%02x" % tuple(int(255 * c) for c in rgb)
    return color_map


def load_group_map(metadata_csv: str, sample_col: str, group_col: str):
    """Load optional sample->group metadata mapping."""
    if metadata_csv is None:
        return None
    if not os.path.exists(metadata_csv):
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_csv}")

    meta_df = pd.read_csv(metadata_csv)
    if sample_col not in meta_df.columns or group_col not in meta_df.columns:
        raise ValueError(
            f"Metadata CSV must contain columns '{sample_col}' and '{group_col}'"
        )

    mapping = {}
    for _, row in meta_df[[sample_col, group_col]].dropna().iterrows():
        sample_key = str(row[sample_col]).strip()
        mapping[sample_key] = str(row[group_col]).strip().lower()
    return mapping


def assign_groups(samples: pd.Series, group_map=None):
    """Assign groups from metadata map if available; otherwise infer from names."""
    if group_map is not None:
        assigned = samples.apply(lambda s: group_map.get(str(s).strip(), "unknown"))
        if assigned.nunique() > 1:
            return assigned.astype(str), "metadata"
    return samples.apply(infer_group_from_sample_name).astype(str), "inferred"


def p_to_stars(p_value: float) -> str:
    """Map p-values to significance stars."""
    if pd.isna(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def format_p_value(p_value: float) -> str:
    """Readable p-value text."""
    if pd.isna(p_value):
        return "nan"
    p_value = float(p_value)
    if p_value < 1e-5:
        return f"{p_value:.2e}"
    return f"{p_value:.5f}".rstrip("0").rstrip(".")


def draw_group_circle(ax, x_values, y_values, color):
    """Draw an enclosing circle around a group on a PCA scatter plot."""
    if len(x_values) == 0:
        return

    center_x = float(np.mean(x_values))
    center_y = float(np.mean(y_values))
    distances = np.sqrt((x_values - center_x) ** 2 + (y_values - center_y) ** 2)

    base_radius = float(np.max(distances)) if len(distances) > 0 else 0.0
    padding = max(0.12, base_radius * 0.15)
    radius = base_radius + padding

    circle = Circle(
        (center_x, center_y),
        radius,
        fill=False,
        edgecolor=color,
        linewidth=1.8,
        linestyle="--",
        alpha=0.85,
    )
    ax.add_patch(circle)


def filter_iqr_outliers(values: np.ndarray, iqr_scale: float = 1.5) -> np.ndarray:
    """Return values with IQR-based outliers removed (1.5×IQR rule)."""
    values = np.asarray(values, dtype=float)
    if len(values) < 4:
        return values
    q1, q3 = np.nanpercentile(values, [25, 75])
    iqr = q3 - q1
    lower = q1 - iqr_scale * iqr
    upper = q3 + iqr_scale * iqr
    return values[(values >= lower) & (values <= upper)]


def permutation_test_p_value(
    values_a: np.ndarray,
    values_b: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Two-sided permutation test using absolute difference in group means."""
    values_a = np.asarray(values_a, dtype=float)
    values_b = np.asarray(values_b, dtype=float)

    if len(values_a) == 0 or len(values_b) == 0:
        return np.nan, np.nan

    observed = float(abs(np.mean(values_a) - np.mean(values_b)))
    combined = np.concatenate([values_a, values_b])
    n_a = len(values_a)
    exceed_count = 0

    for _ in range(max(1, int(n_permutations))):
        shuffled = rng.permutation(combined)
        permuted = float(abs(np.mean(shuffled[:n_a]) - np.mean(shuffled[n_a:])))
        if permuted >= observed - 1e-15:
            exceed_count += 1

    p_value = (exceed_count + 1) / (max(1, int(n_permutations)) + 1)
    return observed, float(p_value)


def load_input_tables(input_csv: str, group_map=None, keep_unknown: bool = False):
    """Load raw table and convert to long + sample-feature matrix forms."""
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    raw_df = pd.read_csv(input_csv)
    if raw_df.shape[1] < 3:
        raise ValueError("Input CSV must contain at least one tissue column and two sample columns")

    tissue_col = raw_df.columns[0]
    raw_df = raw_df.rename(columns={tissue_col: "tissue"}).copy()
    raw_df["tissue"] = raw_df["tissue"].astype(str)

    sample_cols = [c for c in raw_df.columns if c != "tissue"]
    for col in sample_cols:
        raw_df[col] = pd.to_numeric(raw_df[col], errors="coerce")

    long_df = raw_df.melt(
        id_vars=["tissue"],
        value_vars=sample_cols,
        var_name="sample",
        value_name="fraction",
    )
    long_df["group"], group_source = assign_groups(long_df["sample"], group_map=group_map)
    long_df = long_df.dropna(subset=["fraction"]).copy()
    if not keep_unknown:
        long_df = long_df[long_df["group"] != "unknown"].copy()

    sample_feature_df = long_df.pivot_table(
        index=["sample", "group"],
        columns="tissue",
        values="fraction",
        aggfunc="mean",
        fill_value=0.0,
    )

    if sample_feature_df.empty:
        raise ValueError("No valid sample/group rows after parsing input CSV")

    sample_feature_df = sample_feature_df.sort_index()
    return raw_df, long_df, sample_feature_df, group_source


def compute_pairwise_permutation_tests(
    values: np.ndarray,
    groups: np.ndarray,
    metric_name: str,
    analysis_groups,
    n_permutations: int,
    rng: np.random.Generator,
):
    """Compute pairwise permutation tests across requested groups."""
    records = []
    for g1, g2 in combinations(analysis_groups, 2):
        x = np.asarray(values[groups == g1], dtype=float)
        y = np.asarray(values[groups == g2], dtype=float)
        if len(x) == 0 or len(y) == 0:
            stat, p_val = np.nan, np.nan
        else:
            stat, p_val = permutation_test_p_value(x, y, n_permutations=n_permutations, rng=rng)

        records.append(
            {
                "metric": metric_name,
                "group1": g1,
                "group2": g2,
                "n1": int(len(x)),
                "n2": int(len(y)),
                "statistic": float(stat) if not pd.isna(stat) else np.nan,
                "p_value": float(p_val) if not pd.isna(p_val) else np.nan,
                "stars": p_to_stars(p_val),
            }
        )
    return records


def perform_pca(sample_feature_df: pd.DataFrame, analysis_groups, n_permutations: int):
    """Run PCA on sample-by-tissue deconvolution matrix."""
    sample_names = np.array(sample_feature_df.index.get_level_values("sample"), dtype=object)
    groups = np.array(sample_feature_df.index.get_level_values("group"), dtype=object)
    tissue_cols = list(sample_feature_df.columns)
    X = sample_feature_df.to_numpy(dtype=float)
    sample_signal = X.sum(axis=1)

    if X.shape[0] < 4:
        raise ValueError(f"Need at least 4 samples for PCA, found {X.shape[0]}")

    column_std = np.std(X, axis=0)
    valid_columns = column_std > 1e-10
    if np.sum(valid_columns) < 2:
        raise ValueError("Fewer than 2 tissue features with variance")

    X_filtered = X[:, valid_columns]
    tissue_cols_filtered = [tissue_cols[i] for i in range(len(tissue_cols)) if valid_columns[i]]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_filtered)

    n_components = min(X_scaled.shape[0] - 1, X_scaled.shape[1], 10)
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_scaled)

    explained_var = pca.explained_variance_ratio_
    cumulative_var = np.cumsum(explained_var)

    pc1_sig_corr, pc1_sig_p = pearsonr(X_pca[:, 0], sample_signal)
    pc2_sig_corr, pc2_sig_p = pearsonr(X_pca[:, 1], sample_signal)

    if len(np.unique(groups)) > 1 and X_pca.shape[1] >= 2:
        sil_score = silhouette_score(X_pca[:, :2], groups)
    else:
        sil_score = None

    feature_contributions = {}
    for pc_idx in range(min(2, len(pca.components_))):
        pc_loadings = pca.components_[pc_idx]
        top_indices = np.argsort(np.abs(pc_loadings))[-8:][::-1]

        feature_contributions[f"PC{pc_idx + 1}"] = []
        for idx in top_indices:
            feature_contributions[f"PC{pc_idx + 1}"].append(
                {
                    "feature": tissue_cols_filtered[idx],
                    "loading": float(pc_loadings[idx]),
                    "abs_loading": float(abs(pc_loadings[idx])),
                }
            )

    pc_pairwise_stats = []
    rng = np.random.default_rng(42)
    for pc_idx in range(min(2, X_pca.shape[1])):
        pc_pairwise_stats.extend(
            compute_pairwise_permutation_tests(
                X_pca[:, pc_idx],
                groups,
                metric_name=f"PC{pc_idx + 1}",
                analysis_groups=analysis_groups,
                n_permutations=n_permutations,
                rng=rng,
            )
        )

    return {
        "pca": pca,
        "X_pca": X_pca,
        "X_scaled": X_scaled,
        "explained_variance": explained_var,
        "cumulative_variance": cumulative_var,
        "groups": groups,
        "sample_names": sample_names,
        "tissue_cols": tissue_cols_filtered,
        "feature_contributions": feature_contributions,
        "silhouette_score": sil_score,
        "sample_signal": sample_signal,
        "pc1_signal_correlation": float(pc1_sig_corr),
        "pc2_signal_correlation": float(pc2_sig_corr),
        "pc1_signal_p_value": float(pc1_sig_p),
        "pc2_signal_p_value": float(pc2_sig_p),
        "pc_pairwise_stats": pd.DataFrame(pc_pairwise_stats),
        "analysis_groups": list(analysis_groups),
        "n_permutations": int(n_permutations),
    }


def add_significance_bars(ax, p_rows: pd.DataFrame, y_values: np.ndarray, analysis_groups):
    """Draw pairwise MWU brackets labeled with p-values, matching cfDNA style."""
    entries = []
    for _, row in p_rows.iterrows():
        p_value = row["p_value"]
        if pd.isna(p_value):
            label = "ns (p=nan)"
        elif float(p_value) < 0.05:
            stars = p_to_stars(float(p_value)).replace("ns", "")
            p_txt = format_p_value(float(p_value))
            label = f"{stars} (p={p_txt})" if stars else f"p={p_txt}"
        else:
            label = f"ns (p={format_p_value(float(p_value))})"
        entries.append((row["group1"], row["group2"], label))

    if not entries:
        return

    finite_vals = np.asarray(y_values, dtype=float)
    finite_vals = finite_vals[np.isfinite(finite_vals)]
    if finite_vals.size == 0:
        return

    y_min = float(np.min(finite_vals))
    y_max = float(np.max(finite_vals))
    value_range = y_max - y_min
    if value_range <= 0:
        value_range = max(abs(y_max), 1.0) * 0.1

    base_y = y_max + 0.06 * value_range
    step = 0.11 * value_range
    bracket_h = 0.28 * step
    text_pad = 0.016 * value_range

    group_x = {group: idx for idx, group in enumerate(analysis_groups)}
    top_y = base_y

    for level, (g1, g2, label) in enumerate(entries):
        if g1 not in group_x or g2 not in group_x:
            continue

        x1 = group_x[g1]
        x2 = group_x[g2]
        if x1 > x2:
            x1, x2 = x2, x1

        y = base_y + level * step
        ax.plot([x1, x1, x2, x2], [y, y + bracket_h, y + bracket_h, y], lw=1.0, c="#333333")
        ax.text(
            (x1 + x2) / 2,
            y + bracket_h + text_pad,
            label,
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            color="#111111",
        )
        top_y = max(top_y, y + bracket_h + 2 * text_pad)

    current_top = ax.get_ylim()[1]
    ax.set_ylim(top=max(current_top, top_y + 0.05 * value_range))


def style_boxplot_axis(
    ax,
    analysis_groups,
    xlabel: str,
    ylabel: str,
    title: str,
    rotation: int = 0,
    group_counts=None,
):
    """Apply cfDNA-style bold titles, labels, and ticks to boxplot axes."""
    ax.set_xticks(range(len(analysis_groups)))
    labels = []
    for group in analysis_groups:
        label = format_group_label(group)
        if group_counts is not None:
            label = f"{label}\n(n={int(group_counts.get(group, 0))})"
        labels.append(label)
    ax.set_xticklabels(labels, rotation=rotation)
    ax.set_xlabel(xlabel, fontsize=13, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")


def create_comprehensive_plot(results, output_dir: str):
    """Create 3x3 PCA panel mirroring the existing detailed cfDNA layout."""
    sns.set_theme(style="whitegrid")

    fig = plt.figure(figsize=(20, 16))
    unique_groups = get_ordered_groups(results["groups"])
    analysis_groups = results["analysis_groups"]
    color_map = build_group_color_map(unique_groups)

    # 1. PCA Scatter Plot (PC1 vs PC2)
    ax1 = plt.subplot(3, 3, 1)
    for group in unique_groups:
        mask = results["groups"] == group
        x_group = results["X_pca"][mask, 0]
        y_group = results["X_pca"][mask, 1]
        color = color_map.get(group, "#9467bd")

        plt.scatter(
            x_group,
            y_group,
            c=color,
            label=f"{format_group_label(group)} (n={np.sum(mask)})",
            alpha=0.72,
            s=60,
        )
        draw_group_circle(ax1, x_group, y_group, color)

    plt.xlabel(f"PC1 ({results['explained_variance'][0]:.1%})")
    plt.ylabel(f"PC2 ({results['explained_variance'][1]:.1%})")
    plt.title("MethAtlas Deconvolution PCA - PC1 vs PC2")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # 2. Scree Plot
    ax2 = plt.subplot(3, 3, 2)
    pc_numbers = range(1, len(results["explained_variance"]) + 1)
    plt.plot(pc_numbers, results["explained_variance"] * 100, "bo-", linewidth=2, markersize=6)
    plt.xlabel("Principal Component")
    plt.ylabel("Explained Variance (%)")
    plt.title("Scree Plot")
    plt.xticks(list(pc_numbers))
    plt.grid(True, alpha=0.3)

    ax2_twin = ax2.twinx()
    ax2_twin.plot(
        range(1, len(results["cumulative_variance"]) + 1),
        results["cumulative_variance"] * 100,
        marker="o",
        linestyle="-",
        color="red",
        linewidth=2,
    )
    ax2_twin.set_ylabel("Cumulative Variance (%)", color="red")
    ax2_twin.tick_params(axis="y", labelcolor="red")

    # 3. PC1 Box Plot
    ax3 = plt.subplot(3, 3, 3)
    pc1_df = pd.DataFrame({"group": results["groups"], "value": results["X_pca"][:, 0]})
    pc1_counts = pc1_df["group"].value_counts().to_dict()
    sns.boxplot(
        data=pc1_df,
        x="group",
        y="value",
        hue="group",
        order=analysis_groups,
        hue_order=analysis_groups,
        palette=color_map,
        ax=ax3,
        showfliers=True,
        flierprops={"marker": "o", "markersize": 4.5, "markerfacecolor": "white", "markeredgecolor": "#222222", "alpha": 1.0},
        dodge=False,
        legend=False,
    )
    sns.stripplot(data=pc1_df, x="group", y="value", order=analysis_groups, color="black", alpha=0.5, size=3, jitter=0.15, ax=ax3)
    style_boxplot_axis(ax3, analysis_groups, xlabel="Group", ylabel="PC1 Score", title="PC1 Distribution", rotation=10, group_counts=pc1_counts)
    pc1_stats = results["pc_pairwise_stats"][results["pc_pairwise_stats"]["metric"] == "PC1"]
    add_significance_bars(ax3, pc1_stats, pc1_df["value"].to_numpy(dtype=float), analysis_groups)

    # 4. PC2 Box Plot
    ax4 = plt.subplot(3, 3, 4)
    pc2_df = pd.DataFrame({"group": results["groups"], "value": results["X_pca"][:, 1]})
    pc2_counts = pc2_df["group"].value_counts().to_dict()
    sns.boxplot(
        data=pc2_df,
        x="group",
        y="value",
        hue="group",
        order=analysis_groups,
        hue_order=analysis_groups,
        palette=color_map,
        ax=ax4,
        showfliers=True,
        flierprops={"marker": "o", "markersize": 4.5, "markerfacecolor": "white", "markeredgecolor": "#222222", "alpha": 1.0},
        dodge=False,
        legend=False,
    )
    sns.stripplot(data=pc2_df, x="group", y="value", order=analysis_groups, color="black", alpha=0.5, size=3, jitter=0.15, ax=ax4)
    style_boxplot_axis(ax4, analysis_groups, xlabel="Group", ylabel="PC2 Score", title="PC2 Distribution", rotation=10, group_counts=pc2_counts)
    pc2_stats = results["pc_pairwise_stats"][results["pc_pairwise_stats"]["metric"] == "PC2"]
    add_significance_bars(ax4, pc2_stats, pc2_df["value"].to_numpy(dtype=float), analysis_groups)

    # 5. PC1 vs PC2 Feature Loadings
    ax5 = plt.subplot(3, 3, 5)
    pc1_loadings = results["pca"].components_[0, :]
    pc2_loadings = results["pca"].components_[1, :]
    tissue_codes = np.arange(len(results["tissue_cols"]))

    scatter = ax5.scatter(
        pc1_loadings,
        pc2_loadings,
        c=tissue_codes,
        cmap="coolwarm",
        s=55,
        alpha=0.8,
        edgecolors="black",
        linewidth=0.4,
    )
    cbar = plt.colorbar(scatter, ax=ax5, shrink=0.8)
    cbar.set_label("Tissue Index", rotation=270, labelpad=18)

    plt.xlabel("PC1 Loading")
    plt.ylabel("PC2 Loading")
    plt.title("PC Loadings (Tissue Features)")
    plt.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    plt.axvline(x=0, color="black", linestyle="-", alpha=0.3)
    plt.grid(True, alpha=0.3)

    # 6. PC1 vs PC3 (if available)
    if results["X_pca"].shape[1] > 2:
        ax6 = plt.subplot(3, 3, 6)
        for group in unique_groups:
            mask = results["groups"] == group
            x_group = results["X_pca"][mask, 0]
            y_group = results["X_pca"][mask, 2]
            color = color_map.get(group, "#9467bd")
            plt.scatter(
                x_group,
                y_group,
                c=color,
                label=format_group_label(group),
                alpha=0.72,
                s=60,
            )
            draw_group_circle(ax6, x_group, y_group, color)

        plt.xlabel(f"PC1 ({results['explained_variance'][0]:.1%})")
        plt.ylabel(f"PC3 ({results['explained_variance'][2]:.1%})")
        plt.title("PC1 vs PC3")
        plt.legend()
        plt.grid(True, alpha=0.3)

    # 7. Statistics Summary
    ax7 = plt.subplot(3, 3, 7)
    ax7.axis("off")

    stats_text = "MethAtlas Deconvolution PCA Summary:\n\n"
    stats_text += f"Samples: {results['X_pca'].shape[0]}\n"
    stats_text += f"Tissue Features: {len(results['tissue_cols'])}\n"
    stats_text += f"PC1 Variance: {results['explained_variance'][0]:.1%}\n"
    stats_text += f"PC2 Variance: {results['explained_variance'][1]:.1%}\n"
    stats_text += f"PC1+PC2 Total: {results['cumulative_variance'][1]:.1%}\n\n"

    if results["silhouette_score"] is not None:
        stats_text += f"Silhouette Score: {results['silhouette_score']:.3f}\n\n"

    stats_text += "Signal Correlations:\n"
    stats_text += (
        f"PC1-Signal: r={results['pc1_signal_correlation']:.3f} "
        f"{p_to_stars(results['pc1_signal_p_value'])}\n"
    )
    stats_text += (
        f"PC2-Signal: r={results['pc2_signal_correlation']:.3f} "
        f"{p_to_stars(results['pc2_signal_p_value'])}\n\n"
    )

    stats_text += f"PC Permutation ({results['n_permutations']:,} permutations):\n"
    for metric in ["PC1", "PC2"]:
        sub = results["pc_pairwise_stats"][results["pc_pairwise_stats"]["metric"] == metric]
        for _, row in sub.iterrows():
            stats_text += (
                f"{metric} {format_group_label(row['group1'])} vs "
                f"{format_group_label(row['group2'])}: "
                f"p={format_p_value(row['p_value'])} ({row['stars']})\n"
            )

    ax7.text(
        0.05,
        0.95,
        stats_text,
        transform=ax7.transAxes,
        fontsize=10,
        verticalalignment="top",
        fontfamily="monospace",
    )

    # 8. Feature Contributions Bar Plot
    ax8 = plt.subplot(3, 3, 8)
    if "PC1" in results["feature_contributions"]:
        contrib_data = results["feature_contributions"]["PC1"]
        features = [c["feature"] for c in contrib_data]
        loadings = [c["loading"] for c in contrib_data]

        plt.bar(
            range(len(features)),
            loadings,
            color=["red" if x < 0 else "blue" for x in loadings],
            alpha=0.75,
        )
        plt.xlabel("Top Tissue Features")
        plt.ylabel("PC1 Loading")
        plt.title("PC1 Top Contributors")
        plt.xticks(range(len(features)), features, rotation=50, ha="right", fontsize=8)
        plt.axhline(y=0, color="black", linestyle="-", alpha=0.3)

    # 9. PCA Heatmap-style panel by total deconvolution signal
    ax9 = plt.subplot(3, 3, 9)
    scatter = plt.scatter(
        results["X_pca"][:, 0],
        results["X_pca"][:, 1],
        c=results["sample_signal"],
        cmap="plasma",
        s=80,
        alpha=0.8,
        edgecolors="black",
        linewidth=0.5,
    )
    cbar = plt.colorbar(scatter, ax=ax9, shrink=0.8)
    cbar.set_label("Sample Total Deconvolution Signal", rotation=270, labelpad=20)

    plt.xlabel(f"PC1 ({results['explained_variance'][0]:.1%})")
    plt.ylabel(f"PC2 ({results['explained_variance'][1]:.1%})")
    plt.title("PCA Colored by Sample Signal")
    plt.grid(True, alpha=0.3)

    corr_text = (
        f"PC1-Signal: r={results['pc1_signal_correlation']:.3f}\n"
        f"PC2-Signal: r={results['pc2_signal_correlation']:.3f}"
    )
    plt.text(
        0.02,
        0.98,
        corr_text,
        transform=ax9.transAxes,
        verticalalignment="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "black"},
    )

    plt.tight_layout(pad=2.0)

    ensure_output_directory(output_dir)
    plot_file = os.path.join(output_dir, "methatlas_deconv_detailed_PCA_analysis.png")
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved detailed PCA plot: {plot_file}")


def create_tissue_boxplots(
    long_df: pd.DataFrame,
    output_dir: str,
    analysis_groups,
    top_n_tissues: int = 12,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    remove_iqr_outliers: bool = True,
):
    """Create tissue-level group boxplots with permutation p-values."""
    sns.set_theme(style="whitegrid")

    available_tissues = long_df["tissue"].unique()
    selected_tissues = []
    for target in SELECTED_TISSUES:
        target_lower = target.lower()
        matched = next(
            (t for t in available_tissues if target_lower in t.lower()),
            None,
        )
        if matched is not None:
            selected_tissues.append(matched)
        else:
            print(f"Warning: tissue '{target}' not found in data (skipped)")
    if not selected_tissues:
        print("No selected tissues found in data — falling back to top-variance selection")
        tissue_var = (
            long_df.groupby("tissue", as_index=False)["fraction"]
            .var()
            .sort_values("fraction", ascending=False)
        )
        selected_tissues = tissue_var["tissue"].head(top_n_tissues).tolist()
    plot_df = long_df[long_df["tissue"].isin(selected_tissues)].copy()
    color_map = build_group_color_map(analysis_groups)

    permutation_records = []
    rng = np.random.default_rng(42)
    for tissue in selected_tissues:
        sub = plot_df[plot_df["tissue"] == tissue]
        for g1, g2 in combinations(analysis_groups, 2):
            x = sub.loc[sub["group"] == g1, "fraction"].to_numpy(dtype=float)
            y = sub.loc[sub["group"] == g2, "fraction"].to_numpy(dtype=float)
            if len(x) == 0 or len(y) == 0:
                stat, p_val = np.nan, np.nan
            else:
                stat, p_val = permutation_test_p_value(x, y, n_permutations=n_permutations, rng=rng)

            permutation_records.append(
                {
                    "tissue": tissue,
                    "group1": g1,
                    "group2": g2,
                    "n1": int(len(x)),
                    "n2": int(len(y)),
                    "statistic": float(stat) if not pd.isna(stat) else np.nan,
                    "p_value": float(p_val) if not pd.isna(p_val) else np.nan,
                    "stars": p_to_stars(p_val),
                }
            )

    permutation_df = pd.DataFrame(permutation_records)

    n_panels = len(selected_tissues)
    n_cols = 3
    n_rows = int(np.ceil(n_panels / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5.5 * n_rows), squeeze=False)

    for idx, tissue in enumerate(selected_tissues):
        ax = axes[idx // n_cols, idx % n_cols]
        sub_raw = plot_df[plot_df["tissue"] == tissue].copy()

        if remove_iqr_outliers:
            filtered_rows = []
            for g in analysis_groups:
                g_df = sub_raw[sub_raw["group"] == g].copy()
                vals = g_df["fraction"].to_numpy(dtype=float)
                if len(vals) >= 4:
                    q1, q3 = np.nanpercentile(vals, [25, 75])
                    iqr = q3 - q1
                    keep = (vals >= q1 - 1.5 * iqr) & (vals <= q3 + 1.5 * iqr)
                    filtered_rows.append(g_df[keep])
                else:
                    filtered_rows.append(g_df)
            sub = pd.concat(filtered_rows, ignore_index=True)
        else:
            sub = sub_raw.copy()

        group_counts = sub.groupby("group")["sample"].nunique().to_dict()

        sns.boxplot(
            data=sub,
            x="group",
            y="fraction",
            hue="group",
            order=analysis_groups,
            hue_order=analysis_groups,
            palette=color_map,
            showfliers=False,
            dodge=False,
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x="group",
            y="fraction",
            order=analysis_groups,
            color="black",
            alpha=0.55,
            size=3,
            jitter=0.15,
            ax=ax,
        )

        style_boxplot_axis(
            ax,
            analysis_groups,
            xlabel="Group",
            ylabel="Deconvolution Fraction",
            title=tissue,
            rotation=10,
            group_counts=group_counts,
        )

        p_rows = permutation_df[permutation_df["tissue"] == tissue]
        add_significance_bars(ax, p_rows, sub["fraction"].to_numpy(dtype=float), analysis_groups)

    # Hide unused axes if top_n_tissues is not divisible by n_cols.
    total_axes = n_rows * n_cols
    for idx in range(n_panels, total_axes):
        axes[idx // n_cols, idx % n_cols].axis("off")

    fig.suptitle(
        "MethAtlas Tissue Fractions by Group (Selected Tissues)\n"
        + (
            f"Permutation test p-values ({n_permutations:,} permutations), display excludes IQR outliers"
            if remove_iqr_outliers
            else f"Permutation test p-values ({n_permutations:,} permutations), raw values shown without IQR filtering"
        ),
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    ensure_output_directory(output_dir)
    plot_file = os.path.join(
        output_dir,
        "methatlas_tissue_group_boxplots_permutation_iqr_filtered.png"
        if remove_iqr_outliers
        else "methatlas_tissue_group_boxplots_permutation_raw.png",
    )
    fig.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    stats_file = os.path.join(output_dir, "methatlas_tissue_permutation_results.csv")
    permutation_df.to_csv(stats_file, index=False)

    print(f"Saved tissue boxplots: {plot_file}")
    print(f"Saved tissue permutation table: {stats_file}")


def save_summary_and_tables(results, sample_feature_df: pd.DataFrame, output_dir: str):
    """Save summary text and PCA tables."""
    ensure_output_directory(output_dir)

    summary_file = os.path.join(output_dir, "methatlas_deconv_PCA_summary.txt")
    scores_file = os.path.join(output_dir, "methatlas_deconv_PCA_scores.csv")
    loadings_file = os.path.join(output_dir, "methatlas_deconv_PCA_loadings.csv")
    pc_stats_file = os.path.join(output_dir, "methatlas_deconv_PC_permutation_results.csv")

    scores_df = pd.DataFrame(
        {
            "sample": results["sample_names"],
            "group": results["groups"],
            "sample_signal": results["sample_signal"],
        }
    )
    for i in range(results["X_pca"].shape[1]):
        scores_df[f"PC{i + 1}"] = results["X_pca"][:, i]
    scores_df.to_csv(scores_file, index=False)

    loadings_df = pd.DataFrame(
        {
            "tissue": results["tissue_cols"],
            "PC1_loading": results["pca"].components_[0, :],
            "PC2_loading": results["pca"].components_[1, :],
        }
    )
    if results["pca"].components_.shape[0] > 2:
        loadings_df["PC3_loading"] = results["pca"].components_[2, :]
    loadings_df.to_csv(loadings_file, index=False)

    results["pc_pairwise_stats"].to_csv(pc_stats_file, index=False)

    with open(summary_file, "w") as f:
        f.write("=== MethAtlas Deconvolution PCA Summary ===\n\n")
        f.write(f"Samples analyzed: {results['X_pca'].shape[0]}\n")
        f.write(f"Tissue features used: {len(results['tissue_cols'])}\n")
        f.write(f"PC1 explained variance: {results['explained_variance'][0]:.1%}\n")
        f.write(f"PC2 explained variance: {results['explained_variance'][1]:.1%}\n")
        f.write(f"PC1+PC2 combined: {results['cumulative_variance'][1]:.1%}\n")
        if results["silhouette_score"] is not None:
            f.write(f"Silhouette score: {results['silhouette_score']:.3f}\n")

        f.write("\nGroup counts:\n")
        group_counts = pd.Series(results["groups"]).value_counts()
        for group, count in group_counts.items():
            f.write(f"  {format_group_label(group)}: {count}\n")

        f.write("\nSignal correlations:\n")
        f.write(
            f"  PC1 vs sample signal: r={results['pc1_signal_correlation']:.3f}, "
            f"p={results['pc1_signal_p_value']:.4g} {p_to_stars(results['pc1_signal_p_value'])}\n"
        )
        f.write(
            f"  PC2 vs sample signal: r={results['pc2_signal_correlation']:.3f}, "
            f"p={results['pc2_signal_p_value']:.4g} {p_to_stars(results['pc2_signal_p_value'])}\n"
        )

        f.write(f"\nPC-level permutation tests ({results['n_permutations']:,} permutations):\n")
        for _, row in results["pc_pairwise_stats"].iterrows():
            f.write(
                f"  {row['metric']} {format_group_label(row['group1'])} vs "
                f"{format_group_label(row['group2'])}: "
                f"p={format_p_value(row['p_value'])} ({row['stars']})\n"
            )

        f.write("\nTop PC1 contributors:\n")
        for contrib in results["feature_contributions"].get("PC1", []):
            f.write(f"  {contrib['feature']}: {contrib['loading']:+.3f}\n")

        f.write("\nTop PC2 contributors:\n")
        for contrib in results["feature_contributions"].get("PC2", []):
            f.write(f"  {contrib['feature']}: {contrib['loading']:+.3f}\n")

    sample_feature_file = os.path.join(output_dir, "methatlas_sample_feature_matrix.csv")
    sample_feature_df.reset_index().to_csv(sample_feature_file, index=False)

    print(f"Saved summary: {summary_file}")
    print(f"Saved PCA scores: {scores_file}")
    print(f"Saved PCA loadings: {loadings_file}")
    print(f"Saved PC permutation table: {pc_stats_file}")
    print(f"Saved sample-feature matrix: {sample_feature_file}")


def parse_args():
    """CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run PCA on MethAtlas deconvolution output CSV and generate matching "
            "detailed PCA plots plus MWU boxplot summaries"
        ),
    )
    parser.add_argument(
        "--input-csv",
        default=DEFAULT_INPUT_CSV,
        help="Path to input CSV (Windows or WSL path accepted)",
    )
    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory (Windows or WSL path accepted)",
    )
    parser.add_argument(
        "--sample-metadata-csv",
        default=DEFAULT_SAMPLE_METADATA_CSV,
        help=(
            "Optional CSV with sample->group mapping. "
            "If omitted, groups are inferred from sample names using GROUP_RULES"
        ),
    )
    parser.add_argument(
        "--metadata-sample-col",
        default=DEFAULT_METADATA_SAMPLE_COL,
        help="Sample ID column in metadata CSV",
    )
    parser.add_argument(
        "--metadata-group-col",
        default=DEFAULT_METADATA_GROUP_COL,
        help="Group column in metadata CSV",
    )
    parser.add_argument(
        "--keep-unknown",
        action="store_true",
        default=KEEP_UNKNOWN_GROUPS,
        help="Keep samples assigned to unknown group (default is to drop)",
    )
    parser.add_argument(
        "--top-n-tissues",
        type=int,
        default=12,
        help="Number of top-variance tissues to include in MWU boxplot panel",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=DEFAULT_N_PERMUTATIONS,
        help="Number of permutations for pairwise group tests",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_csv = normalize_runtime_path(args.input_csv)
    out_dir = normalize_runtime_path(args.out_dir)
    metadata_csv = (
        normalize_runtime_path(args.sample_metadata_csv)
        if args.sample_metadata_csv
        else None
    )

    group_map = load_group_map(
        metadata_csv,
        sample_col=args.metadata_sample_col,
        group_col=args.metadata_group_col,
    )

    print("=== MethAtlas Deconvolution PCA + MWU Boxplot Analysis ===")
    print(f"Input CSV: {input_csv}")
    print(f"Output directory: {out_dir}")
    print(f"Metadata CSV: {metadata_csv}")

    raw_df, long_df, sample_feature_df, group_source = load_input_tables(
        input_csv,
        group_map=group_map,
        keep_unknown=bool(args.keep_unknown),
    )

    analysis_groups = [g for g in get_ordered_groups(long_df["group"].unique()) if g != "unknown"]
    if len(analysis_groups) < 2:
        raise ValueError(
            "Need at least 2 non-unknown groups for comparisons. "
            "Adjust GROUP_RULES or provide --sample-metadata-csv"
        )

    sample_feature_df = sample_feature_df[
        sample_feature_df.index.get_level_values("group").isin(analysis_groups)
    ].copy()
    long_df = long_df[long_df["group"].isin(analysis_groups)].copy()

    print(f"Raw table shape: {raw_df.shape}")
    print(f"Long table shape: {long_df.shape}")
    print(f"Sample-feature matrix shape: {sample_feature_df.shape}")
    print(f"Group source: {group_source}")
    print(f"Analysis groups: {analysis_groups}")
    print(
        "Detected groups: "
        f"{sample_feature_df.index.get_level_values('group').value_counts().to_dict()}"
    )

    results = perform_pca(
        sample_feature_df,
        analysis_groups=analysis_groups,
        n_permutations=max(1, int(args.n_permutations)),
    )
    print(f"PC1 explains {results['explained_variance'][0]:.1%} of variance")
    print(f"PC2 explains {results['explained_variance'][1]:.1%} of variance")
    print(f"PC1+PC2 explains {results['cumulative_variance'][1]:.1%} of variance")

    create_comprehensive_plot(results, out_dir)
    create_tissue_boxplots(
        long_df,
        out_dir,
        analysis_groups=analysis_groups,
        top_n_tissues=max(1, int(args.top_n_tissues)),
        n_permutations=max(1, int(args.n_permutations)),
        remove_iqr_outliers=True,
    )
    create_tissue_boxplots(
        long_df,
        out_dir,
        analysis_groups=analysis_groups,
        top_n_tissues=max(1, int(args.top_n_tissues)),
        n_permutations=max(1, int(args.n_permutations)),
        remove_iqr_outliers=False,
    )
    save_summary_and_tables(results, sample_feature_df, out_dir)

    print("\n=== Analysis Complete ===")


if __name__ == "__main__":
    main()
