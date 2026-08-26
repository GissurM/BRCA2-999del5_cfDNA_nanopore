#!/usr/bin/env python3
"""
Group-wise cfDNA size analysis from detailed 5 bp bin files.

This script generates:
1) A KS-based line plot of group-average size distributions (25-1000 bp)
   and KS p-values comparing peak-size distributions between groups.
2) Boxplots comparing group distributions of percent cfDNA in each
   fragment class (mono/di/tri/HMW), including a combined panel.
3) Boxplots comparing group distributions of mean fragment size in each
   fragment class, including a combined panel.
4) Additional percent/mean boxplots using an alternate mononucleosomal
    definition (100-200 bp instead of 25-200 bp), including combined panels.
5) CSV summaries for KS and Mann-Whitney U tests.
"""

import argparse
import os
import re
import warnings
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ks_2samp, mannwhitneyu


FRAGMENT_TYPES = ["mononucleosomal", "dinucleosomal", "trinucleosomal", "hmw"]
GROUPS = ["brca2", "control"]
GROUP_COLOR_MAP = {
    "brca2": "red",
    "control": "blue",
}
GROUP_DISPLAY_MAP = {
    "brca2": r"$\it{BRCA2\!\!-\!\!999del5}$",
    "control": "control",
}


def normalize_runtime_path(path_str: str) -> str:
    """Normalize path style between Windows and WSL-style paths."""
    path_str = str(path_str)

    # Convert /mnt/d/... -> D:\... when running on Windows.
    if os.name == "nt" and path_str.startswith("/mnt/"):
        parts = path_str.split("/")
        if len(parts) >= 4 and len(parts[2]) == 1:
            drive = parts[2].upper() + ":"
            tail = parts[3:]
            return os.path.join(drive + os.sep, *tail)

    # Convert D:\... -> /mnt/d/... when running on POSIX/WSL.
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", path_str):
        drive = path_str[0].lower()
        rest = path_str[2:].replace("\\", "/").lstrip("/")
        return f"/mnt/{drive}/{rest}"

    return path_str


def infer_group(sample_name: str) -> str:
    """Infer group from sample naming convention."""
    normalized = os.path.splitext(os.path.basename(str(sample_name)))[0].lower()

    if re.match(r"^brca2([_-]|$)", normalized):
        return "brca2"
    if re.match(r"^ctrl([_-]|$)", normalized):
        return "control"
    return "unknown"


def parse_size_columns(df: pd.DataFrame, fragment_type: str, size_min=None, size_max=None):
    """Return ordered size columns and their integer bp centers for a fragment type."""
    pattern = re.compile(rf"^{re.escape(fragment_type)}_(\d+)bp$")
    col_to_size = {}

    for col in df.columns:
        match = pattern.match(col)
        if not match:
            continue

        size_bp = int(match.group(1))
        if size_min is not None and size_bp < size_min:
            continue
        if size_max is not None and size_bp > size_max:
            continue
        col_to_size[col] = size_bp

    if not col_to_size:
        raise ValueError(f"No size columns found for {fragment_type}")

    ordered_cols = sorted(col_to_size, key=lambda col: col_to_size[col])
    ordered_sizes = np.array([col_to_size[col] for col in ordered_cols], dtype=float)
    return ordered_cols, ordered_sizes


def load_fragment_tables(data_dir: str, suffix: str):
    """Load fragment CSV tables (distributions or counts) with inferred groups."""
    tables = {}

    for fragment_type in FRAGMENT_TYPES:
        file_name = f"{fragment_type}_size_{suffix}.csv"
        file_path = os.path.join(data_dir, file_name)

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing required file: {file_path}")

        df = pd.read_csv(file_path)
        if "bam_file" not in df.columns:
            raise ValueError(f"Expected bam_file column in {file_path}")

        df = df.copy()
        df["group"] = df["bam_file"].apply(infer_group)
        df = df[df["group"].isin(GROUPS)].copy()

        if df.empty:
            raise ValueError(f"No usable samples (BRCA2/ctrl prefixes) found in {file_path}")

        tables[fragment_type] = df

    return tables


def filter_mononucleosomal_min_bp(tables, min_bp: int):
    """Copy tables and drop mononucleosomal size columns below min_bp."""
    filtered = {}
    mono_pattern = re.compile(r"^mononucleosomal_(\d+)bp$")

    for fragment_type, df in tables.items():
        sub = df.copy()

        if fragment_type == "mononucleosomal":
            keep_cols = []
            for col in sub.columns:
                match = mono_pattern.match(col)
                if match and int(match.group(1)) < min_bp:
                    continue
                keep_cols.append(col)

            sub = sub[keep_cols].copy()

            has_mono_bins = any(mono_pattern.match(col) for col in sub.columns)
            if not has_mono_bins:
                raise ValueError(
                    f"No mononucleosomal size columns remain after filtering <{min_bp} bp"
                )

        filtered[fragment_type] = sub

    return filtered


def filter_hmw_min_bp(tables, min_bp: int):
    """Copy tables and drop HMW size columns <= min_bp threshold."""
    filtered = {}
    hmw_pattern = re.compile(r"^hmw_(\d+)bp$")

    for fragment_type, df in tables.items():
        sub = df.copy()

        if fragment_type == "hmw":
            keep_cols = []
            for col in sub.columns:
                match = hmw_pattern.match(col)
                if match and int(match.group(1)) <= min_bp:
                    continue
                keep_cols.append(col)

            sub = sub[keep_cols].copy()

            has_hmw_bins = any(hmw_pattern.match(col) for col in sub.columns)
            if not has_hmw_bins:
                raise ValueError(
                    f"No HMW size columns remain after filtering <= {min_bp} bp"
                )

        filtered[fragment_type] = sub

    return filtered


def build_distribution_long_for_ks(dist_tables, size_min: int, size_max: int):
    """Build long table of per-size percentage values for KS/line plotting."""
    parts = []

    for fragment_type, df in dist_tables.items():
        try:
            size_cols, size_vals = parse_size_columns(
                df, fragment_type, size_min=size_min, size_max=size_max
            )
        except ValueError:
            # Some classes may not contribute bins within the requested KS size window.
            continue

        sub = df[["bam_file", "group"] + size_cols].copy()
        col_to_size = {col: int(bp) for col, bp in zip(size_cols, size_vals)}

        long_df = sub.melt(
            id_vars=["bam_file", "group"],
            value_vars=size_cols,
            var_name="size_col",
            value_name="percent_cfDNA",
        )
        long_df["size_bp"] = long_df["size_col"].map(col_to_size).astype(int)
        long_df = long_df.drop(columns=["size_col"])
        parts.append(long_df)

    if not parts:
        raise ValueError(
            f"No size columns found in requested KS range {size_min}-{size_max} bp"
        )

    combined = pd.concat(parts, ignore_index=True)
    combined = (
        combined.groupby(["bam_file", "group", "size_bp"], as_index=False)["percent_cfDNA"]
        .sum()
        .sort_values(["group", "bam_file", "size_bp"])
    )
    return combined


def build_percent_by_class(dist_tables, rebase_total_min_bp=None):
    """Compute per-sample class percentages, with optional denominator rebasing."""
    rows = []
    totals_for_rebase = []

    for fragment_type, df in dist_tables.items():
        size_cols, size_vals = parse_size_columns(df, fragment_type)
        percent_in_class = df[size_cols].sum(axis=1).to_numpy(dtype=float)

        if rebase_total_min_bp is not None:
            rebase_cols = [
                col for col, size_bp in zip(size_cols, size_vals) if size_bp >= rebase_total_min_bp
            ]
            if rebase_cols:
                included_total = df[rebase_cols].sum(axis=1).to_numpy(dtype=float)
            else:
                included_total = np.zeros(len(df), dtype=float)

            totals_for_rebase.append(
                pd.DataFrame(
                    {
                        "bam_file": df["bam_file"].values,
                        "group": df["group"].values,
                        "included_percent_total": included_total,
                    }
                )
            )

        rows.append(
            pd.DataFrame(
                {
                    "bam_file": df["bam_file"].values,
                    "group": df["group"].values,
                    "fragment_type": fragment_type,
                    "percent_cfDNA": percent_in_class,
                }
            )
        )

    out_df = pd.concat(rows, ignore_index=True)

    if rebase_total_min_bp is None:
        return out_df

    total_df = (
        pd.concat(totals_for_rebase, ignore_index=True)
        .groupby(["bam_file", "group"], as_index=False)["included_percent_total"]
        .sum()
    )
    out_df = out_df.merge(total_df, on=["bam_file", "group"], how="left")

    numer = 100.0 * out_df["percent_cfDNA"].to_numpy(dtype=float)
    denom = out_df["included_percent_total"].to_numpy(dtype=float)
    out_df["percent_cfDNA"] = np.divide(
        numer,
        denom,
        out=np.full(numer.shape, np.nan, dtype=float),
        where=denom > 0,
    )

    return out_df.drop(columns=["included_percent_total"])


def build_mean_size_by_class(count_tables):
    """Compute per-sample weighted mean fragment size (bp) in each class."""
    rows = []

    for fragment_type, df in count_tables.items():
        size_cols, size_vals = parse_size_columns(df, fragment_type)

        counts = df[size_cols].to_numpy(dtype=float)
        totals = counts.sum(axis=1)
        weighted_sum = counts @ size_vals

        mean_size = np.divide(
            weighted_sum,
            totals,
            out=np.full(weighted_sum.shape, np.nan, dtype=float),
            where=totals > 0,
        )

        rows.append(
            pd.DataFrame(
                {
                    "bam_file": df["bam_file"].values,
                    "group": df["group"].values,
                    "fragment_type": fragment_type,
                    "mean_size_bp": mean_size,
                    "fragments_in_class": totals,
                }
            )
        )

    return pd.concat(rows, ignore_index=True)


def run_ks_test(x, y):
    """Run KS test with exact-first strategy and stable p-value fallback."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ks_stat, p_value = ks_2samp(x, y, alternative="two-sided", mode="exact")

    exact_warning = any("Exact calculation unsuccessful" in str(w.message) for w in caught)
    if exact_warning:
        ks_stat, p_value = ks_2samp(x, y, alternative="two-sided", mode="asymp")

    if p_value == 0.0:
        p_value = np.finfo(float).tiny

    return ks_stat, p_value


def compute_peak_ks_tests(dist_long):
    """Find sample-level peak bins and run KS tests between groups."""
    peak_rows = []

    for (bam_file, group), sub in dist_long.groupby(["bam_file", "group"]):
        idx = sub["percent_cfDNA"].idxmax()
        peak_rows.append(
            {
                "bam_file": bam_file,
                "group": group,
                "peak_size_bp": int(sub.loc[idx, "size_bp"]),
                "peak_percent_cfDNA": float(sub.loc[idx, "percent_cfDNA"]),
            }
        )

    peak_df = pd.DataFrame(peak_rows)

    ks_records = []
    for group1, group2 in combinations(GROUPS, 2):
        x = peak_df.loc[peak_df["group"] == group1, "peak_size_bp"].dropna()
        y = peak_df.loc[peak_df["group"] == group2, "peak_size_bp"].dropna()
        ks_stat, p_value = run_ks_test(x, y)

        ks_records.append(
            {
                "test": "KS_on_peak_size_distribution",
                "group1": group1,
                "group2": group2,
                "n1": int(len(x)),
                "n2": int(len(y)),
                "statistic": ks_stat,
                "p_value": p_value,
            }
        )

    ks_df = pd.DataFrame(ks_records)
    return peak_df, ks_df


def _pick_local_maxima(size_vals, pct_vals, min_peak_bp, top_n, min_sep_bp):
    """Pick top local maxima from one curve, enforcing minimal separation."""
    candidates = []
    for i in range(1, len(pct_vals) - 1):
        if (
            size_vals[i] >= min_peak_bp
            and pct_vals[i] > pct_vals[i - 1]
            and pct_vals[i] >= pct_vals[i + 1]
        ):
            candidates.append((int(size_vals[i]), float(pct_vals[i])))

    if not candidates and len(pct_vals) > 0:
        valid_idx = [i for i in range(len(pct_vals)) if size_vals[i] >= min_peak_bp]
        if not valid_idx:
            valid_idx = list(range(len(pct_vals)))
        top_idx = [valid_idx[i] for i in np.argsort(pct_vals[valid_idx])[-max(1, top_n):]]
        candidates = [(int(size_vals[i]), float(pct_vals[i])) for i in top_idx]

    candidates = sorted(candidates, key=lambda t: t[1], reverse=True)
    selected = []
    for size_bp, pct in candidates:
        if all(abs(size_bp - s) >= min_sep_bp for s, _ in selected):
            selected.append((size_bp, pct))
        if len(selected) >= top_n:
            break

    return selected


def _find_peak_bound_left(sizes: np.ndarray, values: np.ndarray, peak_idx: int) -> int:
    """Walk left from peak_idx to find the nearest valley (start of the peak)."""
    for i in range(peak_idx - 1, 0, -1):
        # Once the next step leftward would rise, we have found the valley floor.
        if values[i - 1] >= values[i]:
            return i
    return 0


def _find_peak_bound_right(sizes: np.ndarray, values: np.ndarray, peak_idx: int) -> int:
    """Walk right from peak_idx to find the nearest valley (end of the peak)."""
    n = len(values)
    for i in range(peak_idx + 1, n - 1):
        if values[i + 1] >= values[i]:
            return i
    return n - 1


def _peak_auc_above_baseline(
    sizes: np.ndarray,
    values: np.ndarray,
    start_idx: int,
    end_idx: int,
) -> float:
    """Compute baseline-corrected peak area (%*bp) between local start/end valleys."""
    if end_idx <= start_idx:
        return 0.0

    x = sizes[start_idx : end_idx + 1]
    y = values[start_idx : end_idx + 1]
    if x.size < 2:
        return 0.0

    x0 = float(x[0])
    x1 = float(x[-1])
    y0 = float(y[0])
    y1 = float(y[-1])
    if x1 == x0:
        return 0.0

    baseline = y0 + (y1 - y0) * ((x - x0) / (x1 - x0))
    y_excess = np.maximum(y - baseline, 0.0)
    return float(np.trapezoid(y_excess, x))


def compute_major_peaks(
    dist_long,
    peaks_per_group=3,
    max_total_peaks=4,
    min_sep_bp=70,
    min_peak_bp=25,
    merge_tol_bp=25,
    source_prominence_ratio=0.25,
):
    """
    Detect biologically relevant peaks by combining top peaks from each group.

    This keeps group-specific peaks (including coronary-only early peaks) while
    merging near-duplicate centers across groups.
    """
    mean_curve = (
        dist_long.groupby(["group", "size_bp"], as_index=False)["percent_cfDNA"]
        .mean()
        .sort_values(["group", "size_bp"])
    )

    candidates = []
    for group in GROUPS:
        sub = mean_curve[mean_curve["group"] == group]
        if sub.empty:
            continue

        sizes = sub["size_bp"].to_numpy(dtype=float)
        pcts = sub["percent_cfDNA"].to_numpy(dtype=float)
        peaks = _pick_local_maxima(
            sizes,
            pcts,
            min_peak_bp=min_peak_bp,
            top_n=peaks_per_group,
            min_sep_bp=min_sep_bp,
        )

        for size_bp, pct in peaks:
            candidates.append(
                {
                    "group": group,
                    "peak_size_bp": int(size_bp),
                    "peak_percent_cfDNA": float(pct),
                }
            )

    if not candidates:
        return pd.DataFrame(columns=["peak_id", "peak_size_bp", "peak_percent_cfDNA", "source_groups"])

    # Merge nearby candidate peaks into shared peak centers.
    clusters = []
    for cand in sorted(candidates, key=lambda r: r["peak_size_bp"]):
        placed = False
        for cluster in clusters:
            if abs(cand["peak_size_bp"] - cluster["center_bp"]) <= merge_tol_bp:
                cluster["members"].append(cand)
                sizes = np.array([m["peak_size_bp"] for m in cluster["members"]], dtype=float)
                weights = np.array([m["peak_percent_cfDNA"] for m in cluster["members"]], dtype=float)
                cluster["center_bp"] = float(np.average(sizes, weights=weights))
                cluster["max_pct"] = max(cluster["max_pct"], cand["peak_percent_cfDNA"])
                cluster["groups"].add(cand["group"])
                cluster["group_max_pct"][cand["group"]] = max(
                    cluster["group_max_pct"].get(cand["group"], 0.0),
                    cand["peak_percent_cfDNA"],
                )
                placed = True
                break

        if not placed:
            clusters.append(
                {
                    "center_bp": float(cand["peak_size_bp"]),
                    "max_pct": float(cand["peak_percent_cfDNA"]),
                    "groups": {cand["group"]},
                    "group_max_pct": {cand["group"]: float(cand["peak_percent_cfDNA"])},
                    "members": [cand],
                }
            )

    # Keep strongest clusters, then sort by position.
    clusters = sorted(clusters, key=lambda c: c["max_pct"], reverse=True)[:max_total_peaks]
    clusters = sorted(clusters, key=lambda c: c["center_bp"])

    rows = []
    for i, cluster in enumerate(clusters, start=1):
        dominant_groups = sorted(
            [
                g
                for g, g_peak in cluster["group_max_pct"].items()
                if g_peak >= source_prominence_ratio * cluster["max_pct"]
            ]
        )
        if not dominant_groups:
            dominant_groups = sorted(cluster["groups"])

        rows.append(
            {
                "peak_id": f"peak_{i}",
                "peak_size_bp": int(round(cluster["center_bp"])),
                "peak_percent_cfDNA": float(cluster["max_pct"]),
                "source_groups": "/".join(dominant_groups),
            }
        )

    return pd.DataFrame(rows)


def compute_peak_window_ks_tests(dist_long, peak_df, window_bp=25):
    """
    Compute peak-window test features and KS tests for each detected peak.

    For each peak center, each sample contributes:
    - local peak position (bp)
    - local peak height (% cfDNA)
    - local peak area (% cfDNA in the window)

    KS tests are computed between group pairs for all three feature types.
    """
    local_peak_rows = []
    ks_rows = []

    for _, peak_row in peak_df.iterrows():
        peak_id = peak_row["peak_id"]
        center = int(peak_row["peak_size_bp"])

        window = dist_long[
            (dist_long["size_bp"] >= center - window_bp)
            & (dist_long["size_bp"] <= center + window_bp)
        ]

        for (bam_file, group), sub in window.groupby(["bam_file", "group"]):
            idx = sub["percent_cfDNA"].idxmax()
            local_area = float(sub["percent_cfDNA"].sum())
            local_peak_rows.append(
                {
                    "peak_id": peak_id,
                    "window_center_bp": center,
                    "window_bp": int(window_bp),
                    "bam_file": bam_file,
                    "group": group,
                    "local_peak_size_bp": int(sub.loc[idx, "size_bp"]),
                    "local_peak_percent_cfDNA": float(sub.loc[idx, "percent_cfDNA"]),
                    "local_peak_area_percent_cfDNA": local_area,
                }
            )

        local_peak_df = pd.DataFrame(
            [
                row
                for row in local_peak_rows
                if row["peak_id"] == peak_id and row["window_center_bp"] == center
            ]
        )

        metrics = {
            "local_peak_height_percent": "local_peak_percent_cfDNA",
            "local_peak_position_bp": "local_peak_size_bp",
            "local_peak_area_percent": "local_peak_area_percent_cfDNA",
        }

        for metric_name, metric_col in metrics.items():
            for group1, group2 in combinations(GROUPS, 2):
                x = local_peak_df.loc[local_peak_df["group"] == group1, metric_col].dropna()
                y = local_peak_df.loc[local_peak_df["group"] == group2, metric_col].dropna()
                ks_stat, p_value = run_ks_test(x, y)

                ks_rows.append(
                    {
                        "test": "KS_on_local_peak_feature_distribution",
                        "metric": metric_name,
                        "peak_id": peak_id,
                        "window_center_bp": center,
                        "window_bp": int(window_bp),
                        "group1": group1,
                        "group2": group2,
                        "n1": int(len(x)),
                        "n2": int(len(y)),
                        "statistic": ks_stat,
                        "p_value": p_value,
                    }
                )

    return pd.DataFrame(local_peak_rows), pd.DataFrame(ks_rows)


def compute_mwu_by_fragment(metric_df, value_col: str, metric_name: str):
    """Run Mann-Whitney U tests for each fragment class and group pair."""
    records = []

    for fragment_type in FRAGMENT_TYPES:
        sub = metric_df[metric_df["fragment_type"] == fragment_type]

        for group1, group2 in combinations(GROUPS, 2):
            x = sub.loc[sub["group"] == group1, value_col].dropna()
            y = sub.loc[sub["group"] == group2, value_col].dropna()

            if len(x) == 0 or len(y) == 0:
                stat, p_value = np.nan, np.nan
            else:
                stat, p_value = mannwhitneyu(x, y, alternative="two-sided")

            records.append(
                {
                    "test": f"MWU_{metric_name}",
                    "metric": metric_name,
                    "fragment_type": fragment_type,
                    "group1": group1,
                    "group2": group2,
                    "n1": int(len(x)),
                    "n2": int(len(y)),
                    "statistic": stat,
                    "p_value": p_value,
                }
            )

    return pd.DataFrame(records)


def summarize_metric_by_fragment(metric_df, value_col: str, metric_name: str):
    """Summarize each fragment class by group with n, median, and sample SD."""
    records = []

    for fragment_type in FRAGMENT_TYPES:
        sub = metric_df[metric_df["fragment_type"] == fragment_type]

        for group in GROUPS:
            vals = sub.loc[sub["group"] == group, value_col].dropna().to_numpy(dtype=float)

            records.append(
                {
                    "metric": metric_name,
                    "fragment_type": fragment_type,
                    "group": group,
                    "n": int(vals.size),
                    "median": float(np.median(vals)) if vals.size > 0 else np.nan,
                    "sd": float(np.std(vals, ddof=1)) if vals.size > 1 else np.nan,
                }
            )

    return pd.DataFrame(records)


def _format_pvalue_text(p):
    if pd.isna(p):
        return "nan"
    return f"{p:.3g}"


def pvalue_to_stars(p):
    """Map p-values to significance stars using user-requested thresholds."""
    if pd.isna(p):
        return "ns"
    if p < 0.0005:
        return "***"
    if p < 0.005:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def add_significance_bar(ax, x1, x2, y, h, text):
    """Draw a significance bar between two category positions."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, c="black")
    ax.text((x1 + x2) / 2.0, y + h, text, ha="center", va="bottom", fontsize=11)


def add_single_pair_significance(ax, sub_df, value_col, p_value):
    """Annotate one pairwise significance bar for 2-group comparisons."""
    if pd.isna(p_value):
        return

    stars = pvalue_to_stars(p_value)
    if stars == "ns":
        return

    y_vals = sub_df[value_col].dropna().to_numpy(dtype=float)
    if y_vals.size == 0:
        return

    y_min = float(np.nanmin(y_vals))
    y_max = float(np.nanmax(y_vals))
    span = max(y_max - y_min, 1e-9)
    y = y_max + 0.08 * span
    h = 0.04 * span
    add_significance_bar(ax, 0, 1, y, h, stars)

    top = y + h + 0.08 * span
    cur_bottom, cur_top = ax.get_ylim()
    ax.set_ylim(cur_bottom, max(cur_top, top))


def add_pairwise_text_box(ax, stats_df):
    """Add pairwise p-values text box to an axis."""
    lines = []
    for _, row in stats_df.iterrows():
        g1 = GROUP_DISPLAY_MAP.get(row["group1"], row["group1"])
        g2 = GROUP_DISPLAY_MAP.get(row["group2"], row["group2"])
        stars = pvalue_to_stars(row["p_value"])
        lines.append(f"{g1} vs {g2}: p={_format_pvalue_text(row['p_value'])} ({stars})")

    if not lines:
        return

    text = "MWU p-values\n" + "\n".join(lines)
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "gray"},
    )


def plot_ks_distribution(
    dist_long,
    major_peak_df,
    peakwise_ks_df,
    out_path,
    size_min,
    size_max,
):
    """Plot group-average size distributions with KS peak p-values."""
    sns.set_theme(style="whitegrid")

    mean_curve = (
        dist_long.groupby(["group", "size_bp"], as_index=False)["percent_cfDNA"]
        .mean()
        .sort_values(["group", "size_bp"])
    )
    std_curve = (
        dist_long.groupby(["group", "size_bp"], as_index=False)["percent_cfDNA"]
        .std()
        .rename(columns={"percent_cfDNA": "std_percent_cfDNA"})
    )
    mean_curve = mean_curve.merge(std_curve, on=["group", "size_bp"], how="left")
    mean_curve["std_percent_cfDNA"] = mean_curve["std_percent_cfDNA"].fillna(0.0)
    sample_counts = (
        dist_long[["bam_file", "group"]]
        .drop_duplicates()
        .groupby("group")
        .size()
        .to_dict()
    )

    fig, ax = plt.subplots(figsize=(13, 7))
    line_styles = {"brca2": "-", "control": "--"}

    for group in GROUPS:
        sub = mean_curve[mean_curve["group"] == group]
        if sub.empty:
            continue

        color = GROUP_COLOR_MAP[group]
        x_vals = sub["size_bp"].to_numpy(dtype=float)
        y_vals = sub["percent_cfDNA"].to_numpy(dtype=float)
        sd_vals = sub["std_percent_cfDNA"].to_numpy(dtype=float)

        ax.fill_between(
            x_vals,
            np.maximum(y_vals - sd_vals, 0.0),
            y_vals + sd_vals,
            color=color,
            alpha=0.12,
            linewidth=0,
            zorder=1,
        )

        line = ax.plot(
            sub["size_bp"],
            sub["percent_cfDNA"],
            color=color,
            linewidth=3.0,
            linestyle=line_styles.get(group, "-"),
            label=f"{GROUP_DISPLAY_MAP.get(group, group)} (n={sample_counts.get(group, 0)})",
            zorder=3,
        )[0]
        line.set_path_effects([pe.Stroke(linewidth=4.2, foreground="white"), pe.Normal()])

        peak_idx = sub["percent_cfDNA"].idxmax()
        peak_size = int(sub.loc[peak_idx, "size_bp"])
        peak_val = float(sub.loc[peak_idx, "percent_cfDNA"])
        ax.scatter([peak_size], [peak_val], color=color, s=50, edgecolor="white", linewidth=0.7, zorder=5)

    ax.set_xlim(size_min, size_max)
    ax.set_xlabel("Fragment size (bp)")
    ax.set_ylabel("% of cfDNA at each size bin")
    ax.set_title(
        "Kolmogorov-Smirnov Group Peak Comparison\n"
        "Average cfDNA size distributions (25-1000 bp) with peak-window KS p-values",
        pad=18,
    )
    ax.legend(loc="upper right")

    y_max = float(mean_curve["percent_cfDNA"].max()) if not mean_curve.empty else 1.0
    y_top = max(y_max * 1.35, y_max + 1.0)
    ax.set_ylim(bottom=0, top=y_top)

    label_positions = []
    min_label_dx = 260
    min_label_dy = max(0.22 * y_max, 0.25)

    for i, peak_row in major_peak_df.reset_index(drop=True).iterrows():
        peak_id = peak_row["peak_id"]
        center = int(peak_row["peak_size_bp"])
        peak_val = float(peak_row["peak_percent_cfDNA"])
        color = "#444444"

        ax.axvline(center, color=color, linestyle="--", linewidth=1.4, alpha=0.45)

        sub_stats = peakwise_ks_df[
            (peakwise_ks_df["peak_id"] == peak_id)
            & (peakwise_ks_df["window_center_bp"] == center)
            & (peakwise_ks_df["metric"] == "local_peak_height_percent")
        ]
        stat_lines = []
        for _, row in sub_stats.iterrows():
            g1 = GROUP_DISPLAY_MAP.get(row["group1"], row["group1"])
            g2 = GROUP_DISPLAY_MAP.get(row["group2"], row["group2"])
            stars = pvalue_to_stars(row["p_value"])
            stat_lines.append(f"{g1} vs {g2}: p={_format_pvalue_text(row['p_value'])} ({stars})")
        stat_text = "\n".join(stat_lines) if stat_lines else "No test data"
        source_groups = peak_row.get("source_groups", "")
        source_text = f" [from {source_groups}]" if source_groups else ""

        x_offset = 14
        y_offset = 0.22 + 0.18 * i
        y_text = min(peak_val + y_offset, y_top - 0.18 * y_max)
        x_text = min(center + x_offset, size_max - 240)

        # Nudge labels vertically to reduce overlap among nearby annotation boxes.
        for _ in range(20):
            has_conflict = any(
                abs(x_text - prev_x) < min_label_dx and abs(y_text - prev_y) < min_label_dy
                for prev_x, prev_y in label_positions
            )
            if not has_conflict:
                break
            y_text += min_label_dy
            if y_text > y_top - 0.1 * y_max:
                y_text = max(peak_val * 0.6, y_top - (0.35 + 0.12 * i) * y_max)
                break

        ax.annotate(
            f"{peak_id} ({center} bp){source_text}\n{stat_text}",
            xy=(center, peak_val),
            xytext=(x_text, y_text),
            fontsize=8.5,
            ha="left",
            va="bottom",
            arrowprops={"arrowstyle": "-", "color": color, "lw": 0.9, "alpha": 0.75},
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.82,
                "edgecolor": color,
            },
        )
        label_positions.append((x_text, y_text))

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_single_group_distribution(
    dist_long: pd.DataFrame,
    group: str,
    out_path: str,
    size_min: int,
    size_max: int,
    peaks_per_plot: int = 4,
    min_sep_bp: int = 100,
    min_peak_bp: int = 25,
) -> None:
    """
    Plot the cfDNA size distribution for a single group.

    Individual sample traces are drawn as faint lines behind a bold group mean ± 1 SD
    band.  For each detected peak the following positions are marked on the mean curve:

        ◄  start – nearest valley to the left of the apex
        ●  apex  – highest point of the peak
        ►  end   – nearest valley to the right of the apex
    """
    sns.set_theme(style="whitegrid")

    group_data = dist_long[dist_long["group"] == group].copy()
    if group_data.empty:
        return

    color = GROUP_COLOR_MAP[group]
    display_name = GROUP_DISPLAY_MAP.get(group, group)

    mean_curve = (
        group_data.groupby("size_bp", as_index=False)["percent_cfDNA"]
        .mean()
        .sort_values("size_bp")
        .reset_index(drop=True)
    )
    std_vals = (
        group_data.groupby("size_bp", as_index=False)["percent_cfDNA"]
        .std()
        .rename(columns={"percent_cfDNA": "std_percent_cfDNA"})
    )
    mean_curve = mean_curve.merge(std_vals, on="size_bp", how="left")
    mean_curve["std_percent_cfDNA"] = mean_curve["std_percent_cfDNA"].fillna(0.0)

    sizes = mean_curve["size_bp"].to_numpy(dtype=float)
    means = mean_curve["percent_cfDNA"].to_numpy(dtype=float)
    sds = mean_curve["std_percent_cfDNA"].to_numpy(dtype=float)

    sample_names = sorted(group_data["bam_file"].unique())
    n_samples = len(sample_names)

    # Detect group-specific peaks and sort by position.
    peaks = _pick_local_maxima(
        sizes, means,
        min_peak_bp=min_peak_bp,
        top_n=peaks_per_plot,
        min_sep_bp=min_sep_bp,
    )
    peaks = sorted(peaks, key=lambda t: t[0])

    fig, ax = plt.subplots(figsize=(13, 7))

    # Individual sample traces (no legend entry).
    for bam_file in sample_names:
        sample_sub = (
            group_data[group_data["bam_file"] == bam_file]
            .sort_values("size_bp")
        )
        ax.plot(
            sample_sub["size_bp"],
            sample_sub["percent_cfDNA"],
            color=color,
            alpha=0.05,
            linewidth=0.5,
            zorder=0,
        )

    # Mean line.
    mean_line = ax.plot(
        sizes, means,
        color=color,
        linewidth=3.0,
        zorder=4,
        label=f"{display_name} mean (n={n_samples})",
    )[0]
    mean_line.set_path_effects([
        pe.Stroke(linewidth=4.5, foreground="white"),
        pe.Normal(),
    ])

    y_max = float(np.nanmax(means)) if means.size > 0 else 1.0
    y_top = max(y_max * 1.50, y_max + 1.0)
    ax.set_ylim(bottom=0, top=y_top)
    ax.set_xlim(size_min, size_max)
    total_curve_auc = float(np.trapezoid(means, sizes)) if sizes.size >= 2 else float("nan")

    label_positions: list = []
    min_label_dx = 200
    min_label_dy = max(0.20 * y_max, 0.20)

    for peak_rank, (peak_size_bp, peak_pct) in enumerate(peaks):
        apex_idx = int(np.argmin(np.abs(sizes - peak_size_bp)))
        start_idx = _find_peak_bound_left(sizes, means, apex_idx)
        end_idx = _find_peak_bound_right(sizes, means, apex_idx)

        start_bp = int(sizes[start_idx])
        end_bp = int(sizes[end_idx])
        start_pct = float(means[start_idx])
        end_pct = float(means[end_idx])

        # Simple zero-baseline shading within each peak window for cleaner visuals.
        mask = (sizes >= start_bp) & (sizes <= end_bp)
        x_seg = sizes[mask]
        y_seg = means[mask]

        shade_mask = y_seg > 0
        ax.fill_between(
            x_seg,
            0.0,
            y_seg,
            where=shade_mask,
            alpha=0.10,
            color=color,
            zorder=1,
        )

        peak_auc = _peak_auc_above_baseline(
            sizes,
            means,
            start_idx,
            end_idx,
        )
        peak_auc_pct_total = (
            (100.0 * peak_auc / total_curve_auc)
            if np.isfinite(total_curve_auc) and total_curve_auc > 0
            else float("nan")
        )

        # Start marker (left-pointing triangle).
        ax.scatter(
            [start_bp], [start_pct],
            marker="<", s=70, color=color,
            edgecolor="white", linewidth=0.9, zorder=6,
        )
        # Apex marker (circle).
        ax.scatter(
            [peak_size_bp], [peak_pct],
            marker="o", s=90, color=color,
            edgecolor="white", linewidth=1.0, zorder=6,
        )
        # End marker (right-pointing triangle).
        ax.scatter(
            [end_bp], [end_pct],
            marker=">", s=70, color=color,
            edgecolor="white", linewidth=0.9, zorder=6,
        )

        ann_text = (
            f"Peak {peak_rank + 1}\n"
            f"\u25c4 Start:  {start_bp} bp\n"
            f"\u25cf Apex:   {peak_size_bp} bp ({peak_pct:.2f}%)\n"
            f"\u25ba End:    {end_bp} bp\n"
            f"AUCex/total: {peak_auc_pct_total:.2f}%"
        )

        x_text = min(float(peak_size_bp) + 18.0, float(size_max) - 220.0)
        y_text = min(peak_pct + 0.35 + 0.18 * peak_rank, y_top - 0.18 * y_max)

        for _ in range(20):
            conflict = any(
                abs(x_text - px) < min_label_dx and abs(y_text - py) < min_label_dy
                for px, py in label_positions
            )
            if not conflict:
                break
            y_text += min_label_dy
            if y_text > y_top - 0.10 * y_max:
                y_text = max(peak_pct * 0.65, 0.15)
                break

        ax.annotate(
            ann_text,
            xy=(peak_size_bp, peak_pct),
            xytext=(x_text, y_text),
            fontsize=8.5,
            ha="left",
            va="bottom",
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.9, "alpha": 0.75},
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.85,
                "edgecolor": color,
            },
            zorder=7,
        )
        label_positions.append((x_text, y_text))

    # Dummy scatter calls so marker types appear in the legend.
    ax.scatter([], [], marker="<", s=60, color="#555555", label="Peak start")
    ax.scatter([], [], marker="o", s=75, color="#555555", label="Peak apex")
    ax.scatter([], [], marker=">", s=60, color="#555555", label="Peak end")

    ax.set_xlabel("Fragment size (bp)")
    ax.set_ylabel("% of cfDNA at each size bin")
    ax.set_title(
        f"{display_name} \u2014 cfDNA fragment size distribution\n"
        "Individual samples (faint) \u00b7 Group mean (bold) \u00b7 Peak markers",
        pad=16,
    )
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_metric_grid(metric_df, value_col, ylabel, title, out_path, stats_df):
    """Plot 2x2 boxplot grid for one metric across fragment classes."""
    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(2, 2, figsize=(16, 11))

    for ax, fragment_type in zip(axes.flatten(), FRAGMENT_TYPES):
        sub = metric_df[metric_df["fragment_type"] == fragment_type]
        frag_stats = stats_df[stats_df["fragment_type"] == fragment_type]

        sns.boxplot(
            data=sub,
            x="group",
            y=value_col,
            hue="group",
            order=GROUPS,
            hue_order=GROUPS,
            palette=GROUP_COLOR_MAP,
            showfliers=False,
            showmeans=False,
            medianprops={"color": "black", "linewidth": 2.0},
            width=0.62,
            ax=ax,
            legend=False,
        )
        sns.stripplot(
            data=sub,
            x="group",
            y=value_col,
            order=GROUPS,
            color="black",
            alpha=0.55,
            jitter=0.15,
            size=4,
            ax=ax,
        )

        ax.set_title(fragment_type)
        ax.set_xlabel("Group")
        ax.set_ylabel(ylabel)

        # Replace internal group keys with publication labels.
        ax.set_xticks(np.arange(len(GROUPS)))
        ax.set_xticklabels([GROUP_DISPLAY_MAP[g] for g in GROUPS])

        if len(GROUPS) == 2 and not frag_stats.empty:
            p_val = float(frag_stats.iloc[0]["p_value"])
            add_single_pair_significance(ax, sub, value_col, p_val)

        add_pairwise_text_box(ax, frag_stats)

    fig.suptitle(title, fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_metric_combined(metric_df, value_col, ylabel, title, out_path):
    """Plot one combined grouped boxplot across all fragment classes."""
    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(14, 7))

    sns.boxplot(
        data=metric_df,
        x="fragment_type",
        y=value_col,
        hue="group",
        order=FRAGMENT_TYPES,
        hue_order=GROUPS,
        palette=GROUP_COLOR_MAP,
        showfliers=False,
        showmeans=False,
        medianprops={"color": "black", "linewidth": 2.0},
        ax=ax,
    )
    sns.stripplot(
        data=metric_df,
        x="fragment_type",
        y=value_col,
        hue="group",
        order=FRAGMENT_TYPES,
        hue_order=GROUPS,
        dodge=True,
        alpha=0.3,
        size=2.8,
        palette={group: "black" for group in GROUPS},
        ax=ax,
        legend=False,
    )

    # Keep only one legend (from the boxplot handles).
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) >= len(GROUPS):
        display_labels = [GROUP_DISPLAY_MAP.get(g, g) for g in GROUPS]
        ax.legend(handles[: len(GROUPS)], display_labels, title="Group", loc="best")

    ax.set_xlabel("Fragment class")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Add one significance annotation per fragment class (2 groups expected).
    for i, fragment_type in enumerate(FRAGMENT_TYPES):
        sub = metric_df[metric_df["fragment_type"] == fragment_type]
        x_vals = sub[value_col].dropna().to_numpy(dtype=float)
        if x_vals.size == 0:
            continue

        g1_vals = sub.loc[sub["group"] == GROUPS[0], value_col].dropna()
        g2_vals = sub.loc[sub["group"] == GROUPS[1], value_col].dropna()
        if len(g1_vals) == 0 or len(g2_vals) == 0:
            continue

        _, p_val = mannwhitneyu(g1_vals, g2_vals, alternative="two-sided")
        stars = pvalue_to_stars(p_val)
        if stars == "ns":
            continue

        y_min = float(np.nanmin(x_vals))
        y_max = float(np.nanmax(x_vals))
        span = max(y_max - y_min, 1e-9)
        y = y_max + 0.08 * span
        h = 0.04 * span
        add_significance_bar(ax, i - 0.2, i + 0.2, y, h, stars)

    y_all = metric_df[value_col].dropna().to_numpy(dtype=float)
    if y_all.size > 0:
        y_min = float(np.nanmin(y_all))
        y_max = float(np.nanmax(y_all))
        span = max(y_max - y_min, 1e-9)
        _, cur_top = ax.get_ylim()
        ax.set_ylim(top=max(cur_top, y_max + 0.3 * span))

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Generate KS and MWU cfDNA group comparison plots from size-bin CSV files"
    )
    parser.add_argument(
        "--data-dir",
        default="/mnt/d/BRCA2-misc_files/cfDNA_detailed_size_analysis",
        help="Directory containing *_size_distributions.csv and *_size_counts.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="/mnt/d/BRCA2-misc_files/group_comparison_plots",
        help="Directory for output plots and test summaries",
    )
    parser.add_argument(
        "--size-min",
        type=int,
        default=25,
        help="Minimum size (bp) for KS line plot",
    )
    parser.add_argument(
        "--size-max",
        type=int,
        default=1000,
        help="Maximum size (bp) for KS line plot",
    )
    parser.add_argument(
        "--mono-alt-min-bp",
        type=int,
        default=100,
        help=(
            "Lower bound (bp) for alternate mononucleosomal definition used in "
            "additional boxplots (upper bound remains 200 bp)"
        ),
    )
    parser.add_argument(
        "--hmw-min-bp",
        type=int,
        default=1000,
        help="Lower cutoff for HMW bins; only HMW fragments > this bp are included",
    )
    args = parser.parse_args()

    if args.mono_alt_min_bp < 0:
        raise ValueError("--mono-alt-min-bp must be >= 0")
    if args.hmw_min_bp < 0:
        raise ValueError("--hmw-min-bp must be >= 0")

    args.data_dir = normalize_runtime_path(args.data_dir)
    args.out_dir = normalize_runtime_path(args.out_dir)

    os.makedirs(args.out_dir, exist_ok=True)

    print("=== cfDNA Group KS + MWU Analysis ===")
    print(f"Data directory: {args.data_dir}")
    print(f"Output directory: {args.out_dir}")

    # Load raw tables first; keep a copy for the KS plot so all bins ≤ size_max are used.
    dist_tables_raw = load_fragment_tables(args.data_dir, "distributions")
    count_tables = load_fragment_tables(args.data_dir, "counts")

    # Enforce HMW definition as fragments strictly above threshold (default >1000 bp).
    # Only the boxplot/MWU analyses use these filtered tables.
    hmw_min_bp = int(args.hmw_min_bp)
    dist_tables = filter_hmw_min_bp(dist_tables_raw, hmw_min_bp)
    count_tables = filter_hmw_min_bp(count_tables, hmw_min_bp)

    mono_alt_min_bp = int(args.mono_alt_min_bp)
    dist_tables_mono_alt = filter_mononucleosomal_min_bp(dist_tables, mono_alt_min_bp)
    count_tables_mono_alt = filter_mononucleosomal_min_bp(count_tables, mono_alt_min_bp)

    print(f"HMW minimum (strictly greater than): {hmw_min_bp} bp")
    print(f"Alternate mononucleosomal minimum: {mono_alt_min_bp} bp")

    # Use raw (unfiltered) distributions for the KS curve so HMW bins within 25-1000 bp are included.
    dist_long = build_distribution_long_for_ks(dist_tables_raw, args.size_min, args.size_max)
    peak_df, ks_df = compute_peak_ks_tests(dist_long)
    major_peak_df = compute_major_peaks(
        dist_long,
        peaks_per_group=3,
        max_total_peaks=4,
        min_sep_bp=100,   # raised from 70 to suppress sub-nucleosomal artefact peaks
        min_peak_bp=25,
        merge_tol_bp=25,
    )
    local_peak_df, peakwise_ks_df = compute_peak_window_ks_tests(
        dist_long, major_peak_df, window_bp=25
    )

    percent_df = build_percent_by_class(dist_tables)
    mean_size_df = build_mean_size_by_class(count_tables)

    percent_df_mono_alt = build_percent_by_class(
        dist_tables_mono_alt,
        rebase_total_min_bp=mono_alt_min_bp,
    )
    mean_size_df_mono_alt = build_mean_size_by_class(count_tables_mono_alt)

    mwu_percent_df = compute_mwu_by_fragment(percent_df, "percent_cfDNA", "percent_cfDNA")
    mwu_size_df = compute_mwu_by_fragment(mean_size_df, "mean_size_bp", "mean_size_bp")
    mwu_percent_df_mono_alt = compute_mwu_by_fragment(
        percent_df_mono_alt,
        "percent_cfDNA",
        f"percent_cfDNA_mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}",
    )
    mwu_size_df_mono_alt = compute_mwu_by_fragment(
        mean_size_df_mono_alt,
        "mean_size_bp",
        f"mean_size_bp_mono_{mono_alt_min_bp}_200bp",
    )

    summary_percent_df = summarize_metric_by_fragment(
        percent_df,
        "percent_cfDNA",
        "percent_cfDNA",
    )
    summary_size_df = summarize_metric_by_fragment(
        mean_size_df,
        "mean_size_bp",
        "mean_size_bp",
    )
    summary_percent_df_mono_alt = summarize_metric_by_fragment(
        percent_df_mono_alt,
        "percent_cfDNA",
        f"percent_cfDNA_mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}",
    )
    summary_size_df_mono_alt = summarize_metric_by_fragment(
        mean_size_df_mono_alt,
        "mean_size_bp",
        f"mean_size_bp_mono_{mono_alt_min_bp}_200bp",
    )

    ks_plot_path = os.path.join(args.out_dir, "ks_group_size_distribution_25_1000.png")
    brca2_solo_plot_path = os.path.join(args.out_dir, "size_distribution_brca2_solo.png")
    ctrl_solo_plot_path = os.path.join(args.out_dir, "size_distribution_control_solo.png")
    pct_grid_path = os.path.join(
        args.out_dir,
        "median_boxplot_percent_cfDNA_by_fragment_class.png",
    )
    pct_combined_path = os.path.join(
        args.out_dir,
        "median_boxplot_percent_cfDNA_all_fragment_classes.png",
    )
    size_grid_path = os.path.join(
        args.out_dir,
        "median_boxplot_mean_size_by_fragment_class.png",
    )
    size_combined_path = os.path.join(
        args.out_dir,
        "median_boxplot_mean_size_all_fragment_classes.png",
    )
    pct_grid_mono_alt_path = os.path.join(
        args.out_dir,
        (
            "median_boxplot_percent_cfDNA_by_fragment_class_"
            f"mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}.png"
        ),
    )
    pct_combined_mono_alt_path = os.path.join(
        args.out_dir,
        (
            "median_boxplot_percent_cfDNA_all_fragment_classes_"
            f"mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}.png"
        ),
    )
    size_grid_mono_alt_path = os.path.join(
        args.out_dir,
        f"median_boxplot_mean_size_by_fragment_class_mono_{mono_alt_min_bp}_200bp.png",
    )
    size_combined_mono_alt_path = os.path.join(
        args.out_dir,
        f"median_boxplot_mean_size_all_fragment_classes_mono_{mono_alt_min_bp}_200bp.png",
    )
    alt_percent_ylabel = f"% of cfDNA (total >= {mono_alt_min_bp} bp)"

    plot_ks_distribution(
        dist_long,
        major_peak_df,
        peakwise_ks_df,
        ks_plot_path,
        args.size_min,
        args.size_max,
    )

    for _group, _out in [
        ("brca2", brca2_solo_plot_path),
        ("control", ctrl_solo_plot_path),
    ]:
        plot_single_group_distribution(
            dist_long,
            group=_group,
            out_path=_out,
            size_min=args.size_min,
            size_max=args.size_max,
        )

    plot_metric_grid(
        percent_df,
        value_col="percent_cfDNA",
        ylabel="% of total cfDNA",
        title="Group comparison (median-based boxplots): % of fragments by class",
        out_path=pct_grid_path,
        stats_df=mwu_percent_df,
    )
    plot_metric_combined(
        percent_df,
        value_col="percent_cfDNA",
        ylabel="% of total cfDNA",
        title="All fragment classes together (median-based boxplots): % of total cfDNA",
        out_path=pct_combined_path,
    )

    plot_metric_grid(
        mean_size_df,
        value_col="mean_size_bp",
        ylabel="Mean fragment size (bp)",
        title="Group comparison (median-based boxplots): mean fragment size by class",
        out_path=size_grid_path,
        stats_df=mwu_size_df,
    )
    plot_metric_combined(
        mean_size_df,
        value_col="mean_size_bp",
        ylabel="Mean fragment size (bp)",
        title="All fragment classes together (median-based boxplots): mean fragment size",
        out_path=size_combined_path,
    )

    plot_metric_grid(
        percent_df_mono_alt,
        value_col="percent_cfDNA",
        ylabel=alt_percent_ylabel,
        title=(
            "Group comparison (median-based boxplots): % of fragments by class "
            f"(mono = {mono_alt_min_bp}-200 bp; 100% = fragments >= {mono_alt_min_bp} bp)"
        ),
        out_path=pct_grid_mono_alt_path,
        stats_df=mwu_percent_df_mono_alt,
    )
    plot_metric_combined(
        percent_df_mono_alt,
        value_col="percent_cfDNA",
        ylabel=alt_percent_ylabel,
        title=(
            "All fragment classes together (median-based boxplots): % of total cfDNA "
            f"(mono = {mono_alt_min_bp}-200 bp; 100% = fragments >= {mono_alt_min_bp} bp)"
        ),
        out_path=pct_combined_mono_alt_path,
    )

    plot_metric_grid(
        mean_size_df_mono_alt,
        value_col="mean_size_bp",
        ylabel="Mean fragment size (bp)",
        title=(
            "Group comparison (median-based boxplots): mean fragment size by class "
            f"(mono = {mono_alt_min_bp}-200 bp)"
        ),
        out_path=size_grid_mono_alt_path,
        stats_df=mwu_size_df_mono_alt,
    )
    plot_metric_combined(
        mean_size_df_mono_alt,
        value_col="mean_size_bp",
        ylabel="Mean fragment size (bp)",
        title=(
            "All fragment classes together (median-based boxplots): mean fragment size "
            f"(mono = {mono_alt_min_bp}-200 bp)"
        ),
        out_path=size_combined_mono_alt_path,
    )

    ks_df.to_csv(os.path.join(args.out_dir, "ks_peak_size_tests.csv"), index=False)
    peakwise_ks_df.to_csv(os.path.join(args.out_dir, "ks_peak_window_tests.csv"), index=False)
    major_peak_df.to_csv(os.path.join(args.out_dir, "detected_peak_centers.csv"), index=False)
    peak_df.to_csv(os.path.join(args.out_dir, "peak_sizes_by_sample.csv"), index=False)
    local_peak_df.to_csv(os.path.join(args.out_dir, "local_peak_sizes_by_sample.csv"), index=False)
    mwu_percent_df.to_csv(os.path.join(args.out_dir, "mwu_percent_fragment_class_tests.csv"), index=False)
    mwu_size_df.to_csv(os.path.join(args.out_dir, "mwu_mean_size_fragment_class_tests.csv"), index=False)
    mwu_percent_df_mono_alt.to_csv(
        os.path.join(
            args.out_dir,
            (
                "mwu_percent_fragment_class_tests_"
                f"mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}.csv"
            ),
        ),
        index=False,
    )
    mwu_size_df_mono_alt.to_csv(
        os.path.join(
            args.out_dir,
            f"mwu_mean_size_fragment_class_tests_mono_{mono_alt_min_bp}_200bp.csv",
        ),
        index=False,
    )
    summary_percent_df.to_csv(
        os.path.join(args.out_dir, "summary_percent_fragment_class_stats.csv"),
        index=False,
    )
    summary_size_df.to_csv(
        os.path.join(args.out_dir, "summary_mean_size_fragment_class_stats.csv"),
        index=False,
    )
    summary_percent_df_mono_alt.to_csv(
        os.path.join(
            args.out_dir,
            (
                "summary_percent_fragment_class_stats_"
                f"mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}.csv"
            ),
        ),
        index=False,
    )
    summary_size_df_mono_alt.to_csv(
        os.path.join(
            args.out_dir,
            f"summary_mean_size_fragment_class_stats_mono_{mono_alt_min_bp}_200bp.csv",
        ),
        index=False,
    )

    print("\nGenerated files:")
    print(f"  - {ks_plot_path}")
    print(f"  - {brca2_solo_plot_path}")
    print(f"  - {ctrl_solo_plot_path}")
    print(f"  - {pct_grid_path}")
    print(f"  - {pct_combined_path}")
    print(f"  - {size_grid_path}")
    print(f"  - {size_combined_path}")
    print(f"  - {pct_grid_mono_alt_path}")
    print(f"  - {pct_combined_mono_alt_path}")
    print(f"  - {size_grid_mono_alt_path}")
    print(f"  - {size_combined_mono_alt_path}")
    print(f"  - {os.path.join(args.out_dir, 'ks_peak_size_tests.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'ks_peak_window_tests.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'detected_peak_centers.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'peak_sizes_by_sample.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'local_peak_sizes_by_sample.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'mwu_percent_fragment_class_tests.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'mwu_mean_size_fragment_class_tests.csv')}")
    print(
        "  - "
        + os.path.join(
            args.out_dir,
            (
                "mwu_percent_fragment_class_tests_"
                f"mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}.csv"
            ),
        )
    )
    print(
        "  - "
        + os.path.join(
            args.out_dir,
            f"mwu_mean_size_fragment_class_tests_mono_{mono_alt_min_bp}_200bp.csv",
        )
    )
    print(f"  - {os.path.join(args.out_dir, 'summary_percent_fragment_class_stats.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'summary_mean_size_fragment_class_stats.csv')}")
    print(
        "  - "
        + os.path.join(
            args.out_dir,
            (
                "summary_percent_fragment_class_stats_"
                f"mono_{mono_alt_min_bp}_200bp_rebased_total_ge{mono_alt_min_bp}.csv"
            ),
        )
    )
    print(
        "  - "
        + os.path.join(
            args.out_dir,
            f"summary_mean_size_fragment_class_stats_mono_{mono_alt_min_bp}_200bp.csv",
        )
    )


if __name__ == "__main__":
    main()
