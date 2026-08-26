#!/usr/bin/env python3
"""
Genome-wide filler-region methylation PERMANOVA/PERMDISP for BRCA2/control cfDNA.

Groups are inferred from BED file name prefixes:
- ctrl-*  -> control
- brca2-* -> brca2

Pipeline
1) Build non-coding filler regions as chromosome-wise gaps between merged GENCODE CDS intervals.
2) Aggregate methylated/total CpG counts per filler region per sample via bedtools intersect.
3) Build sample x filler-region methylation matrix (percent methylation).
4) Filter filler regions by per-group minimum observed samples, median-impute remaining sparse NAs.
5) Run global and pairwise PERMANOVA (Euclidean on arcsin-sqrt transformed methylation).
6) Run global and pairwise PERMDISP (distance-to-centroid permutation ANOVA).
7) Export tables + simple QC plots.

Notes
- Filler-region coordinates are built from CDS gaps within each chromosome's GTF-annotated span.
- CDS intervals are merged first; each uncovered gap becomes one filler region.
- Regions are labeled chrN:filler_1, chrN:filler_2, ... in genomic order.
- The SVD-based coordinate compression is exact for Euclidean distances between samples,
    and makes permutation testing much faster when regions >> samples.
- PERMANOVA significance should be interpreted alongside PERMDISP.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.stats import kruskal, mannwhitneyu


DEFAULT_GENCODE_GTF = "/mnt/c/Users/gissu/Downloads/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
DEFAULT_BED_DIR = "/mnt/d/BRCA2-misc-files/BRCA2-ctrl-mod"
DEFAULT_OUTPUT_DIR = "/mnt/d/BRCA2-misc-files/mod_BED_coronary/permanova_filler_output"

DEFAULT_MIN_CPGS = 3
DEFAULT_MIN_SAMPLES = 4
DEFAULT_PERMUTATIONS = 10000
DEFAULT_RANDOM_SEED = 42

# bedMethyl field indices (0-based)
MOD_BASE_CODE_FIELD = 3
NVALID_COV_FIELD = 9
NMOD_5MC_FIELD = 11

DEFAULT_GROUP_ALIASES = {
    "ctrl": "control",
    "brca2": "brca2",
}
DEFAULT_GROUP_LABELS = {
    "control": "Control",
    "brca2": "BRCA2-999del5",
}
DEFAULT_GROUP_ORDER = ["brca2", "control"]


def format_comparison_label(comp: str, group_labels: Dict[str, str]) -> str:
    parts = comp.split("_vs_")
    if len(parts) < 2:
        return comp
    return " vs ".join(group_labels.get(p, p) for p in parts)


def parse_key_value_arg(raw: str, lower_values: bool = True) -> Dict[str, str]:
    """Parse comma-separated key:value entries into a dict."""
    out: Dict[str, str] = {}
    if not raw:
        return out
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid mapping item '{item}'. Expected key:value format.")
        key, value = item.split(":", 1)
        key = key.strip().lower()
        value = value.strip().lower() if lower_values else value.strip()
        if not key or not value:
            raise ValueError(f"Invalid mapping item '{item}'. Empty key or value.")
        out[key] = value
    return out


def parse_group_order_arg(raw: str) -> List[str]:
    if not raw:
        return []
    return [x.strip().lower() for x in raw.split(",") if x.strip()]


@dataclass
class SampleInfo:
    sample_id: str
    barcode: str
    group: str
    bed_file: str


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_field.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


class GencodeFillerRegionCatalog:
    """Build non-coding filler regions from gaps between merged GENCODE CDS intervals."""

    def __init__(self, gtf_path: str):
        self.gtf_path = Path(gtf_path)
        self._catalog: Optional[List[Dict[str, object]]] = None

    def build_catalog(self) -> List[Dict[str, object]]:
        if self._catalog is not None:
            return self._catalog

        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GENCODE GTF not found: {self.gtf_path}")

        opener = gzip.open if self.gtf_path.suffix == ".gz" else open
        cds_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
        chrom_bounds: Dict[str, Tuple[int, int]] = {}

        with opener(self.gtf_path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                chrom, _, feature, start, end, _, strand, _, attrs = parts
                try:
                    start_i = int(start)
                    end_i = int(end)
                except ValueError:
                    continue

                chrom_name = chrom if chrom.startswith("chr") else f"chr{chrom}"

                if chrom_name not in chrom_bounds:
                    chrom_bounds[chrom_name] = (start_i, end_i)
                else:
                    lo, hi = chrom_bounds[chrom_name]
                    chrom_bounds[chrom_name] = (min(lo, start_i), max(hi, end_i))

                if feature != "CDS":
                    continue

                cds_by_chrom.setdefault(chrom_name, []).append((start_i, end_i))

        regions: List[Dict[str, object]] = []
        for chrom_name in sorted(chrom_bounds.keys()):
            chrom_min, chrom_max = chrom_bounds[chrom_name]
            cds_intervals = sorted(set(cds_by_chrom.get(chrom_name, [])), key=lambda x: (x[0], x[1]))

            merged: List[Tuple[int, int]] = []
            for start_i, end_i in cds_intervals:
                if not merged or start_i > merged[-1][1] + 1:
                    merged.append((start_i, end_i))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end_i))

            filler_idx = 1
            cursor = chrom_min
            for start_i, end_i in merged:
                if start_i > cursor:
                    regions.append(
                        {
                            "region_id": f"{chrom_name}:filler_{filler_idx}",
                            "gene_name": "FILLER",
                            "region_index": filler_idx,
                            "gene_id": "FILLER",
                            "chromosome": chrom_name,
                            "strand": ".",
                            "region_start": cursor,
                            "region_end": start_i - 1,
                        }
                    )
                    filler_idx += 1
                cursor = max(cursor, end_i + 1)

            if cursor <= chrom_max:
                regions.append(
                    {
                        "region_id": f"{chrom_name}:filler_{filler_idx}",
                        "gene_name": "FILLER",
                        "region_index": filler_idx,
                        "gene_id": "FILLER",
                        "chromosome": chrom_name,
                        "strand": ".",
                        "region_start": cursor,
                        "region_end": chrom_max,
                    }
                )

        if not regions:
            raise RuntimeError("No filler regions could be derived from the provided GTF.")

        regions.sort(
            key=lambda r: (
                str(r["chromosome"]),
                int(r["region_start"]),
                int(r["region_end"]),
                str(r["gene_name"]),
                int(r["region_index"]),
            )
        )

        self._catalog = regions
        return self._catalog


def aggregate_sample_with_bedtools(task: Tuple[str, str, int]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate per-region 5mC/5hmC and total CpG counts for one sample using bedtools."""
    bed_file, region_bed, gene_count = task
    meth_5mc_counts = np.zeros(gene_count, dtype=np.int64)
    meth_5hmc_counts = np.zeros(gene_count, dtype=np.int64)
    total_counts = np.zeros(gene_count, dtype=np.int64)

    cmd = [
        "bedtools",
        "intersect",
        "-a",
        bed_file,
        "-b",
        region_bed,
        "-wa",
        "-wb",
        "-f",
        "1.0",
    ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 17:
            continue
        try:
            mod_code = str(parts[MOD_BASE_CODE_FIELD]).strip().lower()
            total_cpg = int(float(parts[NVALID_COV_FIELD]))
            nmod = int(float(parts[NMOD_5MC_FIELD]))
            gene_idx = int(parts[-1])
        except ValueError:
            continue

        if 0 <= gene_idx < gene_count:
            if mod_code == "m":
                meth_5mc_counts[gene_idx] += nmod
                total_counts[gene_idx] += total_cpg
            elif mod_code == "h":
                meth_5hmc_counts[gene_idx] += nmod

    stderr_text = proc.stderr.read() if proc.stderr is not None else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"bedtools intersect failed for {bed_file}: {stderr_text.strip()}")

    return meth_5mc_counts, meth_5hmc_counts, total_counts


def discover_samples(bed_dir: Path, group_aliases: Dict[str, str]) -> List[SampleInfo]:
    samples: List[SampleInfo] = []
    bed_files = sorted(bed_dir.glob("*barcode*.bed"))
    if not bed_files:
        raise FileNotFoundError(f"No .bed files found in {bed_dir}")

    for bed_file in bed_files:
        m = re.match(r"([A-Za-z0-9]+)[_\-]barcode(\d+)(?:[_\-]\d+)?$", bed_file.stem, re.IGNORECASE)
        if not m:
            continue
        prefix = m.group(1).lower()
        barcode = m.group(2)
        group = group_aliases.get(prefix, prefix)

        sample_id = f"{group}-barcode{barcode}"
        samples.append(SampleInfo(sample_id=sample_id, barcode=barcode, group=group, bed_file=str(bed_file)))

    if not samples:
        raise RuntimeError("No files matched expected prefixes ctrl/brca2 with barcode naming.")

    return samples


def write_region_bed(outdir: Path, region_items: List[Dict[str, object]]) -> Path:
    out_path = outdir / "_tmp_regions_permanova.bed"
    with open(out_path, "w", encoding="utf-8") as handle:
        for idx, info in enumerate(region_items):
            chrom = str(info["chromosome"])
            # Convert 1-based closed GTF coordinates to 0-based half-open BED.
            start = max(0, int(info["region_start"]) - 1)
            end = int(info["region_end"])
            if end <= start:
                continue
            handle.write(f"{chrom}\t{start}\t{end}\tregion_{idx}\t{idx}\n")
    return out_path


def aggregate_all_samples(
    samples: Sequence[SampleInfo],
    region_bed: Path,
    gene_count: int,
    max_workers: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    tasks = [(s.bed_file, str(region_bed), gene_count) for s in samples]
    meth_5mc_matrix = np.zeros((len(samples), gene_count), dtype=np.int64)
    meth_5hmc_matrix = np.zeros((len(samples), gene_count), dtype=np.int64)
    total_matrix = np.zeros((len(samples), gene_count), dtype=np.int64)

    workers = max(1, min(max_workers, len(samples)))
    print(f"Aggregating {len(samples)} samples with {workers} workers")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_idx = {
            executor.submit(aggregate_sample_with_bedtools, task): i for i, task in enumerate(tasks)
        }
        done = 0
        for future in as_completed(future_to_idx):
            i = future_to_idx[future]
            meth_5mc_counts, meth_5hmc_counts, total_counts = future.result()
            meth_5mc_matrix[i, :] = meth_5mc_counts
            meth_5hmc_matrix[i, :] = meth_5hmc_counts
            total_matrix[i, :] = total_counts
            done += 1
            print(f"  aggregated {done}/{len(samples)}: {Path(samples[i].bed_file).name}")

    return meth_5mc_matrix, meth_5hmc_matrix, total_matrix


def compute_percent_matrix(
    mod_counts: np.ndarray,
    total_matrix: np.ndarray,
    min_cpgs: int,
    keep_mask: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    pct = np.full(mod_counts.shape, np.nan, dtype=float)
    valid = total_matrix >= min_cpgs
    pct[valid] = (mod_counts[valid] / total_matrix[valid]) * 100.0
    return pct[:, keep_mask], total_matrix[:, keep_mask]


def run_metric_outputs(
    metric_key: str,
    metric_title: str,
    metric_pct: np.ndarray,
    mod_counts_kept: np.ndarray,
    total_kept: np.ndarray,
    groups: Sequence[str],
    sample_ids: Sequence[str],
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    pairwise: Sequence[Tuple[str, str]],
    outdir: Path,
    permutations: int,
    rng: np.random.Generator,
    skip_spearman: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metric_imputed, na_fraction = median_impute_by_region(metric_pct)

    eps = 1e-6
    beta = np.clip(metric_imputed / 100.0, eps, 1.0 - eps)
    transformed = np.arcsin(np.sqrt(beta))
    stats_available = True
    try:
        coords = to_euclidean_preserving_coords(transformed)
    except RuntimeError as exc:
        if str(exc) != "All samples are numerically identical after preprocessing.":
            raise
        stats_available = False
        coords = np.zeros((len(sample_ids), 1), dtype=float)
        print(
            f"WARNING: {metric_title} values are identical across all samples after preprocessing; "
            "writing placeholder statistical outputs."
        )

    total_all = total_kept.sum(axis=1).astype(float)
    mod_all = mod_counts_kept.sum(axis=1).astype(float)
    global_pct = np.divide(mod_all, total_all, out=np.zeros_like(mod_all), where=total_all > 0) * 100.0

    embedding_df = pd.DataFrame({"sample_id": sample_ids, "group": groups, "axis1": coords[:, 0]})
    if coords.shape[1] > 1:
        embedding_df["axis2"] = coords[:, 1]
    else:
        embedding_df["axis2"] = 0.0
    embedding_df.to_csv(outdir / f"sample_embedding_{metric_key}.csv", index=False)

    plot_sample_scatter(
        coords,
        groups,
        sample_ids,
        group_order,
        group_labels,
        metric_title,
        outdir / f"sample_embedding_{metric_key}.png",
    )

    global_comparison_name = "_vs_".join(group_order)
    if stats_available:
        global_perm = permanova(coords, groups, permutations, rng)
        global_disp = permdisp(coords, groups, permutations, rng)

        permanova_rows: List[Dict[str, object]] = [
            {
                "comparison": global_comparison_name,
                "test_scope": "global",
                **global_perm,
            }
        ]
        permdisp_rows: List[Dict[str, object]] = [
            {
                "comparison": global_comparison_name,
                "test_scope": "global",
                **global_disp,
            }
        ]

        for g1, g2 in pairwise:
            idx = [i for i, g in enumerate(groups) if g in (g1, g2)]
            sub_coords = coords[idx, :]
            sub_groups = [groups[i] for i in idx]

            p_res = permanova(sub_coords, sub_groups, permutations, rng)
            d_res = permdisp(sub_coords, sub_groups, permutations, rng)

            permanova_rows.append(
                {
                    "comparison": f"{g1}_vs_{g2}",
                    "test_scope": "pairwise",
                    **p_res,
                }
            )
            permdisp_rows.append(
                {
                    "comparison": f"{g1}_vs_{g2}",
                    "test_scope": "pairwise",
                    **d_res,
                }
            )
    else:
        permanova_rows = [
            {
                "comparison": global_comparison_name,
                "test_scope": "global",
                "pseudo_f": np.nan,
                "r_squared": np.nan,
                "p_value": np.nan,
                "n_samples": int(len(groups)),
                "n_groups": int(len(np.unique(groups))),
                "permutations": int(permutations),
            }
        ]
        permdisp_rows = [
            {
                "comparison": global_comparison_name,
                "test_scope": "global",
                "f_stat": np.nan,
                "p_value": np.nan,
                "n_samples": int(len(groups)),
                "n_groups": int(len(np.unique(groups))),
                "permutations": int(permutations),
            }
        ]
        for g1, g2 in pairwise:
            sub_groups = [g for g in groups if g in (g1, g2)]
            permanova_rows.append(
                {
                    "comparison": f"{g1}_vs_{g2}",
                    "test_scope": "pairwise",
                    "pseudo_f": np.nan,
                    "r_squared": np.nan,
                    "p_value": np.nan,
                    "n_samples": int(len(sub_groups)),
                    "n_groups": int(len(np.unique(sub_groups))),
                    "permutations": int(permutations),
                }
            )
            permdisp_rows.append(
                {
                    "comparison": f"{g1}_vs_{g2}",
                    "test_scope": "pairwise",
                    "f_stat": np.nan,
                    "p_value": np.nan,
                    "n_samples": int(len(sub_groups)),
                    "n_groups": int(len(np.unique(sub_groups))),
                    "permutations": int(permutations),
                }
            )

    permanova_df = pd.DataFrame(permanova_rows)
    permdisp_df = pd.DataFrame(permdisp_rows)

    pair_mask = permanova_df["test_scope"] == "pairwise"
    permanova_df.loc[pair_mask, "p_adj_bh"] = bh_fdr(permanova_df.loc[pair_mask, "p_value"].tolist())
    permanova_df.loc[~pair_mask, "p_adj_bh"] = np.nan

    permanova_df.to_csv(outdir / f"genomewide_{metric_key}_permanova_results.csv", index=False)
    permdisp_df.to_csv(outdir / f"genomewide_{metric_key}_permdisp_results.csv", index=False)

    readable_table = build_readable_summary_table(permanova_df, permdisp_df, group_labels)
    readable_table.to_csv(outdir / f"genomewide_{metric_key}_permanova_readable_table.csv", index=False)
    plot_readable_summary_table(readable_table, outdir / f"genomewide_{metric_key}_permanova_readable_table.png")
    write_readable_summary_html(readable_table, outdir / f"genomewide_{metric_key}_permanova_readable_table.html")

    sample_metric_df = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "group": groups,
            f"global_{metric_key}_pct": global_pct,
        }
    )
    sample_metric_df.to_csv(outdir / f"sample_global_{metric_key}.csv", index=False)

    plot_global_boxplot(
        sample_metric_df.rename(columns={f"global_{metric_key}_pct": "global_methylation_pct"}),
        group_order,
        group_labels,
        metric_title,
        outdir / f"global_{metric_key}_boxplot.png",
    )

    group_arrays = {
        g: sample_metric_df.loc[sample_metric_df["group"] == g, f"global_{metric_key}_pct"].to_numpy(dtype=float)
        for g in group_order
    }
    kw_stat, kw_p = kruskal(*(group_arrays[g] for g in group_order))
    uw_rows: List[Dict[str, object]] = [
        {
            "test": f"kruskal_global_{metric_key}_{len(group_order)}group",
            "statistic": float(kw_stat),
            "p_value": float(kw_p),
        }
    ]
    for g1, g2 in pairwise:
        stat, p = mannwhitneyu(group_arrays[g1], group_arrays[g2], alternative="two-sided")
        uw_rows.append({"test": f"mannwhitney_{metric_key}_{g1}_vs_{g2}", "statistic": float(stat), "p_value": float(p)})

    uw_df = pd.DataFrame(uw_rows)
    mw_mask = uw_df["test"].str.startswith("mannwhitney_")
    uw_df.loc[mw_mask, "p_adj_bh"] = bh_fdr(uw_df.loc[mw_mask, "p_value"].tolist())
    uw_df.to_csv(outdir / f"global_{metric_key}_univariate_tests.csv", index=False)

    if not skip_spearman:
        write_spearman_outputs(
            meth_matrix=mod_counts_kept,
            total_matrix=total_kept,
            groups=groups,
            group_order=group_order,
            group_labels=group_labels,
            min_cpgs=1,
            outdir=outdir,
            output_prefix=metric_key,
            metric_label=metric_title,
        )

    qc_df = pd.DataFrame(
        {
            "metric": [metric_title],
            "mean_missing_fraction_before_impute": [float(np.mean(na_fraction))],
            "median_missing_fraction_before_impute": [float(np.median(na_fraction))],
        }
    )
    qc_df.to_csv(outdir / f"filler_region_matrix_qc_{metric_key}.csv", index=False)

    return permanova_df, permdisp_df, sample_metric_df, metric_imputed


def region_level_kruskal_tests(
    metric_pct: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    region_ids: Sequence[str],
    region_genes: Sequence[str],
    region_indices: Sequence[int],
) -> pd.DataFrame:
    groups_arr = np.asarray(groups)
    rows: List[Dict[str, object]] = []

    for j in range(metric_pct.shape[1]):
        arrays = []
        valid_gene = True
        for g in group_order:
            vals = metric_pct[groups_arr == g, j]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                valid_gene = False
                break
            arrays.append(vals)
        if not valid_gene:
            rows.append(
                {
                    "region_id": region_ids[j],
                    "gene": region_genes[j],
                    "region_index": int(region_indices[j]),
                    "region_matrix_index": j,
                    "statistic": np.nan,
                    "p_value": np.nan,
                }
            )
            continue
        try:
            stat, p_value = kruskal(*arrays)
        except ValueError:
            stat, p_value = np.nan, np.nan
        rows.append(
            {
                "region_id": region_ids[j],
                "gene": region_genes[j],
                "region_index": int(region_indices[j]),
                "region_matrix_index": j,
                "statistic": float(stat),
                "p_value": float(p_value),
            }
        )

    out = pd.DataFrame(rows)
    valid_mask = out["p_value"].notna()
    out["q_value_bh"] = np.nan
    if valid_mask.any():
        out.loc[valid_mask, "q_value_bh"] = bh_fdr(out.loc[valid_mask, "p_value"].to_numpy(dtype=float))
    return out


def region_pairwise_mannwhitney_tests(
    metric_pct: np.ndarray,
    groups: Sequence[str],
    pairwise: Sequence[Tuple[str, str]],
    region_ids: Sequence[str],
    region_genes: Sequence[str],
    region_indices: Sequence[int],
) -> pd.DataFrame:
    groups_arr = np.asarray(groups)
    rows: List[Dict[str, object]] = []

    for g1, g2 in pairwise:
        comp = f"{g1}_vs_{g2}"
        idx1 = np.where(groups_arr == g1)[0]
        idx2 = np.where(groups_arr == g2)[0]
        for j in range(metric_pct.shape[1]):
            vals1 = metric_pct[idx1, j]
            vals2 = metric_pct[idx2, j]
            vals1 = vals1[np.isfinite(vals1)]
            vals2 = vals2[np.isfinite(vals2)]

            if vals1.size == 0 or vals2.size == 0:
                stat, p_value = np.nan, np.nan
            else:
                try:
                    stat, p_value = mannwhitneyu(vals1, vals2, alternative="two-sided")
                    stat, p_value = float(stat), float(p_value)
                except ValueError:
                    stat, p_value = np.nan, np.nan

            rows.append(
                {
                    "comparison": comp,
                    "region_id": region_ids[j],
                    "gene": region_genes[j],
                    "region_index": int(region_indices[j]),
                    "region_matrix_index": j,
                    "statistic": stat,
                    "p_value": p_value,
                }
            )

    out = pd.DataFrame(rows)
    out["q_value_bh"] = np.nan
    for comp in out["comparison"].unique().tolist():
        mask = (out["comparison"] == comp) & out["p_value"].notna()
        if mask.any():
            out.loc[mask, "q_value_bh"] = bh_fdr(out.loc[mask, "p_value"].to_numpy(dtype=float))
    return out


def gene_calls_from_regions(
    region_pairwise_df: pd.DataFrame,
    q_threshold: float = 0.05,
    fraction_threshold: float = 0.5,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    if region_pairwise_df.empty:
        return pd.DataFrame(
            columns=[
                "comparison",
                "gene",
                "n_regions_total",
                "n_regions_significant",
                "fraction_significant",
                "is_gene_differential",
                "significant_regions",
            ]
        )

    for comp in sorted(region_pairwise_df["comparison"].dropna().unique().tolist()):
        sub = region_pairwise_df.loc[region_pairwise_df["comparison"] == comp].copy()
        for gene, gene_sub in sub.groupby("gene", sort=True):
            total_regions = int(gene_sub["region_id"].nunique())
            sig_sub = gene_sub.loc[gene_sub["q_value_bh"] < q_threshold].copy()
            n_sig = int(sig_sub["region_id"].nunique())
            frac = float(n_sig / total_regions) if total_regions > 0 else 0.0
            sig_regions = sig_sub.sort_values("region_index")["region_id"].drop_duplicates().tolist()
            rows.append(
                {
                    "comparison": comp,
                    "gene": gene,
                    "n_regions_total": total_regions,
                    "n_regions_significant": n_sig,
                    "fraction_significant": frac,
                    "is_gene_differential": bool(frac > fraction_threshold),
                    "significant_regions": ";".join(sig_regions),
                }
            )

    return pd.DataFrame(rows)


def export_significant_gene_table(
    test_df: pd.DataFrame,
    out_html: Path,
    out_png: Path,
    title: str,
    q_threshold: float = 0.05,
) -> pd.DataFrame:
    """Export significant-gene summaries as HTML and PNG tables."""
    has_comparison = "comparison" in test_df.columns
    has_region = "region_id" in test_df.columns

    required_cols = ["gene", "p_value", "q_value_bh"]
    for col in required_cols:
        if col not in test_df.columns:
            raise ValueError(f"Missing required column '{col}' in test dataframe.")

    sig = test_df.loc[test_df["q_value_bh"].notna() & (test_df["q_value_bh"] < q_threshold)].copy()

    if sig.empty:
        summary_cols = ["gene", "n_significant_tests", "min_p_value", "min_q_value_bh"]
        if has_comparison:
            summary_cols.extend(["n_significant_comparisons", "comparisons"])
        if has_region:
            summary_cols.extend(["n_significant_regions", "significant_regions"])
        summary = pd.DataFrame(columns=summary_cols)
    else:
        rows: List[Dict[str, object]] = []
        for gene, sub in sig.groupby("gene", sort=True):
            row: Dict[str, object] = {
                "gene": str(gene),
                "n_significant_tests": int(sub.shape[0]),
                "min_p_value": float(sub["p_value"].min()),
                "min_q_value_bh": float(sub["q_value_bh"].min()),
            }
            if has_comparison:
                comps = sorted(set(sub["comparison"].dropna().astype(str).tolist()))
                row["n_significant_comparisons"] = int(len(comps))
                row["comparisons"] = "; ".join(comps)
            if has_region:
                regions = sorted(set(sub["region_id"].dropna().astype(str).tolist()))
                row["n_significant_regions"] = int(len(regions))
                row["significant_regions"] = "; ".join(regions)
            rows.append(row)

        summary = pd.DataFrame(rows).sort_values(["min_q_value_bh", "gene"]).reset_index(drop=True)

    fmt_map = {
        "min_p_value": lambda v: "" if pd.isna(v) else f"{float(v):.3e}",
        "min_q_value_bh": lambda v: "" if pd.isna(v) else f"{float(v):.3e}",
    }
    html_df = summary.copy()
    for col, formatter in fmt_map.items():
        if col in html_df.columns:
            html_df[col] = html_df[col].map(formatter)

    html_parts = [
        "<html><head><meta charset='utf-8'>",
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 24px; }",
        "h2 { margin-bottom: 8px; }",
        "p { margin-top: 0; color: #444; }",
        "table { border-collapse: collapse; width: 100%; font-size: 12px; }",
        "th, td { border: 1px solid #ddd; padding: 6px 8px; text-align: left; vertical-align: top; }",
        "th { background: #f2f2f2; }",
        "tr:nth-child(even) { background: #fafafa; }",
        "</style></head><body>",
        f"<h2>{title}</h2>",
        f"<p>q-value threshold: {q_threshold:.3g} | Significant genes: {int(summary.shape[0])}</p>",
        html_df.to_html(index=False, escape=True),
        "</body></html>",
    ]
    out_html.write_text("\n".join(html_parts), encoding="utf-8")

    fig_w = min(22.0, max(10.0, 1.6 * max(1, summary.shape[1])))
    fig_h = min(24.0, max(2.8, 0.38 * (max(1, summary.shape[0]) + 2)))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")

    if summary.empty:
        ax.text(0.5, 0.5, f"No significant genes (q < {q_threshold:.3g})", ha="center", va="center", fontsize=12)
    else:
        display_df = summary.copy()
        for col, formatter in fmt_map.items():
            if col in display_df.columns:
                display_df[col] = display_df[col].map(formatter)
        if display_df.shape[0] > 120:
            display_df = display_df.head(120)

        table = ax.table(
            cellText=display_df.values,
            colLabels=display_df.columns,
            cellLoc="left",
            colLoc="left",
            loc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.15)

    ax.set_title(title, fontsize=12, pad=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250, bbox_inches="tight")
    plt.close(fig)

    return summary


def plot_significant_gene_stacked_composition(
    meth_5mc: np.ndarray,
    meth_5hmc: np.ndarray,
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    significant_idx: np.ndarray,
    out_png: Path,
    out_csv: Path,
) -> None:
    if significant_idx.size == 0:
        pd.DataFrame(
            columns=["group", "label", "pct_5mc", "pct_5hmc", "pct_unmodified", "n_significant_regions"]
        ).to_csv(out_csv, index=False)
        return

    groups_arr = np.asarray(groups)
    rows: List[Dict[str, object]] = []
    for g in group_order:
        idx = np.where(groups_arr == g)[0]
        if idx.size == 0:
            continue
        total = float(total_matrix[idx][:, significant_idx].sum())
        m5 = float(meth_5mc[idx][:, significant_idx].sum())
        hm5 = float(meth_5hmc[idx][:, significant_idx].sum())

        if total <= 0:
            pct_5mc = 0.0
            pct_5hmc = 0.0
            pct_unmod = 0.0
        else:
            pct_5mc = 100.0 * m5 / total
            pct_5hmc = 100.0 * hm5 / total
            pct_unmod = max(0.0, 100.0 - pct_5mc - pct_5hmc)

        rows.append(
            {
                "group": g,
                "label": group_labels.get(g, g),
                "pct_5mc": pct_5mc,
                "pct_5hmc": pct_5hmc,
                "pct_unmodified": pct_unmod,
                "n_significant_regions": int(significant_idx.size),
            }
        )

    comp_df = pd.DataFrame(rows)
    comp_df.to_csv(out_csv, index=False)
    if comp_df.empty:
        return

    x = np.arange(len(comp_df))
    plt.figure(figsize=(10, 6))
    plt.bar(x, comp_df["pct_5mc"], label="5mC", color="#d62728", edgecolor="black", linewidth=0.4)
    plt.bar(
        x,
        comp_df["pct_5hmc"],
        bottom=comp_df["pct_5mc"],
        label="5hmC",
        color="#1f77b4",
        edgecolor="black",
        linewidth=0.4,
    )
    plt.bar(
        x,
        comp_df["pct_unmodified"],
        bottom=comp_df["pct_5mc"] + comp_df["pct_5hmc"],
        label="Unmodified CpG",
        color="#bdbdbd",
        edgecolor="black",
        linewidth=0.4,
    )

    plt.xticks(x, comp_df["label"].tolist())
    plt.ylabel("CpG composition (%)")
    plt.title(
        "Significant coding-region CpG composition by group\n"
        f"(regions significant in >=1 group comparison for 5mC or 5hmC; n={int(significant_idx.size)})"
    )
    plt.ylim(0, 100)
    plt.grid(axis="y", alpha=0.25)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def filter_regions_by_group_coverage(
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    min_cpgs: int,
    min_samples: int,
) -> np.ndarray:
    groups_arr = np.asarray(groups)
    keep = np.ones(total_matrix.shape[1], dtype=bool)
    for group in group_order:
        idx = np.where(groups_arr == group)[0]
        if len(idx) == 0:
            continue
        observed = (total_matrix[idx, :] >= min_cpgs).sum(axis=0)
        keep &= observed >= min_samples
    return keep


def median_impute_by_region(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Median-impute NaNs per region; return imputed matrix and per-region NA fraction."""
    x = values.copy()
    na_mask = np.isnan(x)
    na_fraction = na_mask.mean(axis=0)

    medians = np.nanmedian(x, axis=0)
    medians = np.where(np.isnan(medians), 0.0, medians)
    rows, cols = np.where(na_mask)
    x[rows, cols] = medians[cols]
    return x, na_fraction


def to_euclidean_preserving_coords(x: np.ndarray) -> np.ndarray:
    """Project samples to low dimension that preserves pairwise Euclidean distances exactly."""
    xc = x - x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(xc, full_matrices=False)
    tol = 1e-12
    r = int(np.sum(s > tol))
    if r == 0:
        raise RuntimeError("All samples are numerically identical after preprocessing.")
    return u[:, :r] * s[:r]


def _ss_from_labels(coords: np.ndarray, labels: np.ndarray) -> Tuple[float, float, int, int]:
    unique = np.unique(labels)
    n = coords.shape[0]
    k = len(unique)
    if k < 2:
        raise ValueError("Need at least two groups.")
    if n <= k:
        raise ValueError("Need more samples than groups for pseudo-F denominator.")

    grand_centroid = coords.mean(axis=0)
    ss_between = 0.0
    ss_within = 0.0
    for g in unique:
        idx = labels == g
        group_coords = coords[idx]
        centroid = group_coords.mean(axis=0)
        ss_between += group_coords.shape[0] * float(np.sum((centroid - grand_centroid) ** 2))
        ss_within += float(np.sum((group_coords - centroid) ** 2))

    return ss_between, ss_within, n, k


def permanova(coords: np.ndarray, labels: Sequence[str], permutations: int, rng: np.random.Generator) -> Dict[str, float]:
    labels_arr = np.asarray(labels)
    ss_between, ss_within, n, k = _ss_from_labels(coords, labels_arr)
    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n - k)
    pseudo_f = np.inf if ms_within == 0 else ms_between / ms_within
    r2 = ss_between / (ss_between + ss_within) if (ss_between + ss_within) > 0 else 0.0

    if permutations <= 0:
        p_value = 1.0
    else:
        count = 0
        for _ in range(permutations):
            perm_labels = rng.permutation(labels_arr)
            p_ss_between, p_ss_within, p_n, p_k = _ss_from_labels(coords, perm_labels)
            p_ms_between = p_ss_between / (p_k - 1)
            p_ms_within = p_ss_within / (p_n - p_k)
            p_f = np.inf if p_ms_within == 0 else p_ms_between / p_ms_within
            if p_f >= pseudo_f:
                count += 1
        p_value = (count + 1) / (permutations + 1)

    return {
        "pseudo_f": float(pseudo_f),
        "r_squared": float(r2),
        "p_value": float(p_value),
        "n_samples": int(n),
        "n_groups": int(k),
        "permutations": int(permutations),
    }


def _one_way_anova_f(values: np.ndarray, labels: np.ndarray) -> Tuple[float, int, int]:
    unique = np.unique(labels)
    n = values.shape[0]
    k = len(unique)
    if k < 2:
        raise ValueError("Need at least two groups.")
    if n <= k:
        raise ValueError("Need more samples than groups for ANOVA denominator.")

    grand_mean = float(np.mean(values))
    ss_between = 0.0
    ss_within = 0.0
    for g in unique:
        group_vals = values[labels == g]
        group_mean = float(np.mean(group_vals))
        ss_between += len(group_vals) * (group_mean - grand_mean) ** 2
        ss_within += float(np.sum((group_vals - group_mean) ** 2))

    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n - k)
    f_stat = np.inf if ms_within == 0 else ms_between / ms_within
    return float(f_stat), n, k


def permdisp(coords: np.ndarray, labels: Sequence[str], permutations: int, rng: np.random.Generator) -> Dict[str, float]:
    labels_arr = np.asarray(labels)
    unique = np.unique(labels_arr)

    distances = np.zeros(coords.shape[0], dtype=float)
    mean_distance_by_group: Dict[str, float] = {}
    for g in unique:
        idx = labels_arr == g
        centroid = coords[idx].mean(axis=0)
        d = np.linalg.norm(coords[idx] - centroid, axis=1)
        distances[idx] = d
        mean_distance_by_group[str(g)] = float(np.mean(d))

    observed_f, n, k = _one_way_anova_f(distances, labels_arr)

    if permutations <= 0:
        p_value = 1.0
    else:
        count = 0
        for _ in range(permutations):
            perm_labels = rng.permutation(labels_arr)
            perm_f, _, _ = _one_way_anova_f(distances, perm_labels)
            if perm_f >= observed_f:
                count += 1
        p_value = (count + 1) / (permutations + 1)

    result: Dict[str, float] = {
        "f_stat": float(observed_f),
        "p_value": float(p_value),
        "n_samples": int(n),
        "n_groups": int(k),
        "permutations": int(permutations),
    }
    for g in unique:
        result[f"mean_distance_{str(g)}"] = mean_distance_by_group[str(g)]
    return result


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    sorted_p = p[order]
    adj = sorted_p * n / np.arange(1, n + 1)
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    adj = np.minimum(adj, 1.0)
    out = np.zeros(n, dtype=float)
    out[order] = adj
    return out


def _group_palette(group_order: Sequence[str], group_labels: Dict[str, str]) -> Dict[str, str]:
    fixed = {
        "control": "#1f77b4",
        "brca2": "#d62728",
    }
    cmap = plt.get_cmap("tab20")
    palette: Dict[str, str] = {}
    for i, g in enumerate(group_order):
        palette[g] = fixed.get(g, cmap(i % 20))
    return palette


def plot_sample_scatter(
    coords: np.ndarray,
    groups: Sequence[str],
    sample_ids: Sequence[str],
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    metric_label: str,
    out_png: Path,
) -> None:
    groups_arr = np.asarray(groups)
    if coords.shape[1] == 1:
        xs = coords[:, 0]
        ys = np.zeros(coords.shape[0], dtype=float)
        xlabel = "Axis 1"
        ylabel = "Axis 2 (constant)"
    else:
        xs = coords[:, 0]
        ys = coords[:, 1]
        xlabel = "Axis 1"
        ylabel = "Axis 2"

    palette = _group_palette(group_order, group_labels)
    plt.figure(figsize=(8, 6))
    for g in group_order:
        idx = np.where(groups_arr == g)[0]
        if len(idx) == 0:
            continue
        plt.scatter(
            xs[idx],
            ys[idx],
            s=80,
            alpha=0.9,
            label=group_labels.get(g, g),
            color=palette[g],
            edgecolors="black",
            linewidths=0.4,
        )
        for i in idx:
            plt.text(xs[i], ys[i], sample_ids[i], fontsize=7, alpha=0.8)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"Genome-wide coding-region {metric_label} sample embedding")
    plt.grid(alpha=0.2)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def plot_global_boxplot(
    sample_df: pd.DataFrame,
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    metric_label: str,
    out_png: Path,
) -> None:
    palette = _group_palette(group_order, group_labels)
    order = [g for g in group_order if g in set(sample_df["group"])]

    plt.figure(figsize=(7, 5))
    positions = np.arange(len(order))
    data = [sample_df.loc[sample_df["group"] == g, "global_methylation_pct"].to_numpy() for g in order]
    labels = [group_labels.get(g, g) for g in order]

    bp = plt.boxplot(data, positions=positions, widths=0.6, patch_artist=True)
    for patch, g in zip(bp["boxes"], order):
        patch.set_facecolor(palette[g])
        patch.set_alpha(0.35)

    for i, g in enumerate(order):
        vals = sample_df.loc[sample_df["group"] == g, "global_methylation_pct"].to_numpy()
        jitter = np.random.default_rng(42 + i).normal(0, 0.04, size=len(vals))
        plt.scatter(np.full(len(vals), i) + jitter, vals, s=45, color=palette[g], edgecolors="black", linewidths=0.3)

    plt.xticks(positions, labels)
    plt.ylabel(f"Coverage-weighted global coding-region {metric_label} (%)")
    plt.title(f"Global coding-region {metric_label} by cohort")
    plt.grid(axis="y", alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def build_readable_summary_table(
    permanova_df: pd.DataFrame,
    permdisp_df: pd.DataFrame,
    group_labels: Dict[str, str],
) -> pd.DataFrame:
    merged = permanova_df.merge(
        permdisp_df[["comparison", "test_scope", "f_stat", "p_value"]],
        on=["comparison", "test_scope"],
        how="left",
        suffixes=("_permanova", "_permdisp"),
    )

    rows: List[Dict[str, object]] = []
    for _, row in merged.iterrows():
        p_perm = float(row["p_value_permanova"])
        p_disp = float(row["p_value_permdisp"]) if pd.notna(row["p_value_permdisp"]) else np.nan
        q_perm = float(row["p_adj_bh"]) if pd.notna(row.get("p_adj_bh", np.nan)) else np.nan

        if p_perm < 0.05 and pd.notna(p_disp) and p_disp < 0.05:
            interpretation = "PERMANOVA sig.; dispersion also differs (interpret centroid shift cautiously)"
        elif p_perm < 0.05 and (pd.isna(p_disp) or p_disp >= 0.05):
            interpretation = "PERMANOVA sig.; no strong dispersion evidence"
        else:
            interpretation = "No PERMANOVA significance"

        rows.append(
            {
                "comparison": format_comparison_label(str(row["comparison"]), group_labels),
                "scope": "Global" if row["test_scope"] == "global" else "Pairwise",
                "n_samples": int(row["n_samples"]),
                "pseudo_f": float(row["pseudo_f"]),
                "r_squared": float(row["r_squared"]),
                "p_permanova": p_perm,
                "q_permanova_bh_pairwise_only": q_perm,
                "f_permdisp": float(row["f_stat"]) if pd.notna(row["f_stat"]) else np.nan,
                "p_permdisp": p_disp,
                "interpretation": interpretation,
            }
        )

    out = pd.DataFrame(rows)
    scope_order = {"Global": 0, "Pairwise": 1}
    out["_scope_order"] = out["scope"].map(scope_order).fillna(9)
    out = out.sort_values(["_scope_order", "comparison"]).drop(columns=["_scope_order"]).reset_index(drop=True)
    return out


def plot_readable_summary_table(table_df: pd.DataFrame, out_png: Path) -> None:
    display_df = table_df.copy()
    for col in [
        "pseudo_f",
        "r_squared",
        "p_permanova",
        "q_permanova_bh_pairwise_only",
        "f_permdisp",
        "p_permdisp",
    ]:
        if col in display_df.columns:
            display_df[col] = display_df[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.4g}")

    fig_h = max(2.5, 0.55 * (len(display_df) + 1))
    fig, ax = plt.subplots(figsize=(18, fig_h))
    ax.axis("off")

    col_widths = [0.16, 0.07, 0.06, 0.07, 0.07, 0.08, 0.11, 0.07, 0.07, 0.24]
    tbl = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        colWidths=col_widths,
        cellLoc="center",
        loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.35)

    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#e9ecef")
            cell.set_text_props(weight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f8f9fa")

    ax.set_title("Genome-wide Coding-region Methylation PERMANOVA/PERMDISP Summary", fontsize=13, fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def write_readable_summary_html(table_df: pd.DataFrame, out_html: Path) -> None:
    html = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'><title>PERMANOVA/PERMDISP Summary</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;margin:24px;background:#f6f8fa;color:#222;}",
        "h2{margin-top:0;}",
        "table{border-collapse:collapse;width:100%;background:white;}",
        "th,td{border:1px solid #d0d7de;padding:8px 10px;font-size:13px;}",
        "th{background:#e9ecef;text-align:center;}",
        "tr:nth-child(even){background:#f8f9fa;}",
        "td:last-child{text-align:left;}",
        "</style></head><body>",
        "<h2>Genome-wide Coding-region Methylation PERMANOVA/PERMDISP Summary</h2>",
        table_df.to_html(index=False, border=0, escape=False),
        "</body></html>",
    ]
    out_html.write_text("\n".join(html), encoding="utf-8")


def compute_group_weighted_region_means(
    meth_matrix: np.ndarray,
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    min_cpgs: int,
) -> Dict[str, np.ndarray]:
    """Compute per-region coverage-weighted mean methylation per group."""
    group_means: Dict[str, np.ndarray] = {}
    groups_arr = np.asarray(groups)
    for group in group_order:
        idx = np.where(groups_arr == group)[0]
        n_genes = meth_matrix.shape[1]
        means = np.full(n_genes, np.nan, dtype=float)
        if len(idx) == 0:
            group_means[group] = means
            continue

        g_meth = meth_matrix[idx, :].astype(float)
        g_total = total_matrix[idx, :].astype(float)
        valid = g_total >= float(min_cpgs)

        weighted_meth = np.where(valid, g_meth, 0.0).sum(axis=0)
        weighted_total = np.where(valid, g_total, 0.0).sum(axis=0)
        np.divide(weighted_meth, weighted_total, out=means, where=weighted_total > 0)
        means *= 100.0
        group_means[group] = means

    return group_means


def plot_spearman_pair(
    mean_a: np.ndarray,
    mean_b: np.ndarray,
    label_a: str,
    label_b: str,
    metric_label: str,
    out_png: Path,
    out_summary: Path,
) -> None:
    valid = np.isfinite(mean_a) & np.isfinite(mean_b)
    x = mean_a[valid]
    y = mean_b[valid]

    if x.size > 1:
        rho, pvalue = scipy_stats.spearmanr(x, y)
        rho = float(rho)
        pvalue = float(pvalue)
    else:
        rho = float("nan")
        pvalue = float("nan")

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, alpha=0.25, s=10, color="#1f77b4", edgecolors="none")

    if x.size > 0 and y.size > 0:
        min_val = float(np.nanmin([np.nanmin(x), np.nanmin(y)]))
        max_val = float(np.nanmax([np.nanmax(x), np.nanmax(y)]))
        plt.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black", linewidth=1)

    plt.xlabel(
        f"{label_a} weighted mean filler-region {metric_label} (%)",
        fontsize=15,
        fontweight="bold",
    )
    plt.ylabel(
        f"{label_b} weighted mean filler-region {metric_label} (%)",
        fontsize=15,
        fontweight="bold",
    )
    plt.title(
        f"Filler-region {metric_label} concordance: {label_a} vs {label_b}\n"
        f"Spearman rho = {rho:.4f}, N = {int(valid.sum())}",
        fontsize=16,
        fontweight="bold",
    )
    ax = plt.gca()
    ax.tick_params(axis="both", labelsize=13)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")
    plt.text(
        0.05,
        0.95,
        f"Spearman rho = {rho:.4f}\nN = {int(valid.sum())}",
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    with open(out_summary, "w", encoding="utf-8") as f:
        f.write(f"Spearman analysis ({metric_label}): {label_a} vs {label_b}\n")
        f.write(f"N regions: {int(valid.sum())}\n")
        f.write(f"Spearman rho: {rho:.8f}\n")
        f.write(f"Spearman p-value: {pvalue:.8e}\n")


def write_spearman_outputs(
    meth_matrix: np.ndarray,
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    min_cpgs: int,
    outdir: Path,
    output_prefix: str,
    metric_label: str,
) -> None:
    group_means = compute_group_weighted_region_means(
        meth_matrix=meth_matrix,
        total_matrix=total_matrix,
        groups=groups,
        group_order=group_order,
        min_cpgs=min_cpgs,
    )

    for g1, g2 in itertools.combinations(group_order, 2):
        label_a = group_labels.get(g1, g1)
        label_b = group_labels.get(g2, g2)
        slug = f"{g1}_vs_{g2}"
        out_png = outdir / f"spearman_{output_prefix}_{slug}.png"
        out_txt = outdir / f"spearman_{output_prefix}_{slug}_summary.txt"
        plot_spearman_pair(group_means[g1], group_means[g2], label_a, label_b, metric_label, out_png, out_txt)


def run_analysis(args: argparse.Namespace) -> None:
    np.random.seed(args.random_seed)
    rng = np.random.default_rng(args.random_seed)

    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found in PATH. Install bedtools first.")

    bed_dir = Path(args.bed_dir)
    if not bed_dir.exists():
        raise FileNotFoundError(f"Input BED directory not found: {bed_dir}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    group_aliases = parse_key_value_arg(args.group_aliases, lower_values=True)
    if not group_aliases:
        group_aliases = DEFAULT_GROUP_ALIASES.copy()

    group_labels = parse_key_value_arg(args.group_labels, lower_values=False)
    if not group_labels:
        group_labels = DEFAULT_GROUP_LABELS.copy()

    samples = discover_samples(bed_dir, group_aliases=group_aliases)
    sample_df = pd.DataFrame([s.__dict__ for s in samples])
    sample_df = sample_df.sort_values(["group", "sample_id"]).reset_index(drop=True)

    samples = [SampleInfo(**row) for row in sample_df.to_dict(orient="records")]

    group_counts = sample_df["group"].value_counts().to_dict()
    discovered_groups = sorted(group_counts.keys())
    requested_order = parse_group_order_arg(args.group_order) if args.group_order else []

    group_order: List[str] = []
    for g in requested_order if requested_order else DEFAULT_GROUP_ORDER:
        if g in discovered_groups and g not in group_order:
            group_order.append(g)
    for g in discovered_groups:
        if g not in group_order:
            group_order.append(g)

    for g in discovered_groups:
        if g not in group_labels:
            group_labels[g] = g.upper() if len(g) <= 5 else g.title()

    print("Discovered samples:")
    for g in group_order:
        print(f"  {group_labels.get(g, g):<12}: {group_counts.get(g, 0)}")

    if len(group_order) < 2:
        raise RuntimeError("Need at least 2 groups for between-group tests.")

    min_group_size = min(group_counts.get(g, 0) for g in group_order)
    if min_group_size < 2:
        raise RuntimeError("Need at least 2 samples per group for stable group-level tests.")

    min_samples = args.min_samples
    if min_samples > min_group_size:
        print(
            f"WARNING: min-samples={min_samples} exceeds smallest group size={min_group_size}. "
            f"Using {min_group_size}."
        )
        min_samples = min_group_size

    region_items = GencodeFillerRegionCatalog(args.gencode_gtf).build_catalog()
    region_count = len(region_items)
    print(f"Loaded filler regions: {region_count:,}")

    region_bed = write_region_bed(outdir, region_items)
    try:
        meth_5mc_matrix, meth_5hmc_matrix, total_matrix = aggregate_all_samples(
            samples=samples,
            region_bed=region_bed,
            gene_count=region_count,
            max_workers=args.workers if args.workers is not None else max(1, os.cpu_count() or 1),
        )
    finally:
        try:
            region_bed.unlink(missing_ok=True)
        except Exception:
            pass

    groups = [s.group for s in samples]
    sample_ids = [s.sample_id for s in samples]

    keep = filter_regions_by_group_coverage(total_matrix, groups, group_order, args.min_cpgs, min_samples)
    kept_regions = np.where(keep)[0]
    if len(kept_regions) < 20:
        raise RuntimeError(
            "Too few filler regions left after filtering. Try lowering --min-cpgs or --min-samples."
        )

    meth_5mc_pct, total_kept = compute_percent_matrix(meth_5mc_matrix, total_matrix, args.min_cpgs, keep)
    meth_5hmc_pct, _ = compute_percent_matrix(meth_5hmc_matrix, total_matrix, args.min_cpgs, keep)
    meth_5mc_kept = meth_5mc_matrix[:, keep]
    meth_5hmc_kept = meth_5hmc_matrix[:, keep]
    kept_region_ids = [str(region_items[i]["region_id"]) for i in kept_regions]
    kept_region_genes = [str(region_items[i]["gene_name"]) for i in kept_regions]
    kept_region_indices = [int(region_items[i]["region_index"]) for i in kept_regions]

    pairwise = list(itertools.combinations(group_order, 2))

    _, _, sample_5mc_df, _ = run_metric_outputs(
        metric_key="5mc",
        metric_title="5mC",
        metric_pct=meth_5mc_pct,
        mod_counts_kept=meth_5mc_kept,
        total_kept=total_kept,
        groups=groups,
        sample_ids=sample_ids,
        group_order=group_order,
        group_labels=group_labels,
        pairwise=pairwise,
        outdir=outdir,
        permutations=args.permutations,
        rng=rng,
        skip_spearman=args.skip_spearman,
    )

    _, _, sample_5hmc_df, _ = run_metric_outputs(
        metric_key="5hmc",
        metric_title="5hmC",
        metric_pct=meth_5hmc_pct,
        mod_counts_kept=meth_5hmc_kept,
        total_kept=total_kept,
        groups=groups,
        sample_ids=sample_ids,
        group_order=group_order,
        group_labels=group_labels,
        pairwise=pairwise,
        outdir=outdir,
        permutations=args.permutations,
        rng=rng,
        skip_spearman=args.skip_spearman,
    )

    sample_out = sample_df.copy()
    sample_out = sample_out.merge(sample_5mc_df, on=["sample_id", "group"], how="left")
    sample_out = sample_out.merge(sample_5hmc_df, on=["sample_id", "group"], how="left")
    sample_out.to_csv(outdir / "sample_global_5mc_5hmc.csv", index=False)

    region_tests_5mc = region_level_kruskal_tests(
        metric_pct=meth_5mc_pct,
        groups=groups,
        group_order=group_order,
        region_ids=kept_region_ids,
        region_genes=kept_region_genes,
        region_indices=kept_region_indices,
    )
    region_tests_5mc.to_csv(outdir / "region_level_5mc_group_kruskal.csv", index=False)

    region_tests_5hmc = region_level_kruskal_tests(
        metric_pct=meth_5hmc_pct,
        groups=groups,
        group_order=group_order,
        region_ids=kept_region_ids,
        region_genes=kept_region_genes,
        region_indices=kept_region_indices,
    )
    region_tests_5hmc.to_csv(outdir / "region_level_5hmc_group_kruskal.csv", index=False)

    region_pairwise_5mc = region_pairwise_mannwhitney_tests(
        metric_pct=meth_5mc_pct,
        groups=groups,
        pairwise=pairwise,
        region_ids=kept_region_ids,
        region_genes=kept_region_genes,
        region_indices=kept_region_indices,
    )
    region_pairwise_5mc.to_csv(outdir / "region_level_5mc_pairwise_mannwhitney.csv", index=False)
    export_significant_gene_table(
        test_df=region_pairwise_5mc,
        out_html=outdir / "region_level_5mc_pairwise_mannwhitney_significant_genes.html",
        out_png=outdir / "region_level_5mc_pairwise_mannwhitney_significant_genes.png",
        title="Significant genes from pairwise Mann-Whitney tests (5mC)",
        q_threshold=0.05,
    )

    region_pairwise_5hmc = region_pairwise_mannwhitney_tests(
        metric_pct=meth_5hmc_pct,
        groups=groups,
        pairwise=pairwise,
        region_ids=kept_region_ids,
        region_genes=kept_region_genes,
        region_indices=kept_region_indices,
    )
    region_pairwise_5hmc.to_csv(outdir / "region_level_5hmc_pairwise_mannwhitney.csv", index=False)
    export_significant_gene_table(
        test_df=region_pairwise_5hmc,
        out_html=outdir / "region_level_5hmc_pairwise_mannwhitney_significant_genes.html",
        out_png=outdir / "region_level_5hmc_pairwise_mannwhitney_significant_genes.png",
        title="Significant genes from pairwise Mann-Whitney tests (5hmC)",
        q_threshold=0.05,
    )

    gene_calls_5mc = gene_calls_from_regions(region_pairwise_5mc, q_threshold=0.05, fraction_threshold=0.5)
    gene_calls_5mc.to_csv(outdir / "gene_level_5mc_from_regions_majority_rule.csv", index=False)

    gene_calls_5hmc = gene_calls_from_regions(region_pairwise_5hmc, q_threshold=0.05, fraction_threshold=0.5)
    gene_calls_5hmc.to_csv(outdir / "gene_level_5hmc_from_regions_majority_rule.csv", index=False)

    sig_5mc = set(
        region_tests_5mc.loc[region_tests_5mc["q_value_bh"] < 0.05, "region_matrix_index"].dropna().astype(int).tolist()
    )
    sig_5hmc = set(
        region_tests_5hmc.loc[region_tests_5hmc["q_value_bh"] < 0.05, "region_matrix_index"].dropna().astype(int).tolist()
    )
    sig_union = np.array(sorted(sig_5mc | sig_5hmc), dtype=int)

    sig_gene_majority_5mc = set(gene_calls_5mc.loc[gene_calls_5mc["is_gene_differential"], "gene"].tolist())
    sig_gene_majority_5hmc = set(gene_calls_5hmc.loc[gene_calls_5hmc["is_gene_differential"], "gene"].tolist())
    sig_gene_majority_union = sig_gene_majority_5mc | sig_gene_majority_5hmc

    plot_significant_gene_stacked_composition(
        meth_5mc=meth_5mc_kept,
        meth_5hmc=meth_5hmc_kept,
        total_matrix=total_kept,
        groups=groups,
        group_order=group_order,
        group_labels=group_labels,
        significant_idx=sig_union,
        out_png=outdir / "stacked_sig_region_cpg_composition_5mc_5hmc.png",
        out_csv=outdir / "stacked_sig_region_cpg_composition_5mc_5hmc.csv",
    )

    qc_df = pd.DataFrame(
        {
            "n_regions_total": [region_count],
            "n_regions_kept": [len(kept_regions)],
            "fraction_regions_kept": [len(kept_regions) / region_count],
            "mean_missing_fraction_5mc_before_impute": [float(np.mean(np.isnan(meth_5mc_pct)))],
            "mean_missing_fraction_5hmc_before_impute": [float(np.mean(np.isnan(meth_5hmc_pct)))],
            "n_sig_regions_union_5mc_5hmc": [int(sig_union.size)],
            "n_sig_genes_majority_rule_union_5mc_5hmc": [int(len(sig_gene_majority_union))],
        }
    )
    qc_df.to_csv(outdir / "filler_region_matrix_qc.csv", index=False)

    lines = [
        "# Genome-wide filler-region methylation PERMANOVA/PERMDISP",
        "",
        "## Inputs",
        f"- BED directory: {bed_dir}",
        f"- GENCODE GTF: {args.gencode_gtf}",
        f"- Group aliases: {group_aliases}",
        f"- Group order used: {group_order}",
        "",
        "## Matrix",
        f"- Samples: {len(samples)}",
        f"- Filler regions total: {region_count:,}",
        f"- Filler regions kept: {len(kept_regions):,}",
        f"- Significant filler regions (union of 5mC/5hmC, q<0.05): {int(sig_union.size):,}",
        f"- Majority-rule differential genes (union of 5mC/5hmC): {int(len(sig_gene_majority_union)):,}",
        "",
        "## Key outputs",
        "- genomewide_5mc_permanova_results.csv",
        "- genomewide_5hmc_permanova_results.csv",
        "- genomewide_5mc_permdisp_results.csv",
        "- genomewide_5hmc_permdisp_results.csv",
        "- genomewide_5mc_permanova_readable_table.csv",
        "- genomewide_5hmc_permanova_readable_table.csv",
        "- global_5mc_univariate_tests.csv",
        "- global_5hmc_univariate_tests.csv",
        "- sample_global_5mc.csv",
        "- sample_global_5hmc.csv",
        "- sample_global_5mc_5hmc.csv",
        "- sample_embedding_5mc.csv",
        "- sample_embedding_5hmc.csv",
        "- sample_embedding_5mc.png",
        "- sample_embedding_5hmc.png",
        "- global_5mc_boxplot.png",
        "- global_5hmc_boxplot.png",
        "- region_level_5mc_group_kruskal.csv",
        "- region_level_5hmc_group_kruskal.csv",
        "- region_level_5mc_pairwise_mannwhitney.csv",
        "- region_level_5hmc_pairwise_mannwhitney.csv",
        "- region_level_5mc_pairwise_mannwhitney_significant_genes.html",
        "- region_level_5mc_pairwise_mannwhitney_significant_genes.png",
        "- region_level_5hmc_pairwise_mannwhitney_significant_genes.html",
        "- region_level_5hmc_pairwise_mannwhitney_significant_genes.png",
        "- gene_level_5mc_from_regions_majority_rule.csv",
        "- gene_level_5hmc_from_regions_majority_rule.csv",
        "- stacked_sig_region_cpg_composition_5mc_5hmc.csv",
        "- stacked_sig_region_cpg_composition_5mc_5hmc.png",
        "- spearman_5mc_<group1>_vs_<group2>.png",
        "- spearman_5hmc_<group1>_vs_<group2>.png",
    ]
    (outdir / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")

    print("Analysis complete")
    print(f"Results directory: {outdir}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Genome-wide filler-region methylation PERMANOVA/PERMDISP for BRCA2/control cfDNA"
    )
    parser.add_argument("--bed-dir", default=DEFAULT_BED_DIR, help="Directory with bedMethyl BED files")
    parser.add_argument("--gencode-gtf", default=DEFAULT_GENCODE_GTF, help="GENCODE GTF path (.gz supported)")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Output directory")
    parser.add_argument("--min-cpgs", type=int, default=DEFAULT_MIN_CPGS, help="Minimum CpGs required per region/sample")
    parser.add_argument(
        "--min-samples",
        type=int,
        default=DEFAULT_MIN_SAMPLES,
        help="Minimum observed samples per group for a region to be retained",
    )
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS, help="Permutation count")
    parser.add_argument("--workers", type=int, default=None, help="Max worker processes for bedtools aggregation")
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED, help="Random seed")
    parser.add_argument(
        "--skip-spearman",
        action="store_true",
        help="Skip Spearman concordance plots/tables to reduce runtime.",
    )
    parser.add_argument(
        "--group-aliases",
        type=str,
        default="",
        help=(
            "Optional prefix:group mapping (comma-separated), e.g. "
            "'ctrl:control,brca2:brca2'. Unmapped prefixes use themselves."
        ),
    )
    parser.add_argument(
        "--group-labels",
        type=str,
        default="",
        help=(
            "Optional group:display_label mapping (comma-separated), e.g. "
            "'control:Control,brca2:BRCA2-999del5'."
        ),
    )
    parser.add_argument(
        "--group-order",
        type=str,
        default=",".join(DEFAULT_GROUP_ORDER),
        help=(
            "Optional comma-separated group order for reporting/plotting, e.g. "
            "'brca2,control'. Discovered groups not listed are appended."
        ),
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    run_analysis(args)


if __name__ == "__main__":
    main()

