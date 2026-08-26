#!/usr/bin/env python3
"""
Convert BED Files to 450K Array Format for MethAtlas Package (Tight Filtering + Count-based Beta)

This script combines the strict filtering logic of *_tight (Nvalid_cov >= 3, exact/approximate match, etc.) with the count-based beta calculation of *_counts (sum Nmod / sum Nvalid_cov per probe).
"""

import os
import sys
import glob
import csv
import time
import bisect
import traceback
from collections import Counter, defaultdict
import pandas as pd
import numpy as np

# Made for modkit version 0.2.0+ with Python 3.8+ and pandas 1.0+

# Configure paths (update as needed)
BED_MANIFEST_FILE = "/set/to/your/manifest"
BED_DIR = "/set/to/your/bed/file/dir"
OUTPUT_DIR = "/set/to/your/desired/output_dir"
LOG_FILE = os.path.join(OUTPUT_DIR, "mapping_log.txt")

# Parameters
MAX_DISTANCE = 5  # Only allow exact matches or within 5bp
BETA_VALUE_FIELD = 10  # Fraction modified (percentage)
NMOD_FIELD = 11        # Nmod (number of methylated reads)
NVALID_COV_FIELD = 9   # Nvalid_cov (total coverage)
USE_BEST_MATCH = True  # Allow best match within MAX_DISTANCE
NVALID_COV_THRESHOLD = 1  # Minimum coverage required

# Create output directory
try:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Created output directory: {OUTPUT_DIR}")
except Exception as e:
    print(f"Error creating output directory: {e}")
    sys.exit(1)

def log(message):
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    full_message = f"[{timestamp}] {message}"
    print(full_message, flush=True)
    try:
        with open(LOG_FILE, 'a') as f:
            f.write(full_message + '\n')
    except Exception as e:
        print(f"ERROR writing to log file: {e}", flush=True)

def normalize_chromosome(chrom):
    norm_chrom = str(chrom).replace('chr', '').upper()
    if '_' in norm_chrom:
        norm_chrom = norm_chrom.split('_')[0]
    if norm_chrom == 'X':
        return 'X'
    elif norm_chrom == 'Y':
        return 'Y'
    elif norm_chrom in ['M', 'MT']:
        return 'MT'
    try:
        return str(int(norm_chrom))
    except ValueError:
        if any(x in chrom.upper() for x in ['_ALT', '_RANDOM', 'HSCHR', 'UN_', 'UNLOCALIZED']):
            return None
        return norm_chrom

def format_illumina_id(probe_id):
    if probe_id.isdigit() or not probe_id.startswith('cg'):
        return f"cg{probe_id}"
    return probe_id

def load_bed_manifest():
    log(f"Loading BED Illumina 450K manifest: {BED_MANIFEST_FILE}")
    manifest_dict = {}
    manifest_positions = {}
    manifest_count = 0
    try:
        for chr_num in range(1, 23):
            manifest_dict[str(chr_num)] = {}
            manifest_positions[str(chr_num)] = []
        manifest_dict['X'] = {}
        manifest_dict['Y'] = {}
        manifest_dict['MT'] = {}
        manifest_positions['X'] = []
        manifest_positions['Y'] = []
        manifest_positions['MT'] = []
        with open(BED_MANIFEST_FILE, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track') or line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    continue
                chr_val = fields[0]
                start_pos = int(fields[1])
                probe_id = fields[3]
                chr_norm = normalize_chromosome(chr_val)
                if not chr_norm or chr_norm not in manifest_dict:
                    continue
                formatted_probe_id = format_illumina_id(probe_id)
                manifest_dict[chr_norm][start_pos] = formatted_probe_id
                manifest_positions[chr_norm].append(start_pos)
                manifest_count += 1
                if manifest_count % 100000 == 0:
                    log(f"Processed {manifest_count:,} BED manifest lines...")
        # Sort manifest positions once
        for chr_key in manifest_positions:
            manifest_positions[chr_key].sort()
        log(f"Successfully loaded {manifest_count:,} positions from BED manifest")
        for chr_key, pos_dict in manifest_dict.items():
            if pos_dict:
                log(f"  BED Chr {chr_key}: {len(pos_dict):,} positions")
        return manifest_dict, manifest_positions
    except Exception as e:
        log(f"ERROR loading BED manifest: {e}")
        log(traceback.format_exc())
        return None

def process_bed_file(bed_file, manifest_tuple, max_distance=MAX_DISTANCE):
    manifest_dict, manifest_positions = manifest_tuple
    manifest_sets = {k: set(v) for k, v in manifest_positions.items()}
    sample_name = os.path.basename(bed_file).replace('.bed', '')
    log(f"Processing {sample_name}...")
    results = {
        'sample_name': sample_name,
        'mapped_positions': {},  # {probe_id: beta_value}
        'stats': {
            'total_positions': 0,
            'mapped_positions': 0,
            'exact_matches': 0,
            'approximate_matches': 0,
            'chromosomes': Counter(),
            'mapped_by_chr': Counter(),
            'distances': Counter(),
        }
    }
    # For each probe, sum Nmod and Nvalid_cov
    probe_counts = defaultdict(lambda: {'nmod': 0, 'nvalid_cov': 0})

    # --- Coordinate system detection ---
    offset_detected = 0
    try:
        with open(bed_file, 'r') as f:
            test_lines = []
            for _ in range(1000):
                line = f.readline()
                if not line:
                    break
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                try:
                    pos = int(fields[1])
                except (IndexError, ValueError):
                    continue
                # Heuristic: if any position is 0, likely 0-based; if all >=1, likely 1-based
                test_lines.append(pos)
            if test_lines:
                if any(p == 0 for p in test_lines):
                    offset_detected = 0
                    log("Detected 0-based coordinates in BED file (positions start at 0). No adjustment will be made.")
                elif all(p >= 1 for p in test_lines):
                    offset_detected = -1
                    log("Detected 1-based coordinates in BED file (positions start at 1). Will subtract 1 from positions for manifest matching.")
                else:
                    log("Warning: Unable to confidently determine coordinate system. Assuming 0-based.")
            else:
                log("Warning: No valid positions found in first 1000 lines for coordinate check. Assuming 0-based.")
    except Exception as e:
        log(f"ERROR during coordinate system detection: {e}")
        log(traceback.format_exc())

    try:
        with open(bed_file, 'r') as f:
            for line_idx, line in enumerate(f, 1):
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                try:
                    nvalid_cov = float(fields[NVALID_COV_FIELD])
                except (ValueError, IndexError):
                    nvalid_cov = None
                if nvalid_cov is None or nvalid_cov < NVALID_COV_THRESHOLD:
                    continue
                try:
                    nmod = float(fields[NMOD_FIELD])
                except (ValueError, IndexError):
                    nmod = None
                if nmod is None:
                    continue
                chr_val = fields[0]
                start_pos = int(fields[1]) + offset_detected
                chr_norm = normalize_chromosome(chr_val)
                results['stats']['chromosomes'][chr_norm] += 1
                results['stats']['total_positions'] += 1
                # Use local variables for manifest lookups
                chr_manifest_dict = manifest_dict.get(chr_norm)
                chr_manifest_set = manifest_sets.get(chr_norm)
                chr_positions = manifest_positions.get(chr_norm)
                if not chr_manifest_dict or not chr_manifest_set or not chr_positions:
                    continue
                best_match_probe = None
                best_match_dist = float('inf')
                # Fast exact match
                if start_pos in chr_manifest_set:
                    best_match_probe = chr_manifest_dict[start_pos]
                    best_match_dist = 0
                    results['stats']['exact_matches'] += 1
                elif USE_BEST_MATCH:
                    idx = bisect.bisect_left(chr_positions, start_pos)
                    if idx > 0:
                        left_pos = chr_positions[idx-1]
                        left_dist = abs(start_pos - left_pos)
                        if left_dist <= max_distance and left_dist < best_match_dist:
                            best_match_probe = chr_manifest_dict[left_pos]
                            best_match_dist = left_dist
                            results['stats']['approximate_matches'] += 1
                    if idx < len(chr_positions):
                        if chr_positions[idx] == start_pos:
                            best_match_probe = chr_manifest_dict[start_pos]
                            best_match_dist = 0
                            results['stats']['exact_matches'] += 1
                        else:
                            right_pos = chr_positions[idx]
                            right_dist = abs(start_pos - right_pos)
                            if right_dist <= max_distance and right_dist < best_match_dist:
                                best_match_probe = chr_manifest_dict[right_pos]
                                best_match_dist = right_dist
                                results['stats']['approximate_matches'] += 1
                if best_match_probe is not None:
                    probe_counts[best_match_probe]['nmod'] += nmod
                    probe_counts[best_match_probe]['nvalid_cov'] += nvalid_cov
                    results['stats']['mapped_positions'] += 1
                    results['stats']['distances'][best_match_dist] += 1
                    results['stats']['mapped_by_chr'][chr_norm] += 1
                if line_idx % 1000000 == 0:
                    log(f"  Processed {results['stats']['total_positions']:,} positions...")
        # Calculate beta for each probe
        for probe_id, counts in probe_counts.items():
            nmod = counts['nmod']
            nvalid_cov = counts['nvalid_cov']
            if nvalid_cov > 0:
                beta = nmod / nvalid_cov
                results['mapped_positions'][probe_id] = beta
        mapping_rate = (len(results['mapped_positions']) / results['stats']['total_positions'] * 100) if results['stats']['total_positions'] > 0 else 0
        log(f"  {sample_name}: {len(results['mapped_positions']):,} mapped of {results['stats']['total_positions']:,} positions ({mapping_rate:.2f}%)")
        log(f"  Exact matches: {results['stats']['exact_matches']:,}, Approximate matches: {results['stats']['approximate_matches']:,}")
        return results
    except Exception as e:
        log(f"ERROR processing {sample_name}: {e}")
        log(traceback.format_exc())
        return None

def export_to_csv(results, output_file, methatlas_format=True):
    log(f"Exporting data to {output_file}")
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            if methatlas_format:
                writer.writerow(['IlmnID', 'Beta_Value'])
            else:
                writer.writerow(['ProbeID', 'Beta_Value'])
            for probe_id, beta_value in results['mapped_positions'].items():
                writer.writerow([probe_id, f"{beta_value:.6f}"])
        log(f"Successfully exported {len(results['mapped_positions']):,} positions to {output_file}")
        return True
    except Exception as e:
        log(f"ERROR exporting data: {e}")
        log(traceback.format_exc())
        return False

def export_combined_methatlas_format(all_results, output_file):
    log(f"Creating combined MethAtlas format output: {output_file}")
    try:
        all_probe_ids = set()
        for result in all_results.values():
            all_probe_ids.update(result['mapped_positions'].keys())
        log(f"Total unique probe IDs across all samples: {len(all_probe_ids):,}")
        matrix_df = pd.DataFrame(index=sorted(all_probe_ids))
        matrix_df.index.name = 'IlmnID'
        for sample_name, result in all_results.items():
            beta_values = pd.Series(result['mapped_positions'], name=sample_name)
            matrix_df[sample_name] = beta_values
        matrix_df = matrix_df.fillna('NA')
        matrix_df.to_csv(output_file)
        log(f"Successfully exported combined matrix with {len(matrix_df.index)} rows and {len(matrix_df.columns)} columns")
        return True
    except Exception as e:
        log(f"ERROR creating combined output: {e}")
        log(traceback.format_exc())
        return False

def main():
    start_time = time.time()
    with open(LOG_FILE, 'w') as f:
        f.write("BED to 450K Array Transformation for MethAtlas (Tight Filtering + Count-based)\n")
        f.write("===========================================\n\n")
    log("Starting BED to 450K array transformation using BED manifest approach (tight filtering + count-based beta)")
    log(f"BED Manifest: {BED_MANIFEST_FILE}")
    log(f"BED Directory: {BED_DIR}")
    log(f"Maximum distance allowed: {MAX_DISTANCE}bp")
    manifest = load_bed_manifest()
    if not manifest or (isinstance(manifest, tuple) and not manifest[0]):
        log("Failed to load BED manifest. Exiting.")
        return
    manifest_dict, manifest_positions = manifest
    bed_files = glob.glob(os.path.join(BED_DIR, "*.bed"))
    if not bed_files:
        log(f"No BED files found in {BED_DIR}")
        return
    log(f"Found {len(bed_files)} BED files")
    all_results = {}
    for file_idx, bed_file in enumerate(bed_files, 1):
        log(f"[{file_idx}/{len(bed_files)}] Processing file: {os.path.basename(bed_file)}...")
        result = process_bed_file(bed_file, (manifest_dict, manifest_positions))
        if result:
            all_results[result['sample_name']] = result
            output_file = os.path.join(OUTPUT_DIR, f"{result['sample_name']}_450k.csv")
            export_to_csv(result, output_file, methatlas_format=True)
    if all_results:
        combined_output = os.path.join(OUTPUT_DIR, "methatlas_combined_450k.csv")
        export_combined_methatlas_format(all_results, combined_output)
    elapsed_time = time.time() - start_time
    log(f"Transformation completed in {elapsed_time:.2f} seconds")
    log(f"Results saved to: {OUTPUT_DIR}")
    log(f"The files are now ready to be used with the MethAtlas package")

if __name__ == "__main__":
    main()
