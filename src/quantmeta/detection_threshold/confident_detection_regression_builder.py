#!/usr/bin/env python3
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "argparse>=1.4.0",
#     "importlib>=1.0.4",
#     "matplotlib>=3.10.8",
#     "numpy>=2.4.4",
#     "pandas>=3.0.2",
#     "path>=17.1.1",
#     "scikit-learn>=1.8.0",
# ]
# ///
"""Purpose: This script develops a regression to determine E_detect (minimum E_rel) 
for confident detection with respect to a target's length. The regression from 
Langenfeld et al. (2025) may be adopted, but the minimum E_rel is based on the 
standards proposed by FastViromeExplorer (Lithi et al. 2018) (10% read coverage and 
0.3 observed/expected read distribution). If these parameters are deemed inappropriate 
for specific applications, it is recommended that new E_detect regressions be created 
for specific research needs.

Adjust the minimum coverage (min_coverage) and minimum read distribution 
(min_distribution) parameters to fit your needs.

Usage:
python confident_detection_regression_builder.py --sample-file Config/sample_list.txt --std_file Spike-ins/sequins_Mix_A.txt --ssDNA-std-file Spike-ins/ssDNA_standards.txt --min-coverage 0.1 --min-distribution 0.3 --test-database RefSeq_Viruses

By default, the script uses the example layout:
- sample:          example_1, example_2
- downsample:      1, 10, 100
- fail_set:        fail_standards_r{1..5}

Outputs:
- Regressions/detection/detection_threshold_custom.json
- Regressions/detection/optimal_thresholds_plot.png
"""

import argparse
import json
from pathlib import Path
import importlib.resources as resources
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from sklearn.linear_model import LogisticRegression
except ImportError:
    raise ImportError("scikit-learn is required: pip install scikit-learn")

try:
    from sklearn.metrics import roc_curve
except ImportError:
    raise ImportError("scikit-learn is required for ROC curve calculation: pip install scikit-learn")

def get_sample_file_path(*parts) -> Path:
    """Get path to package data files."""
    with resources.as_file(resources.files('quantmeta.data') / 'Config' / 'sample_list.txt') as p:
        base = Path(p).parent.parent.parent
    return base / 'data' / Path(*parts)

def get_std_file_path(*parts) -> Path:
    """Get path to package data files."""
    with resources.as_file(resources.files('quantmeta.data') / 'Spike-ins' / 'sequins_Mix_A.txt') as p:
        base = Path(p).parent.parent.parent
    return base / 'data' / Path(*parts)

def get_ssdna_file_path(*parts) -> Path:
    """Get path to package data files."""
    with resources.as_file(resources.files('quantmeta.data') / 'Spike-ins' / 'ssDNA_stds.txt') as p:
        base = Path(p).parent.parent.parent
    return base / 'data' / Path(*parts)

def mapping_analysis(mapping_results, target_lengths, min_coverage, min_distribution):
    # mapping_results: DataFrame with ID, read_depth
    # target_lengths: DataFrame with ID, length

    # Start with unique IDs and merge lengths
    targets = pd.DataFrame({'ID': mapping_results['ID'].unique()})
    targets = pd.merge(targets, target_lengths, on='ID', how='left')
    
    # Group mapping_results by ID for efficient aggregation
    grouped = mapping_results.groupby('ID')
    
    # Compute basic aggregations
    cover_count = grouped['read_depth'].agg(lambda x: (x > 0).sum())
    B_G = grouped['read_depth'].sum()
    gene_copies = grouped['read_depth'].mean()
    
    # Compute I_G using vectorized operations within groups
    def compute_I_G(group):
        B_G_val = group['read_depth'].sum()
        if B_G_val > 0:
            ratios = group['read_depth'] / B_G_val
            log_ratios = np.log(ratios.replace(0, np.nan))  # Avoid log(0)
            I_G_x = ratios * log_ratios
            return -np.nansum(I_G_x.replace([-np.inf, np.inf], np.nan))
        return 0.0
    
    I_G = grouped.apply(compute_I_G)

    # Merge aggregated results into targets
    targets = pd.merge(targets,
        pd.DataFrame({
            'ID': B_G.index,
            'cover_count': cover_count.values,
            'B_G': B_G.values,
            'gene_copies': gene_copies.values,
            'I_G': I_G.values
        }),
        on='ID',
        how='left'
    ).fillna({'cover_count': 0, 'B_G': 0.0, 'gene_copies': 0.0, 'I_G': 0.0})
    
    # Compute remaining metrics using vectorized operations
    targets = targets.assign(
        E_rel=lambda df: df['I_G'] / np.log(df['length']).where(df['length'] > 1, 1),
        C_o=lambda df: df['cover_count'] / df['length'].where(df['length'] > 0, 1),
        C_e=lambda df: (1 - np.exp(-df['B_G'] / df['length'])).where(df['length'] > 0, 0),
        R_FVE=lambda df: (df['C_o'] / df['C_e']).where(df['C_e'] > 0, 0),
        LOD=lambda df: ((df['C_o'] > min_coverage) & (df['R_FVE'] > min_distribution)).astype(int)
    )
    
    return targets

def find_optimal_cutpoint(y_true, y_scores, n_boot=1000):
    """
    Find the optimal cutpoint by maximizing the bootstrapped Youden's J statistic.
    
    This replicates the behavior of cutpointr::maximize_boot_metric in R.
    
    Args:
        y_true (array-like): True binary labels (0 or 1).
        y_scores (array-like): Predicted scores (probabilities or continuous values).
        n_boot (int): Number of bootstrap resamples (default: 1000).
    
    Returns:
        float: The optimal cutpoint (threshold) value.
    """
    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)
    
    if len(np.unique(y_true)) != 2:
        raise ValueError("y_true must be binary (0 and 1).")
    
    # Use unique scores as candidate thresholds (matching cutpointr's approach)
    thresholds = np.unique(y_scores)
    n = len(y_true)
    
    # Store bootstrapped metrics for each threshold
    boot_metrics = np.zeros((n_boot, len(thresholds)))
    
    for b in range(n_boot):
        # Bootstrap resample with replacement
        indices = np.random.choice(n, n, replace=True)
        y_true_boot = y_true[indices]
        y_scores_boot = y_scores[indices]
        
        for i, thresh in enumerate(thresholds):
            # Predict based on threshold
            y_pred = (y_scores_boot >= thresh).astype(int)
            
            # Compute confusion matrix
            tp = np.sum((y_pred == 1) & (y_true_boot == 1))
            tn = np.sum((y_pred == 0) & (y_true_boot == 0))
            fp = np.sum((y_pred == 1) & (y_true_boot == 0))
            fn = np.sum((y_pred == 0) & (y_true_boot == 1))
            
            # Compute Youden's J: sensitivity + specificity - 1
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            youden_j = sensitivity + specificity - 1
            
            boot_metrics[b, i] = youden_j
    
    # Average Youden's J across bootstraps for each threshold
    avg_metrics = np.mean(boot_metrics, axis=0)
    
    # Find the threshold with the maximum average metric
    optimal_idx = np.argmax(avg_metrics)
    optimal_cutpoint = thresholds[optimal_idx]
    
    return optimal_cutpoint

def main():
    parser = argparse.ArgumentParser(description='Build confident detection threshold regression')
    parser.add_argument('--sample-file', default='Config/sample_list.txt', help='TSV with columns sample')
    parser.add_argument('--std-file', default='Spike-ins/sequins_Mix_A.txt', help='Table of dsDNA standards (ID, Mass, Rel_Abund, length)')
    parser.add_argument('--ssDNA-std-file', default=None, help='Optional table of ssDNA standards (ID, Mass, Rel_Abund, length)')
    parser.add_argument('--min-coverage', type=float, default=0.1, help='Minimum read coverage threshold for detection')
    parser.add_argument('--min-distribution', type=float, default=0.3, help='Minimum observed/Poisson distributed read distribution threshold for detection')
    parser.add_argument('--test-database', default=None, help='Read mapping to a test dataset with a range of relevant target lengths')
    parser.add_argument('--output-dir', default='QuantMeta/', help='Output directory for the project')
    
    args = parser.parse_args()

    out_dir = args.output_dir

    sample_file = args.sample_file
    if sample_file == 'Config/sample_list.txt':  # default
        sample_file = get_sample_file_path('Config', 'sample_list.txt')
    sample_names = pd.read_csv(sample_file, sep='\t', header = None)
    sample_names.columns = ['sample']

    samples = []
    
    for sample in sample_names['sample']:
        for ds in [1, 10, 100]:
            samples.append({'sample': sample, 'downsample': ds})
    samples = pd.DataFrame(samples)

    fail_samples = []

    for sample in sample_names['sample']:
        for fail_set in ['fail_standards_r1', 'fail_standards_r2', 'fail_standards_r3', 'fail_standards_r4', 'fail_standards_r5']:
            fail_samples.append({'sample': sample, 'fail_set': fail_set})
    fail_samples = pd.DataFrame(fail_samples)

    std_file = args.std_file
    if std_file == 'Spike-ins/sequins_Mix_A.txt':  # default
        std_file = get_std_file_path('Spike-ins', 'sequins_Mix_A.txt')
    stds = pd.read_csv(std_file, sep='\t')
    stds_list = stds.loc[~stds['Rel_Abund'].isna(), ['ID', 'length']].copy()
    if args.ssDNA_std_file:
        ssDNA_file = args.ssDNA_std_file
        if ssDNA_file == 'Spike-ins/ssDNA_stds.txt':  #from Langenfeld et al. 2025
            ssDNA_file = get_ssdna_file_path('Spike-ins', 'ssDNA_stds.txt')
        ssdna = pd.read_csv(ssDNA_file, sep='\t')
        ssdna = ssdna.loc[~ssdna['Rel_Abund'].isna(), ['ID', 'length']].copy()
        stds_list = pd.concat([stds_list, ssdna], ignore_index=True)

    min_coverage = args.min_coverage
    min_distribution = args.min_distribution

    # standard sample mapping
    logit = []
    mapping_results = []

    for _, row in samples.iterrows():
        sample_name = str(row['sample'])
        downsample = str(row['downsample'])
        mapping_path = Path(f'{out_dir}/Mapping/{sample_name}_{downsample}/standards_mapping.txt')
        if not mapping_path.exists():
            print(f"WARNING: mapping file not found: {mapping_path}")
            continue

        mapping_out = pd.read_csv(mapping_path, sep='\t')

        if mapping_out.empty == False:
            out_logit = mapping_analysis(mapping_out, stds_list, min_coverage, min_distribution)
            out_map = out_logit.copy()
            out_logit['downsample'] = downsample
            out_logit['sample'] = sample_name
            out_map['downsample'] = downsample
            out_map['sample'] = sample_name
            logit.append(out_logit)
            mapping_results.append(out_map)

    logit = pd.concat(logit, ignore_index=True) if logit else pd.DataFrame()
    mapping_results = pd.concat(mapping_results, ignore_index=True) if mapping_results else pd.DataFrame()

    # fail standards mapping
    fail_logit = []
    fail_mapping_results = []

    for _, row in fail_samples.iterrows():
        sample_name = str(row['sample'])
        fail_set = str(row['fail_set'])
        mapping_path = Path(f'{out_dir}/Mapping/{sample_name}/{fail_set}_mapping.txt')
        if not mapping_path.exists():
            print(f"WARNING: fail mapping file not found: {mapping_path}")
            continue

        mapping_out = pd.read_csv(mapping_path, sep='\t')

        if mapping_out.empty == False:
            fail_lengths = pd.read_csv(f"{out_dir}/Map_Indexes/{fail_set}/fail_standards_lengths.txt", sep='\t')
            fail_lengths.columns = ['ID', 'length']
            out_logit = mapping_analysis(mapping_out, fail_lengths, min_coverage, min_distribution)
            out_map = out_logit.copy()
            out_logit['ID'] = list(map(lambda x: f"{x}_{fail_set}", out_logit['ID']))
            out_logit['sample'] = sample_name
            out_map['ID'] = list(map(lambda x: f"{x}_{fail_set}", out_map['ID']))
            out_map['sample'] = sample_name
            fail_logit.append(out_logit)
            fail_mapping_results.append(out_map)

    fail_logit = pd.concat(fail_logit, ignore_index=True) if fail_logit else pd.DataFrame()
    fail_mapping_results = pd.concat(fail_mapping_results, ignore_index=True) if fail_mapping_results else pd.DataFrame()

    # filter for fails
    fail_logit = fail_logit.loc[fail_logit['LOD'] == 0].copy()
    fail_mapping_results = fail_mapping_results.loc[fail_mapping_results['sample'].isin(fail_logit['sample']) &
                                                    fail_mapping_results['ID'].isin(fail_logit['ID'])].copy()

    if not logit.empty:
        logit['fail_set'] = 'Not_applicable'
    if not mapping_results.empty:
        mapping_results['fail_set'] = 'Not_applicable'

    fail_logit['downsample'] = 100
    fail_mapping_results['downsample'] = 100

    combined_logit = pd.concat([logit, fail_logit], ignore_index=True) if not logit.empty or not fail_logit.empty else pd.DataFrame()
    combined_mapping = pd.concat([mapping_results, fail_mapping_results], ignore_index=True) if not mapping_results.empty or not fail_mapping_results.empty else pd.DataFrame()

    if combined_logit.empty:
        raise RuntimeError('No combined logit data available for regression')

    # logistic regression LOD ~ E_rel
    X = combined_logit[['E_rel']].fillna(0.0)
    y = combined_logit['LOD']

    model = LogisticRegression(solver='lbfgs', max_iter=10000)
    model.fit(X, y)

    # test dataset sample mapping
    test_logit = []
    test_mapping_results = []

    for _, row in sample_names.iterrows():
        sample_name = str(row['sample'])
        test_db = str(args.test_database)
        mapping_path = Path(f'{out_dir}/Mapping/{sample_name}/{test_db}_mapping.txt')
        if not mapping_path.exists():
            print(f"WARNING: mapping file not found: {mapping_path}")
            continue

        mapping_out = pd.read_csv(mapping_path, sep='\t')
        test_lengths = pd.read_csv(f"{out_dir}/Map_Indexes/{test_db}_lengths.txt", sep='\t')
        test_lengths.columns = ['ID', 'length']
        out_logit = mapping_analysis(mapping_out, test_lengths, min_coverage, min_distribution)
        out_map = out_logit.copy()
        out_logit['downsample'] = 100
        out_logit['sample'] = sample_name
        out_map['downsample'] = 100
        out_map['sample'] = sample_name
        test_logit.append(out_logit)
        test_mapping_results.append(out_map)

    test_logit = pd.concat(test_logit, ignore_index=True) if test_logit else pd.DataFrame()
    test_mapping_results = pd.concat(test_mapping_results, ignore_index=True) if test_mapping_results else pd.DataFrame()

    # create 50 length bins for test dataset and predict LOD probabilities
    if not test_logit.empty:
        test_logit['length_bin'] = pd.qcut(test_logit['length'], q=50, duplicates='drop')
        test_logit['LOD_prob'] = model.predict_proba(test_logit[['E_rel']].fillna(0.0))[:, 1]
        
        # Calculate optimal cutpoint (threshold on LOD_prob) for each length bin
        optimal_thresholds = {}
        for bin_name, group in test_logit.groupby('length_bin'):
            # Ensure the group has at least 2 samples and variability in LOD (to avoid errors)
            if len(group) >= 2 and group['LOD'].nunique() > 1:
                threshold = find_optimal_cutpoint(group['LOD'], group['LOD_prob'])
                optimal_thresholds[bin_name] = threshold
            else:
                optimal_thresholds[bin_name] = None  # No threshold if insufficient data
        
        # Plot the thresholds
        lengths = []
        thresholds = []
        for bin_name, thresh in optimal_thresholds.items():
            if thresh is not None:
                # Calculate midpoint of the bin interval
                mid_length = np.log10((bin_name.left + bin_name.right) / 2)
                lengths.append(mid_length)
                thresholds.append(thresh)
        
        if lengths and thresholds:
            from scipy.stats import linregress
            
            # Fit linear regression: y = slope * x + intercept
            slope, intercept, r_value, p_value, std_err = linregress(lengths, thresholds)
            
            # Generate points for the regression line
            lengths_sorted = np.sort(lengths)
            regression_line = slope * lengths_sorted + intercept
            
            plt.figure(figsize=(10, 6))
            plt.plot(lengths, thresholds, 'o-', color='blue', markersize=4, label='Optimal Thresholds')
            plt.plot(lengths_sorted, regression_line, 'r--', linewidth=2, label=f'Linear Fit: y = {slope:.4f}x + {intercept:.4f}\nR² = {r_value**2:.4f}')
            plt.xlabel('Target Length (midpoint of bin, log scale)')
            plt.ylabel('Optimal Threshold (on LOD Probability)')
            plt.title('Optimal Threshold vs. Target Length')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            output_path = Path(f'{out_dir}/Regressions/detection/optimal_thresholds_plot.png')
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_path, dpi=300)
            print("Plot saved to Regressions/detection/optimal_thresholds_plot.png")
            print(f"Linear regression: slope={slope:.4f}, intercept={intercept:.4f}, R²={r_value**2:.4f}")
        else:
            print("No valid thresholds to plot.")
        
        # Add the threshold back to test_logit if needed (per row)
        test_logit['optimal_threshold'] = test_logit['length_bin'].map(optimal_thresholds)

    model_out = {
        'intercept': intercept,
        'slope': slope,
        'equation': f"E_detect = {intercept:.4f} + {slope:.4f} * log10(length)"
    }

    res_file = Path(f'{out_dir}/Regressions/detection/detection_threshold_custom.json')
    res_file.write_text(json.dumps(model_out, indent=2))

    print('✅ Detection threshold regression custom build complete')
    print(f"Model equation: {model_out['equation']}")


if __name__ == '__main__':
    main()
