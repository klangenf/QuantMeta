#!/usr/bin/env python
"""
Quantification Correction Regression Builder

Purpose: Detecting and correcting quantification errors caused by non-specific 
mapping or assembly errors requires limitations on the acceptable read depth 
variability across target sequences. The acceptable read depth variability may 
differ depending on the library preparation, if PCR amplification was performed, 
and the sequencing technology used as each may introduce different bias and 
increase or decrease how much read depth may be an intrinsic result of sequencing. 
Langenfeld et al. (2025) used Swift 1S Plus library prep for simultaneous sequencing 
of dsDNA and ssDNA with Illumina NovaSeq on SP flowcells to produce 251-bp paired-end 
reads. It is recommended that specific read depth variability thresholds are developed 
for each sequencing protocol.
"""

import pandas as pd
import numpy as np
import json
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import linregress
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')


def Shannon_entropy(x):
    """Calculate Shannon entropy for a sequence."""
    if len(x) == 0:
        return 0
    counts = pd.Series(list(x)).value_counts()
    proportions = counts / len(x)
    entropy = -np.sum(proportions * np.log2(proportions))
    return entropy


def mapping_analysis(df):
    """
    Calculate mapping metrics including I_G (GC information content).
    Similar to R's mapping_analysis function.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with 'sequence' column
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with added mapping metrics
    """
    # Vectorized GC content calculation
    sequences = df['sequence'].fillna('')
    gc_content = (sequences.str.count('G') + sequences.str.count('C')) / sequences.str.len()
    gc_content = gc_content.replace([np.inf, -np.inf], 0).fillna(0)
    
    # Vectorized entropy calculation
    def vectorized_entropy(seq_series):
        """Calculate Shannon entropy for each sequence in a Series."""
        results = []
        for seq in seq_series:
            if len(seq) == 0:
                results.append(0)
                continue
            counts = pd.Series(list(seq)).value_counts()
            proportions = counts / len(seq)
            entropy = -np.sum(proportions * np.log2(proportions))
            results.append(entropy)
        return pd.Series(results, index=seq_series.index)
    
    entropy = vectorized_entropy(sequences)
    
    # Vectorized I_G calculation (GC information content)
    I_G = -2 * (gc_content * np.log2(gc_content + 1e-10) + 
                (1 - gc_content) * np.log2(1 - gc_content + 1e-10))
    
    # Create result DataFrame
    result_df = pd.DataFrame({
        'gc_content': gc_content,
        'entropy': entropy,
        'I_G': I_G
    })
    
    return pd.concat([df.reset_index(drop=True), result_df], axis=1)


def load_detection_model(model_path):
    """
    Load detection model (log10-linear regression for detection threshold).
    
    Parameters:
    -----------
    model_path : str
        Path to detection model (JSON format with 'intercept', 'slope', 'type')
    
    Returns:
    --------
    dict : Model parameters
    """
    try:
        with open(model_path, 'r') as f:
            model = json.load(f)
        return model
    except FileNotFoundError:
        print(f"Warning: Detection model not found at {model_path}")
        return None


def predict_E_detect(E_detect_model, lengths):
    """
    Predict detection threshold based on sequence length.
    Uses log10 linear formula: E_detect = intercept + slope * log10(length)
    
    Parameters:
    -----------
    E_detect_model : dict
        Model with 'intercept' and 'slope'
    lengths : array-like
        Sequence lengths
    
    Returns:
    --------
    np.array : Predicted detection thresholds
    """
    if E_detect_model is None:
        return np.zeros_like(lengths)
    
    return E_detect_model['intercept'] + E_detect_model['slope'] * np.log10(np.asarray(lengths))


def fit_read_var_regression(X, y, bin_name):
    """
    Fit polynomial regression with appropriate terms for each depth bin.
    
    Parameters:
    -----------
    X : pd.DataFrame
        Feature matrix with columns: avg_GC, total_avg_depth, E_rel (if applicable)
    y : np.array
        Response variable (avg_depth or log-transformed)
    bin_name : str
        Bin identifier ('bin1', 'bin2', 'bin3', or 'bin4')
    
    Returns:
    --------
    dict : Model with coefficients, type, and bin information
    """
    
    if bin_name == 'bin4':  # >= 1000 reads/bp
        # avg_depth ~ poly(avg_GC, 2) + log(total_avg_depth)
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            np.log(X['total_avg_depth'])
        ])
        response_type = 'linear'
        
    elif bin_name == 'bin3':  # 100-1000 reads/bp
        # log(avg_depth) ~ poly(avg_GC, 2) + log(total_avg_depth)
        y = np.log(y)
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            np.log(X['total_avg_depth'])
        ])
        response_type = 'log'
        
    elif bin_name == 'bin2':  # 10-100 reads/bp
        # log(avg_depth+1) ~ poly(avg_GC, 2) + poly(log(total_avg_depth), 2)
        y = np.log(y + 1)
        log_depth = np.log(X['total_avg_depth'])
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            log_depth**2, log_depth  # poly(log(total_avg_depth), 2)
        ])
        response_type = 'log'
        
    elif bin_name == 'bin1':  # 0-10 reads/bp
        # log(avg_depth+1) ~ poly(avg_GC, 2) + poly(log(total_avg_depth), 2) + poly(E_rel, 2)
        y = np.log(y + 1)
        log_depth = np.log(X['total_avg_depth'])
        X_fit = np.column_stack([
            X['avg_GC']**2, X['avg_GC'], np.ones(len(X)),  # poly(avg_GC, 2)
            log_depth**2, log_depth,  # poly(log(total_avg_depth), 2)
            X['E_rel']**2, X['E_rel']  # poly(E_rel, 2)
        ])
        response_type = 'log'
    
    # Fit OLS regression
    coeffs = np.linalg.lstsq(X_fit, y, rcond=None)[0]
    residuals = y - X_fit @ coeffs
    sse = np.sum(residuals**2)
    mse = sse / (len(y) - len(coeffs))
    rmse = np.sqrt(mse)
    
    return {
        'coefficients': coeffs,
        'response_type': response_type,
        'bin': bin_name,
        'X_test': X,
        'y_original': y,
        'residuals': residuals,
        'rmse': rmse
    }


def predict_read_var(model, X_new):
    """
    Predict using fitted read depth variability model.
    
    Parameters:
    -----------
    model : dict
        Fitted model from fit_read_var_regression
    X_new : pd.DataFrame
        New features with same structure
    
    Returns:
    --------
    np.array : Predictions
    """
    bin_name = model['bin']
    coeffs = model['coefficients']
    
    if bin_name == 'bin4':
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            np.log(X_new['total_avg_depth'])
        ])
        preds = X_fit @ coeffs
        
    elif bin_name == 'bin3':
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            np.log(X_new['total_avg_depth'])
        ])
        preds = X_fit @ coeffs
        preds = np.exp(preds)
        
    elif bin_name == 'bin2':
        log_depth = np.log(X_new['total_avg_depth'])
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            log_depth**2, log_depth
        ])
        preds = X_fit @ coeffs
        preds = np.exp(preds) - 1
        
    elif bin_name == 'bin1':
        log_depth = np.log(X_new['total_avg_depth'])
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            log_depth**2, log_depth,
            X_new['E_rel']**2, X_new['E_rel']
        ])
        preds = X_fit @ coeffs
        preds = np.exp(preds) - 1
    
    return preds

def calculate_rmse(standards_filtered, reg1, reg2, reg3, reg4):
    """
    Vectorized RMSE calculation using groupby operations.
    """
    # Add bin assignment
    def assign_bin(depth):
        if depth >= 1000:
            return "≥1,000 reads/bp"
        elif 100 <= depth < 1000:
            return "100-1,000 reads/bp"
        elif 10 <= depth < 100:
            return "10-100 reads/bp"
        elif depth < 10:
            return "0-10 reads/bp"
        else:
            return None
    
    standards_filtered = standards_filtered.copy()
    standards_filtered['bin'] = standards_filtered['total_avg_depth'].apply(assign_bin)
    
    # Group by unique_ID and calculate predictions
    def predict_for_group(group):
        depth = group['total_avg_depth'].iloc[0]
        bin_label = group['bin'].iloc[0]
        
        # Select appropriate regression
        if depth >= 1000 and reg4 is not None:
            pred = predict_read_var(reg4, group[['avg_GC', 'total_avg_depth', 'E_rel']])
        elif 100 <= depth < 1000 and reg3 is not None:
            pred = predict_read_var(reg3, group[['avg_GC', 'total_avg_depth', 'E_rel']])
        elif 10 <= depth < 100 and reg2 is not None:
            pred = predict_read_var(reg2, group[['avg_GC', 'total_avg_depth', 'E_rel']])
        elif depth < 10 and reg1 is not None:
            pred = predict_read_var(reg1, group[['avg_GC', 'total_avg_depth', 'E_rel']])
        else:
            return pd.Series(dtype=float)
        
        # Calculate RMSE for this group
        rmse = np.sqrt(np.sum((pred - group['avg_depth'].values)**2) / len(group))
        norm_rmse = rmse / depth
        
        return pd.Series({
            'RMSE': rmse,
            'norm_RMSE': norm_rmse,
            'total_avg_depth': depth,
            'bin': bin_label
        })
    
    # Apply predictions group-wise
    rmse_results = standards_filtered.groupby('unique_ID').apply(predict_for_group).reset_index()
    
    # Add sample info
    sample_info = standards_filtered[['unique_ID', 'sample', 'downsample', 'ID']].drop_duplicates('unique_ID')
    rmse_results = rmse_results.merge(sample_info, on='unique_ID', how='left')
    
    # Set categorical bin order
    rmse_results['bin'] = pd.Categorical(
        rmse_results['bin'],
        categories=["0-10 reads/bp", "10-100 reads/bp", "100-1,000 reads/bp", "≥1,000 reads/bp"],
        ordered=True
    )
    
    return rmse_results

def calculate_standardized_residuals(model):
    """
    Calculate standardized residuals from fitted model.
    Approximation of R's rstandard() function.
    
    Parameters:
    -----------
    model : dict
        Fitted model
    
    Returns:
    --------
    np.array : Standardized residuals
    """
    residuals = model['residuals']
    rmse = model['rmse']
    return residuals / rmse


def main(standards_file, output_dir="Regressions", E_detect_model_path=None):
    """
    Main workflow to build read depth variability regressions.
    
    Parameters:
    -----------
    standards_file : str
        Path to standards CSV/TSV with required columns
    output_dir : str
        Output directory for models and plots
    E_detect_model_path : str
        Path to detection threshold model (JSON)
    """
    
    # Create output directories
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(f"{output_dir}/read_depth_variability").mkdir(parents=True, exist_ok=True)
    Path(f"{output_dir}/threshold_read_depth_variability").mkdir(parents=True, exist_ok=True)
    
    # Load standards data
    print("Loading standards data...")
    if standards_file.endswith('.csv'):
        standards = pd.read_csv(standards_file)
    else:
        standards = pd.read_csv(standards_file, sep='\t')
    
    # Ensure numeric columns
    standards['total_avg_depth'] = pd.to_numeric(standards['total_avg_depth'], errors='coerce')
    standards['E_rel'] = pd.to_numeric(standards['E_rel'], errors='coerce')
    standards['avg_depth'] = pd.to_numeric(standards['avg_depth'], errors='coerce')
    standards['avg_GC'] = pd.to_numeric(standards['avg_GC'], errors='coerce')
    
    # Load and apply detection threshold filter
    print("Filtering by detection threshold...")
    if E_detect_model_path and Path(E_detect_model_path).exists():
        E_detect_model = load_detection_model(E_detect_model_path)
        standards['detect_threshold'] = predict_E_detect(E_detect_model, standards['length'])
        standards_filtered = standards[standards['E_rel'] > standards['detect_threshold']].copy()
    else:
        print("Warning: Detection model not provided; using all standards")
        standards_filtered = standards.copy()
    
    print(f"Retained {len(standards_filtered)} standards after detection filtering")
    
    # Fit regressions for each depth bin
    print("Fitting read depth variability regressions...")
    
    # Bin 4: >= 1000 reads/bp
    bin4_data = standards_filtered[standards_filtered['total_avg_depth'] >= 1000].copy()
    if len(bin4_data) > 5:
        reg4 = fit_read_var_regression(
            bin4_data[['avg_GC', 'total_avg_depth', 'E_rel']],
            bin4_data['avg_depth'].values,
            'bin4'
        )
    else:
        print("Warning: Insufficient data for bin4 (>= 1000)")
        reg4 = None
    
    # Bin 3: 100-1000 reads/bp
    bin3_data = standards_filtered[
        (standards_filtered['total_avg_depth'] >= 100) & 
        (standards_filtered['total_avg_depth'] < 1000)
    ].copy()
    if len(bin3_data) > 5:
        reg3 = fit_read_var_regression(
            bin3_data[['avg_GC', 'total_avg_depth', 'E_rel']],
            bin3_data['avg_depth'].values,
            'bin3'
        )
    else:
        print("Warning: Insufficient data for bin3 (100-1000)")
        reg3 = None
    
    # Bin 2: 10-100 reads/bp
    bin2_data = standards_filtered[
        (standards_filtered['total_avg_depth'] >= 10) & 
        (standards_filtered['total_avg_depth'] < 100)
    ].copy()
    if len(bin2_data) > 5:
        reg2 = fit_read_var_regression(
            bin2_data[['avg_GC', 'total_avg_depth', 'E_rel']],
            bin2_data['avg_depth'].values,
            'bin2'
        )
    else:
        print("Warning: Insufficient data for bin2 (10-100)")
        reg2 = None
    
    # Bin 1: 0-10 reads/bp
    bin1_data = standards_filtered[standards_filtered['total_avg_depth'] < 10].copy()
    if len(bin1_data) > 5:
        reg1 = fit_read_var_regression(
            bin1_data[['avg_GC', 'total_avg_depth', 'E_rel']],
            bin1_data['avg_depth'].values,
            'bin1'
        )
    else:
        print("Warning: Insufficient data for bin1 (0-10)")
        reg1 = None
    
    # RMSE analysis by unique standard ID
    print("Calculating RMSE across samples...")
    standards_filtered['unique_ID'] = standards_filtered['sample'] + "_" + standards_filtered['ID'].astype(str)
    
    RMSE_results = calculate_rmse(standards_filtered, reg1, reg2, reg3, reg4)
    
    RMSE_results['bin'] = pd.Categorical(
        RMSE_results['bin'],
        categories=["0-10 reads/bp", "10-100 reads/bp", "100-1,000 reads/bp", "≥1,000 reads/bp"],
        ordered=True
    )
    
    # Boxplot of normalized RMSE
    print("Creating boxplot of normalized RMSE...")
    fig, ax = plt.subplots(figsize=(10, 6))
    RMSE_results.boxplot(column='norm_RMSE', by='bin', ax=ax)
    ax.set_xlabel('Read Depth Bin')
    ax.set_ylabel('Normalized RMSE')
    ax.set_title('Normalized RMSE by Read Depth Bin')
    plt.suptitle('')  # Remove automatic title
    plt.tight_layout()
    plt.savefig(f"{output_dir}/read_depth_variability/boxplot_regressions.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    # RMSE scatter plot with thresholds
    print("Creating RMSE vs depth scatter plot...")
    fig, ax = plt.subplots(figsize=(10, 7))
    for downsample in RMSE_results['downsample'].unique():
        subset = RMSE_results[RMSE_results['downsample'] == downsample]
        ax.scatter(subset['total_avg_depth'], subset['RMSE'], alpha=0.6, label=f"D={downsample}")
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Average Read Depth (reads/bp)')
    ax.set_ylabel('RMSE')
    ax.legend()
    ax.axvline(x=10, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=100, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=1000, linestyle='--', color='gray', alpha=0.5)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/read_depth_variability/RMSE_scatter.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    # Q-Q plots for standardized residuals
    print("Creating Q-Q plots...")
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    if reg1 is not None:
        stdres1 = calculate_standardized_residuals(reg1)
        stats.probplot(stdres1, dist="norm", plot=axes[0])
        axes[0].set_title('0-10 reads/bp')
    
    if reg2 is not None:
        stdres2 = calculate_standardized_residuals(reg2)
        stats.probplot(stdres2, dist="norm", plot=axes[1])
        axes[1].set_title('10-100 reads/bp')
    
    if reg3 is not None:
        stdres3 = calculate_standardized_residuals(reg3)
        stats.probplot(stdres3, dist="norm", plot=axes[2])
        axes[2].set_title('100-1,000 reads/bp')
    
    if reg4 is not None:
        stdres4 = calculate_standardized_residuals(reg4)
        stats.probplot(stdres4, dist="norm", plot=axes[3])
        axes[3].set_title('≥1,000 reads/bp')
    
    for ax in axes:
        ax.set_xlabel('Normal Scores')
        ax.set_ylabel('Standardized Residuals')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/read_depth_variability/qq_plots_regressions.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    # Save models
    print("Saving models...")
    models = {
        'reg1': reg1,
        'reg2': reg2,
        'reg3': reg3,
        'reg4': reg4
    }
    
    for name, model in models.items():
        if model is not None:
            with open(f"{output_dir}/read_depth_variability/{name}.pkl", 'wb') as f:
                pickle.dump(model, f)
    
    # Build RMSE threshold functions
    print("Building RMSE threshold functions...")
    
    # Aggregate RMSE by rounded depth
    RMSE_limits = RMSE_results.copy()
    RMSE_limits['total_avg_depth_rounded'] = np.round(RMSE_limits['total_avg_depth'], 1)
    RMSE_limits = RMSE_limits.groupby('total_avg_depth_rounded').agg({
        'RMSE': 'max'
    }).reset_index()
    RMSE_limits.columns = ['total_avg_depth', 'RMSE']
    
    # Fit RMSE limit functions for each bin
    rmse_limit_funcs = {}
    
    # Bin 1: 0-10
    bin1_limits = RMSE_limits[RMSE_limits['total_avg_depth'] < 10]
    if len(bin1_limits) > 2:
        log_rmse = np.log(bin1_limits['RMSE'].values) + 0.25
        log_depth = np.log(bin1_limits['total_avg_depth'].values)
        slope, intercept, _, _, _ = linregress(log_depth, log_rmse)
        rmse_limit_funcs['func1'] = {'intercept': intercept, 'slope': slope, 'type': 'log_linear'}
    
    # Bin 2: 10-100
    bin2_limits = RMSE_limits[(RMSE_limits['total_avg_depth'] >= 10) & 
                               (RMSE_limits['total_avg_depth'] < 100)]
    if len(bin2_limits) > 2:
        log_rmse = np.log(bin2_limits['RMSE'].values) + 0.25
        log_depth = np.log(bin2_limits['total_avg_depth'].values)
        slope, intercept, _, _, _ = linregress(log_depth, log_rmse)
        rmse_limit_funcs['func2'] = {'intercept': intercept, 'slope': slope, 'type': 'log_linear'}
    
    # Bin 3: 100-1000
    bin3_limits = RMSE_limits[(RMSE_limits['total_avg_depth'] >= 100) & 
                               (RMSE_limits['total_avg_depth'] < 1000)]
    if len(bin3_limits) > 2:
        log_rmse = np.log(bin3_limits['RMSE'].values) + 0.25
        log_depth = np.log(bin3_limits['total_avg_depth'].values)
        slope, intercept, _, _, _ = linregress(log_depth, log_rmse)
        rmse_limit_funcs['func3'] = {'intercept': intercept, 'slope': slope, 'type': 'log_linear'}
    
    # Bin 4: >= 1000 (use constant)
    bin4_limits = RMSE_limits[RMSE_limits['total_avg_depth'] >= 1000]
    if len(bin4_limits) > 0:
        rmse_limit_funcs['func4'] = {'value': float(bin4_limits['RMSE'].max() + np.exp(0.25)), 'type': 'constant'}
    
    # Save RMSE limit functions
    with open(f"{output_dir}/threshold_read_depth_variability/rmse_limit_functions.json", 'w') as f:
        json.dump(rmse_limit_funcs, f, indent=2)
    
    # Visualize RMSE thresholds
    print("Creating RMSE threshold visualization...")
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for downsample in RMSE_results['downsample'].unique():
        subset = RMSE_results[RMSE_results['downsample'] == downsample]
        ax.scatter(subset['total_avg_depth'], subset['RMSE'], alpha=0.6, label=f"D={downsample}")
    
    # Plot threshold lines
    depths_log = np.logspace(0, 4, 100)
    
    # Bin 1
    if 'func1' in rmse_limit_funcs:
        func = rmse_limit_funcs['func1']
        preds = np.exp(func['intercept'] + func['slope'] * np.log(np.minimum(depths_log, 10)))
        ax.plot(depths_log[depths_log < 10], preds[depths_log < 10], 'r-', linewidth=2)
    
    # Bin 2
    if 'func2' in rmse_limit_funcs:
        func = rmse_limit_funcs['func2']
        mask = (depths_log >= 10) & (depths_log < 100)
        preds = np.exp(func['intercept'] + func['slope'] * np.log(depths_log[mask]))
        ax.plot(depths_log[mask], preds, 'r-', linewidth=2)
    
    # Bin 3
    if 'func3' in rmse_limit_funcs:
        func = rmse_limit_funcs['func3']
        mask = (depths_log >= 100) & (depths_log < 1000)
        preds = np.exp(func['intercept'] + func['slope'] * np.log(depths_log[mask]))
        ax.plot(depths_log[mask], preds, 'r-', linewidth=2)
    
    # Bin 4
    if 'func4' in rmse_limit_funcs:
        func = rmse_limit_funcs['func4']
        ax.axhline(y=func['value'], xmin=0.75, color='r', linestyle='-', linewidth=2)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Average Read Depth (reads/bp)')
    ax.set_ylabel('RMSE')
    ax.axvline(x=10, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=100, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(x=1000, linestyle='--', color='gray', alpha=0.5)
    ax.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/threshold_read_depth_variability/RMSE_thresholds.png", dpi=400, bbox_inches='tight')
    plt.close()
    
    print("Done! Models and visualizations saved.")
    
    return {
        'regressions': models,
        'RMSE_results': RMSE_results,
        'RMSE_limits': RMSE_limits,
        'RMSE_limit_functions': rmse_limit_funcs
    }


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python quant_correct_regression_builder.py <standards_file> [output_dir] [E_detect_model]")
        sys.exit(1)
    
    standards_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "Regressions"
    E_detect_model = sys.argv[3] if len(sys.argv) > 3 else None
    
    main(standards_file, output_dir, E_detect_model)
