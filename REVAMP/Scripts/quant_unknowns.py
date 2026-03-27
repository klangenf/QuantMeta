#!/usr/bin/env python
"""
Quantification of Unknown Targets

This script quantifies unknown targets in metagenomic data by:
1. Assessing detection thresholds for standards
2. Correcting for non-specific mapping artifacts
3. Converting relative abundances to absolute concentrations

Based on the original quant_unknowns.R script.
"""

import pandas as pd
import numpy as np
import pickle
import json
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')


def load_detection_model(model_path):
    """
    Load detection threshold model.

    Parameters:
    -----------
    model_path : str
        Path to detection model (JSON format)

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


def load_pickle_model(model_path):
    """
    Load a pickle model file.

    Parameters:
    -----------
    model_path : str
        Path to pickle file

    Returns:
    --------
    dict : Model data
    """
    try:
        with open(model_path, 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        print(f"Warning: Model not found at {model_path}")
        return None


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


def predict_rmse_limit(rmse_model, depth):
    """
    Predict RMSE limit for a given depth.

    Parameters:
    -----------
    rmse_model : dict
        RMSE limit model
    depth : float
        Average read depth

    Returns:
    --------
    float : RMSE limit
    """
    if rmse_model['type'] == 'constant':
        return rmse_model['value']
    elif rmse_model['type'] == 'log_linear':
        return np.exp(rmse_model['intercept'] + rmse_model['slope'] * np.log(depth))
    else:
        return 0


def detection_threshold(mapping_results_file, lengths, detect_thresh):
    """
    Determine coverage, average read depth, and detection threshold parameters.

    Parameters:
    -----------
    mapping_results_file : str
        Path to mapping results TSV
    lengths_file : str
        Path to sequence lengths TSV
    detect_thresh : dict
        Length-dependent entropy detection threshold model

    Returns:
    --------
    pd.DataFrame : Detection analysis results
    """
    # Load data
    mapping_results = pd.read_csv(mapping_results_file, sep='\t', header=0)
    mapping_results = mapping_results.astype({
        'ID': str,
        'read_depth': float
    })

    # mapping_results: DataFrame with ID, read_depth
    # lengths: DataFrame with ID, length

    # Start with unique IDs and merge lengths
    targets = pd.DataFrame({'ID': mapping_results['ID'].unique()})
    targets = targets.merge(lengths[['ID', 'length']], on='ID', how='left')
    
    # Group mapping_results by ID for efficient aggregation
    grouped = mapping_results.groupby('ID')
    
    # Compute basic aggregations
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
    targets = targets.merge(
        pd.DataFrame({
            'ID': B_G.index,
            'B_G': B_G.values,
            'gene_copies': gene_copies.values,
            'I_G': I_G.values
        }),
        on='ID',
        how='left'
    ).fillna({'B_G': 0.0, 'gene_copies': 0.0, 'I_G': 0.0})
    
    # Compute remaining metrics using vectorized operations
    targets = targets.assign(
        E_rel=lambda df: df['I_G'] / np.log(df['length']).where(df['length'] > 1, 1),
        E_detect=predict_E_detect(detect_thresh, df['length'].astype(float).values),
        detection_status=lambda df: np.where(df['E_rel'] >= df['E_detect'], 'detected', 'not_detected')
    )

    out_table = targets[['ID', 'gene_copies', 'E_rel', 'E_detect', 'detection_status']].copy()

    return out_table


def quant_correction_v5(input_info, sliding_window, quad_reg1, quad_reg2, quad_reg3, quad_reg4,
                       cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4,
                       output_file=None):
    """
    Detect and correct non-specific mapping for database targets.

    Parameters:
    -----------
    input_info : pd.DataFrame
        Target information with RMSE and depth data
    sliding_window : pd.DataFrame
        Sliding window read depth data
    quad_reg1-4 : dict
        Read depth variability regression models
    cutoff_function1-4 : dict
        RMSE threshold functions
    output_file : str, optional
        Path to save corrected mapping data

    Returns:
    --------
    pd.DataFrame : Corrected target information
    """
    # Calculate RMSE limits
    input_info = input_info.copy()
    input_info['RMSE_limit'] = 0.0

    # Apply appropriate RMSE limits based on depth
    mask_1000_plus = input_info['total_avg_depth'] >= 1000
    mask_100_1000 = (input_info['total_avg_depth'] >= 100) & (input_info['total_avg_depth'] < 1000)
    mask_10_100 = (input_info['total_avg_depth'] >= 10) & (input_info['total_avg_depth'] < 100)
    mask_under_10 = input_info['total_avg_depth'] < 10

    input_info.loc[mask_1000_plus, 'RMSE_limit'] = predict_rmse_limit(cutoff_function4, input_info.loc[mask_1000_plus, 'total_avg_depth'])
    input_info.loc[mask_100_1000, 'RMSE_limit'] = predict_rmse_limit(cutoff_function3, input_info.loc[mask_100_1000, 'total_avg_depth'])
    input_info.loc[mask_10_100, 'RMSE_limit'] = predict_rmse_limit(cutoff_function2, input_info.loc[mask_10_100, 'total_avg_depth'])
    input_info.loc[mask_under_10, 'RMSE_limit'] = predict_rmse_limit(cutoff_function1, input_info.loc[mask_under_10, 'total_avg_depth'])

    # For log-linear models, exponentiate the predictions
    input_info.loc[input_info['total_avg_depth'] < 1000, 'RMSE_limit'] = \
        np.exp(input_info.loc[input_info['total_avg_depth'] < 1000, 'RMSE_limit'])

    print(f"Targets requiring correction: {len(input_info[input_info['RMSE'] > input_info['RMSE_limit']])}")

    # Save initial RMSE
    input_info['prev_RMSE'] = input_info['RMSE']
    input_info['frac_corrected'] = 0.0
    input_info['correction_status'] = 'correctable'

    # Get IDs requiring correction
    correction_ids = input_info[input_info['RMSE'] > input_info['RMSE_limit']]['ID'].unique()

    corrected_mapping = []

    if len(correction_ids) > 0:
        for target_id in correction_ids:
            print(f"Correcting target: {target_id}")

            # Get sliding window data for this target
            wind_clean = sliding_window[sliding_window['ID'] == target_id].copy()
            wind_clean['initial_avg_depth'] = wind_clean['avg_depth']

            cycle = 1
            max_cycles = 20

            # Iterative correction loop
            while (input_info.loc[input_info['ID'] == target_id, 'RMSE'].iloc[0] >
                   input_info.loc[input_info['ID'] == target_id, 'RMSE_limit'].iloc[0] and
                   cycle <= max_cycles and
                   input_info.loc[input_info['ID'] == target_id, 'total_avg_depth'].iloc[0] > 0):

                # Select appropriate regression model
                depth = wind_clean['total_avg_depth'].iloc[0]

                if depth < 10 and quad_reg1 is not None:
                    pred = predict_read_var(quad_reg1, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 10 <= depth < 100 and quad_reg2 is not None:
                    pred = predict_read_var(quad_reg2, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 100 <= depth < 1000 and quad_reg3 is not None:
                    pred = predict_read_var(quad_reg3, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred)
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif depth >= 1000 and quad_reg4 is not None:
                    pred = predict_read_var(quad_reg4, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                else:
                    break

                # Ensure non-negative bounds
                lower = np.maximum(lower, 0)
                upper = np.maximum(upper, 0)

                # Identify non-specific regions
                wind_clean['nonspec_region'] = 'NO'
                nonspec_mask = ((wind_clean['avg_depth'] > upper) |
                               (wind_clean['avg_depth'] < lower)) & (wind_clean['avg_depth'] != 0)
                wind_clean.loc[nonspec_mask, 'nonspec_region'] = 'YES'

                # Recalculate average depth excluding non-specific regions
                specific_regions = wind_clean[wind_clean['nonspec_region'] == 'NO']
                if len(specific_regions) > 0:
                    wind_clean['total_avg_depth'] = specific_regions['avg_depth'].mean()
                else:
                    wind_clean['total_avg_depth'] = 0
                    break

                # Recalculate predictions with updated depth
                depth = wind_clean['total_avg_depth'].iloc[0]

                if depth < 10 and quad_reg1 is not None:
                    pred = predict_read_var(quad_reg1, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 10 <= depth < 100 and quad_reg2 is not None:
                    pred = predict_read_var(quad_reg2, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 100 <= depth < 1000 and quad_reg3 is not None:
                    pred = predict_read_var(quad_reg3, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred)
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif depth >= 1000 and quad_reg4 is not None:
                    pred = predict_read_var(quad_reg4, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                # Ensure non-negative bounds
                lower = np.maximum(lower, 0)
                upper = np.maximum(upper, 0)

                # Apply corrections
                wind_clean.loc[wind_clean['nonspec_region'] == 'YES', 'avg_depth'] = \
                    np.clip(wind_clean.loc[wind_clean['nonspec_region'] == 'YES', 'avg_depth'], lower, upper)

                # Update total average depth
                wind_clean['total_avg_depth'] = wind_clean['avg_depth'].mean()
                input_info.loc[input_info['ID'] == target_id, 'total_avg_depth'] = wind_clean['total_avg_depth'].iloc[0]

                # Recalculate RMSE
                depth = wind_clean['total_avg_depth'].iloc[0]

                if depth < 10 and quad_reg1 is not None:
                    pred = predict_read_var(quad_reg1, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = predict_rmse_limit(cutoff_function1, depth)
                    rmse_limit = np.exp(rmse_limit)

                elif 10 <= depth < 100 and quad_reg2 is not None:
                    pred = predict_read_var(quad_reg2, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = predict_rmse_limit(cutoff_function2, depth)
                    rmse_limit = np.exp(rmse_limit)

                elif 100 <= depth < 1000 and quad_reg3 is not None:
                    pred = predict_read_var(quad_reg3, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred)
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = predict_rmse_limit(cutoff_function3, depth)
                    rmse_limit = np.exp(rmse_limit)

                elif depth >= 1000 and quad_reg4 is not None:
                    pred = predict_read_var(quad_reg4, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = cutoff_function4['value']

                input_info.loc[input_info['ID'] == target_id, 'RMSE'] = rmse
                input_info.loc[input_info['ID'] == target_id, 'RMSE_limit'] = rmse_limit

                print(f"Cycle {cycle}: RMSE = {rmse:.4f}, Limit = {rmse_limit:.4f}")
                cycle += 1

            # Calculate fraction corrected
            frac_corrected = np.sum(wind_clean['avg_depth'] != wind_clean['initial_avg_depth']) / len(wind_clean)
            input_info.loc[input_info['ID'] == target_id, 'frac_corrected'] = frac_corrected

            if frac_corrected > 0.2:
                input_info.loc[input_info['ID'] == target_id, 'correction_status'] = 'error'

            # Keep only relevant columns
            wind_clean = wind_clean[[
                'ID', 'avg_GC', 'avg_depth', 'E_rel', 'total_avg_depth',
                'E_detect', 'detection_status', 'length', 'initial_avg_depth'
            ]]

            corrected_mapping.append(wind_clean)

        # Combine all corrected mapping data
        if corrected_mapping:
            corrected_mapping_df = pd.concat(corrected_mapping, ignore_index=True)
            if output_file:
                corrected_mapping_df.to_csv(output_file, sep='\t', index=False)
        else:
            corrected_mapping_df = pd.DataFrame({"message": ["No targets required correction!"]})
            if output_file:
                corrected_mapping_df.to_csv(output_file, index=False, header=False)

    return input_info


def quant_correction_contigs(input_info, sliding_window, bin_assignment, quad_reg1, quad_reg2, quad_reg3, quad_reg4,
                           cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4,
                           output_file=None):
    """
    Detect and correct non-specific mapping for contig targets.

    Parameters:
    -----------
    input_info : pd.DataFrame
        Target information with RMSE and depth data
    sliding_window : pd.DataFrame
        Sliding window read depth data
    bin_assignment : str or None
        Path to bin assignment file
    quad_reg1-4 : dict
        Read depth variability regression models
    cutoff_function1-4 : dict
        RMSE threshold functions
    output_file : str, optional
        Path to save corrected mapping data

    Returns:
    --------
    pd.DataFrame : Corrected target information
    """
    # Calculate RMSE limits (same as database version)
    input_info = input_info.copy()
    input_info['RMSE_limit'] = 0.0

    mask_1000_plus = input_info['total_avg_depth'] >= 1000
    mask_100_1000 = (input_info['total_avg_depth'] >= 100) & (input_info['total_avg_depth'] < 1000)
    mask_10_100 = (input_info['total_avg_depth'] >= 10) & (input_info['total_avg_depth'] < 100)
    mask_under_10 = input_info['total_avg_depth'] < 10

    input_info.loc[mask_1000_plus, 'RMSE_limit'] = predict_rmse_limit(cutoff_function4, input_info.loc[mask_1000_plus, 'total_avg_depth'])
    input_info.loc[mask_100_1000, 'RMSE_limit'] = predict_rmse_limit(cutoff_function3, input_info.loc[mask_100_1000, 'total_avg_depth'])
    input_info.loc[mask_10_100, 'RMSE_limit'] = predict_rmse_limit(cutoff_function2, input_info.loc[mask_10_100, 'total_avg_depth'])
    input_info.loc[mask_under_10, 'RMSE_limit'] = predict_rmse_limit(cutoff_function1, input_info.loc[mask_under_10, 'total_avg_depth'])

    input_info.loc[input_info['total_avg_depth'] < 1000, 'RMSE_limit'] = \
        np.exp(input_info.loc[input_info['total_avg_depth'] < 1000, 'RMSE_limit'])

    print(f"Targets requiring correction: {len(input_info[input_info['RMSE'] > input_info['RMSE_limit']])}")

    # Save initial RMSE
    input_info['prev_RMSE'] = input_info['RMSE']
    input_info['frac_corrected'] = 0.0
    input_info['correction_status'] = 'correctable'

    # Get IDs requiring correction
    correction_ids = input_info[input_info['RMSE'] > input_info['RMSE_limit']]['ID'].unique()

    corrected_mapping = []

    if len(correction_ids) > 0:
        for target_id in correction_ids:
            print(f"Correcting contig: {target_id}")

            # Get sliding window data for this target
            wind_clean = sliding_window[sliding_window['ID'] == target_id].copy()
            wind_clean['initial_avg_depth'] = wind_clean['avg_depth']

            cycle = 1
            max_cycles = 20

            # Iterative correction loop
            while (input_info.loc[input_info['ID'] == target_id, 'RMSE'].iloc[0] >
                   input_info.loc[input_info['ID'] == target_id, 'RMSE_limit'].iloc[0] and
                   cycle <= max_cycles and
                   input_info.loc[input_info['ID'] == target_id, 'total_avg_depth'].iloc[0] > 0):

                # Select appropriate regression model
                depth = wind_clean['total_avg_depth'].iloc[0]

                if depth < 10 and quad_reg1 is not None:
                    pred = predict_read_var(quad_reg1, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 10 <= depth < 100 and quad_reg2 is not None:
                    pred = predict_read_var(quad_reg2, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 100 <= depth < 1000 and quad_reg3 is not None:
                    pred = predict_read_var(quad_reg3, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred)
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif depth >= 1000 and quad_reg4 is not None:
                    pred = predict_read_var(quad_reg4, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    sd_mult = np.std(wind_clean['avg_depth'])
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                else:
                    break

                # Ensure non-negative bounds
                lower = np.maximum(lower, 0)
                upper = np.maximum(upper, 0)

                # Identify non-specific regions
                wind_clean['nonspec_region'] = 'NO'
                nonspec_mask = ((wind_clean['avg_depth'] > upper) |
                               (wind_clean['avg_depth'] < lower)) & (wind_clean['avg_depth'] != 0)
                wind_clean.loc[nonspec_mask, 'nonspec_region'] = 'YES'

                # Check for convergence issues
                specific_regions = wind_clean[wind_clean['nonspec_region'] == 'NO']
                all_regions = wind_clean
                if len(specific_regions) == 0 or len(specific_regions) == len(all_regions):
                    cycle = 21  # Force exit
                    input_info.loc[input_info['ID'] == target_id, 'correction_status'] = 'error_no convergence'
                    print("No convergence to a solution, evaluate contig(s) quality")
                    break

                # Recalculate average depth excluding non-specific regions
                wind_clean['total_avg_depth'] = specific_regions['avg_depth'].mean()

                # Recalculate predictions with updated depth
                depth = wind_clean['total_avg_depth'].iloc[0]

                if depth < 10 and quad_reg1 is not None:
                    pred = predict_read_var(quad_reg1, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 10 <= depth < 100 and quad_reg2 is not None:
                    pred = predict_read_var(quad_reg2, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif 100 <= depth < 1000 and quad_reg3 is not None:
                    pred = predict_read_var(quad_reg3, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred)
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                elif depth >= 1000 and quad_reg4 is not None:
                    pred = predict_read_var(quad_reg4, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    specific_depths = wind_clean.loc[wind_clean['nonspec_region'] == 'NO', 'avg_depth']
                    sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
                    upper = pred + 1.5 * sd_mult
                    lower = pred - 1.5 * sd_mult

                # Ensure non-negative bounds
                lower = np.maximum(lower, 0)
                upper = np.maximum(upper, 0)

                # Apply corrections
                wind_clean.loc[wind_clean['nonspec_region'] == 'YES', 'avg_depth'] = \
                    np.clip(wind_clean.loc[wind_clean['nonspec_region'] == 'YES', 'avg_depth'], lower, upper)

                # Update total average depth
                wind_clean['total_avg_depth'] = wind_clean['avg_depth'].mean()
                input_info.loc[input_info['ID'] == target_id, 'total_avg_depth'] = wind_clean['total_avg_depth'].iloc[0]

                # Recalculate RMSE
                depth = wind_clean['total_avg_depth'].iloc[0]

                if depth < 10 and quad_reg1 is not None:
                    pred = predict_read_var(quad_reg1, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = predict_rmse_limit(cutoff_function1, depth)
                    rmse_limit = np.exp(rmse_limit)

                elif 10 <= depth < 100 and quad_reg2 is not None:
                    pred = predict_read_var(quad_reg2, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred) - 1
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = predict_rmse_limit(cutoff_function2, depth)
                    rmse_limit = np.exp(rmse_limit)

                elif 100 <= depth < 1000 and quad_reg3 is not None:
                    pred = predict_read_var(quad_reg3, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    pred = np.exp(pred)
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = predict_rmse_limit(cutoff_function3, depth)
                    rmse_limit = np.exp(rmse_limit)

                elif depth >= 1000 and quad_reg4 is not None:
                    pred = predict_read_var(quad_reg4, wind_clean[['avg_GC', 'total_avg_depth', 'E_rel']])
                    rmse = np.sqrt(np.sum((pred - wind_clean['avg_depth'])**2) / len(wind_clean))
                    rmse_limit = cutoff_function4['value']

                input_info.loc[input_info['ID'] == target_id, 'RMSE'] = rmse
                input_info.loc[input_info['ID'] == target_id, 'RMSE_limit'] = rmse_limit

                print(f"Cycle {cycle}: RMSE = {rmse:.4f}, Limit = {rmse_limit:.4f}")
                cycle += 1

            # Calculate fraction corrected
            frac_corrected = np.sum(wind_clean['avg_depth'] != wind_clean['initial_avg_depth']) / len(wind_clean)
            input_info.loc[input_info['ID'] == target_id, 'frac_corrected'] = frac_corrected

            if frac_corrected > 0.2:
                input_info.loc[input_info['ID'] == target_id, 'correction_status'] = 'error'

            # Keep only relevant columns (include contig_ID for contigs)
            if bin_assignment is not None:
                wind_clean = wind_clean[[
                    'contig_ID', 'ID', 'avg_GC', 'avg_depth', 'E_rel', 'total_avg_depth',
                    'E_detect', 'detection_status', 'length', 'initial_avg_depth'
                ]]
            else:
                wind_clean = wind_clean[[
                    'ID', 'avg_GC', 'avg_depth', 'E_rel', 'total_avg_depth',
                    'E_detect', 'detection_status', 'length', 'initial_avg_depth'
                ]]

            corrected_mapping.append(wind_clean)

        # Combine all corrected mapping data
        if corrected_mapping:
            corrected_mapping_df = pd.concat(corrected_mapping, ignore_index=True)
            if output_file:
                corrected_mapping_df.to_csv(output_file, sep='\t', index=False)
        else:
            corrected_mapping_df = pd.DataFrame({"message": ["No targets required correction!"]})
            if output_file:
                corrected_mapping_df.to_csv(output_file, index=False, header=False)

    return input_info


def quant_correction(sample_name, descript, sliding_window_file, lengths, results, bin_assignment,
                    target_type, quad_reg_paths, cutoff_function_paths, output_file=None):
    """
    Main correction function that orchestrates non-specific mapping detection and correction.

    Parameters:
    -----------
    sample_name : str
        Sample name
    descript : str
        Description for filenames
    sliding_window_file : str
        Path to sliding window data
    lengths : pd.DataFrame
        Sequence lengths
    results : pd.DataFrame
        Detection results
    bin_assignment : str or None
        Path to bin assignment file (for contigs)
    target_type : str
        'database' or 'contigs'
    quad_reg_paths : list
        Paths to read depth variability models
    cutoff_function_paths : list
        Paths to RMSE threshold functions
    output_file : str, optional
        Path to save corrected results

    Returns:
    --------
    pd.DataFrame : Corrected results
    """
    # Load sliding window data
    sliding_window = pd.read_csv(sliding_window_file, sep='\t')

    # Load models
    quad_reg1 = load_pickle_model(quad_reg_paths[0]) if len(quad_reg_paths) > 0 else None
    quad_reg2 = load_pickle_model(quad_reg_paths[1]) if len(quad_reg_paths) > 1 else None
    quad_reg3 = load_pickle_model(quad_reg_paths[2]) if len(quad_reg_paths) > 2 else None
    quad_reg4 = load_pickle_model(quad_reg_paths[3]) if len(quad_reg_paths) > 3 else None

    cutoff_function1 = load_pickle_model(cutoff_function_paths[0]) if len(cutoff_function_paths) > 0 else None
    cutoff_function2 = load_pickle_model(cutoff_function_paths[1]) if len(cutoff_function_paths) > 1 else None
    cutoff_function3 = load_pickle_model(cutoff_function_paths[2]) if len(cutoff_function_paths) > 2 else None
    cutoff_function4 = load_pickle_model(cutoff_function_paths[3]) if len(cutoff_function_paths) > 3 else None

    # Filter targets present in sliding window
    targets_list = lengths[lengths['ID'].isin(sliding_window['ID'])].copy()

    # Handle bin assignment for contigs
    if bin_assignment is not None:
        bins = pd.read_csv(bin_assignment, sep='\t', header=None, names=['contig_ID', 'ID'])
        results.columns = ['contig_ID', 'E_rel', 'total_avg_depth', 'E_detect', 'detection_status']
        results = results.merge(bins, on='contig_ID', how='left')
        targets_list.columns = ['contig_ID', 'length']
        results = results.merge(targets_list, on='contig_ID', how='left')
    else:
        results = results.merge(targets_list, on='ID', how='left')

    # Merge sliding window with results
    if bin_assignment is not None:
        sliding_window.columns = ['contig_ID', 'avg_GC', 'avg_depth']
        sliding_window = sliding_window.merge(results, on='contig_ID', how='left')
        sliding_window.columns = ['contig_ID', 'avg_GC', 'avg_depth', 'ID', 'E_rel', 'total_avg_depth',
                                 'E_detect', 'detection_status', 'length']
        results.columns = ['contig_ID', 'ID', 'E_rel', 'total_avg_depth', 'E_detect', 'detection_status', 'length']

        # Update total_avg_depth for sliding window
        temp = sliding_window.groupby('ID')['avg_depth'].mean().reset_index()
        temp.columns = ['ID', 'total_avg_depth']
        sliding_window = sliding_window.merge(temp, on='ID', how='left')
        sliding_window = sliding_window[[
            'contig_ID', 'avg_GC', 'avg_depth', 'ID', 'E_rel', 'total_avg_depth_y',
            'E_detect', 'detection_status', 'length'
        ]]
        sliding_window.columns = ['contig_ID', 'avg_GC', 'avg_depth', 'ID', 'E_rel', 'total_avg_depth',
                                 'E_detect', 'detection_status', 'length']
        results = results.merge(temp, on='ID', how='left')
        results = results[[
            'contig_ID', 'ID', 'E_rel', 'total_avg_depth_y', 'E_detect', 'detection_status', 'length'
        ]]
        results.columns = ['contig_ID', 'ID', 'E_rel', 'total_avg_depth', 'E_detect', 'detection_status', 'length']
    else:
        sliding_window = sliding_window.merge(results, on='ID', how='left')
        sliding_window.columns = ['ID', 'avg_GC', 'avg_depth', 'E_rel', 'total_avg_depth',
                                'E_detect', 'detection_status', 'length']
        results.columns = ['ID', 'E_rel', 'total_avg_depth', 'E_detect', 'detection_status', 'length']

    # Calculate RMSE for each target
    unique_ids = results['ID'].unique()
    for target_id in unique_ids:
        test_set = sliding_window[sliding_window['ID'] == target_id]
        depth = test_set['total_avg_depth'].iloc[0]

        if depth < 10 and quad_reg1 is not None:
            pred = predict_read_var(quad_reg1, test_set[['avg_GC', 'total_avg_depth', 'E_rel']])
            pred = np.exp(pred) - 1
        elif 10 <= depth < 100 and quad_reg2 is not None:
            pred = predict_read_var(quad_reg2, test_set[['avg_GC', 'total_avg_depth', 'E_rel']])
            pred = np.exp(pred) - 1
        elif 100 <= depth < 1000 and quad_reg3 is not None:
            pred = predict_read_var(quad_reg3, test_set[['avg_GC', 'total_avg_depth', 'E_rel']])
            pred = np.exp(pred)
        elif depth >= 1000 and quad_reg4 is not None:
            pred = predict_read_var(quad_reg4, test_set[['avg_GC', 'total_avg_depth', 'E_rel']])
        else:
            continue

        rmse = np.sqrt(np.sum((pred - test_set['avg_depth'])**2) / len(test_set))
        results.loc[results['ID'] == target_id, 'RMSE'] = rmse

    # Apply correction based on target type
    if target_type == "database":
        results_v2 = quant_correction_v5(results, sliding_window, quad_reg1, quad_reg2, quad_reg3, quad_reg4,
                                       cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4)
    else:
        results_v2 = quant_correction_contigs(results, sliding_window, bin_assignment, quad_reg1, quad_reg2, quad_reg3, quad_reg4,
                                            cutoff_function1, cutoff_function2, cutoff_function3, cutoff_function4)

    # Add reliability assessment
    results_v2['reliability'] = 'low conf/high error'
    results_v2.loc[results_v2['RMSE'] <= results_v2['RMSE_limit'], 'reliability'] = 'high conf/low error'

    # Standardize column names
    if bin_assignment is not None:
        results_v2.columns = ['ID', 'contig_ID', 'E_rel', 'gene_copies', 'E_detect',
                             'detection_status', 'length', 'RMSE', 'RMSE_limit', 'prev_RMSE',
                             'frac_corrected', 'correction_status', 'reliability']
    else:
        results_v2.columns = ['ID', 'E_rel', 'gene_copies', 'E_detect',
                             'detection_status', 'length', 'RMSE', 'RMSE_limit', 'prev_RMSE',
                             'frac_corrected', 'correction_status', 'reliability']

    # Save results
    if output_file:
        results_v2.to_csv(output_file, sep='\t', index=False)

    return results_v2


def prediction(mapping, DNA_input, DNA_conc, target_type):
    """
    Calculate predicted concentration (gene copies/µL DNA extract).

    Parameters:
    -----------
    mapping : pd.DataFrame
        Mapping results with gene_copies
    DNA_input : float
        Mass of DNA used in library prep (ng)
    DNA_conc : float
        DNA concentration in extract (ng/µL)
    target_type : str
        'database' or 'contigs'

    Returns:
    --------
    pd.DataFrame : Predicted concentrations
    """
    def pred_conc_calc(gc, mass, conc):
        """Convert gene copies to concentration in DNA extract."""
        return gc * conc / mass

    if target_type == "database":
        result = pd.DataFrame({
            'ID': mapping['ID'],
            'predicted_conc': pred_conc_calc(mapping['gene_copies'], DNA_input, DNA_conc)
        })
    else:
        result = pd.DataFrame({
            'ID': mapping['ID'],
            'contig_ID': mapping['contig_ID'],
            'predicted_conc': pred_conc_calc(mapping['gene_copies'], DNA_input, DNA_conc)
        })

    return result

def quant_unknown(sample_name, target_name, database_lengths, quant_regression_path, mapping_results, sliding_window_file,
                 bin_assignment, target_type, DNA_input, DNA_conc,
                 E_detect_model_path, quad_reg_paths, cutoff_function_paths,
                 output_files=None):
    """
    Main function to quantify unknown targets.

    Parameters:
    -----------
    sample_name : str
        Sample name for output files
    target_name : str
        Target name for output files
    database_lengths : str
        Path to sequence lengths file
    quant_regression_path : str
        Path to quantification regression model
    mapping_results : str
        Path to mapping results file
    sliding_window_file : str
        Path to sliding window data
    bin_assignment : str or None
        Path to bin assignment file (for contigs)
    target_type : str
        'database' or 'contigs'
    DNA_input : float
        DNA mass used in library prep (ng)
    DNA_conc : float
        DNA concentration (ng/µL)
    E_detect_model_path : str
        Path to detection model
    quad_reg_paths : list
        Paths to read depth variability models
    cutoff_function_paths : list
        Paths to RMSE threshold functions

    Returns:
    --------
    pd.DataFrame : Final quantification results
    """
    # Load sequence lengths
    lengths = pd.read_csv(database_lengths, sep='\t', header=None, names=['ID', 'length'])

    # Load quantification regression model
    quant_rel = load_pickle_model(quant_regression_path)

    # Assess detection thresholds
    detect_thresh = load_detection_model(E_detect_model_path)

    mapping = detection_threshold(mapping_results, lengths, detect_thresh)

    output_path = Path('Mapping')/ sample_name / target_name '_mapping_analysis.txt'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(output_path, sep='\t', index=False)

    # Filter targets above detection threshold
    mapping = mapping[mapping['E_rel'] >= mapping['E_detect']]

    # Detect and correct non-specific mapping
    if len(mapping) > 0:
        mapping = quant_correction(
            "sample", "descript", sliding_window_file, lengths, mapping, bin_assignment, target_type,
            quad_reg_paths, cutoff_function_paths,
            output_files.get('corrected_mapping') if output_files else None
        )

        # Remove targets that are not quantifiable
        mapping = mapping[~mapping['correction_status'].isin([
            'error_no convergence', 'error_high fraction corrected'
        ])]
    else:
        print("No targets above detection threshold!")
        if output_files:
            pd.DataFrame({"message": ["No targets above detection!"]}).to_csv(
                output_files.get('corrected_mapping', ''), index=False, header=False
            )
            pd.DataFrame({"message": ["No targets above detection!"]}).to_csv(
                output_files.get('corrected_results', ''), index=False, header=False
            )

    # Convert relative abundances to absolute concentrations
    results = prediction(mapping, DNA_input, DNA_conc, target_type)

    # Convert to absolute abundance using regression model
    pred = 10 ** (quant_rel['intercept'] + quant_rel['slope'] * np.log10(results['predicted_conc']))
    pred_se = quant_rel['slope'] * np.log(10) * pred  # Standard error approximation

    if target_type == "database":
        results_final = pd.DataFrame({
            'ID': results['ID'],
            'concentration (gc/µL)': pred,
            'std deviation (gc/µL)': pred_se
        })
    else:
        results_final = pd.DataFrame({
            'ID': results['ID'],
            'contig_ID': results['contig_ID'],
            'concentration (gc/µL)': pred,
            'std deviation (gc/µL)': pred_se
        })

    # Save final results
    if output_files and 'final' in output_files:
        results_final.to_csv(output_files['final'], sep='\t', index=False)

    return results_final


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Quantify unknown targets in metagenomic data')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--target_name', required=True, help='Target name')
    parser.add_argument('--database-lengths', required=True, help='Path to sequence lengths file')
    parser.add_argument('--quant-regression', required=True, help='Path to quantification regression model')
    parser.add_argument('--mapping-results', required=True, help='Path to mapping results file')
    parser.add_argument('--sliding-window', required=True, help='Path to sliding window data')
    parser.add_argument('--bin-assignment', help='Path to bin assignment file (for contigs)')
    parser.add_argument('--target-type', required=True, choices=['database', 'contigs'], help='Target type')
    parser.add_argument('--dna-input', type=float, required=True, help='DNA mass used in library prep (ng)')
    parser.add_argument('--dna-conc', type=float, required=True, help='DNA concentration (ng/µL)')
    parser.add_argument('--detect-thresh', type=float, required=True, help='Detection threshold')
    parser.add_argument('--quad-regs', nargs=4, required=True, help='Paths to read depth variability models')
    parser.add_argument('--cutoff-functions', nargs=4, required=True, help='Paths to RMSE threshold functions')

    args = parser.parse_args()

    

    results = quant_unknown(
        args.sample_name, args.database_lengths, args.quant_regression, args.mapping_results, args.sliding_window,
        args.bin_assignment, args.target_type, args.dna_input, args.dna_conc,
        args.recovery, args.conc_factor, args.e_detect_model, args.quad_regs, args.cutoff_functions,
        output_files
    )

    print("Quantification complete!")
    print(f"Quantified {len(results)} targets")