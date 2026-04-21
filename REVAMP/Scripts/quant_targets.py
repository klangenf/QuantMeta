#!/usr/bin/env python
"""
Quantification of Unknown Targets

This script quantifies unknown targets in metagenomic data by:
1. Assessing detection thresholds for standards
2. Correcting for non-specific mapping artifacts
3. Converting relative abundances to absolute concentrations
"""

import pandas as pd
import numpy as np
import pickle
import json
import collections
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
            np.log(X_new['gene_copies'])
        ])
        preds = X_fit @ coeffs
        preds = np.vstack(preds)[:,0]

    elif bin_name == 'bin3':
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            np.log(X_new['gene_copies'])
        ])
        preds = X_fit @ coeffs
        preds = np.vstack(preds)[:,0]

    elif bin_name == 'bin2':
        log_depth = np.log(X_new['gene_copies'])
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            log_depth**2, log_depth
        ])
        preds = X_fit @ coeffs
        preds = np.vstack(preds)[:,0]

    elif bin_name == 'bin1':
        log_depth = np.log(X_new['gene_copies'])
        X_fit = np.column_stack([
            X_new['avg_GC']**2, X_new['avg_GC'], np.ones(len(X_new)),
            log_depth**2, log_depth,
            X_new['E_rel']**2, X_new['E_rel']
        ])
        preds = X_fit @ coeffs
        preds = np.vstack(preds)[:,0]

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


def detection_threshold(mapping_results, lengths, detect_thresh):
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
        E_detect=predict_E_detect(detect_thresh, targets['length'].astype(float).values),
        detection_status=lambda df: np.where(df['E_rel'] >= df['E_detect'], 'detected', 'not_detected')
    )

    out_table = targets[['ID', 'length', 'gene_copies', 'E_rel', 'E_detect', 'detection_status']].copy()

    return out_table


def quant_correction(sample_name, input_info, sliding_window, quad_reg1, quad_reg2, 
                     quad_reg3, quad_reg4, cutoff_function1, cutoff_function2, 
                     cutoff_function3, cutoff_function4):
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

    Returns:
    --------
    pd.DataFrame : Corrected target information
    """
    # Save initial RMSE
    input_info['initial_gene_copies'] = input_info['gene_copies']
    input_info['initial_RMSE'] = input_info['RMSE']
    input_info['frac_corrected'] = 0.0
    input_info['cycle'] = 0

    max_cycles = 20

    sliding_window['initial_avg_depth'] = sliding_window['avg_depth']

    # Select appropriate models based on depth
    if input_info['gene_copies'].iloc[0] < 10:
        read_depth_var_model = quad_reg1
    elif input_info['gene_copies'].iloc[0] < 100:
        read_depth_var_model = quad_reg2
    elif input_info['gene_copies'].iloc[0] < 1000:
        read_depth_var_model = quad_reg3
    else:
        read_depth_var_model = quad_reg4


    # Iterative correction loop
    while (input_info['RMSE'].iloc[0] > input_info['RMSE_limit'].iloc[0] and
           input_info['cycle'].iloc[0] <= max_cycles and
           input_info['frac_corrected'].iloc[0] <= 0.2 and
           input_info['gene_copies'].iloc[0] > 0):
           
           pred = predict_read_var(read_depth_var_model, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
           
           sd_mult = np.std(sliding_window['avg_depth'])
           upper = pred + 1.5 * sd_mult
           lower = pred - 1.5 * sd_mult
           
           # Ensure non-negative bounds
           lower = np.maximum(lower, 0)
           upper = np.maximum(upper, 0)
           
           # Identify mapping error regions
           sliding_window['error_region'] = False
           sliding_window.loc[((sliding_window['avg_depth'] != 0.0) & ((sliding_window['avg_depth'] > upper) | (sliding_window['avg_depth'] < lower))), 'error_region'] = True
           
           # Recalculate average depth excluding error regions
           specific_regions = sliding_window.loc[sliding_window['error_region'] == False]
           if len(specific_regions) > 0: 
               sliding_window['gene_copies'] = specific_regions['avg_depth'].mean()
           else:
               input_info['frac_corrected'].iloc[0] = 1
               break
           
           # Select appropriate models based on depth
           if sliding_window['gene_copies'].iloc[0] < 10:
               read_depth_var_model = quad_reg1
           elif sliding_window['gene_copies'].iloc[0] < 100:
               read_depth_var_model = quad_reg2
           elif sliding_window['gene_copies'].iloc[0] < 1000:
               read_depth_var_model = quad_reg3
           else:
               read_depth_var_model = quad_reg4

           pred = predict_read_var(read_depth_var_model, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
           
           specific_depths = sliding_window.loc[sliding_window['error_region'] == False, 'avg_depth']
           sd_mult = np.std(specific_depths) if len(specific_depths) > 0 else 0
           upper = pred + 1.5 * sd_mult
           lower = pred - 1.5 * sd_mult
           
           # Ensure non-negative bounds
           lower = np.maximum(lower, 0)
           upper = np.maximum(upper, 0)

           temp = pd.concat([sliding_window, pd.DataFrame({'lower': lower, 'upper': upper})], axis=1)

           # Apply corrections
           sliding_window.loc[sliding_window['error_region'] == True, 'avg_depth'] = np.clip(sliding_window.loc[sliding_window['error_region'] == True, 'avg_depth'], temp.loc[temp['error_region'] == True, 'lower'], temp.loc[temp['error_region'] == True, 'upper'])
           
           # Update total average depth
           sliding_window['gene_copies'] = sliding_window['avg_depth'].mean()
           input_info['gene_copies'].iloc[0] = sliding_window['gene_copies'].iloc[0]
           
           # Select appropriate models based on depth
           if input_info['gene_copies'].iloc[0] < 10:
               read_depth_var_model = quad_reg1
               rmse_threshold_model = cutoff_function1
           elif input_info['gene_copies'].iloc[0] < 100:
               read_depth_var_model = quad_reg2
               rmse_threshold_model = cutoff_function2
           elif input_info['gene_copies'].iloc[0] < 1000:
               read_depth_var_model = quad_reg3
               rmse_threshold_model = cutoff_function3
           else:
               read_depth_var_model = quad_reg4
               rmse_threshold_model = cutoff_function4

           # Recalculate RMSE and limits
           pred = predict_read_var(read_depth_var_model, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
           
           rmse = np.sqrt(np.sum((pred - sliding_window['avg_depth'])**2) / len(sliding_window))
           RMSE_limit = predict_rmse_limit(rmse_threshold_model, input_info['gene_copies'].iloc[0])
           
           input_info['RMSE'].iloc[0] = rmse
           input_info['RMSE_limit'].iloc[0] = RMSE_limit
           input_info['cycle'].iloc[0] += 1
           
           # Calculate fraction corrected
           input_info['frac_corrected'].iloc[0] = np.sum(sliding_window['avg_depth'] != sliding_window['initial_avg_depth']) / len(sliding_window)
           
    # Save corrected mapping sliding windows for this target
    output_dir = f'Results/{sample_name}/Ind_Correction_Results/{input_info["ID"].iloc[0]}_corrected_mapping.tsv'
    output_dir = Path(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    sliding_window.to_csv(output_dir, sep='\t', index=False)

    return input_info

def compute_stats(seq, window_size=49, contig=None):
    stats = []
    seq_array = seq['nucleic_acid']
    depth_array = seq['read_depth']
    seq_len = seq.shape[0]

    # check if seq is greater than sliding window
    if seq_len < window_size:
        print("does not meet length req")
        #return -1

    # get remaining windows and stats
    for i in range(seq_len - window_size):
        end = i + window_size
        seq_window = seq_array[i:end]
        depth_window = depth_array[i:end]
        avg_depth = compute_avg_depth(depth_window)
        gc = compute_gc(seq_window)
        stats.append({"ID":contig,"avg_GC":gc,"avg_depth":avg_depth})
    
    return(stats)
    
def compute_gc(seq_array):
    return((seq_array.str.count('C').sum() + seq_array.str.count('G').sum()) / len(seq_array))

def compute_avg_depth(depth_array):
    return(sum(depth_array)/len(depth_array))

def quant_correct_analysis(sample_name, descript, mapping_results, results, quad_reg1_path, 
                           quad_reg2_path, quad_reg3_path, quad_reg4_path, cutoff_function1_path, 
                           cutoff_function2_path, cutoff_function3_path, cutoff_function4_path, window_size=49):
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
    # Load models
    quad_reg1 = load_pickle_model(quad_reg1_path)
    quad_reg2 = load_pickle_model(quad_reg2_path)
    quad_reg3 = load_pickle_model(quad_reg3_path)
    quad_reg4 = load_pickle_model(quad_reg4_path)

    cutoff_function1 = load_pickle_model(cutoff_function1_path)
    cutoff_function2 = load_pickle_model(cutoff_function2_path)
    cutoff_function3 = load_pickle_model(cutoff_function3_path)
    cutoff_function4 = load_pickle_model(cutoff_function4_path)

    results['RMSE'] = 0.0
    results['RMSE_limit'] = 0.0
    results['status'] = 'no_correction'
    corrected = pd.DataFrame()

    for ID in results['ID'].unique():
        map_ID = mapping_results[mapping_results['ID'] == ID]
        if len(map_ID) > window_size:
            sliding_window = pd.DataFrame(compute_stats(map_ID, window_size=window_size, contig=ID), columns=['ID', 'avg_GC', 'avg_depth'])
            sliding_window = pd.merge(sliding_window, results, on='ID', how='left')
        else:
            print(f"ID {ID} does not meet length requirement for sliding window analysis.")
            continue
        
        # Calculate RMSE and RMSE limit for the target
        depth = results[results['ID'] == ID]['gene_copies'].iloc[0]

        if depth < 10 and quad_reg1 is not None:
            pred = predict_read_var(quad_reg1, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
            RMSE_limit = predict_rmse_limit(cutoff_function1, depth)
        elif 10 <= depth < 100 and quad_reg2 is not None:
            pred = predict_read_var(quad_reg2, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
            RMSE_limit = predict_rmse_limit(cutoff_function2, depth)
        elif 100 <= depth < 1000 and quad_reg3 is not None:
            pred = predict_read_var(quad_reg3, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
            RMSE_limit = predict_rmse_limit(cutoff_function3, depth)
        elif depth >= 1000 and quad_reg4 is not None:
            pred = predict_read_var(quad_reg4, sliding_window[['avg_GC', 'gene_copies', 'E_rel']])
            RMSE_limit = predict_rmse_limit(cutoff_function4, depth)
        else:
            continue

        rmse = np.sqrt(np.sum((pred - sliding_window['avg_depth'])**2) / sliding_window.shape[0])

        results.loc[(results['ID'] == ID), 'RMSE'] = rmse
        results.loc[(results['ID'] == ID), 'RMSE_limit'] = RMSE_limit

        # Apply correction if RMSE exceeds limit (based on average read depth)
        if results[results['ID'] == ID]['RMSE'].iloc[0] > results[results['ID'] == ID]['RMSE_limit'].iloc[0]:
            results.loc[results['ID'] == ID, 'status'] = 'needs_correction'
            sliding_window['status'] = 'needs_correction'
            correct = quant_correction(sample_name, results[results['ID'] == ID], sliding_window, 
                                          quad_reg1, quad_reg2, quad_reg3, quad_reg4, cutoff_function1, 
                                          cutoff_function2, cutoff_function3, cutoff_function4)
            corrected = pd.concat([corrected, correct], ignore_index=True)

    if len(corrected) > 0:
        results = results[~results['ID'].isin(corrected['ID'])]

    results = results[['ID', 'RMSE', 'RMSE_limit', 'gene_copies', 'status']]
    results['initial_gene_copies'] = results['gene_copies']
    results['frac_corrected'] = 0.0
    results['cycle'] = 0
    
    if len(corrected) > 0:
        corrected = corrected[['ID', 'RMSE', 'RMSE_limit', 'gene_copies', 'status', 'initial_gene_copies', 'frac_corrected', 'cycle']]
        results = pd.concat([results, corrected], ignore_index=True)

    # Add reliability assessment
    results['reliability'] = 'High RMSE'
    results.loc[(results['RMSE'] <= results['RMSE_limit']), 'reliability'] = 'Low RMSE'

    results['correction_results'] = 'Accurate'
    mask = "status == 'needs_correction' and (frac_corrected > 0.2 or cycle > 20 or RMSE > RMSE_limit or gene_copies == 0)"
    results.loc[results.query(mask).index, 'correction_results'] = 'Review'
    
    # Save results
    output_dir = f'Results/{sample_name}/Ind_Correction_Results/{descript}_corrected_results.tsv'
    output_dir = Path(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir, sep='\t', index=False)

    return results


def prediction(mapping, DNA_input, DNA_conc):
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

    Returns:
    --------
    pd.DataFrame : Predicted concentrations
    """
    def pred_conc_calc(copies, mass, conc):
        """Convert gene copies to concentration in DNA extract."""
        return copies * conc / mass

    result = pd.DataFrame({
        'ID': mapping['ID'],
        'predicted_conc': pred_conc_calc(mapping['gene_copies'], DNA_input, DNA_conc)
    })

    return result

def quant_unknown(sample_name, target_name, database_lengths, mapping_results_path, E_detect_model_path,
                  quad_reg1_path, quad_reg2_path, quad_reg3_path, quad_reg4_path,
                  cutoff_function1_path, cutoff_function2_path, cutoff_function3_path, cutoff_function4_path,
                  quant_regression_path, DNA_input, DNA_conc, window_size=49):
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
    mapping_results_path : str
        Path to mapping results file
    E_detect_model_path : str
        Path to detection model
    window_size : int
        Size of sliding window for correction analysis, default: 49
    quad_reg1..4_path : str
        Path to each read depth variability model
    cutoff_function1..4_path : str
        Path to each RMSE threshold function
    quant_regression_path : str
        Path to quantification regression model
    DNA_input : float
        DNA mass used in library prep (ng)
    DNA_conc : float
        DNA concentration (ng/µL)
    

    Returns:
    --------
    pd.DataFrame : Final quantification results
    """
    output_dir = f'Results/{sample_name}/{target_name}_concentrations.tsv'
    output_dir = Path(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    # Load sequence lengths
    lengths = pd.read_csv(database_lengths, sep='\t', header=None, names=['ID', 'length'])

    # Load mapping data
    mapping_results = pd.read_csv(mapping_results_path, sep='\t', header=0)

    # Assess detection thresholds
    detect_thresh = load_detection_model(E_detect_model_path)

    mapping = detection_threshold(mapping_results, lengths, detect_thresh)

    output_path = Path('Mapping')/ sample_name / f'{target_name}_mapping_analysis.txt'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(output_path, sep='\t', index=False)

    # Filter targets above detection threshold
    mapping = mapping[mapping['E_rel'] >= mapping['E_detect']]

    # Detect and correct read mapping errors
    if len(mapping) > 0:
        mapping = quant_correct_analysis(sample_name, target_name, mapping_results, mapping, quad_reg1_path, 
                                         quad_reg2_path, quad_reg3_path, quad_reg4_path, cutoff_function1_path, 
                                         cutoff_function2_path, cutoff_function3_path, cutoff_function4_path, window_size)

        # Remove targets that are not quantifiable
        mapping = mapping[mapping['correction_results'].isin(['Accurate'])]
    else:
        print("No targets above detection and quantifiable!")
        pd.DataFrame({"message": ["No targets above detection and quantifiable!"]}).to_csv(output_dir, index=False, header=False)

    # Convert relative abundances from fractions to copies/µL DNA extract
    results = prediction(mapping, DNA_input, DNA_conc)

    # Convert to absolute abundance using regression model
    quant_regression = load_pickle_model(quant_regression_path)
    pred = quant_regression['intercept'] + quant_regression['slope'] * np.log10(results['predicted_conc'])
    pred_se = np.sqrt((4*quant_regression['r_sq'] * (1 - quant_regression['r_sq'])**2 * (quant_regression['n'] - 2)**2)/((quant_regression['n']**2 - 1) * (quant_regression['n'] + 3)))
    pred_se_up = quant_regression['r_sq'] + pred_se
    pred_se_low = quant_regression['r_sq'] - pred_se

    results_final = pd.DataFrame({
        'ID': results['ID'],
        'concentration (copies/µL DNA extract)': 10**pred,
        'upper 95% CI (copies/µL DNA extract)': 10**pred_se_up,
        'lower 95% CI (copies/µL DNA extract)': 10**pred_se_low
    })

    # Save final results
    results_final.to_csv(output_dir, sep='\t', index=False)

    return results_final


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Quantify unknown targets in metagenomic data')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--target_name', required=True, help='Target name')
    parser.add_argument('--sample_info', required=True, help='Path to sample information file')
    parser.add_argument('--database_lengths', required=True, help='Path to sequence lengths file')
    parser.add_argument('--mapping_results', required=True, help='Path to mapping results file')
    parser.add_argument('--window_size', type=int, default=49, help='Size of sliding window for correction analysis')
    parser.add_argument('--quant_regression', required=True, help='Path to quantification regression model')
    parser.add_argument('--detect_thresh', required=True, help='Path to detection threshold')
    parser.add_argument('--quad_reg1', required=True, help='Path to read depth variability model for 0-10 reads/bp')
    parser.add_argument('--quad_reg2', required=True, help='Path to read depth variability model for 10-100 reads/bp')
    parser.add_argument('--quad_reg3', required=True, help='Path to read depth variability model for 100-1000 reads/bp')
    parser.add_argument('--quad_reg4', required=True, help='Path to read depth variability model for 1000+ reads/bp')
    parser.add_argument('--cutoff_function1', required=True, help='Paths to RMSE threshold function for 0-10 reads/bp')
    parser.add_argument('--cutoff_function2', required=True, help='Paths to RMSE threshold function for 10-100 reads/bp')
    parser.add_argument('--cutoff_function3', required=True, help='Paths to RMSE threshold function for 100-1000 reads/bp')
    parser.add_argument('--cutoff_function4', required=True, help='Paths to RMSE threshold function for 1000+ reads/bp')

    args = parser.parse_args()

    sample_name = args.sample_name
    target_name = args.target_name
    sample_info = pd.read_csv(args.sample_info, sep='\t', header=0)
    DNA_input = sample_info.loc[sample_info['Sample'] == sample_name, 'Library_Mass']
    DNA_conc = sample_info.loc[sample_info['Sample'] == sample_name, 'DNA_Extract_Conc']

    results = quant_unknown(sample_name, target_name, args.database_lengths, args.mapping_results, 
                            args.detect_thresh, args.quad_reg1, args.quad_reg2, args.quad_reg3, 
                            args.quad_reg4, args.cutoff_function1, args.cutoff_function2, 
                            args.cutoff_function3, args.cutoff_function4, args.quant_regression, 
                            DNA_input, DNA_conc, args.window_size)

    print("Quantification complete!")
    print(f"Quantified {len(results)} targets.")