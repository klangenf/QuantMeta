#!/usr/bin/env python3
"""QuantMeta 

Purpose: generating standard curves relating relative 
abundance to absolute abundance for each sample.
"""

import argparse
import json
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress


def load_detection_model(detect_thresh_path):
    detect_thresh_path = Path(detect_thresh_path)
    if detect_thresh_path.suffix.lower() in ['.json']:
        with open(detect_thresh_path, 'r') as f:
            model = json.load(f)
        # model expected to have intercept, coef_E_rel, and optionally y_thresholds.P0.95
        return model

    # fallback: plain table with length and E_detect
    if detect_thresh_path.exists():
        df = pd.read_csv(detect_thresh_path, sep='\t', comment='#')
        if {'length', 'E_detect'}.issubset(df.columns):
            df = df.sort_values('length').reset_index(drop=True)
            return df

    raise ValueError(f"Cannot parse detect_thresh model from {detect_thresh_path}")


def predict_E_detect_by_length(lengths, detect_model):
    if isinstance(detect_model, dict):
        if detect_model.get('type') == 'log10_length':
            intercept = float(detect_model['intercept'])
            slope = float(detect_model['slope'])
            with np.errstate(divide='ignore', invalid='ignore'):
                l = np.array(lengths, dtype=float)
                result = intercept + slope * np.log10(l)
            return result

    raise ValueError('Unsupported detect_model type')


def detection_threshold(sample_name, mapping_results_path, target_length_df, detect_thresh, out_dir):
    mapping_results = pd.read_csv(mapping_results_path, sep='\t', header=0)
    mapping_results = mapping_results.astype({
        'ID': str,
        'read_depth': float
    })

    # mapping_results: DataFrame with ID, read_depth
    # target_lengths: DataFrame with ID, length

    # Start with unique IDs and merge lengths
    targets = pd.DataFrame({'ID': mapping_results['ID'].unique()})
    targets = targets.merge(target_length_df[['ID', 'length']], on='ID', how='left')
    
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
        E_detect=predict_E_detect_by_length(targets['length'].astype(float).values, detect_thresh),
        detection_status=lambda df: np.where(df['E_rel'] >= df['E_detect'], 'detected', 'not_detected')
    )

    out_table = targets[['ID', 'gene_copies', 'E_rel', 'E_detect', 'detection_status']].copy()

    output_path = Path(out_dir / 'Mapping')/ sample_name /'standards_mapping_analysis.txt'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_table.to_csv(output_path, sep='\t', index=False)

    return out_table


def expected_conc(MIX, spike, DNA_conc):
    result = MIX[['ID']].copy()
    result['known_conc'] = MIX['Rel_Abund'] * spike * DNA_conc / MIX['Mass']
    return result.rename(columns={'ID': 'ID', 'known_conc': 'known_conc'})


def prediction(mapping, DNA_input, DNA_conc):
    result = pd.DataFrame({'ID': mapping['ID']})
    result['predicted_conc'] = mapping['gene_copies'] * DNA_input / DNA_conc
    return result


def visualize(sample_name, results, out_dir):
    results = results.copy().dropna(subset=['predicted_conc', 'known_conc'])
    results = results[(results['predicted_conc'] > 0) & (results['known_conc'] > 0)]

    if results.empty:
        print(f'No data to visualize for {sample_name}')
        return

    x = np.log10(results['predicted_conc'].values)
    y = np.log10(results['known_conc'].values)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    plt.figure(figsize=(8, 6))
    plt.scatter(results['predicted_conc'], results['known_conc'], edgecolor='k', facecolor='none', alpha=0.7)

    x_line = np.logspace(np.log10(results['predicted_conc'].min()), np.log10(results['predicted_conc'].max()), 100)
    y_line = 10 ** (intercept + slope * np.log10(x_line))

    plt.plot(x_line, y_line, color='red', linestyle='--', linewidth=1.5,
             label=f'log10(y) = {intercept:.2f} + {slope:.2f}log10(x)  (R^2={r_value**2:.3f})')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Predicted concentration (gc/µL)')
    plt.ylabel('Known concentration (gc/µL)')
    plt.title(f'{sample_name}: Absolute vs. Rela')
    plt.legend()
    plt.grid(True, which='both', ls=':', alpha=0.3)

    output_path = Path(out_dir / 'Regressions/quantification') / f'{sample_name}_standards_rel_to_abs.png'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=320, bbox_inches='tight')
    plt.close()

    return {'slope': slope, 'intercept': intercept, 'r2': r_value**2}


def quantmeta(sample_name,
              dsDNA_stds,
              ssDNA_stds,
              mapping_results_path,
              DNA_input,
              DNA_conc,
              sequins_spike,
              ssDNA_spike,
              detect_thresh,
              out_dir):

    dsDNA_STD = dsDNA_stds
    dsDNA_STD['Type'] = 'dsDNA'
    STD_MIX = dsDNA_STD.dropna().reset_index(drop=True)

    if ssDNA_spike > 0:
        ssDNA_STD = ssDNA_stds
        ssDNA_STD = ssDNA_STD.dropna().reset_index(drop=True)
        ssDNA_STD['Type'] = 'ssDNA'
        STD_MIX = pd.concat([STD_MIX, ssDNA_STD], ignore_index=True)

    std_length = STD_MIX[['ID', 'length']].copy()

    mapping = detection_threshold(sample_name, mapping_results_path, std_length, detect_thresh, out_dir)
    mapping = mapping.astype({'E_rel': float, 'gene_copies': float, 'E_detect': float})
    mapping = mapping[mapping['E_rel'] >= mapping['E_detect']].copy()

    dsDNA_mapping = mapping[mapping['ID'].isin(dsDNA_STD['ID'])].copy()

    results = expected_conc(dsDNA_STD.copy(), sequins_spike, DNA_conc)
    predict_df = prediction(dsDNA_mapping, DNA_input, DNA_conc)
    results = results.merge(predict_df, on='ID', how='inner')

    if ssDNA_spike > 0:
        ssDNA_mapping = mapping[mapping['ID'].isin(ssDNA_STD['ID'])].copy()

        ssDNA_results = expected_conc(ssDNA_STD.copy(), ssDNA_spike, DNA_conc)
        ssDNA_predict = prediction(ssDNA_mapping, DNA_input, DNA_conc)
        ssDNA_predict['predicted_conc'] = ssDNA_predict['predicted_conc']*2 # adjust for ssDNA vs dsDNA
        ssDNA_results = ssDNA_results.merge(ssDNA_predict, on='ID', how='inner')
        results = pd.concat([results, ssDNA_results], ignore_index=True)

    # linear regression on log10 values
    valid = results[(results['predicted_conc'] > 0) & (results['known_conc'] > 0)].copy()
    if valid.empty:
        raise ValueError('No valid rows for quantification regression')

    log_x = np.log10(valid['predicted_conc'].values)
    log_y = np.log10(valid['known_conc'].values)
    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)

    model = {
        'slope': float(slope),
        'intercept': float(intercept),
        'r_sq': float(r_value**2),
        'n': int(len(results))
    }

    reg_name = Path(out_dir / 'Regressions/quantification') / f'{sample_name}_rel_to_abs.pkl'
    reg_name.parent.mkdir(parents=True, exist_ok=True)
    with open(reg_name, 'wb') as f:
        pickle.dump(model, f)

    visualize(sample_name, results, out_dir)

    return results


def main():
    parser = argparse.ArgumentParser(description='QuantMeta: Standard curve generation for relating relative abundance to absolute abundance.')
    parser.add_argument('--sample-info', required=True, default='Config/spike_in_info.txt', help='Tab-separated sample info table')
    parser.add_argument('--dsDNA-std-mixes', required=True, default='Spike-ins/sequins_Mix_A.txt', help='Table of dsDNA standards (ID, Mass, Rel_Abund, length)')
    parser.add_argument('--ssDNA-std-mixes', required=False, default=None, help='Optional table of ssDNA standards (ID, Mass, Rel_Abund, length)')
    parser.add_argument('--detect-threshold', required=True, default='Regressions/detect/Langenfeld_2025_E_detect.json', help='Length-dependent entropy detection threshold')
    parser.add_argument('--output-dir', required=True, default='QuantMeta/', help='Output directory for project')

    args = parser.parse_args()

    sample_info = pd.read_csv(args.sample_info, sep='\t', header=0)
    dsDNA_stds = pd.read_csv(args.dsDNA_std_mixes, sep='\t', header=0)
    ssDNA_stds = pd.read_csv(args.ssDNA_std_mixes, sep='\t', header=0) if args.ssDNA_std_mixes else None
    detect_model = load_detection_model(args.detect_threshold)
    out_dir = args.output_dir

    for _, row in sample_info.iterrows():
        sample_name = row['Sample']
        mapping_results_path = f'{out_dir}/Mapping/{sample_name}/standards_mapping.txt'
        DNA_input = row['Library_Mass']
        DNA_conc = row['DNA_Extract_Conc']
        sequins_spike = row['Spike_Frac']
        ssDNA_spike = row['ssDNA'] if 'ssDNA' in row else 0

        results = quantmeta(
            sample_name=sample_name,
            dsDNA_stds=dsDNA_stds,
            ssDNA_stds=ssDNA_stds,
            mapping_results_path=mapping_results_path,
            DNA_input=DNA_input,
            DNA_conc=DNA_conc,
            sequins_spike=sequins_spike,
            ssDNA_spike=ssDNA_spike,
            detect_thresh=detect_model,
            out_dir=out_dir
        )

        output = Path(out_dir / 'Regressions/quantification') / f'{sample_name}_quantmeta_results.txt'
        output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(output, sep='\t', index=False)
        
        print('✅ Quantification regression custom build complete for sample:', sample_name)

if __name__ == '__main__':
    main()
