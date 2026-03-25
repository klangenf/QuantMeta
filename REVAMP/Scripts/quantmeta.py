#!/usr/bin/env python3
"""QuantMeta 

Purpose: generating standard curves relating relative 
abundance to absolute abundance for each sample.
"""

import argparse
import json
import os
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


def detection_threshold(sample_name, mapping_results_path, target_length_df, detect_thresh_path):
    detect_model = load_detection_model(detect_thresh_path)

    mapping_results = pd.read_csv(mapping_results_path, sep='\t', header=0)
    mapping_results = mapping_results.astype({
        'ID': str,
        'read_depth': float
    })

    targets = pd.DataFrame({'ID': mapping_results['ID'].unique()})
    targets = targets.merge(target_length_df[['ID', 'length']], on='ID', how='left')

    targets['B_G'] = 0.0
    targets['gene_copies'] = 0.0
    targets['I_G'] = 0.0
    targets['E_rel'] = 0.0

    for i, row in targets.iterrows():
        target_id = row['ID']
        target = mapping_results[mapping_results['ID'] == target_id]

        B_G = target['read_depth'].sum()
        gene_copies = target['read_depth'].mean() if len(target) > 0 else 0.0

        if B_G > 0 and len(target) > 0:
            ratios = target['read_depth'] / B_G
            ratio_log = np.log(ratios.replace({0: np.nan}))
            I_G_x = ratios * ratio_log
            I_G = -np.nansum(I_G_x.replace([np.inf, -np.inf], np.nan))
        else:
            I_G = 0.0

        length_val = float(row['length']) if not pd.isna(row['length']) else np.nan
        E_rel = I_G / np.log(length_val) if length_val > 1 else 0.0

        targets.at[i, 'B_G'] = B_G
        targets.at[i, 'gene_copies'] = gene_copies
        targets.at[i, 'I_G'] = I_G
        targets.at[i, 'E_rel'] = E_rel

    # compute E_detect for each length from model (either constant or length-based table)
    targets['E_detect'] = predict_E_detect_by_length(targets['length'].astype(float).values, detect_model)

    targets['detection_status'] = np.where(targets['E_rel'] >= targets['E_detect'], 'detected', 'not_detected')

    out_table = targets[['ID', 'E_rel', 'gene_copies', 'E_detect', 'detection_status']].copy()

    output_path = Path('Mapping') / sample_name / 'standards_mapping_analysis.txt'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_table.to_csv(output_path, sep='\t', index=False)

    return out_table


def expected_conc(MIX, spike, DNA_conc):
    result = MIX[['ID']].copy()
    result['known_conc'] = MIX['MIX'] * spike * DNA_conc / MIX['Mass']
    return result.rename(columns={'ID': 'ID', 'known_conc': 'known_conc'})


def prediction(mapping, DNA_input, DNA_conc):
    result = pd.DataFrame({'ID': mapping['ID']})
    result['predicted_conc'] = mapping['gene_copies'] * DNA_input / DNA_conc
    return result


def visualize(sample_name, results):
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
             label=f'y = {10**intercept:.2f} x^{slope:.2f}  (R^2={r_value**2:.3f})')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Predicted concentration (gc/µL)')
    plt.ylabel('Known concentration (gc/µL)')
    plt.title(f'{sample_name}: Known vs predicted concentration')
    plt.legend()
    plt.grid(True, which='both', ls=':', alpha=0.3)

    output_path = Path('Results') / sample_name / 'standards_rel_to_abs.png'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=320, bbox_inches='tight')
    plt.close()

    return {'slope': slope, 'intercept': intercept, 'r2': r_value**2}


def quantmeta(sample_name,
              STD_MIXES,
              mapping_results_path,
              DNA_input,
              DNA_conc,
              MIX,
              sequins_spike,
              PCR_quant_file,
              ssDNA_stds,
              ssDNA_spike,
              detect_thresh):

    STD_MIXES = pd.read_csv(STD_MIXES, sep='\t', header=0)
    keep = ['ID', 'length', 'Mass', MIX]
    STD_MIX = STD_MIXES[keep].copy()
    STD_MIX.columns = ['ID', 'length', 'Mass', 'MIX']
    STD_MIX = STD_MIX.dropna().reset_index(drop=True)

    if ssDNA_stds.upper() == 'YES':
        ssDNA_STD_MIX = STD_MIXES[['ID', 'length', 'Mass', 'ssDNA']].copy()
        ssDNA_STD_MIX.columns = ['ID', 'length', 'Mass', 'MIX']
        ssDNA_STD_MIX = ssDNA_STD_MIX.dropna().reset_index(drop=True)
        STD_MIX = pd.concat([STD_MIX, ssDNA_STD_MIX], ignore_index=True)
    else:
        ssDNA_STD_MIX = pd.DataFrame(columns=['ID', 'length', 'Mass', 'MIX'])

    std_length = STD_MIX[['ID', 'length']].copy()

    mapping = detection_threshold(sample_name, mapping_results_path, std_length, detect_thresh)
    mapping = mapping.astype({'E_rel': float, 'gene_copies': float, 'E_detect': float})
    mapping = mapping[mapping['E_rel'] >= mapping['E_detect']].copy()

    dsDNA_mapping = mapping[~mapping['ID'].isin(ssDNA_STD_MIX['ID'])].copy()

    if os.path.getsize(PCR_quant_file) < 500:
        results = expected_conc(STD_MIX[~STD_MIX['ID'].isin(ssDNA_STD_MIX['ID'])].copy(), sequins_spike, DNA_conc)
    else:
        results = pd.read_csv(PCR_quant_file, sep='\t', header=0)
        if sample_name in results.columns:
            results = results[['ID', sample_name]].copy()
            results.columns = ['ID', 'known_conc']
        elif 'known_conc' not in results.columns:
            raise ValueError('PCR quant file must contain known_conc or sample column with known concentration')
        results = results[~results['ID'].isin(ssDNA_STD_MIX['ID'])].copy()

    predict_df = prediction(dsDNA_mapping, DNA_input, DNA_conc)
    results = results.merge(predict_df, on='ID', how='inner')

    if ssDNA_stds.upper() == 'YES':
        ssDNA_mapping = mapping[mapping['ID'].isin(ssDNA_STD_MIX['ID'])].copy()

        if os.path.getsize(PCR_quant_file) < 500:
            ssDNA_results = expected_conc(ssDNA_STD_MIX, ssDNA_spike, DNA_conc)
        else:
            ssDNA_results = pd.read_csv(PCR_quant_file, sep='\t', header=0)
            if sample_name in ssDNA_results.columns:
                ssDNA_results = ssDNA_results[['ID', sample_name]].copy()
                ssDNA_results.columns = ['ID', 'known_conc']
            ssDNA_results = ssDNA_results[ssDNA_results['ID'].isin(ssDNA_STD_MIX['ID'])].copy()

        ssDNA_predict = prediction(ssDNA_mapping, DNA_input, DNA_conc)
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
        'r_sq': float(r_value**2)
    }

    reg_name = Path('Regressions/quantification') / f'{sample_name}_rel_to_abs.pkl'
    reg_name.parent.mkdir(parents=True, exist_ok=True)
    with open(reg_name, 'wb') as f:
        pickle.dump(model, f)

    visualize(sample_name, results)

    return results


def main():
    parser = argparse.ArgumentParser(description='Convert quantmeta.R workflow into python.')
    parser.add_argument('--sample-info', required=False, help='Tab-separated sample info table')
    parser.add_argument('--sample-name', required=False, help='Sample name to process')
    parser.add_argument('--mapping-results', required=False, help='Mapping results table path')
    parser.add_argument('--std-mixes', required=False, help='STD_MIXES file path')
    parser.add_argument('--dna-input', type=float, required=False, help='DNA input mass (ng)')
    parser.add_argument('--dna-conc', type=float, required=False, help='DNA concentration (ng/µL)')
    parser.add_argument('--mix', required=False, help='MIX name column in STD_MIXES')
    parser.add_argument('--dsdna-spike', type=float, required=False, help='dsDNA spike ratio')
    parser.add_argument('--pcr-quant-file', required=False, help='PCR quantification file path')
    parser.add_argument('--ssdna-stds', default='NO', choices=['YES', 'NO'], help='Use ssDNA standards')
    parser.add_argument('--ssdna-spike', type=float, default=0.0, help='ssDNA spike ratio')
    parser.add_argument('--detect-thresh', required=False, help='Path to detection threshold model (json or table)')
    parser.add_argument('--output', required=False, help='Output table file path')

    args = parser.parse_args()

    if 'snakemake' in globals():
        snakemake = globals()['snakemake']

        sample_info = pd.read_csv(snakemake.input[1], sep='\t', header=0)
        sample_info = sample_info[sample_info['Sample'] == snakemake.params[0]].iloc[0]

        results = quantmeta(
            sample_name=sample_info['Sample'],
            STD_MIXES=snakemake.input[2],
            mapping_results_path=snakemake.input[0],
            DNA_input=float(sample_info['lib_mass']),
            DNA_conc=float(sample_info['DNA_conc']),
            MIX=sample_info['mix'],
            sequins_spike=float(sample_info['dsDNA_spike']),
            PCR_quant_file=snakemake.input[3],
            ssDNA_stds=sample_info['ssDNA'],
            ssDNA_spike=float(sample_info['ssDNA_spike']),
            detect_thresh=snakemake.input[4]
        )

        results.to_csv(snakemake.output[0], sep='\t', index=False, quoting=1)

    else:
        if not args.sample_name or not args.mapping_results or not args.std_mixes or not args.pcr_quant_file or not args.detect_thresh:
            parser.error('When running outside Snakemake, you must provide --sample-name, --mapping-results, --std-mixes, --pcr-quant-file, --detect-thresh')

        results = quantmeta(
            sample_name=args.sample_name,
            STD_MIXES=args.std_mixes,
            mapping_results_path=args.mapping_results,
            DNA_input=args.dna_input,
            DNA_conc=args.dna_conc,
            MIX=args.mix,
            sequins_spike=args.dsdna_spike,
            PCR_quant_file=args.pcr_quant_file,
            ssDNA_stds=args.ssdna_stds,
            ssDNA_spike=args.ssdna_spike,
            detect_thresh=args.detect_thresh
        )

        if args.output:
            Path(args.output).parent.mkdir(exist_ok=True, parents=True)
            results.to_csv(args.output, sep='\t', index=False, quoting=1)
        else:
            print(results.head())


if __name__ == '__main__':
    main()
