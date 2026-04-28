#!/usr/bin/env python3
"""Wrapper script to run quant_targets.sh from the package."""

import argparse
import subprocess
import sys
from pathlib import Path
import importlib.resources as resources


def get_script_path(script_name: str) -> Path:
    """Get the path to a shell script in the package."""
    try:
        with resources.as_file(resources.files('quantmeta.quant_targets') / script_name) as p:
            return Path(p)
    except Exception:
        return Path(__file__).parent / script_name


def main():
    parser = argparse.ArgumentParser(
        description='Quantify target sequences using standard curves and apply detection thresholds and assess/correct for read mapping errors'
    )
    parser.add_argument('--config', '-c', default='Config/sample_list.txt', help='Config file of samples (default: Config/sample_list.txt)')
    parser.add_argument('--spike-info', '-s', default='Config/spike_in_info.txt', help='Table of spike-in information (Sample, Library Mass, DNA Concentration, Spike-in Fraction, ssDNA Fraction) (default: Config/spike_in_info.txt)')
    parser.add_argument('--targets', '-T', required=True, help='Fasta file with target sequences (default: none, required)')
    parser.add_argument('--targets-bam-dir', '-Tb', default='Mapping/', help='Directory containing sorted bam files of mapping reads to target database named as [sample_name]_[target_name]_sorted.bam (default: Mapping/)')
    parser.add_argument('--target-name', '-N', required=True, help='Name for target database (default: none, required, used for naming output files)')
    parser.add_argument('--window-size', '-w', type=int, default=49, help='Window size for sliding window analysis, must be the same as used in standard_curve_generator.sh (default: 49)')
    parser.add_argument('--detection-threshold', '-detect', default='Regressions/detection/Langenfeld_2025_E_detect.json', help='Detection threshold json file (default: Regressions/detect/Langenfeld_2025_E_detect.json)')
    parser.add_argument('--read-depth-variability-model1', '-rdv1', default='Regressions/read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl', help='Read depth variability model for 0-10 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl)')
    parser.add_argument('--read-depth-variability-model2', '-rdv2', default='Regressions/read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl', help='Read depth variability model for 10-100 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl)')
    parser.add_argument('--read-depth-variability-model3', '-rdv3', default='Regressions/read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl', help='Read depth variability model for 100-1000 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl)')
    parser.add_argument('--read-depth-variability-model4', '-rdv4', default='Regressions/read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl', help='Read depth variability model for >=1000 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl)')
    parser.add_argument('--rmse-cutoff-function1', '-rmse1', default='Regressions/threshold_read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl', help='RMSE cutoff function for 0-10 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl)')
    parser.add_argument('--rmse-cutoff-function2', '-rmse2', default='Regressions/threshold_read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl', help='RMSE threshold function model for 10-100 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl)')
    parser.add_argument('--rmse-cutoff-function3', '-rmse3', default='Regressions/threshold_read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl', help='RMSE threshold function model for 100-1000 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl)')
    parser.add_argument('--rmse-cutoff-function4', '-rmse4', default='Regressions/threshold_read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl', help='RMSE threshold function model for >=1000 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl)')
    parser.add_argument('--output-dir', '-o', default='QuantMeta/', help='Directory to project outputs (default: QuantMeta/)')
    parser.add_argument('--cores', '-j', type=int, default=4, help='Number of threads (default: 4)')

    args, unknown = parser.parse_known_args()
    script_path = get_script_path('quant_targets.sh')

    if not script_path.exists():
        print(f"Error: Script not found at {script_path}", file=sys.stderr)
        sys.exit(1)

    cmd_args = [
        'bash', str(script_path),
        '--config', args.config,
        '--spike-info', args.spike_info,
        '--targets', args.targets,
        '--targets-bam-dir', args.targets_bam_dir,
        '--target-name', args.target_name,
        '--window-size', str(args.window_size),
        '--detection-threshold', args.detection_threshold,
        '--read-depth-variability-model1', args.read_depth_variability_model1,
        '--read-depth-variability-model2', args.read_depth_variability_model2,
        '--read-depth-variability-model3', args.read_depth_variability_model3,
        '--read-depth-variability-model4', args.read_depth_variability_model4,
        '--rmse-cutoff-function1', args.rmse_cutoff_function1,
        '--rmse-cutoff-function2', args.rmse_cutoff_function2,
        '--rmse-cutoff-function3', args.rmse_cutoff_function3,
        '--rmse-cutoff-function4', args.rmse_cutoff_function4,
        '--output-dir', args.output_dir,
        '--cores', str(args.cores),
    ]

    cmd_args.extend(unknown)
    result = subprocess.run(cmd_args, check=False)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
