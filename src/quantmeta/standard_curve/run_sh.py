#!/usr/bin/env python3
"""Wrapper script to run standard_curve_generator.sh from the package."""

import argparse
import subprocess
import sys
from pathlib import Path
import importlib.resources as resources


def get_script_path(script_name: str) -> Path:
    """Get the path to a shell script in the package."""
    try:
        with resources.as_file(resources.files('quantmeta.standard_curve') / script_name) as p:
            return Path(p)
    except Exception:
        return Path(__file__).parent / script_name


def main():
    parser = argparse.ArgumentParser(
        description='Generate standard curves for each metagenome and establish quantification correction regressions for the dataset'
    )
    parser.add_argument('--config', '-c', default='Config/sample_list.txt', help='Config file of samples (default: Config/sample_list.txt)')
    parser.add_argument('--bam-dir', '-b', default='Mapping/', help='Directory containing sorted bam files of mapping reads to standards named as [sample_name]_standards_sorted.bam (default: Mapping/)')
    parser.add_argument('--standards', '-std', default='Map_Indexes/Langenfeld_2025_standards.fasta', help='Fasta file with standard sequences (default: Map_Indexes/Langenfeld_2025_standards.fasta)')
    parser.add_argument('--dsDNA-std-file', '-mix', default='Spike-ins/sequins_Mix_A.txt', help='Table of dsDNA standards (ID, Mass, Rel_Abund, length) (default: Spike-ins/sequins_Mix_A.txt)')
    parser.add_argument('--ssDNA-std-file', '-ssmix', default=None, help='Optional table of ssDNA standards (ID, Mass, Rel_Abund, length) (default: none, optional for Langenfeld et al. 2025 ssDNA standards specify Spike-ins/ssDNA_stds.txt)')
    parser.add_argument('--spike-in-info', '-spike', default='Config/spike_in_info.txt', help='Table of sample-specific spike-in information (Sample, Library_Mass (ng), DNA_Extract_Conc (ng/µL), Spike_Frac, ssDNA (0/Spike_Frac) (default: Config/spike_in_info.txt)')
    parser.add_argument('--detection-threshold', '-detect', default='Regressions/detection/Langenfeld_2025_E_detect.json', help='Detection threshold json file (default: Regressions/detection/Langenfeld_2025_E_detect.json)')
    parser.add_argument('--window-size', '-w', type=int, default=49, help='Window size for sliding window analysis (default: 49)')
    parser.add_argument('--output-dir', '-o', default='QuantMeta/', help='Directory for output files (default: QuantMeta/)')
    parser.add_argument('--cores', '-j', type=int, default=4, help='Number of threads (default: 4)')

    args, unknown = parser.parse_known_args()
    script_path = get_script_path('standard_curve_generator.sh')

    if not script_path.exists():
        print(f"Error: Script not found at {script_path}", file=sys.stderr)
        sys.exit(1)

    cmd_args = [
        'bash', str(script_path),
        '--config', args.config,
        '--bam-dir', args.bam_dir,
        '--standards', args.standards,
        '--dsDNA-std-file', args.dsDNA_std_file,
        '--spike-in-info', args.spike_in_info,
        '--detection-threshold', args.detection_threshold,
        '--window-size', str(args.window_size),
        '--output-dir', args.output_dir,
        '--cores', str(args.cores),
    ]

    if args.ssDNA_std_file:
        cmd_args.extend(['--ssDNA-std-file', args.ssDNA_std_file])

    cmd_args.extend(unknown)
    result = subprocess.run(cmd_args, check=False)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
