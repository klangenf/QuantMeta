#!/usr/bin/env python3
"""Wrapper script to run detection_threshold_generator.sh from the package."""

import argparse
import subprocess
import sys
from pathlib import Path
import importlib.resources as resources


def get_script_path(script_name: str) -> Path:
    """Get the path to a shell script in the package."""
    # Try to find the script in the package
    try:
        # When installed as package, look in package resources
        with resources.as_file(resources.files('quantmeta.detection_threshold') / script_name) as p:
            return Path(p)
    except Exception:
        # Fallback: look relative to this file (development mode)
        return Path(__file__).parent / script_name


def main():
    parser = argparse.ArgumentParser(
        description='Generate entropy-based detection thresholds with user-specified minimum coverage and read distribution thresholds. This step is optional if users need different thresholds than those established by Langenfeld et al. 2025.'
    )
    parser.add_argument('--config', '-c', default='Config/sample_list.txt', help='Config file of samples (default: Config/sample_list.txt)')
    parser.add_argument('--fastq-dir', '-fq', default='Reads/', help='Directory containing deinterleaved fastq files named as [sample_name]_R1.fastq.gz and [sample_name]_R2.fastq.gz (default: Reads/)')
    parser.add_argument('--standards', '-std', default='Map_Indexes/Langenfeld_2025_standards.fasta', help='Fasta file with standard sequences (default: Map_Indexes/Langenfeld_2025_standards.fasta)')
    parser.add_argument('--dsDNA-std-file', '-mix', default='Spike-ins/sequins_Mix_A.txt', help='Table of dsDNA standards (ID, Mass, Rel_Abund, length) (default: Spike-ins/sequins_Mix_A.txt)')
    parser.add_argument('--ssDNA-std-file', '-ssmix', default=None, help='Optional table of ssDNA standards (ID, Mass, Rel_Abund, length) (default: none, optional for Langenfeld et al. 2025 ssDNA standards specify Spike-ins/ssDNA_stds.txt)')
    parser.add_argument('--test-database', '-test', required=True, help='Fasta file with test sequences (default: none, required for regression builder)')
    parser.add_argument('--test-bam-dir', '-tb', default='Mapping/', help='Directory containing sorted bam files of mapping reads to test database named as [sample_name]_[test_name]_sorted.bam (default: Mapping/)')
    parser.add_argument('--test-name', '-tn', required=True, help='Name for test database (default: none, required, used for naming output files)')
    parser.add_argument('--min-coverage', '-min_cov', type=float, default=0.1, help='Minimum read coverage threshold for detection (default: 0.1)')
    parser.add_argument('--min-distribution', '-min_dist', type=float, default=0.3, help='Minimum read distribution threshold for detection (default: 0.3)')
    parser.add_argument('--output-dir', '-o', default='QuantMeta/', help='Directory for output files (default: QuantMeta/)')
    parser.add_argument('--cores', '-j', type=int, default=4, help='Number of threads (default: 4)')

    args, unknown = parser.parse_known_args()

    # Get the path to the shell script
    script_path = get_script_path('detection_threshold_generator.sh')

    if not script_path.exists():
        print(f"Error: Script not found at {script_path}", file=sys.stderr)
        sys.exit(1)

    # Build command arguments
    cmd_args = [
        'bash',
        str(script_path),
        '--config', args.config,
        '--fastq-dir', args.fastq_dir,
        '--standards', args.standards,
        '--dsDNA-std-file', args.dsDNA_std_file,
        '--test-database', args.test_database,
        '--test-bam-dir', args.test_bam_dir,
        '--test-name', args.test_name,
        '--min-coverage', str(args.min_coverage),
        '--min-distribution', str(args.min_distribution),
        '--output-dir', args.output_dir,
        '--cores', str(args.cores),
    ]

    # Add optional arguments
    if args.ssDNA_std_file:
        cmd_args.extend(['--ssDNA-std-file', args.ssDNA_std_file])

    # Add any unknown arguments (for forward compatibility)
    cmd_args.extend(unknown)

    # Run the shell script
    result = subprocess.run(cmd_args, check=False)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
