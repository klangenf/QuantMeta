#!/bin/bash
set -euo pipefail

### detection_threshold_generator.sh
### Generate length-dependent detection threshold regression from workflow outputs.
### Works towards:
### - generating failure standards (negative control dataset)
### - downsampling reads (expand dynamic range of standards)
### - mapping reads to standards (bowtie2 + samtools)
### - mapping reads to test database (bowtie2 + samtools) 
### - confident_detection_regression_builder.py (regression model builder)

# Base directory for QuantMeta project
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"

usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -c, --config FILE               Config file of samples (default: Config/sample_list.txt)
  -fq, --fastq-dir DIR            Directory containing deinterleaved fastq files named as {sample_name}_R1.fastq.gz and {sample_name}_R2.fastq.gz (default: Reads/)
  -std, --standards FILE          Fasta file with standard sequences (default: Map_Indexes/Langenfeld_2025_standards.fasta)
  -mix, --dsDNA-std-file FILE     Table of dsDNA standards (ID, Mass, Rel_Abund, length) (default: Spike-ins/sequins_Mix_A.txt)
  -ssmix, --ssDNA-std-file FILE   Optional table of ssDNA standards (ID, Mass, Rel_Abund, length)
  -test, --test-database FILE     Fasta file with test sequences (default: Map_Indexes/RefSeq_viruses.fasta)
  -tb, --test-bam-dir DIR         Directory containing sorted bam files of mapping reads to test database named as {sample_name}_{test_name}_sorted.bam (default: Mapping/example_1_RefSeq_viruses_sorted.bam)
  -tn, --test-name NAME           Name for test database (default: RefSeq_viruses)
  -min_cov, --min-coverage FLOAT  Minimum read coverage threshold for detection (default: 0.1)
  -min_dist, --min-distribution FLOAT Minimum read distribution threshold for detection (default: 0.3)
  -j, --cores N                   Number of cores (default: 4)
  -mem, --memory N                Memory per CPU (default: 10gb)
  -t, --time TIME                 Time limit (default: 2-00:00:00)
  -h, --help                      Show this help message

Example:
  $0 --config Config/config_example_samples.yaml --fastq-dir Reads/ --std Map_Indexes/Langenfeld_2025_standards.fasta --mix Spike-ins/sequins_Mix_A.txt --test Map_Indexes/RefSeq_viruses.fasta --test-name RefSeq_viruses --min-coverage 0.1 --min-distribution 0.3
EOF
  exit 1
}

# Resources
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$cores
#SBATCH --mem-per-cpu=$memory
#SBATCH --time=$time

# Environment
##SBATCH --export=ALL
#SBATCH --array=0-$(($(wc -l < $config) - 1))

echo ${SLURM_ARRAY_TASK_ID}

bash Scripts/detection_threshold_analysis.sh ${SLURM_ARRAY_TASK_ID} $config $standards $tn $test $fq $tb

### once all of the array jobs are complete, run the regression builder
if [[ $(squeue -u $USER -n $(basename $0) -h | wc -l) -eq 0 ]]; then
    python3 Scripts/confident_detection_regression_builder.py --sample-file $config --std-file $mix --ssdna-std-file $ssmix --min_coverage $min_coverage --min_distribution $min_distribution --test-database $tn