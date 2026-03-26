#!/bin/bash
set -euo pipefail

### standard_curve_generator.sh
### Generate standard curves relating relative abundance to absolute abundance.
### Works towards:
### - mapping reads to standards (bowtie2 + samtools)
### - applying detection thresholds (default: Langenfeld_2025_E_detect)
### - quantmeta.py (standard curve builder)

# Base directory for QuantMeta project
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"

usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -c, --config FILE               Config file of samples (default: Config/sample_list.txt)
  -fq, --fastq-dir DIR            Directory containing deinterleaved fastq files named as {sample_name}_R1.fastq.gz and {sample_name}_R2.fastq.gz (default: Reads/)
  -T, --targets FILE              Fasta file with target sequences (default: Map_Indexes/RefSeq_Viruses.fasta)
  -N, --target-name NAME          Name for target database (default: RefSeq_Viruses)
  -C, --contigs True/False        Whether target database are contigs (default: False)
  -spike, --spike-in-info FILE    Table of sample-specific spike-in information (Sample, Library_Mass (ng), DNA_Extract_Conc (ng/µL), Spike_Frac, ssDNA (0/Spike_Frac) (default: Config/spike_in_info.txt)
  -detect, --detection-threshold FILE  Detection threshold json file (default: Regressions/detect/Langenfeld_2025_E_detect.json)
  -j, --cores N                   Number of cores (default: 2)
  -mem, --memory N                Memory per CPU (default: 10gb)
  -t, --time TIME                 Time limit (default: 02:00:00)
  -h, --help                      Show this help message

Example:
  $0 --config Config/config_example_samples.yaml --fastq-dir Reads/ --std Map_Indexes/Langenfeld_2025_standards.fasta --mix Spike-ins/sequins_Mix_A.txt --spike-in-info Config/spike_in_info.txt --detection-threshold Regressions/detect/Langenfeld_2025_E_detect.json
EOF
exit 1
}

# Resources
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$cores
#SBATCH --mem-per-cpu=$memory
#SBATCH --time=$time

bioawk -c fastx '{{ print $name, length($seq) }}' < $targets > Map_Indexes/${target_name}_lengths.txt
bowtie2-build -f $targets Map_Indexes/${target_name}

# Environment
##SBATCH --export=ALL
#SBATCH --array=0-$(($(wc -l < $config) - 1))

echo ${SLURM_ARRAY_TASK_ID}

bash Scripts/quant_targets.sh ${SLURM_ARRAY_TASK_ID} $config $fq $targets $target_name