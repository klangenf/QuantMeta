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
  -std, --standards FILE          Fasta file with standard sequences (default: Map_Indexes/Langenfeld_2025_standards.fasta)
  -mix, --dsDNA-std-file FILE     Table of dsDNA standards (ID, Mass, Rel_Abund, length) (default: Spike-ins/sequins_Mix_A.txt)
  -ssmix, --ssDNA-std-file FILE   Optional table of ssDNA standards (ID, Mass, Rel_Abund, length)
  -spike, --spike-in-info FILE    Table of sample-specific spike-in information (Sample, Library_Mass (ng), DNA_Extract_Conc (ng/µL), Spike_Frac, ssDNA (0/Spike_Frac) (default: Config/spike_in_info.txt)
  -detect, --detection-threshold FILE  Detection threshold json file (default: Regressions/detect/Langenfeld_2025_E_detect.json)
  -w, --window-size N             Window size for sliding window analysis (default: 49)
  -j, --cores N                   Number of cores (default: 2)
  -mem, --memory N                Memory per CPU (default: 10gb)
  -t, --time TIME                 Time limit (default: 02:00:00)
  -h, --help                      Show this help message

Example:
  $0 --config Config/config_example_samples.yaml --fastq-dir Reads/ --std Map_Indexes/Langenfeld_2025_standards.fasta --mix Spike-ins/sequins_Mix_A.txt --spike-in-info Config/spike_in_info.txt --detection-threshold Regressions/detect/Langenfeld_2025_E_detect.json --window-size 49
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

bash Scripts/standard_curve_analysis.sh ${SLURM_ARRAY_TASK_ID} $config $fq $standards $window_size

### once all of the array jobs are complete, run the regression builder
if [[ $(squeue -u $USER -n $(basename $0) -h | wc -l) -eq 0 ]]; then
    python3 Scripts/quantmeta.py --sample-info $spike --dsDNA-std-mixes $mix --ssDNA-std-mixes $ssmix --detect-threshold $detect
    python3 Scripts/quant_correct_regression_builder.py --sample-names $config