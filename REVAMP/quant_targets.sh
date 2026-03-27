#!/bin/bash
set -euo pipefail

### standard_curve_generator.sh
### Generate standard curves relating relative abundance to absolute abundance.
### Works towards:
### - mapping reads to selected targets (bowtie2 + samtools)
### - applying detection thresholds (default: Langenfeld_2025_E_detect)
### - applying quantification correction
### - quant_targets.py (determine concentrations of targets)

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
  -w, --window-size N             Window size for sliding window analysis, must be the same as used in standard_curve_generator.sh (default: 49)
  -detect, --detection-threshold FILE  Detection threshold json file (default: Regressions/detect/Langenfeld_2025_E_detect.json)
  -rdv1, --read-depth-variability-model1 FILE  Read depth variability model for 0-10 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl)
  -rdv2, --read-depth-variability-model2 FILE  Read depth variability model for 10-100 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl)
  -rdv3, --read-depth-variability-model3 FILE  Read depth variability model for 100-1000 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl)
  -rdv4, --read-depth-variability-model4 FILE  Read depth variability model for >1000 reads/bp (default: Regressions/read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl)
  -rmse1, --rmse-cutoff-function1 FILE  RMSE threshold function model for 0-10 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl)
  -rmse2, --rmse-cutoff-function2 FILE  RMSE threshold function model for 10-100 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl)
  -rmse3, --rmse-cutoff-function3 FILE  RMSE threshold function model for 100-1000 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl)
  -rmse4, --rmse-cutoff-function4 FILE  RMSE threshold function model for >1000 reads/bp (default: Regressions/threshold_read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl)
  -j, --cores N                   Number of cores (default: 4)
  -mem, --memory N                Memory per CPU (default: 10gb)
  -t, --time TIME                 Time limit (default: 2-00:00:00)
  -h, --help                      Show this help message

Example:
  $0 --config Config/sample_list.txt --fastq-dir Reads/ --targets Map_Indexes/RefSeq_Viruses.fasta --target-name RefSeq_Viruses --window-size 49 --detection-threshold Regressions/detect/Langenfeld_2025_E_detect.json --read-depth-variability-model1 Regressions/read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl --read-depth-variability-model2 Regressions/read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl --read-depth-variability-model3 Regressions/read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl --read-depth-variability-model4 Regressions/read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl --rmse-cutoff-function1 Regressions/threshold_read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl --rmse-cutoff-function2 Regressions/threshold_read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl --rmse-cutoff-function3 Regressions/threshold_read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl --rmse-cutoff-function4 Regressions/threshold_read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl
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

bash Scripts/quant_targets_unknown.sh ${SLURM_ARRAY_TASK_ID} $config $fq $targets $target_name $window_size $detect $rdv1 $rdv2 $rdv3 $rdv4 $rmse1 $rmse2 $rmse3 $rmse4