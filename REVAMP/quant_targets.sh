#!/bin/bash
set -euo pipefail

### quant_targets.sh
### Generate quantification results for selected targets.
### Works towards:
### - applying detection thresholds (default: Langenfeld_2025_E_detect)
### - applying quantification correction
### - quant_targets.py (determine concentrations of targets)

# Base directory for QuantMeta project
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"

# Default values
config="Config/sample_list.txt"
spike_info="Config/spike_in_info.txt"
targets=""
targets_bam_dir="Mapping/"
target_name=""
window_size=49
detect="Regressions/detection/Langenfeld_2025_E_detect.json"
read_depth_variability_model1="Regressions/read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl"
read_depth_variability_model2="Regressions/read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl"
read_depth_variability_model3="Regressions/read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl"
read_depth_variability_model4="Regressions/read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl"
rmse_cutoff_function1="Regressions/threshold_read_depth_variability/Langenfeld_2025_0to10readsperbp.pkl"
rmse_cutoff_function2="Regressions/threshold_read_depth_variability/Langenfeld_2025_10to100readsperbp.pkl"
rmse_cutoff_function3="Regressions/threshold_read_depth_variability/Langenfeld_2025_100to1000readsperbp.pkl"
rmse_cutoff_function4="Regressions/threshold_read_depth_variability/Langenfeld_2025_gte1000readsperbp.pkl"
cores=4
memory="10gb"
time="2-00:00:00"

usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -c, --config FILE               Config file of samples (default: Config/sample_list.txt)
  -s, --spike-info FILE           Table of spike-in information (Sample, Library Mass, DNA Concentration, Spike-in Fraction, ssDNA Fraction) (default: Config/spike_in_info.txt)
  -T, --targets FILE              Fasta file with target sequences (default: none, required)
  -Tb, --targets-bam-dir DIR      Directory containing sorted bam files of mapping reads to target database named as {sample_name}_{target_name}_sorted.bam (default: Mapping/)
  -N, --target-name NAME          Name for target database (default: none, required, used for naming output files)
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

# Option parsing
while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--config) config="$2"; shift 2 ;;
    -s|--spike-info) spike_info="$2"; shift 2 ;;
    -T|--targets) targets="$2"; shift 2 ;;
    -Tb|--targets-bam-dir) targets_bam_dir="$2"; shift 2 ;;
    -N|--target-name) target_name="$2"; shift 2 ;;
    -w|--window-size) window_size="$2"; shift 2 ;;
    -detect|--detection-threshold) detect="$2"; shift 2 ;;
    -rdv1|--read-depth-variability-model1) rdv1="$2"; shift 2 ;;
    -rdv2|--read-depth-variability-model2) rdv2="$2"; shift 2 ;;
    -rdv3|--read-depth-variability-model3) rdv3="$2"; shift 2 ;;
    -rdv4|--read-depth-variability-model4) rdv4="$2"; shift 2 ;;
    -rmse1|--rmse-cutoff-function1) rmse1="$2"; shift 2 ;;
    -rmse2|--rmse-cutoff-function2) rmse2="$2"; shift 2 ;;
    -rmse3|--rmse-cutoff-function3) rmse3="$2"; shift 2 ;;
    -rmse4|--rmse-cutoff-function4) rmse4="$2"; shift 2 ;;
    -j|--cores) cores="$2"; shift 2 ;;
    -mem|--memory) memory="$2"; shift 2 ;;
    -t|--time) time="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

bioawk -c fastx '{{ print $name, length($seq) }}' < $targets > Map_Indexes/${target_name}_lengths.txt

# Get number of samples
num_samples=$(wc -l < "$config")

# Run analysis for each sample in parallel (up to $cores at a time)
seq 0 $((num_samples - 1)) | xargs -n 1 -P "$cores" -I {} bash -c 'bash Scripts/quant_targets_unknown.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "${12}" "${13}" "${14}" "${15}" "${16}"' _ {} "$config" "$targets_bam_dir" "$targets" "$target_name" "$window_size" "$detect" "$read_depth_variability_model1" "$read_depth_variability_model2" "$read_depth_variability_model3" "$read_depth_variability_model4" "$rmse_cutoff_function1" "$rmse_cutoff_function2" "$rmse_cutoff_function3" "$rmse_cutoff_function4" "$spike_info"