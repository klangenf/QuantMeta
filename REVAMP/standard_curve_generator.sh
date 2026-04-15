#!/bin/bash
set -euo pipefail

### standard_curve_generator.sh
### Generate standard curves relating relative abundance to absolute abundance.
### Modified version of standard_curve_generator.sh without SLURM dependencies.
### Works towards:
### - mapping reads to standards (bowtie2 + samtools)
### - applying detection thresholds (default: Langenfeld_2025_E_detect)
### - quantmeta.py (standard curve builder)

# Base directory for QuantMeta project
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASE_DIR"


# Default values
config="Config/sample_list.txt"
bam_dir="Mapping/"
standards="Map_Indexes/Langenfeld_2025_standards.fasta"
mix="Spike-ins/sequins_Mix_A.txt"
ssmix=""
spike="Config/spike_in_info.txt"
detect="Regressions/detection/Langenfeld_2025_E_detect.json"
window_size=49
cores=2
memory="10gb"
time="02:00:00"

usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -c, --config FILE               Config file of samples (default: Config/sample_list.txt)
  -b, --bam-dir DIR               Directory containing sorted bam files of mapping reads to standards named as {sample_name}_standards_sorted.bam (default: Mapping/)
  -std, --standards FILE          Fasta file with standard sequences (default: Map_Indexes/Langenfeld_2025_standards.fasta)
  -mix, --dsDNA-std-file FILE     Table of dsDNA standards (ID, Mass, Rel_Abund, length) (default: Spike-ins/sequins_Mix_A.txt)
  -ssmix, --ssDNA-std-file FILE   Optional table of ssDNA standards (ID, Mass, Rel_Abund, length)
  -spike, --spike-in-info FILE    Table of sample-specific spike-in information (Sample, Library_Mass (ng), DNA_Extract_Conc (ng/µL), Spike_Frac, ssDNA (0/Spike_Frac) (default: Config/spike_in_info.txt)
  -detect, --detection-threshold FILE  Detection threshold json file (default: Regressions/detection/Langenfeld_2025_E_detect.json)
  -w, --window-size N             Window size for sliding window analysis (default: 49)
  -j, --cores N                   Number of cores (default: 2)
  -mem, --memory N                Memory per CPU (default: 10gb)
  -t, --time TIME                 Time limit (default: 02:00:00)
  -h, --help                      Show this help message

Example:
  $0 --config Config/sample_list.txt --bam-dir Mapping/ --standards Map_Indexes/Langenfeld_2025_standards.fasta --mix Spike-ins/sequins_Mix_A.txt --spike Config/spike_in_info.txt --detect Regressions/detect/Langenfeld_2025_E_detect.json --window-size 49
EOF
  exit 1
}

# Option parsing
while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--config) config="$2"; shift 2 ;;
    -b|--bam-dir) bam_dir="$2"; shift 2 ;;
    -std|--standards) standards="$2"; shift 2 ;;
    -mix|--dsDNA-std-file) mix="$2"; shift 2 ;;
    -ssmix|--ssDNA-std-file) ssmix="$2"; shift 2 ;;
    -spike|--spike-in-info) spike="$2"; shift 2 ;;
    -detect|--detection-threshold) detect="$2"; shift 2 ;;
    -w|--window-size) window_size="$2"; shift 2 ;;
    -j|--cores) cores="$2"; shift 2 ;;
    -mem|--memory) memory="$2"; shift 2 ;;
    -t|--time) time="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done


# Get number of samples
num_samples=$(wc -l < "$config")

# Run analysis for each sample in parallel (up to $cores at a time)
seq 0 $((num_samples - 1)) | xargs -n 1 -P "$cores" -I {} bash -c 'bash Scripts/standard_curve_analysis.sh "$1" "$2" "$3" "$4" "$5"' _ {} "$config" "$bam_dir" "$standards" "$window_size"

# Run the regression builder
echo "Running regression builder"
python3 Scripts/quantmeta.py --sample-info "$spike" --dsDNA-std-mixes "$mix" --ssDNA-std-mixes "$ssmix" --detect-threshold "$detect"
python3 Scripts/quant_correct_regression_builder.py --sample-names "$config"