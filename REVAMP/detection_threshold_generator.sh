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

# Default values
config="Config/sample_list.txt"
fastq_dir="Reads/"
standards="Map_Indexes/Langenfeld_2025_standards.fasta"
mix="Spike-ins/sequins_Mix_A.txt"
ssmix=""
test=""
test_bam_dir="Mapping/"
test_name=""
min_coverage=0.1
min_distribution=0.3
output_dir="QuantMeta/"
cores=4
memory="10gb"
time="2-00:00:00"


usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -c, --config FILE               Config file of samples (default: Config/sample_list.txt)
  -fq, --fastq-dir DIR            Directory containing deinterleaved fastq files named as {sample_name}_R1.fastq.gz and {sample_name}_R2.fastq.gz (default: Reads/)
  -std, --standards FILE          Fasta file with standard sequences (default: Map_Indexes/Langenfeld_2025_standards.fasta)
  -mix, --dsDNA-std-file FILE     Table of dsDNA standards (ID, Mass, Rel_Abund, length) (default: Spike-ins/sequins_Mix_A.txt)
  -ssmix, --ssDNA-std-file FILE   Optional table of ssDNA standards (ID, Mass, Rel_Abund, length)
  -test, --test-database FILE     Fasta file with test sequences (default: none, required for regression builder)
  -tb, --test-bam-dir DIR         Directory containing sorted bam files of mapping reads to test database named as {sample_name}_{test_name}_sorted.bam (default: Mapping/)
  -tn, --test-name NAME           Name for test database (default: none, required, used for naming output files)
  -min_cov, --min-coverage FLOAT  Minimum read coverage threshold for detection (default: 0.1)
  -min_dist, --min-distribution FLOAT Minimum read distribution threshold for detection (default: 0.3)
  -o, --output-dir DIR            Directory for output files (default: QuantMeta/)
  -j, --cores N                   Number of cores (default: 4)
  -mem, --memory N                Memory per CPU (default: 10gb)
  -t, --time TIME                 Time limit (default: 2-00:00:00)
  -h, --help                      Show this help message

Example:
  $0 --config Config/config_example_samples.yaml --fastq-dir Reads/ --std Map_Indexes/Langenfeld_2025_standards.fasta --mix Spike-ins/sequins_Mix_A.txt --test Map_Indexes/RefSeq_viruses.fasta --test-name RefSeq_viruses --min-coverage 0.1 --min-distribution 0.3
EOF
  exit 1
}

# Option parsing
while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--config) config="$2"; shift 2 ;;
    -fq|--fastq-dir) fastq_dir="$2"; shift 2 ;;
    -std|--standards) standards="$2"; shift 2 ;;
    -mix|--dsDNA-std-file) mix="$2"; shift 2 ;;
    -ssmix|--ssDNA-std-file) ssmix="$2"; shift 2 ;;
    -test|--test-database) test="$2"; shift 2 ;;
    -tb|--test-bam-dir) test_bam_dir="$2"; shift 2 ;;
    -tn|--test-name) test_name="$2"; shift 2 ;;
    -min_cov|--min-coverage) min_coverage="$2"; shift 2 ;;
    -min_dist|--min-distribution) min_distribution="$2"; shift 2 ;;
    -o|--output-dir) output_dir="$2"; shift 2 ;;
    -j|--cores) cores="$2"; shift 2 ;;
    -mem|--memory) memory="$2"; shift 2 ;;
    -t|--time) time="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

bowtie2-build -f $standards ${output_dir}/Map_Indexes/standards
bioawk -c fastx '{{ print $name, length($seq) }}' < $test > ${output_dir}/Map_Indexes/${test_name}_lengths.txt

### Curate failure standards (negative controls) and map reads to them
mkdir ${output_dir}/Map_Indexes/fail_standards_r1
mkdir ${output_dir}/Map_Indexes/fail_standards_r2
mkdir ${output_dir}/Map_Indexes/fail_standards_r3
mkdir ${output_dir}/Map_Indexes/fail_standards_r4
mkdir ${output_dir}/Map_Indexes/fail_standards_r5
mutation-simulator -o ${output_dir}/Map_Indexes/fail_standards_r1/fail_standards $standards args -sn 0.1 -in 0.02 -inb 5 -de 0.02 -deb 5 -iv 0.02 -ivb 5 -du 0.02 -dub 5 -tl 0.02 -tlb 5 
mutation-simulator -o ${output_dir}/Map_Indexes/fail_standards_r2/fail_standards ${output_dir}/Map_Indexes/fail_standards_r1/fail_standards_ms.fasta args -sn 0.1 -in 0.02 -inb 5 -de 0.02 -deb 5 -iv 0.02 -ivb 5 -du 0.02 -dub 5 -tl 0.02 -tlb 5 
mutation-simulator -o ${output_dir}/Map_Indexes/fail_standards_r3/fail_standards ${output_dir}/Map_Indexes/fail_standards_r2/fail_standards_ms.fasta args -sn 0.1 -in 0.02 -inb 5 -de 0.02 -deb 5 -iv 0.02 -ivb 5 -du 0.02 -dub 5 -tl 0.02 -tlb 5 
mutation-simulator -o ${output_dir}/Map_Indexes/fail_standards_r4/fail_standards ${output_dir}/Map_Indexes/fail_standards_r3/fail_standards_ms.fasta args -sn 0.1 -in 0.02 -inb 5 -de 0.02 -deb 5 -iv 0.02 -ivb 5 -du 0.02 -dub 5 -tl 0.02 -tlb 5 
mutation-simulator -o ${output_dir}/Map_Indexes/fail_standards_r5/fail_standards ${output_dir}/Map_Indexes/fail_standards_r4/fail_standards_ms.fasta args -sn 0.1 -in 0.02 -inb 5 -de 0.02 -deb 5 -iv 0.02 -ivb 5 -du 0.02 -dub 5 -tl 0.02 -tlb 5

for fail_set in r1 r2 r3 r4 r5; do
    bioawk -c fastx '{{ print $name, length($seq) }}' < ${output_dir}/Map_Indexes/fail_standards_${fail_set}/fail_standards_ms.fasta > ${output_dir}/Map_Indexes/fail_standards_${fail_set}/fail_standards_lengths.txt
    bowtie2-build -f ${output_dir}/Map_Indexes/fail_standards_${fail_set}/fail_standards_ms.fasta ${output_dir}/Map_Indexes/fail_standards_${fail_set}/fail_standards
done

# Get number of samples
num_samples=$(wc -l < "$config")

# Run analysis for each sample in parallel (up to $cores at a time)
seq 0 $((num_samples - 1)) | xargs -n 1 -P "$cores" -I {} bash -c 'bash Scripts/detection_threshold_analysis.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"' _ {} "$config" "$standards" "$test_name" "$test" "$fastq_dir" "$test_bam_dir" "$output_dir"

### once all of the array jobs are complete, run the regression builder
python3 Scripts/confident_detection_regression_builder.py --sample-file $config --std-file $mix --ssdna-std-file $ssmix --min_coverage $min_coverage --min_distribution $min_distribution --test-database $test_name