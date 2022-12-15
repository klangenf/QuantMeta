# QuantMeta
QuantMeta uses the results of sequencing metagenomes spiked with known concentrations of standards to determine the absolute abundance of target genes, genomes, or contigs within metagenomes. The pipeline provides a method for applying rigorous quantitative limitations to potential targets to improve the detection confidence and quantification accuracy.

See more details in the publication. The pipeline was developed using a set of 86 dsDNA standards developed by Hardwick et al. (2018) complimented with a set of 5 ssDNA standards. However, the pipeline allows the spike-in standards to be amended based on the users selected standards.

Please cite our work: Langenfeld K, Hegarty B, Vidaurri S, Crossette E, Duhaime MB, Wigginton K. Evaluating limitations of quantitative metagenomics with synthetic dsDNA and ssDNA standards. bioRxiv. 2022:2022.07.08.499345.

# Installation
QuantMeta has several steps that run with snakemake. Snakemake will automatically install dependencies for each step using conda. However, the user must first set up snakemake. We recommend setting mamba has your base installer. Additional details on snakemake installation may be found at https://snakemake.readthedocs.io/en/stable/getting_started/installation.html. 

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
git clone https://github.com/klangenf/QuantMeta.git
cd QuantMeta
```
# Dependencies
Conda environments for dependencies can be created from the QuantMeta directory with the following commands:
- Bioawk: ```mamba env create --file Envs/bioawk.yaml```
- Bowtie2: ```mamba env create --file Envs/bowtie2.yaml```

# Run Instructions
### 1. Update directories to include sample and spike-in information. 
#### 1.1 Configure Config directory with a sample list and hpc submission information.
- Create a sample list file using the same format as config_example_samples.yaml. Format:
```
samples:
    sample_name_1
    sample_name_2
    ...
    sample_name_x
```
- Update cluster.yaml with your account and email information. If jobs fail due to memory errors or wall time limitations, the run allocations can be updated in the cluster.yaml script.
#### 1.2 Configure Reads directory by adding symbolic links to the quality controlled fastq files.
- Create symbolic links to each sample's quality controlled fastq file. Label each fast file as {sample_id}.fastq.
```
cd Reads
ln -s {location of qc read fastq file} {sample_id}.fastq
```
#### 1.3 Configure Sample_Characteristics directory with concentration factor and recovery results. 
- Edit the example_sample_extraction_info.txt file to include the recovery and concentration factor (CF) of microbiomes from sample collection to DNA extraction. 
    - If this information is unknown or the concentration of targets in DNA extracts is desired, fill in CF and recovery columns with ones.
#### 1.4 Configure Spike-ins directory to provide information on the spike-in standards.
- Provide information on each sample using the format in the example_sample_info.txt template.
    - "Sample": sample names matching the list in the "Config/config_{sample name}_samples.yaml"
    - "mapping_loc": location of bowtie2 mapping results, the results will automatically end up in "Mapping/{sample name}/standards"
    - "lib_mass": DNA mass used for the library preparations (ng)
    - "mix": column of "Spike-ins/STD_MIXES.txt" representing the mix of dsDNA standards spiked in, if a different spike-in was used, then use the column name for the amended "Spike-ins/STD_MIXES.txt" table
    - "dsDNA_spike": dsDNA standards fraction of total DNA mass
    - "ssDNA": "YES" if ssDNA standards were spiked in, otherwise "NO"
    - "ssDNA_spike": ssDNA standards fraction of total DNA mass
    - "downsample": percent downsampling (if no downsampling was performed use 100)
    - "DNA_conc": concentration of DNA in the DNA extract (ng/µL)
```
"Sample"	"mapping_loc"	"lib_mass"	"mix"	"dsDNA_spike"	"ssDNA"	"ssDNA_spike"	"downsample"	"DNA_conc"
"example_1"	"Mapping/example_1/standards"	50	"MIX_A"	0.010229	"YES"	5.7e-05	100	41.4
"example_2"	"Mapping/example_2/standards"	50	"MIX_A"	0.010229	"YES"	5.7e-05	100	128
```
- *Optional:* Provide ddPCR measurements of each spiked-in standard (if applicable) in the same format as example_stds_ddPCR_conc.txt template.
- Check that STD_MIXES.txt includes all relevant spiked-in standard mixes, adjust by adding or removing information as necessary for your particular samples.
#### 1.5 *Optional (only for reference-based quantification):* Configure Map_Indexes directory with information about the target sequences that reads will be mapped onto. 
- For each set of targets, the following steps must be complete. Update {target} in the code to reflect sequences within the fasta files that reads will be mapped to.
  - Download relevant gene or genome fasta files here or create symbolic links to relevant fasta files. Label each fasta file as {target}.fasta (keep target names consistent for snakemake purposes).
  - Create files listing the length of each contig with bioawk.
  - Create bowtie2 indexes in this folder (this bowtie2 line of code may be edited if a fasta file is large, see bowtie2 manual).
```
conda activate bioawk
cd Map_Indexes
bioawk -c fastx '{ print $name, length($seq) }' < {target}.fasta > {target}_lengths.txt
conda deactivate
conda activate bowtie2
bowtie2-build {target}.fasta {target}
```
#### 1.6 *Optional (only for contig-based quantification):* Configure Assembly directory to include contig files and a list of sequences lengths.
- Create symbolic links to each sample's assembly fasta file. Label each fast file as {sample_id}.fasta.
- Create files listing the length of each contig. This can be created by installing bioawk conda environment and loading bioawk.
```
conda activate bioawk
cd Assembly
ln -s {location of assembly fasta file} {sample_id}.fasta
bioawk -c fastx '{ print $name, length($seq) }' < {sample_id}.fasta > {sample_id}_length.txt
```
### 2. Map reads to standard sequences 
#### 2.1 Each step of this process is performed by running “Snakefile-map_standards”.
- Update “Snakefile-map_standards” based on your samples
  - Update configfile for your sample list
    - *Line 6* ```configfile: Config/{your sample list}.yaml```
  - Update SAMPLE wildcard to reflect your sample variable name (if changed from sample)
    - *Line 7* ```SAMPLE = config[“{your variable}”]```
    - OR if running Snakefile-map_standards on a subset of samples,
      update *Line 10* ```SAMPLE = [“{your specific sample name}”, “{another specific sample name}”]```
  - If using a different set of spike-in standards than Langenfeld et al. (2022), update the Snakefile to use the other bowtie indexes and .fasta file in Map-Indexes directory
    - *Line 25* ```“Map_Indexes/{your standards index name}"```
    - *Line 64* ```“Map_Indexes/{your standards fasta name}.fasta"```
#### 2.2 Update “hpc_submission/map_standards_hpc_run.sh” to represent your system configuration
- In the Job Submission Options section, Update the allocation information 
- In the Job Commands section, update the PATH to the user specific storage locations for the anaconda environments that snakemake will create and use 
#### 2.3 Run Snakefile-map_standards on all or selected samples
- From the QuantMeta directory:
  ```
  conda activate snakemake
  sbatch hpc_submission/map_standards_hpc_run.sh
  ```
  Snakemake will automatically submit all of the jobs. Output logs from each individual job will be located in the Logs directory.
- *Optional:* We recommend performing a dry run prior to job submission to the hpc from the QuantMeta directory.
  ```
  conda activate snakemake
  snakemake -prn -s Snakefile-map_standards
  ```
  Check that all expected rules and number of executions are listed in DAG of jobs snakemake prints out.
- *Optional:* If your jobs fail for any reason, snakemake will need to be unlocked prior to rerunning. This can be completed with:
  ```
  snakemake -prn -s Snakefile-map_standards --unlock
  ```
### 3. *Optional:* Determine detection threshold
The provided entropy-based detection threshold was created based on a minimum coverage of 10% and read distribution of 0.3 (ratio of observed read distribution to Poisson distribution of reads). If this threshold is deemed inadequate for your analysis needs, a new detection threshold may be created. The following is the list of steps to create a new detection threshold.
#### 3.1 Perform downsampling of reads and map the downsampled reads to the standard sequences
- Update “Snakefile-downsample” based on your samples
  - Update configfile for your sample list
    - *Line 6* ```configfile: Config/{your sample list}.yaml```
  - Update SAMPLE wildcard to reflect your sample variable name (if changed from sample)
    - *Line 7* ```SAMPLE = config[“{your variable}”]```
    - OR if running Snakefile-map_standards on a subset of samples, update *Line 10* ```SAMPLE = [“{your specific sample name}”, “{another specific sample name}”]```
  - Update FRACTION_LABEL wildcard to reflect the fraction of reads to be downsampled (Langenfeld et al. (2022) performed 20% and 1% downsampling)
    - *Line 12* ```FRACTION_LABEL = ["{integer fraction value}", “{another integer fraction value}”]```
  - If using a different set of spike-in standards than Langenfeld et al. (2022), update the Snakefile to use the other bowtie indexes and .fasta file in Map-Indexes directory
    -  *Line 13* ```STANDARDS = ["Map_Indexes/{your standards index name}”]```
    -  *Line 61* ```“Map_Indexes/{your standards index name}"```
    -  *Line 100* ```“Map_Indexes/{your standards fasta name}.fasta"```
- Update “hpc_submission/downsample_hpc_run.sh” to represent your system configuration
  -  In the Job Submission Options section, Update the allocation information 
  -  In the Job Commands section, update the PATH to the user specific storage locations for the anaconda environments that snakemake will create and use 
##### 3.2 Run Snakefile-downsample on all or selected samples
- From the QuantMeta directory:
  ```
  conda activate snakemake
  sbatch hpc_submission/downsample_hpc_run.sh
  ```
  Snakemake will automatically submit all of the jobs. Output logs from each individual job will be located in the Logs directory.
- *Optional:* We recommend performing a dry run prior to job submission to the hpc from the QuantMeta directory.
  ```
  conda activate snakemake
  snakemake -prn -s Snakefile-downsample
  ```
  Check that all expected rules and number of executions are listed in DAG of jobs snakemake prints out.
- *Optional:* If your jobs fail for any reason, snakemake will need to be unlocked prior to rerunning. This can be completed with:
  ```
  snakemake -prn -s Snakefile-downsample --unlock
  ```
#### 3.3 Create a “failure” set of standards and map reads to the standard sequences
- Update “Snakefile-failure_standards” based on your samples
  - Update configfile for your sample list
    - *Line 6* ```configfile: Config/{your sample list}.yaml```
  - Update SAMPLE wildcard to reflect your sample variable name (if changed from sample)
    - *Line 7* ```SAMPLE = config[“{your variable}”]```
    - OR if running Snakefile-map_standards on a subset of samples, update *Line 10* ```SAMPLE = [“{your specific sample name}”, “{another specific sample name}”]```
  - If using a different set of spike-in standards than Langenfeld et al. (2022), update the Snakefile to use the other bowtie indexes and .fasta file in Map-Indexes directory
    - *Line 12* ```STANDARDS = ["Map_Indexes/{your standards index name}”]```
    - *Line 26* ```“Map_Indexes/{your standards fasta name}.fasta"```
- Update “hpc_submission/failure_standards_hpc_run.sh” to represent your system configuration
  - In the Job Submission Options section, Update the allocation information 
  - In the Job Commands section, update the PATH to the user specific storage locations for the anaconda environments that snakemake will create and use 
#### 3.4 Run Snakefile-failure_standards on all or selected samples
- From the QuantMeta directory:
  ```
  conda activate snakemake
  sbatch hpc_submission/failure_standards_hpc_run.sh
  ```
  Snakemake will automatically submit all of the jobs. Output logs from each individual job will be located in the Logs directory.
- *Optional:* We recommend performing a dry run prior to job submission to the hpc from the QuantMeta directory.
  ```
  conda activate snakemake
  snakemake -prn -s Snakefile-failure_standards
  ```
  Check that all expected rules and number of executions are listed in DAG of jobs snakemake prints out.
- *Optional:* If your jobs fail for any reason, snakemake will need to be unlocked prior to rerunning. This can be completed with:
  ```
  snakemake -prn -s Snakefile-failure_standards --unlock
  ```
### 3.5 Follow the instructions in the **confident_detection_regression_builder.Rmd** R notebook.
### 4. Relate relative to absolute abundances of standards in a linear regression
#### 4.1 Each step of this process is performed by running “Snakefile-quantmeta”
- Update “Snakefile-quantmeta” based on your samples
  - Update configfile for your sample list
    - *Line 6* ```configfile: Config/{your sample list}.yaml```
  - Update SAMPLE wildcard to reflect your sample variable name (if changed from sample)
    - *Line 7* ```SAMPLE = config[“{your variable}”]```
    - OR if running Snakefile-map_standards on a subset of samples, update *Line 10* ```SAMPLE = [“{your specific sample name}”, “{another specific sample name}”]```
  - Update rule quantmeta input files
    - *Line 22* ```spikes="Spike-ins/{your specific info table}.txt"```
    - *Line 39* ```spikes="Spike-ins/{your specific info table}.txt"```
    - *Line 23* ```STD_MIXES="Spike-ins/STD_MIXES.txt"``` (check that this STD_MIXES.txt file includes your specific standards, alter if necessary)
    - *Line 24* ```ddPCR="Spike-ins/{your ddPCR/qPCR standard measurements}.txt"``` (Update to your specific ddPCR/qPCR measurements or a blank file if ddPCR or qPCR measurements were not made of spike-in standards)
    - *Line 25* ```detect_reg="Regressions/detection/{your E_detect regression}"``` (Only update if the optional step 3 was performed)
  - Update “hpc_submission/quantmeta_hpc_run.sh” to represent your system configuration
    - In the Job Submission Options section, Update the allocation information 
    - In the Job Commands section, update the PATH to the user specific storage locations for the anaconda environments that snakemake will create and use 
#### 4.2 Run Snakefile-quantmeta on all or selected samples
- From the QuantMeta directory:
  ```
  conda activate snakemake
  sbatch hpc_submission/quantmeta_hpc_run.sh
  ```
- Snakemake will automatically submit all of the jobs. Output logs from each individual job will be located in the Logs directory.
- *Optional:* We recommend performing a dry run prior to job submission to the hpc from the QuantMeta directory.
  ```
  conda activate snakemake
  snakemake -prn -s Snakefile-quantmeta
  ```
  Check that all expected rules and number of executions are listed in DAG of jobs snakemake prints out
- *Optional:* If your jobs fail for any reason, snakemake will need to be unlocked prior to rerunning. This can be completed with:
  ```
  snakemake -prn -s Snakefile-quantmeta --unlock
  ```
### 5. *Optional:* Determine sequencing specific parameters for quantification correction
The read depth variability regressions and RMSE thresholds may be specific to the sequencing technologies used in Langenfeld et al. (2022) where sequencing was performed with Illumina NovaSeq SP flow cells creating 251-bp paired-end reads with Swift 1S Plus library preparations. If different sequencing technologies were used, we recommend creating new read depth variability regressions and assessing RMSE thresholds.
  - Follow instructions in the **quant_correct_regression_builder.Rmd** R notebook.
### 6. Quantify target genes, genomes, or contigs
#### 6.1 Each step of this process is performed by running “Snakefile-quant_targets”
- Update “Snakefile-quant_targets” based on your samples
  - Update configfile to your sample list
    - *Line 6* ```configfile: Config/{your sample list}.yaml```
  - Update SAMPLE wildcard to reflect your sample variable name (if changed from sample)
    - *Line 7* ```SAMPLE = config[“{your variable}”]```
    - OR if running Snakefile-map_standards on a subset of samples, update *Line 10* ```SAMPLE = [“{your specific sample name}”, “{another specific sample name}”]```
  - Update TARGET wildcard to reflect the target that will be quantified (e.g., “contigs”, “NCBI_viral_database”, “phage_HM1”)
    - *Line 12* ```TARGET = [“{your target}”]```
  - Update rule quantmeta_unknowns input files
    - *Line 99* ```sample_info="Spike-ins/{your specific info table}.txt"```
    - *Line 100* ```char="Sample_Characteristics/{your extraction info}.txt"```
    - *Line 104* ```regression="Regressions/quantification/{output of quantmeta}"``` (If using separate regressions for each sample, update this field to include the {sample} wildcard)
    - *Line 105* ```bin_assign= "NA"``` (Update to a file binning contigs. The file should be a tab separated file with two columns: "contig_ID" and "ID" where "contig_ID" provides unique identifiers for each contig and "ID" is the respective bin name for the contig)
    - *Line 109* ```type="database"``` (type should be database if mapping reads to any sequence from a database of previous sequencing, change to contigs or MAGs if mapping reads to a sequence derived from the samples)
    - *Line 110* ```detect_regression="Regressions/detection/{your E_detect regression}"``` (Only update if the optional step 3 was performed)
    - *Line 111* ```read_var_regression_1="Regressions/read_depth_variability/{your 0-10 reads/bp regression}"``` (Only update if the optional step 5 was performed)
    - *Line 112* ```read_var_regression_2="Regressions/read_depth_variability/{your 10-100 reads/bp regression}"``` (Only update if the optional step 5 was performed)
    - *Line 113* ```read_var_regression_3="Regressions/read_depth_variability/{your 100-1000 reads/bp regression}"``` (Only update if the optional step 5 was performed)
    - *Line 114* ```read_var_regression_4="Regressions/read_depth_variability/{your >1000 reads/bp regression}"``` (Only update if the optional step 5 was performed)
    - *Line 115* ```thresh_var_regression_1="Regressions/threshold_read_depth_variability/{your 0-10 reads/bp RMSE threshold}"``` (Only update if the optional step 5 was performed)
    - *Line 116* ```thresh_var_regression_2="Regressions/threshold_read_depth_variability/{your 10-100 reads/bp RMSE threshold}"``` (Only update if the optional step 5 was performed)
    - *Line 117* ```thresh_var_regression_3="Regressions/threshold_read_depth_variability/{your 100-1000 reads/bp RMSE threshold}"``` (Only update if the optional step 5 was performed)
    - *Line 118* ```thresh_var_regression_4="Regressions/threshold_read_depth_variability/{your >1000 reads/bp RMSE threshold}"``` (Only update if the optional step 5 was performed)
#### 6.2 Update “hpc_submission/quant_targets_hpc_run.sh” to represent your system configuration
- In the Job Submission Options section, Update the allocation information 
- In the Job Commands section, update the PATH to the user specific storage locations for the anaconda environments that snakemake will create and use 
#### 6.3 Run Snakefile-quant_targets on all or selected samples
- From the QuantMeta directory:
  ```
  conda activate snakemake
  sbatch hpc_submission/quant_targets_hpc_run.sh
  ```
  Snakemake will automatically submit all of the jobs. Output logs from each individual job will be located in the Logs directory.
- *Optional:* We recommend performing a dry run prior to job submission to the hpc from the QuantMeta directory.
  ```
  conda activate snakemake
  snakemake -prn -s Snakefile-quant_targets
  ```
  Check that all expected rules and number of executions are listed in DAG of jobs snakemake prints out.
- *Optional:* If your jobs fail for any reason, snakemake will need to be unlocked prior to rerunning. This can be completed with:
  ```
  snakemake -prn -s Snakefile-quant_targets --unlock
  ```
- NOTE: If you performed step 5 to create new quantification correction regression and the >1,000 reads/bp regression is not a constant value as was the case in Langenfeld et al. 2022, you will need to run a different quantification script. From the QuantMeta directory:
  ```
  conda activate snakemake
  sbatch hpc_submission/quant_targets_all_linear_regs_hpc_run.sh
  ```
  *Optional:* If your jobs fail for any reason, snakemake will need to be unlocked prior to rerunning. This can be completed with:
  ```
  snakemake -prn -s Snakefile-quant_targets_all_linear_regs --unlock
  ```


