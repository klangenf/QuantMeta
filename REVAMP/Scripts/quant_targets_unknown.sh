#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

# reads_1="${3}/${sample}_R1.fastq.gz"
# reads_2="${3}/${sample}_R2.fastq.gz"

bam="${3}/${sample}_${5}_sorted.bam"

### Map reads to target fasta
# bowtie2 -x Map_Indexes/${5} -q -1 $reads_1 -2 $reads_2 -S tmp/${sample}/${5}.sam
# samtools view -S -b tmp/${sample}/${5}.sam > tmp/${sample}/${5}.bam
# samtools sort tmp/${sample}/${5}.bam -o tmp/${sample}/${5}_sorted.bam
# samtools index tmp/${sample}/${5}_sorted.bam
samtools depth -a -H $bam --reference $4 > tmp/${sample}/${5}_depth.txt
python3 Scripts/organize_mapping.py -d tmp/${sample}/${5}_depth.txt -f $4 -o Mapping/${sample}/${5}_mapping.txt

### Run quantification script
python3 Scripts/quant_targets.py --sample_name $sample --target_name $5 --sample_info $16 --database_lengths Map_Indexes/${5}_lengths.txt --mapping_results Mapping/${sample}/${5}_mapping.txt --window_size $6 --quant_regression Regressions/quantification/${sample}_rel_to_abs.pkl --detect_thresh $7 --quad_reg1 $8 --quad_reg2 $9 --quad_reg3 $10 --quad_reg4 $11 --cutoff_function1 $12 --cutoff_function2 $13 --cutoff_function3 $14 --cutoff_function4 $15