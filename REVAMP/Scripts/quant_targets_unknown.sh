#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

reads_1="${3}/${sample}_R1.fastq.gz"
reads_2="${3}/${sample}_R2.fastq.gz"

### Map reads to target fasta
bowtie2 -x Map_Indexes/${5} -q -1 $reads_1 -2 $reads_2 -S tmp/${sample}/${5}.sam
samtools view -S -b tmp/${sample}/${5}.sam > tmp/${sample}/${5}.bam
samtools sort tmp/${sample}/${5}.bam -o tmp/${sample}/${5}_sorted.bam
samtools index tmp/${sample}/${5}_sorted.bam
samtools depth -a -H tmp/${sample}/${5}_sorted.bam --reference $4 > tmp/${sample}/${5}_depth.txt
python3 Scripts/organize_mapping.py -d tmp/${sample}/${5}_depth.txt -f $4 -o Mapping/${sample}/${5}_mapping.txt
python3 Scripts/sliding_window.py -m Mapping/${sample}/${5}_mapping.txt -o Mapping/${sample}/${5}_sliding_window.txt
python3 Scripts/quant_targets.py 