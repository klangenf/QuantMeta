#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

# reads_1="${3}/${sample}_R1.fastq.gz"
# reads_2="${3}/${sample}_R2.fastq.gz"

bam="${3}/${sample}_standards_sorted.bam"

### Map reads to standards
# bowtie2-build -f $4 Map_Indexes/standards
# bowtie2 -x Map_Indexes/standards -q -1 $reads_1 -2 $reads_2 -S tmp/${sample}/standards.sam
# samtools view -S -b tmp/${sample}/standards.sam > tmp/${sample}/standards.bam
# samtools sort tmp/${sample}/standards.bam -o tmp/${sample}/standards_sorted.bam
# samtools index tmp/${sample}/standards_sorted.bam
mkdir tmp/${sample}
samtools depth -a -H $bam --reference $4 > tmp/${sample}/standards_depth.txt
python3 Scripts/organize_mapping.py -d tmp/${sample}/standards_depth.txt -f $4 -o Mapping/${sample}/standards_mapping.txt
python3 Scripts/sliding_window.py -m Mapping/${sample}/standards_mapping.txt -o Mapping/${sample}/standards_sliding_window.txt -n $5