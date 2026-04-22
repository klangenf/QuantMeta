#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

bam="${3}/${sample}_standards_sorted.bam"

### Map reads to standards
mkdir ${6}/tmp/${sample}
samtools depth -a -H $bam --reference $4 > ${6}/tmp/${sample}/standards_depth.txt
python3 Scripts/organize_mapping.py -d ${6}/tmp/${sample}/standards_depth.txt -f $4 -o ${6}/Mapping/${sample}/standards_mapping.txt
python3 Scripts/sliding_window.py -m ${6}/Mapping/${sample}/standards_mapping.txt -o ${6}/Mapping/${sample}/standards_sliding_window.txt -n $5