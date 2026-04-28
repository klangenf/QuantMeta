#!/bin/bash

echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

bam="${3}/${sample}_standards_sorted.bam"


if [ $4 == "Map_Indexes/Langenfeld_2025_standards.fasta" ]; then
    fa_input="${7}/../data/${4}"
else
    fa_input=$4
fi

### Map reads to standards
mkdir ${6}/tmp/${sample}
samtools depth -a -H $bam --reference $fa_input > ${6}/tmp/${sample}/standards_depth.txt
mkdir ${6}/Mapping/${sample}
python3 ${7}/../organize_mapping.py -d ${6}/tmp/${sample}/standards_depth.txt -f $fa_input -o ${6}/Mapping/${sample}/standards_mapping.txt
python3 ${7}/sliding_window.py -m ${6}/Mapping/${sample}/standards_mapping.txt -o ${6}/Mapping/${sample}/standards_sliding_window.txt -n $5