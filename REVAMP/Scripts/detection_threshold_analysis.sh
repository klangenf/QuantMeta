#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

reads_1="${6}/${sample}_R1.fastq.gz"
reads_2="${6}/${sample}_R2.fastq.gz"

### Map reads to standards
mkdir ${8}/tmp/${sample}
mkdir ${8}/Mapping/${sample}_100
mkdir ${8}/Mapping/${sample}
bowtie2 -x ${8}/Map_Indexes/standards -q -1 $reads_1 -2 $reads_2 -S ${8}/tmp/${sample}/standards.sam
samtools view -S -b ${8}/tmp/${sample}/standards.sam > ${8}/tmp/${sample}/standards.bam
samtools sort ${8}/tmp/${sample}/standards.bam -o ${8}/tmp/${sample}/standards_sorted.bam
samtools index ${8}/tmp/${sample}/standards_sorted.bam
samtools depth -a -H ${8}/tmp/${sample}/standards_sorted.bam --reference $3 > ${8}/tmp/${sample}/standards_depth.txt
python3 Scripts/organize_mapping.py -d ${8}/tmp/${sample}/standards_depth.txt -f $3 -o ${8}/Mapping/${sample}_100/standards_mapping.txt

### Downsample reads and repeat mapping to standards
read_count=$(($(zcat $reads_1 | wc -l) / 4))

for downsample in 1 10; do
    out_r1="${8}/tmp/${sample}_${downsample}_R1.fastq.gz"
    out_r2="${8}/tmp/${sample}_${downsample}_R2.fastq.gz"
    number=$(($read_count/((100/$downsample))))
    seqtk sample -s100 $reads_1 $number > $out_r1
    seqtk sample -s100 $reads_2 $number > $out_r2
    mkdir ${8}/tmp/${sample}_${downsample}
    mkdir ${8}/Mapping/${sample}_${downsample}
    bowtie2 -x ${8}/Map_Indexes/standards -q -1 $out_r1 -2 $out_r2 -S ${8}/tmp/${sample}_${downsample}/standards.sam
    samtools view -S -b ${8}/tmp/${sample}_${downsample}/standards.sam > ${8}/tmp/${sample}_${downsample}/standards.bam
    samtools sort ${8}/tmp/${sample}_${downsample}/standards.bam -o ${8}/tmp/${sample}_${downsample}/standards_sorted.bam
    samtools index ${8}/tmp/${sample}_${downsample}/standards_sorted.bam
    samtools depth -a -H ${8}/tmp/${sample}_${downsample}/standards_sorted.bam --reference $3 > ${8}/tmp/${sample}_${downsample}/standards_depth.txt
    python3 Scripts/organize_mapping.py -d ${8}/tmp/${sample}_${downsample}/standards_depth.txt -f $3 -o ${8}/Mapping/${sample}_${downsample}/standards_mapping.txt
done

for fail_set in r1 r2 r3 r4 r5; do
    bowtie2 -x ${8}/Map_Indexes/fail_standards_${fail_set}/fail_standards -q -1 $reads_1 -2 $reads_2 -S ${8}/tmp/${sample}/fail_standards_${fail_set}.sam
    samtools view -S -b ${8}/tmp/${sample}/fail_standards_${fail_set}.sam > ${8}/tmp/${sample}/fail_standards_${fail_set}.bam
    samtools sort ${8}/tmp/${sample}/fail_standards_${fail_set}.bam -o ${8}/tmp/${sample}/fail_standards_${fail_set}_sorted.bam
    samtools index ${8}/tmp/${sample}/fail_standards_${fail_set}_sorted.bam
    samtools depth -a -H ${8}/tmp/${sample}/fail_standards_${fail_set}_sorted.bam --reference ${8}/Map_Indexes/fail_standards_${fail_set}/fail_standards_ms.fasta > ${8}/tmp/${sample}/fail_standards_${fail_set}_depth.txt
    python3 Scripts/organize_mapping.py -d ${8}/tmp/${sample}/fail_standards_${fail_set}_depth.txt -f ${8}/Map_Indexes/fail_standards_${fail_set}/fail_standards_ms.fasta -o ${8}/Mapping/${sample}/fail_standards_${fail_set}_mapping.txt
done

### Map reads to test database
samtools depth -a -H ${7}/${sample}_${4}_sorted.bam --reference $5 > ${8}/tmp/${sample}/${4}_depth.txt
python3 Scripts/organize_mapping.py -d ${8}/tmp/${sample}/${4}_depth.txt -f $5 -o ${8}/Mapping/${sample}/${4}_mapping.txt

echo done