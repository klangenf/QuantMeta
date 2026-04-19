#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

reads_1="${6}/${sample}_R1.fastq.gz"
reads_2="${6}/${sample}_R2.fastq.gz"

### Map reads to standards
mkdir tmp/${sample}
mkdir Mapping/${sample}_100
mkdir Mapping/${sample}
bowtie2 -x Map_Indexes/standards -q -1 $reads_1 -2 $reads_2 -S tmp/${sample}/standards.sam
samtools view -S -b tmp/${sample}/standards.sam > tmp/${sample}/standards.bam
samtools sort tmp/${sample}/standards.bam -o tmp/${sample}/standards_sorted.bam
samtools index tmp/${sample}/standards_sorted.bam
samtools depth -a -H tmp/${sample}/standards_sorted.bam --reference $3 > tmp/${sample}/standards_depth.txt
python3 Scripts/organize_mapping.py -d tmp/${sample}/standards_depth.txt -f $3 -o Mapping/${sample}_100/standards_mapping.txt

### Downsample reads and repeat mapping to standards
read_count=$(($(zcat $reads_1 | wc -l) / 4))

for downsample in 1 10; do
    out_r1="tmp/${sample}_${downsample}_R1.fastq.gz"
    out_r2="tmp/${sample}_${downsample}_R2.fastq.gz"
    number=$(($read_count/((100/$downsample))))
    seqtk sample -s100 $reads_1 $number > $out_r1
    seqtk sample -s100 $reads_2 $number > $out_r2
    mkdir tmp/${sample}_${downsample}
    mkdir Mapping/${sample}_${downsample}
    bowtie2 -x Map_Indexes/standards -q -1 $out_r1 -2 $out_r2 -S tmp/${sample}_${downsample}/standards.sam
    samtools view -S -b tmp/${sample}_${downsample}/standards.sam > tmp/${sample}_${downsample}/standards.bam
    samtools sort tmp/${sample}_${downsample}/standards.bam -o tmp/${sample}_${downsample}/standards_sorted.bam
    samtools index tmp/${sample}_${downsample}/standards_sorted.bam
    samtools depth -a -H tmp/${sample}_${downsample}/standards_sorted.bam --reference $3 > tmp/${sample}_${downsample}/standards_depth.txt
    python3 Scripts/organize_mapping.py -d tmp/${sample}_${downsample}/standards_depth.txt -f $3 -o Mapping/${sample}_${downsample}/standards_mapping.txt
done

for fail_set in r1 r2 r3 r4 r5; do
    bowtie2 -x Map_Indexes/fail_standards_${fail_set}/fail_standards -q -1 $reads_1 -2 $reads_2 -S tmp/${sample}/fail_standards_${fail_set}.sam
    samtools view -S -b tmp/${sample}/fail_standards_${fail_set}.sam > tmp/${sample}/fail_standards_${fail_set}.bam
    samtools sort tmp/${sample}/fail_standards_${fail_set}.bam -o tmp/${sample}/fail_standards_${fail_set}_sorted.bam
    samtools index tmp/${sample}/fail_standards_${fail_set}_sorted.bam
    samtools depth -a -H tmp/${sample}/fail_standards_${fail_set}_sorted.bam --reference Map_Indexes/fail_standards_${fail_set}/fail_standards_ms.fasta > tmp/${sample}/fail_standards_${fail_set}_depth.txt
    python3 Scripts/organize_mapping.py -d tmp/${sample}/fail_standards_${fail_set}_depth.txt -f Map_Indexes/fail_standards_${fail_set}/fail_standards_ms.fasta -o Mapping/${sample}/fail_standards_${fail_set}_mapping.txt
done

### Map reads to test database
samtools depth -a -H $7/${sample}_${4}_sorted.bam --reference $5 > tmp/${sample}/${4}_depth.txt
python3 Scripts/organize_mapping.py -d tmp/${sample}/${4}_depth.txt -f $5 -o Mapping/${sample}/${4}_mapping.txt

echo done