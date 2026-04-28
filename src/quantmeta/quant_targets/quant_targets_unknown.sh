#!/bin/bash
echo start
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" $2)
echo $sample

bam="${3}/${sample}_${5}_sorted.bam"

### Map reads to target fasta
samtools depth -a -H $bam --reference $4 > ${17}/tmp/${sample}/${5}_depth.txt
python3 ${18}/../organize_mapping.py -d ${17}/tmp/${sample}/${5}_depth.txt -f $4 -o ${17}/Mapping/${sample}/${5}_mapping.txt

### Run quantification script
python3 ${18}/quant_targets.py --sample_name $sample --target_name $5 --sample_info ${16} --database_lengths ${17}/Map_Indexes/${5}_lengths.txt --mapping_results ${17}/Mapping/${sample}/${5}_mapping.txt --window_size $6 --quant_regression ${17}/Regressions/quantification/${sample}_rel_to_abs.pkl --detect_thresh $7 --quad_reg1 $8 --quad_reg2 $9 --quad_reg3 ${10} --quad_reg4 ${11} --cutoff_function1 ${12} --cutoff_function2 ${13} --cutoff_function3 ${14} --cutoff_function4 ${15} --output_dir ${17}