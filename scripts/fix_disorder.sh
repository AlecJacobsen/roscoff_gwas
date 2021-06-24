#!/bin/bash

#SBATCH --job-name=RosGWAS_main
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=32G
#SBATCH --error=fix_disorder.%J.err
#SBATCH --output=fix_disorder.%J.out
#SBATCH --partition=standard

for file in trimmed/issue_files/*_1P.fastq.gz
do
name=`echo $file | sed 's|_1P.fastq.gz||g'`
repair.sh in1=${name}_1P.fastq.gz in2=${name}_2P.fastq.gz out1=${name}_1P.fixed.fastq.gz out2=${name}_2P.fixed.fastq.gz outs=${name}_singletons.fastq repair
done
