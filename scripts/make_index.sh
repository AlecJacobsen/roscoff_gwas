#!/bin/bash

#SBATCH --job-name=Make_index
#SBATCH --ntasks=6
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=16G
#SBATCH --error=Make_index.%J.err
#SBATCH --output=Make_index.%J.out
#SBATCH --partition=standard

for file in gvcfs/*.vcf
do
gatk --java-options "-Xmx16G" IndexFeatureFile -I $file
done
