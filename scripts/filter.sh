#!/bin/bash

#SBATCH --job-name=filter
#SBATCH --ntasks=6
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=16G
#SBATCH --error=filter_vcf.%J.err
#SBATCH --output=filter_vcf.%J.out
#SBATCH --partition=standard

vcftools --vcf ros_286.vcf --minQ 19 --recode --recode-INFO-all --out ros_286_filt

vcftools --vcf ros_286.vcf --missing-indv --out ros_286_filt_missing
