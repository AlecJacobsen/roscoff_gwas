#!/bin/bash

#SBATCH --job-name=check_missing
#SBATCH --ntasks=6
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=16G
#SBATCH --error=check_missing.%J.err
#SBATCH --output=check_missing.%J.out
#SBATCH --partition=standard

vcftools --vcf ros_286.vcf --missing-indv --out ros_286_missing
