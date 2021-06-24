#!/bin/bash

#SBATCH --job-name=mean_depth
#SBATCH --ntasks=6
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=24G
#SBATCH --error=mean_depth.%J.err
#SBATCH --output=mean_depth.%J.out
#SBATCH --partition=standard

vcftools --vcf ros_286.vcf --depth --out ros_286

