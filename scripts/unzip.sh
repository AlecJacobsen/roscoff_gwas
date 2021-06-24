#!/bin/bash

#SBATCH --job-name=unzip
#SBATCH --ntasks=6
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=24G
#SBATCH --error=unzip.%J.err
#SBATCH --output=unzip.%J.out
#SBATCH --partition=standard

gunzip -c vcf/ros_286.vcf.gz > vcf/ros_286.vcf
