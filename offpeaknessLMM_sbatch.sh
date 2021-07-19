#!/bin/bash

#SBATCH --job-name=Offpeakness_GWAS
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=32G
#SBATCH --error=/mnt/beegfs/alec/GWAS/error/Offpeakness_GWAS.%J.err
#SBATCH --output=/mnt/beegfs/alec/GWAS/error/Offpeakness_GWAS.%J.out
#SBATCH --partition=standard

eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

Rscript --vanilla scripts/BinaryGWAS_LMM.r plink/ros_286_imputed.bed plink/ros_286_imputed.map GWAS/pheno_all.txt off_peakness GWAS/Offpeakness_GWAS.results plots/Offpeakness_GWAS

conda deactivate
