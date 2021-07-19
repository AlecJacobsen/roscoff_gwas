#!/bin/bash

#SBATCH --job-name=lunarCPH
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=60G
#SBATCH --error=/mnt/beegfs/alec/GWAS/error/lunarCHP.%J.err
#SBATCH --output=/mnt/beegfs/alec/GWAS/error/lunarCHP.%J.out
#SBATCH --partition=global

### with PCs added ###


eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

Rscript --vanilla scripts/Circalunar_GWAS.r plink/ros_286_imputed.bed plink/ros_286_imputed.map GWAS/pheno_all.txt GWAS/lunarCHP_PC1_GWAS_results.csv plink/ros_286_imputed.eigenvec

conda deactivate
