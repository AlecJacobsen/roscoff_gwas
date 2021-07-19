#!/bin/bash

#SBATCH --job-name=DielCHP
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=50G
#SBATCH --error=/mnt/beegfs/alec/GWAS/error/dielCHP.%J.err
#SBATCH --output=/mnt/beegfs/alec/GWAS/error/dielCHP.%J.out
#SBATCH --partition=standard

eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

Rscript --vanilla scripts/CircadianCPH_GWAS.r plink/ros_286_imputed.bed plink/ros_286_imputed.map GWAS/pheno_all.txt GWAS/DielCHP_GWAS_results.csv

conda deactivate
