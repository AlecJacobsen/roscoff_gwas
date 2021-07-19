#!/bin/bash

#SBATCH --job-name=strainGLMM_GWAS
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=32G
#SBATCH --error=/mnt/beegfs/alec/GWAS/error/strainGLMM_GWAS.%J.err
#SBATCH --output=/mnt/beegfs/alec/GWAS/error/strainGLMM_GWAS.%J.out
#SBATCH --partition=standard

eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

Rscript --vanilla scripts/BinaryGWAS_GLMM.r plink/ros_286_imputed.bed plink/ros_286_imputed.map GWAS/pheno_all.txt b.strain GWAS/strainGLMM_GWAS.results plots/strainGLMM_GWAS

conda deactivate
