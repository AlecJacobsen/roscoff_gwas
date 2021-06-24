#!/bin/bash

#SBATCH --job-name=RosGWAS_main
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=8G
#SBATCH --error=../RosGWAS_main.%J.err
#SBATCH --output=../RosGWAS_main.%J.out
#SBATCH --partition=standard

eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

snakemake --jobs 560 --latency-wait 60  --cluster "sbatch --job-name="{rule}" --ntasks=8 --nodes=1 --time=48:00:00 --mem=32G --error=../error/"{rule}".%J.err --output=../error/"{rule}".%J.out --partition=standard"

conda deactivate
