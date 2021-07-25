#!/bin/bash

#SBATCH --job-name=SNPS2Genes
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=32G
#SBATCH --error=/mnt/beegfs/alec/GWAS/error/SNPS2Genes.%J.err
#SBATCH --output=/mnt/beegfs/alec/GWAS/error/SNPS2Genes.%J.out
#SBATCH --partition=standard

## sample command: sbatch SNPsToGenes.sh GWAS/strainGLMM.results vcf/ros_286.imputed.vcf.gz

eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

cd /mnt/beegfs/alec/GWAS

gwas_res=$1
vcf=$2
pref=`echo $gwas_res | rev | cut -d'/' -f1 | rev`
res_dir=`echo $gwas_res | sed "s|${pref}||g"`
pref=`echo $pref | cut -d'.' -f1`
vcf_name=`echo $vcf | rev | cut -d'/' -f1 | rev`
vcf_dir=`echo $vcf | sed "s|${vcf_name}||g"`

echo '##### Gene Set Enrichment Pipeline for the Roscoff Strains GWAS #####'

if ! test -f ${res_dir}${pref}.sig_snps; then
  echo '##### Calculating Significant SNPs #####'
  python3 /mnt/beegfs/alec/GWAS/roscoff_gwas/scripts/extract_sigSNPS.py -i $gwas_res
  echo '##### Done #####'
fi

if test -f ${res_dir}${pref}.sig_snps; then
  if ! test -f ${vcf_dir}${pref}_sigSNPs.recode.vcf; then
    echo '##### Extracting Significant SNPs from VCF'
    vcftools --positions ${res_dir}${pref}.sig_snps --gzvcf $vcf --out ${vcf_dir}${pref}_sigSNPs --recode --recode-INFO-all
    echo '##### Done #####'
  fi
fi

if test -f ${vcf_dir}${pref}_sigSNPs.recode.vcf; then
  if ! test -f ${vcf_dir}${pref}_ann.vcf; then
    echo '##### Running snpEff #####'
    java -Xmx8g -jar /mnt/beegfs/alec/snpEff/snpEff.jar -c /mnt/beegfs/alec/snpEff/snpEff.config -v CLUMA2-0_Chromosomes ${vcf_dir}${pref}_sigSNPs.recode.vcf > ${vcf_dir}${pref}_ann.vcf
    echo '##### Done #####'
  fi
fi

if test -f ${vcf_dir}${pref}_ann.vcf; then
  if ! test -f ${res_dir}${pref}_genes.txt; then
    echo '##### Retrieving Genes and GO terms #####'
    python3 /mnt/beegfs/alec/GWAS/roscoff_gwas/scripts/GetGenes.py -v ${vcf_dir}${pref}_ann.vcf -g /mnt/beegfs/alec/reference/Cluma2.0_29102020_one2one_non_elec_GO_annotations_e-10_cov60.txt -o ${res_dir}${pref}_genes.txt -n 5
    echo '##### Done #####'
  fi
fi

echo '##### All Finished! #####'

conda deactivate
