#!/bin/bash

cd /mnt/beegfs/alec/GWAS

gwas_res=$1
vcf=$2
pref=`echo $gwas_res | cut -d'.' -f1`
echo $pref
vcf_name=`echo $vcf | rev | cut -d'/' -f1 | rev`
echo $vcf_name
vcf_dir=`echo $vcf | sed "s|${vcf_name}||g"`
echo $vcf_dir

if test -f ${pref}.sig_snps; then
  echo '##### Calculating Significant SNPs #####'
  python3 /mnt/beegfs/alec/GWAS/roscoff_gwas/scripts/extract_sigSNPS.py $gwas_res
fi

if test -f ${vcf_dir}${pref}_sigSNPs.recode.vcf; then
  echo '##### Extracting Significant SNPs from VCF'
  vcftools --positions ${pref}.sig_snps --vcf $vcf --out ${vcf_dir}${pref}_sigSNPs --recode --recode-INFO-all
fi

if test -f ${vcf_dir}${pref}_ann.vcf; then
  echo '##### Running snpEff #####'
  java -Xmx8g -jar /mnt/beegfs/alec/snpEff/snpEff.jar -c /mnt/beegfs/alec/snpEff/snpEff.config -v CLUMA2-0_Chromosomes ${vcf_dir}${pref}_sigSNPs.recode.vcf > ${vcf_dir}${pref}_ann.vcf
fi
