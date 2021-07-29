#!/bin/bash

#SBATCH --job-name=NPDR_pipe
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --time=76:00:00
#SBATCH --mem=40G
#SBATCH --error=/mnt/beegfs/alec/GWAS/error/NPDR_pipe.%J.err
#SBATCH --output=/mnt/beegfs/alec/GWAS/error/NPDR_pipe.%J.out
#SBATCH --partition=standard

## LD pruning in plink -> Filter SNPs to 10,000 based on Lunar CPH results -> run NPDR in R

## example command: sbatch NPDR_pipe.sh plink/ros_286_imputed.ped GWAS/lunarCHP_GWAS_results.csv

eval "$(/data/modules/python/python-anaconda3-2020.02/bin/conda shell.bash hook)"

conda activate clunio-GWAS

cd /mnt/beegfs/alec/GWAS

ped=$1
ped_name=`echo $ped | rev | cut -d'/' -f1 | rev`
ped_dir=`echo $ped | sed "s|${ped_name}||g"`
ped_name=`echo $ped_name | cut -d'.' -f1`


gwas_res=$2
res_pref=`echo $gwas_res | rev | cut -d'/' -f1 | rev`
res_dir=`echo $gwas_res | sed "s|${res_pref}||g"`
res_pref=`echo $res_pref | cut -d'.' -f1`


echo '##### NPDR pipe #####'

if ! test -f ${ped_dir}${ped_name}.prune.in; then
  echo '##### Calculating LD SNPs #####'
  #plink --file ${ped_dir}${ped_name}
  plink --file ${ped_dir}${ped_name} --indep-pairwise 50 10 0.5 --chr-set 3 no-xy no-mt --nonfounders --set-missing-var-ids @:# --out ${ped_dir}${ped_name}
  echo '##### Done #####'
fi

if test -f ${ped_dir}${ped_name}.prune.in; then
  if ! test -f ${ped_dir}${ped_name}.prune_${res_pref}.in; then
    echo '##### Filtering for Significance #####'
    python3 roscoff_gwas/scripts/PruneForSig.py -r $gwas_res -i ${ped_dir}${ped_name}.prune.in -n 10000 -o ${ped_dir}${ped_name}.prune_${res_pref}.in
    echo '##### Done #####'
  fi
fi

if test -f ${ped_dir}${ped_name}.prune_${res_pref}.in; then
  if ! test -f ${ped_dir}${ped_name}_prune_${res_pref}.bed; then
    echo '##### Pruning BED file #####'
    plink --bfile ${ped_dir}${ped_name} --extract ${ped_dir}${ped_name}.prune_${res_pref}.in --nonfounders --set-missing-var-ids @:# --chr-set 3 no-xy no-mt --make-bed --out ${ped_dir}${ped_name}_prune_${res_pref}
    echo '##### Done #####'
  fi
fi

if test -f ${ped_dir}${ped_name}_prune_${res_pref}.bed; then
  if ! test -f ${res_dir}${res_pref}_NPDR.results; then
    echo '##### Running NPDR #####'
    Rscript --vanilla roscoff_gwas/scripts/NPDR.r ${ped_dir}${ped_name}_prune_${res_pref}.bed ${ped_dir}${ped_name}_prune_${res_pref}.bim GWAS/pheno_all.txt b.strain ${res_dir}${res_pref}_NPDR.results
    echo '##### Done #####'
  fi
fi

echo '##### Pipe complete #####'

conda deactivate
