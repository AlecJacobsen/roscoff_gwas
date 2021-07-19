#Clunio GWAS pipeline
# FastQC/MultiQC -> Trim adapters/Filter out bad reads -> Alignment -> Filter out bad alignments -> sorting -> read groups -> indexing -> haplotype calling -> BQSR -> recalling haplotypes -> joint genotyping ->
import os

SAMPLES = []
directory = '/mnt/beegfs/alec/GWAS/raw-data' #absolute directory to raw-reads
sep = '_'
for file in os.listdir(directory):
    if file.endswith("fastq.gz"):
        name = sep.join(file.split('.')[0].split('_')[0:4])
        SAMPLES.append(name)
SAMPLES = list(set(SAMPLES))
reference = "/mnt/beegfs/alec/reference/CLUMA2-0_Chromosomes.fa"

workdir: "/mnt/beegfs/alec/GWAS"

rule Target:
    input:
        "multiqc_report.html",
        expand("unmapped/{sample}.unmapped.sam",sample = SAMPLES),
        "vcf/ros_286.vcf.gz",
        "plots/ros_286.depth.pdf",
        "plots/ros_286.missing_indv.pdf",
        "plots/ros_286.HWE.pdf",
        "plots/ros_286.Het.pdf",
        "plots/PCA.pdf",
        "GWAS/pheno_all.txt",
        "GWAS/pheno_strain.txt"

        #"output/ros_286_GRM.log.txt",
        #"output/ros_286_GRM.cXX.txt"


rule Trimmomatic:
    input:
        "raw-data/{sample}_R1.fastq.gz",
        "raw-data/{sample}_R2.fastq.gz"
    output:
        temp("trimmed/{sample}_1P.fastq.gz"),
        temp("trimmed/{sample}_1U.fastq.gz"),
        temp("trimmed/{sample}_2P.fastq.gz"),
        temp("trimmed/{sample}_2U.fastq.gz")
    params:
        adapters = '/mnt/beegfs/alec/reference/TruSeq3-PE-2.fa'
    shell:
        "java -jar /data/biosoftware/Trimmomatic/Trimmomatic/trimmomatic-0.38.jar PE {input} {output} ILLUMINACLIP:{params.adapters}:2:30:10:8:true LEADING:20 TRAILING:20 MINLEN:75"

rule FastQC:
    input:
        "trimmed/{sample}_1P.fastq.gz",
        "trimmed/{sample}_2P.fastq.gz"
    output:
        "read-qc/{sample}_1P_fastqc.html",
        "read-qc/{sample}_2P_fastqc.html"
    shell:
        "fastqc -t 4 -o read-qc {input}"

rule MultiQC:
    input:
        expand(["read-qc/{sample}_R1_fastqc.html","read-qc/{sample}_R2_fastqc.html"], sample = SAMPLES)
    output:
        "multiqc_report.html"
    shell:
        "multiqc read-qc/"

rule Mapping:
    input:
        "trimmed/{sample}_1P.fastq.gz",
        "trimmed/{sample}_2P.fastq.gz"
    output:
        temp("mapped/{sample}.sam")
    params:
        ref = reference
    shell:
        "bwa mem -t 16 {params.ref} {input} > {output}"

rule Filter_Sort:
    input:
        "mapped/{sample}.sam"
    output:
        temp("mapped/{sample}.bam")
    params:
        ref = "/mnt/beegfs/alec/reference/CLUMA2-0_Chromosomes.fa.fai"
    shell:
        "samtools view -q 20 -bS -t {params.ref} {input} | samtools sort - -o {output}"

rule Extract_unmapped:
    input:
        "mapped/{sample}.sam"
    output:
        "unmapped/{sample}.unmapped.sam"
    shell:
        "samtools view -f 4 {input} > {output}"

rule ReadGroups:
    input:
        "mapped/{sample}.bam"
    output:
        temp("mapped/{sample}.RG.bam")
    shell:
        "java -jar $PICARD AddOrReplaceReadGroups I={input} O={output} LB=whatever PL=illumina PU=whatever SM={wildcards.sample}"

rule Index:
    input:
        "mapped/{sample}.RG.bam"
    output:
        temp("mapped/{sample}.RG.bam.bai")
    shell:
        "samtools index {input} {output}"

rule HaplotypeCalling:
    input:
        bam = "mapped/{sample}.RG.bam",
        index = "mapped/{sample}.RG.bam.bai"
    output:
        gvcf = temp("gvcfs/{sample}.g.vcf"),
        index = temp("gvcfs/{sample}.g.vcf.idx")
    params:
        ref = reference
    shell:
        'gatk --java-options "-Xmx16G" HaplotypeCaller -ERC GVCF -R {params.ref} -I {input.bam} -O {output.gvcf}'

rule RecalBases:
    input:
        sites = "gvcfs/{sample}.g.vcf",
        bam = "mapped/{sample}.RG.bam",
        index = "gvcfs/{sample}.g.vcf.idx"
    output:
        temp("recalibrated_bams/{sample}.recal_data.table")
    params:
        ref = reference
    shell:
        'gatk --java-options "-Xmx16G" BaseRecalibrator -R {params.ref} -known-sites {input.sites} -I {input.bam} -O {output}'

rule ApplyBQSR:
    input:
        bam = "mapped/{sample}.RG.bam",
        bqsr = "recalibrated_bams/{sample}.recal_data.table",
        index = "gvcfs/{sample}.g.vcf.idx"
    output:
        "recalibrated_bams/{sample}.recalibrated.bam"
    params:
        ref = reference
    shell:
        'gatk --java-options "-Xmx16G" ApplyBQSR -R {params.ref} -I {input.bam} --bqsr-recal-file {input.bqsr} -O {output}'

rule BootstrapHaplotypes:
    input:
        "recalibrated_bams/{sample}.recalibrated.bam"
    output:
        "recalibrated_gvcfs/{sample}.recalibrated.g.vcf"
    params:
        ref = reference
    shell:
        'gatk --java-options "-Xmx16G" HaplotypeCaller -ERC GVCF -R {params.ref} -I {input} -O {output}'

rule gvcfList:
    input:
        expand("recalibrated_gvcfs/{sample}.recalibrated.g.vcf", sample = SAMPLES)
    output:
        temp("gvcfs/gvcf.list")
    shell:
        "ls {input} > {output}"

rule CombineGVCFs:
    input:
        "gvcfs/gvcf.list"
    output:
        temp("gvcfs/ros_286.g.vcf")
    params:
        ref = reference
    shell:
        'gatk --java-options "-Xmx32G" CombineGVCFs -R {params.ref} -V {input} -O {output}'

rule JointGenotype:
    input:
        "gvcfs/ros_286.g.vcf"
    output:
        "vcf/ros_286.vcf"
    params:
        ref = reference
    shell:
        'gatk --java-options "-Xmx32G" GenotypeGVCFs -R {params.ref} -V {input} -O {output}'

rule Compress:
    input:
        "vcf/ros_286.vcf"
    output:
        protected("vcf/ros_286.vcf.gz")
    shell:
        "bgzip -c {input} > {output}"

rule Filter_vcf:
    input:
        "vcf/ros_286.vcf"
    output:
        temp("vcf/ros_286.q20mm005maf005.recode.vcf")
    shell:
        "vcftools --vcf {input} --minQ 20 --maf 0.05 --max-missing 0.05 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out vcf/ros_286.q20mm005maf005"

rule Remove_spanning_deletions:
    input:
        "vcf/ros_286.q20mm005maf005.recode.vcf"
    output:
        temp("vcf/spanning_deletions.txt"),
        temp("vcf/ros_286.q20mm005maf005_nospan.recode.vcf")
    shell:
        """
        awk '$5 ~ /*/' {input} | cut -f1,2 > vcf/spanning_deletions.txt
        vcftools --vcf {input} --exclude-positions vcf/spanning_deletions.txt --recode --recode-INFO-all --out vcf/ros_286.q20mm005maf005_nospan
        """

rule GetDepth:
    input:
        "vcf/ros_286.q20mm005maf005_nospan.recode.vcf"
    output:
        "vcf/ros_286.idepth"
    shell:
        "vcftools --vcf {input} --depth --out vcf/ros_286"

rule Depth_plot:
    input:
        idepth_file = "vcf/ros_286.idepth",
        samples_file = "raw-data/SampleNames"
    output:
        "plots/ros_286.depth.pdf"
    shell:
        "python3 roscoff_gwas/scripts/depthPlot.py -i {input.idepth_file} -s {input.samples_file} -o {output}"

rule Get_Missing_indiv:
    input:
        "vcf/ros_286.q20mm005maf005_nospan.recode.vcf"
    output:
        "vcf/ros_286.imiss"
    shell:
        "vcftools --vcf {input} --missing-indv --out vcf/ros_286"

rule Missing_indiv_plot:
    input:
        imiss = "vcf/ros_286.imiss",
        samples_file = "raw-data/SampleNames"
    output:
        'plots/ros_286.missing_indv.pdf'
    shell:
        "python3 roscoff_gwas/scripts/missingPlot.py -i {input.imiss} -s {input.samples_file} -o {output}"

rule Get_het:
    input:
        "vcf/ros_286.q20mm005maf005_nospan.recode.vcf"
    output:
        "vcf/ros_286.het"
    shell:
        "vcftools --vcf {input} --het --out vcf/ros_286"

rule Het_plot:
    input:
        het = "vcf/ros_286.het",
        samples_file = "raw-data/SampleNames"
    output:
        'plots/ros_286.Het.pdf'
    shell:
        "python3 roscoff_gwas/scripts/HetPlot.py -i {input.het} -s {input.samples_file} -o {output}"

rule Get_hardy:
    input:
        "vcf/ros_286.q20mm005maf005_nospan.recode.vcf"
    output:
        "vcf/ros_286.hwe"
    shell:
        "vcftools --vcf {input} --hardy --out vcf/ros_286"

rule Hardy_plot:
    input:
        "vcf/ros_286.hwe"
    output:
        'plots/ros_286.HWE.pdf'
    shell:
        "python3 roscoff_gwas/scripts/HwePlot.py -i {input} -o {output}"

rule Filter_Poor_Indivs:
    input:
        vcf = "vcf/ros_286.q20mm005maf005.recode.vcf",
        imiss = "vcf/ros_286.imiss"
    output:
        temp("vcf/lowDP.indv"),
        temp("vcf/ros_286_filt.recode.vcf")
    shell:
        """
        awk '$5 > 0.07' vcf/ros_286.imiss | awk 'NR!=1' | cut -f1 > vcf/lowDP.indv
        vcftools --vcf {input.vcf} --remove vcf/lowDP.indv --recode --recode-INFO-all --out vcf/ros_286_filt
        """

rule Fix_lost_gts:
    input:
        "vcf/ros_286_filt.recode.vcf"
    output:
        temp("vcf/ros_286_filt_fixed.vcf")
    shell:
        'cat {input} | perl -pe "s/\s\.:/\t.\/.:/g" > {output}'

rule Impute_gt:
    input:
        "vcf/ros_286_filt_fixed.vcf"
    output:
        "vcf/ros_286.imputed.vcf.gz"
    shell:
        "java -Xmx32g -jar $BEAGLE gt={input} out=vcf/ros_286.imputed"

rule Convert_vcf_to_ped:
    input:
        "vcf/ros_286.imputed.vcf.gz"
    output:
        "plink/ros_286_imputed.ped",
        "plink/ros_286_imputed.map"
    shell:
        "plink --vcf vcf/ros_286.imputed.vcf.gz --const-fid --recode --out plink/ros_286_imputed"

rule PCA:
    input:
        "plink/ros_286_imputed.ped",
        "plink/ros_286_imputed.map"
    output:
        "plink/ros_286_imputed.eigenvec",
        "plink/ros_286_imputed.eigenval"
    shell:
        "plink --file plink/ros_286_imputed --nonfounders --pca var-wts --chr-set 3 no-xy no-mt --out plink/ros_286_imputed"

rule PCA_plot:
    input:
        vec = "plink/ros_286_imputed.eigenvec",
        val = "plink/ros_286_imputed.eigenval"
    output:
        "plots/PCA.pdf"
    shell:
        "python3 roscoff_gwas/scripts/plotPCA.py --vec {input.vec} --val {input.val} --out {output}"

rule Get_phenotypes:
    input:
        "plink/ros_286_imputed.ped"
    output:
        samples = temp('GWAS/samples.txt'),
        pheno_file = "GWAS/pheno_all.txt",
        strain_file = "GWAS/pheno_strain.txt"
    shell:
        """
        cut -d' ' -f2 {input} > {output.samples}
        python3 roscoff_gwas/scripts/GetPheno.py -i {output.samples} -o GWAS/pheno
        """

rule Make_bed:
    input:
        "plink/ros_286_imputed.ped",
        "plink/ros_286_imputed.map",
        "GWAS/pheno_strain.txt"
    output:
        "plink/ros_286_imputed.bed",
        "plink/ros_286_imputed.bim"
    shell:
        "plink --file plink/ros_286_imputed --out plink/ros_286_imputed --make-bed --allow-no-sex --make-pheno GWAS/pheno_strain.txt FM --assoc"

rule Make_raw:
    input:
        "plink/ros_286_imputed.ped",
        "plink/ros_286_imputed.map"
    output:
        "plink/ros_286_imputed.raw"
    shell:
        "plink --file plink/ros_286_imputed --recode A --out plink/ros_286_imputed"

#rule Relatedness_matrix:
#    input:
#        "plink/ros_286_imputed.bed",
#        "plink/ros_286_imputed.bim"
#    output:
#        "output/ros_286_GRM.log.txt",
#        "output/ros_286_GRM.cXX.txt"
#    shell:
#        "gemma -bfile plink/ros_286_imputed -gk 1 -o ros_286_GRM"

rule Strain_GWAS:
    input:
        geno = "plink/ros_286_imputed.bed",
        pheno = "GWAS/pheno_strain.txt",
        map = "plink/ros_286_imputed.map"
    output:
        results = "GWAS/Strain_GWAS.results"
    shell:
        "Rscript --vanilla roscoff_gwas/scripts/BinaryGWAS_LMM.r {input.geno} {input.map} {input.pheno} b_strain {output.results} plots/Strain_GWAS"

#rule Strain_GWAS_GLMM:
#    input:
#        geno = "plink/ros_286_imputed.bed",
#        pheno = "GWAS/pheno_strain.txt",
#        map = "plink/ros_286_imputed.map"
#    output:
#        results = "GWAS/Strain_GWAS_GLMM.txt"
#    shell:
#        "Rscript --vanilla roscoff_gwas/scripts/BinaryGWAS_LMM.r {input.geno} {input.map} {input.pheno} b_strain {output.results} plots/Strain_GWAS"


## use coxmeg for time to event GWAS
