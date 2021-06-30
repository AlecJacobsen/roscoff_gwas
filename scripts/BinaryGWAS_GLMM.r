#################### Binary Trait GWAS with GLMM and Kinship Matrix ##############################

library(qqman)
library(genio)
library(statgenGWAS)
library(GMMAT)
library(data.table)

### arguments are 1) genotype file, 2) SNP map, 3) phenotype file, 4) phenotype, 5) output file, 6) output plots ###
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<6) {
  stop("Missing arguments", call.=FALSE)
}

#### Reading in phenotypes ####
print("Reading in phenotypes")
pheno.file <- args[3] #"GWAS/pheno_strain.txt"
phenotype <- args[4] #"b_strain"

pheno <- read.table(pheno.file, header = TRUE)

Y <- pheno[,c("IID",phenotype)]
rownames(Y) <- pheno$IID

print(paste(nrow(Y),"individuals in phenotype file",sep=" "))
print(paste("Running GWAS for phenotype:", phenotype, sep=" "))


#### Reading in genotypes ####
print("Reading in markers")
geno.file <- args[1] #"plink/ros_286_imputed.bed"
map.file <- args[2] #"plink/ros_286_imputed.map"

SNP_INFO <- read.table(map.file)
names(SNP_INFO) <- c("Chr","SNP","cM","Pos")
SNP_INFO <- SNP_INFO[,c("SNP","Chr","Pos")]

names <- as.vector(pheno$IID)
X <- read_bed(geno.file, names_ind = names, m_loci = nrow(SNP_INFO))

print(paste(ncol(X),"markers loaded", sep=" "))
print(paste(nrow(X),"individuals present", sep=" "))


## Generating Matrix ##
print("Generating Genomic Relationship Matrix")
K <- kinship(t(X), method = "vanRaden")
print("Done!")

## Running the GWAS ##
print("Fitting the GLMM")
model <- glmmkin(b_strain ~ ., data = Y, kins = K, id = "IID", family = binomial(link = "logit"))
print("Running GWAS")
glmm.score(model, infile = args[1], outfile = args[5]) #infile = "plink/ros_286_imputed.bed", outfile = "GWAS/Strain_GWAS_GLMM.txt")
print("Done!")


# plots and saving data
print("Creating figures")

res <- read.csv(args[5],sep = "\t")
gwasResults <- res[,c("CHR","POS","PVAL")]

png(paste(args[6],"manhattan_GMMAT.png",sep="_"))
manhattan(gwasResults, chr = "CHR", bp = "POS", p = "PVAL", suggestiveline = TRUE, col = c("red","blue"),
          genomewideline = FALSE, logp = TRUE, width = 1000, height = 600, res = 100)
dev.off()

## qq-plot
png(paste(args[6],"qqplot_GMMAT.png",sep="_"), width = 600, height = 600)
qq(gwasResults$PVAL)
dev.off()

print("All Done!")
