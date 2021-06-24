## GWAS for a Binary trait
#library(statgenGWAS)
#library(genio)
#library(GMMAT)
#
### arguments are genotype file, number of markers, phenotype file, output file
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args)<4) {
#  stop("Missing arguments", call.=FALSE)
#}
#
#print("reading in phenotypes")
#pheno.file <- args[3]
#pheno <- read.table(pheno.file, header = TRUE)
#print(paste(nrow(pheno),"individuals in pheno file",sep=" "))
#print(paste(ncol(pheno),"phenotypes",sep=" "))
#
#
#print("Generating Genomic Relationship Matrix")
#geno.file <- args[1]
#names <- as.vector(pheno$IID)
#bed.length <- as.numeric(args[2])
#bed <- read_bed(geno.file, names_ind = names, m_loci = bed.length)
#GRM <- kinship(t(bed), method = "vanRaden")
#rm(bed)
#print("Done!")
#
#print("Starting GWAS")
#out.file <- args[4]
#GWAS_result <- glmm.wald(fixed = strain ~ . , data = pheno, kins = GRM, id = "IID", family = binomial(link = "logit"), infile = geno.file, snps = ".")
#write.csv(GWAS_result,out.file)
#print("Done!")

#################### Code from GWAS course ##############################
source("scripts/software/gwas.r")
source("scripts/software/emma.r")

library(qqman)
library(statgenGWAS)

## arguments are genotype file, SNP map, phenotype file, phenotype, output file, output plots
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<4) {
  stop("Missing arguments", call.=FALSE)
}


## Reading in genotypes ##
print("Reading in markers")
geno.file <- args[1]
map.file <- args[2]

snpMatrix <- read.table(geno.file, header = TRUE)
SNP_INFO <- read.table(map.file)
names(SNP_INFO) <- c("Chr","SNP","cM","Pos")

X <- as.matrix(snpMatrix[,-c(1:6)])
colnames(X) <- gsub("\\_[A-Z]{1}$","",colnames(X))
rownames(X) <- snpMatrix$IID

print(paste(ncol(X),"markers loaded", sep=" "))
print(paste(nrow(X),"individuals present", sep=" "))


## Reading in phenotypes ##
print("Reading in phenotypes")
pheno.file <- args[3]
phenotype <- args[4]

pheno <- read.table(pheno.file, header = TRUE)
pheno <- pheno[phenotypes$IID %in% snpMatrix$IID,]

Y <- as.matrix(phenotype)
rownames(Y) <- phenotypes$IID

print(paste(nrow(pheno),"individuals in phenotype file",sep=" "))
print(paste(ncol(pheno),"phenotypes present",sep=" "))
print(paste("Running GWAS for phenotype:", phenotype, sep=" "))


print("Generating Genomic Relationship Matrix")
K <- gVanRaden(X)
print("Done!")


print("Running GWAS")
res <- amm_gwas(Y = Y, X = X, K = K, m = 1, use.SNP_INFO = TRUE)
print("Done!")

print("Exporting results and plots")
out.file <- args[5]
out.plots <- args[6]

gwasResults <- res[,c("SNP","Chr","Pos","Pval")]
names(gwasResults) <- c("SNP","CHR","BP","P")
write.csv(gwasResult,out.file)

pdf(out.plots)
manhattan(gwasResults, suggestiveline = FALSE, col = c("red","blue"))
qq(gwasResults$P)
heatmap(K,col=rev(heat.colors(75)))
dev.off()
print("All done!")
