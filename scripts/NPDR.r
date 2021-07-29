setwd("/mnt/beegfs/alec/GWAS")

library(npdr)
library(genio)
library(statgenGWAS)

### arguments are 1) genotype file, 2) SNP map, 3) phenotype file, 4) phenotype, 5) output file ###
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<5) {
  stop("Missing arguments", call.=FALSE)
}

print("Reading in phenotypes")
pheno.file <- args[3] #"GWAS/pheno_strain.txt"
phenotype <- args[4] #"b_strain"

pheno <- read.table(pheno.file, header = TRUE)

print(phenotype)

Y <- pheno[,c("IID",phenotype)]
rownames(Y) <- pheno$IID

print(paste(nrow(Y),"individuals in phenotype file",sep=" "))
print(paste("Running GWAS for phenotype:", phenotype, sep=" "))


#### Reading in genotypes ####
print("Reading in markers")
geno.file <- args[1]
map.file <- args[2]

SNP_INFO <- read_bim(map.file)

names <- as.vector(pheno$IID)
X <- t(read_bed(geno.file, names_loci = SNP_INFO$id, names_ind = names))

X_df <- as.data.frame(X)
X_df$Y <- as.factor(Y)

print(paste(ncol(X),"markers loaded", sep=" "))
print(paste(nrow(X),"individuals present", sep=" "))

#### Generating Matrix ####
print("Generating Genomic Relationship Matrix")
K <- kinship(t(X), method = "vanRaden")

#### Calculating PCs ####
print('Calculating PCs')
K.pca <- prcomp(K, center = TRUE,scale. = TRUE)
Kpcs <- as.matrix(K.pca$x[,1:5])

#### NDPR GWAS ####
print('Starting Nearest-Neighbor Projected-Distance Regression GWAS')
results <- npdr("Y", X_df,
                 regression.type = "binomial", attr.diff.type = "allele-sharing",
                 nbd.method = "multisurf", nbd.metric = "manhattan", msurf.sd.frac = .5,
                 neighbor.sampling = "none", fast.reg = F, dopar.nn = T, dopar.reg = T,
                 padj.method = "fdr", verbose = T
) #covars = Kpcs,
print('Done!')

#### Writing out results ####
print('Writing out results')
write.csv(results, args[5], row.names = F)
