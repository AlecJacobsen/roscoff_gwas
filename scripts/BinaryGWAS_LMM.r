#################### Bianary Trait GWAS with an LMM and Kinship Matrix ##############################

library(qqman)
library(genio)
library(statgenGWAS)
library(rrBLUP)
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
SNP_INFO <- cbind(SNP_INFO,X)

print(paste(ncol(X),"markers loaded", sep=" "))
print(paste(nrow(X),"individuals present", sep=" "))


## Generating Matrix ##
print("Generating Genomic Relationship Matrix")
K <- kinship(t(X), method = "vanRaden")
print("Done!")

## Running the GWAS ##
print("Running GWAS")
res <- GWAS(
  pheno = Y,
  geno = SNP_INFO,
  K = K,
  plot = FALSE
)
print("Done!")

# plots and saving data
print("writing out results and figures ...")
names(res)[length(res)] <- phenotype
gwasResults <- res[,c("SNP","Chr","Pos",phenotype)]
names(gwasResults) <- c("SNP","CHR","BP","P")

png(paste(args[6],"manhattan_rrBLUP.png",sep="_"))
manhattan(gwasResults, suggestiveline = TRUE, col = c("red","blue"),
          genomewideline = FALSE, logp = FALSE, width = 800, height = 600, res = 100)
dev.off()

# convert -log(p) back to p-values
p <- 10^((-gwasResults$P))

## rename P to log_p (as it is) and add the column with p-values
names(gwasResults)[4] <- "log_p"
gwasResults$P <- p

fwrite(x = gwasResults, file = args[5])

## qq-plot
png(paste(args[6],"qqplot_rrBLUP.png",sep="_"), width = 600, height = 600)
qq(p)
dev.off()

print("All Done!")
