
#################### Bianary Trait GWAS with an LMM and Kinship Matrix
#################### ##############################

#### needs to be edited to run without the PCs

setwd("/mnt/beegfs/alec/GWAS")

library(genio)
library(statgenGWAS)
library(coxme)
library(parallel)
library(doParallel)
library(foreach)

### arguments are 1) genotype file, 2) SNP map, 3) phenotype file, 4) output, 5) PCA file
### file ###
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Missing arguments", call. = FALSE)
}

### Functions ###
extract_coxme_table <- function(mod) {
  beta <- mod$coefficients  #$fixed is not needed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2)
  p <- signif(1 - pchisq((beta/se)^2, 1), 2)
  table = data.frame(cbind(beta, se, z, p))
  return(table)
}

### GWAS Function ###
day_cph <- function(i, geno, snps, day, time, id, k, pc) {
  if ((sum(is.na(geno[, i])) == 0) && (sum(is.nan(geno[, i])) == 0) && (sum(is.infinite(geno[,
    i])) == 0)) {
    if (var(geno[, i]) != 0) {
      row <- cbind(snps[i, ], extract_coxme_table(coxme(day ~ geno[,i] + time + pc + (1 | id), varlist = k))[1, ])
      if ((i%%500) == 0) {
        system(paste("echo '", i, " markers tested'"))
      }
      return(row)
    }
  } else {
    system(paste("echo 'Genotype error:", geno[, i], "'"))
  }
}


### Pheno file ###
print("Reading in Phenotypes")
pheno.file <- args[3]  #'pheno_all.txt'
pheno <- read.table(pheno.file, header = TRUE)
print("Done")

print("Reading in Genotypes")
geno.file <- args[1]  #'ros_286_imputed.bed'
map.file <- args[2]  #'ros_286_imputed.map'

SNP_INFO <- read.table(map.file)
names(SNP_INFO) <- c("Chr", "SNP", "cM", "Pos")
SNP_INFO <- SNP_INFO[, c("SNP", "Chr", "Pos")]

names <- as.vector(pheno$IID)
X <- t(read_bed(geno.file, names_ind = names, m_loci = nrow(SNP_INFO)))
print("Done")
print(paste(nrow(SNP_INFO), "SNPs and", nrow(X), "Individuals in Genotype file",
  sep = " "))

print("Reading in PCA eigenvectors")
pca.file <- args[5] #'plink/ros_286_imputed.eigenvec'
pca <- read.table(pca.file, header = FALSE)
print("Done")

print("Generating Kinship Matrix")
K <- kinship(X, method = "vanRaden")
print("Done")

### Cox proportional hazards model ###
Day <- Surv(pheno$lunar_day, pheno$status)
Time <- pheno$time
Id <- rownames(X)
PC <- pca[,3]


print('Configuring Parallel Processes')
n_cores <- 15
clust <- makeCluster(n_cores, type='FORK')
registerDoParallel(clust) #foreach specific
clusterEvalQ(clust,library(coxme))
clusterExport(clust, c('X','SNP_INFO','Day','Time','Id','K','PC','day_cph','extract_coxme_table'))

print("Beginning GWAS")
#result_list <- parLapply(clust, seq_len(nrow(SNP_INFO)), day_cph, geno = X, snps = SNP_INFO, day = Day, time = Time, id = Id, k = K, pc = PC)
result_df <- foreach(i = 1:nrow(SNP_INFO), .combine = rbind, .errorhandling = 'remove') %dopar%
                    day_cph(i = i, geno = X, snps = SNP_INFO, day = Day, time = Time, id = Id, k = K, pc = PC)
stopCluster(clust)
print("Done")

print("Writing out results")
#result_df <- as.data.frame(do.call(rbind, result_list))
colnames(result_df) <- c("SNP", "CHR", "POS", "beta", "se", "z", "p")
write.csv(result_df, args[4], row.names = F)
print(paste(nrow(result_df)," markers tested"))
print("All Finished!")
