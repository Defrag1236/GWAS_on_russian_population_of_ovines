### make ped file ###

# load pedAB matrix

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/Rdata")

load("pedAB.Rdata")

# load pheno data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d <- read.csv("phenotypes_for_48_sheeps_6d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# make empty pheno file

pheno <- matrix (ncol=3, nrow=48)

pheno[,1] <- seq(1:48)
pheno[,2] <- pheno_6d$sex
pheno[,3] <- -9


# make first 6 columns ped files 

ped6 <- data.frame(
  		FID = 0,
  		IID = pheno[,1],
  		FatherID = 0,
  		MotherID = 0,
  		Sex = pheno[,2],
  		Phenotype = pheno_values[,3])

# make match matrix

match_matrix <- matrix(ncol=2, nrow=48)
match_matrix [,1] <- seq(1:48) 
match_matrix [,2] <- apply(match_matrix, 2, function(x) paste("sample", x, sep="_"))[,1]

to_match <- match(rownames(pedAB), match_matrix[,2])


# make ped file 

ped <- cbind(ped6[to_match,], pedAB)

ped <- ped[order(ped$IID),]

# save ped file

library (data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files")

fwrite(ped, "for_all_samples.ped", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
