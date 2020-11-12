### make ped for 2nd batch ###

# load pedAB matrix

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/Rdata")

load("pedAB.Rdata")

# load pheno for sex and id info

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d <- read.csv("phenotypes_for_48_sheeps_6d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# make first 6 columns ped files 

ped6 <-  data.frame(
  		FID = 0,
  		IID = pheno_6d[,1],
  		FatherID = 0,
  		MotherID = 0,
  		Sex = pheno_6d$sex,
  		Phenotype = -9)

# make match matrix

match_matrix <- matrix(ncol=2, nrow=48)
match_matrix [,1] <- seq(1:48) 
match_matrix [,2] <- apply(match_matrix, 2, function(x) paste("sample", x, sep="_"))[,1]

to_match <- match(rownames(pedAB), match_matrix[,2])

# make ped 

ped_full <- cbind(ped6[to_match,], pedAB)

# save ped 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files")

library(data.table)

fwrite(ped_full, "2nd_batch_48_animals_2020_raw.ped", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")