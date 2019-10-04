### make plink files ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/Rdata")

load("chip_data_pre-processed.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d <- read.csv("phenotypes_for_48_sheeps_6d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")
pheno_42d <- read.csv("phenotypes_for_48_sheeps_42d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")
pheno_3m <- read.csv("phenotypes_for_48_sheeps_3m.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# replace NA in pheno files to -9

pheno_6d[is.na(pheno_6d)] <- -9

# make map file

map <- data.frame(
  CHR = chip_data_list[[1]]$Chr,
  SNP = chip_data_list[[1]]$SNP,
  CM = 0,
  BP = chip_data_list[[1]]$Pos)


# make first 6 columns ped files for each trait for 6d

ped6_1_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,3])
ped6_2_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,4])
ped6_3_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,5])
ped6_4_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,6])
ped6_5_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,7])
ped6_6_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,8])
ped6_7_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,9])
ped6_8_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,10])
ped6_9_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,11])
ped6_10_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,12])
ped6_11_6d <- data.frame(
  FID = 0,
  IID = pheno_6d[,1],
  FatherID = 0,
  MotherID = 0,
  Sex = pheno_6d$sex,
  Phenotype = pheno_6d[,13])

# make A1A2 column for each individual

library(tidyr)

for (n in (1:48)) {

	chip_data_list[[n]] <- unite(data=chip_data_list[[n]], ...=2:3, col=AB, sep=" ")

	}

# make match chip_data_list matrix

match_matrix <- matrix(ncol=2, nrow=48)
match_matrix [,1] <- seq(1:48) 
match_matrix [,2] <- apply(match_matrix, 2, function(x) paste("sample", x, sep="_"))[,1]

chip_data_list[["sample_46"]] <- NULL
match_matrix <- match_matrix[-46,]

to_match <- match(names(chip_data_list), match_matrix[,2])

# make AB matrix for ped files 

pedAB <- matrix(ncol=nrow(chip_data_list[[1]]), nrow=47)
rownames(pedAB) <- names(chip_data_list)

for (n in (c(1:47))) {

	pedAB_n <- t(as.matrix((unlist(chip_data_list[[n]]$AB))))

	pedAB[n,] <- pedAB_n
}

#make ped files

ped6_1_6d <- cbind(ped6_1_6d[to_match,], pedAB_n)
ped6_1_6d <- ped6_1_6d[order(ped6_1_6d$IID),]



