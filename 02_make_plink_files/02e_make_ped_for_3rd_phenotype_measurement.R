### make ped for 1st phenotype measurement ###

# load pedAB matrix

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/Rdata")

load("pedAB.Rdata")

# load pheno data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_3m <- read.csv("phenotypes_for_48_sheeps_3m.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# replace NA in pheno file to -9

pheno_3m[is.na(pheno_3m)] <- -9

# make first 6 columns ped files 

ped6 <- list()

pheno_values <- pheno_3m[,3:13]

for (n in (1:11)) {

	x <- data.frame(
  		FID = 0,
  		IID = pheno_3m[,1],
  		FatherID = 0,
  		MotherID = 0,
  		Sex = pheno_3m$sex,
  		Phenotype = pheno_values[,n])

	ped6[[n]] <- x

	}

# make match matrix

match_matrix <- matrix(ncol=2, nrow=48)
match_matrix [,1] <- seq(1:48) 
match_matrix [,2] <- apply(match_matrix, 2, function(x) paste("sample", x, sep="_"))[,1]


to_match <- match(rownames(pedAB), match_matrix[,2])

# make ped files 

ped_full <- list()

for (n in (1:11)) {

	x <- cbind(ped6[[n]][to_match,], pedAB)
	x <- x[order(x$IID),]

	ped_full [[n]] <- x


	}


# save ped files

ped_names_to_save <- list()

for (n in (1:11)) {

	x <- paste("ped_3m_", n, "_trait.ped", sep="") 
	ped_names_to_save[[n]] <- x
	
	}

library (data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files/ped_3m")

for (n in (1:11)) {

	fwrite(ped_full[[n]], ped_names_to_save[[n]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

	}