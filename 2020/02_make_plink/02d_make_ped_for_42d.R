### make ped for 42d ###

# load data

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

load("pedAB_96_animals.Rdata")

load("chip_data_pre-processed_96_animals_forward_alleles.Rdata")

# load pheno data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

pheno_42d <- read.csv("pheno_index_42d_2020.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# add unknown samples to pheno data

unknown_samples <- matrix(nrow=21, ncol=ncol(pheno_42d))

unknown_samples[,1] <- names(chip_data_list)[!(names(chip_data_list) %in% pheno_42d[,1])]

unknown_samples[,2:4] <- 0

unknown_samples[,5:11] <- NA 

colnames(unknown_samples) <- colnames(pheno_42d)

pheno_42d <- rbind(pheno_42d, unknown_samples)

# replace NA in pheno file to -9

pheno_42d[is.na(pheno_42d)] <- -9

# make first 6 columns ped files 

ped6 <- list()

pheno_values <- pheno_42d[,5:11]

for (n in (1:7)) {

	x <- data.frame(
  		FID = 0,
  		IID = pheno_42d[,1],
  		FatherID = 0,
  		MotherID = 0,
  		Sex = pheno_42d$Sex,
  		Phenotype = pheno_values[,n])

	ped6[[n]] <- x

	}


to_match <- match(rownames(pedAB), ped6[[1]]$IID)	

# make ped files 

ped_full <- list()

for (n in (1:7)) {

	x <- cbind(ped6[[n]][to_match,], pedAB)
	x <- x[order(x$IID),]

	ped_full [[n]] <- x

	print(n)

	}

# save data

ped_names_to_save <- list()

for (n in (1:7)) {

	x <- paste("ped_42d_", n, "_index.ped", sep="") 
	ped_names_to_save[[n]] <- x
	
	}

library (data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_42d")

for (n in (1:7)) {

	fwrite(ped_full[[n]], ped_names_to_save[[n]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

	}