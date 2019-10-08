### pre_processing_chip_data ###

# load data
library(data.table)

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data")

chip_data <- fread ("F1_OAM_RMNV.txt", head=T, stringsAsFactors=F)

# parse the full table for separate tables of samples

chip_data_list <- split(chip_data, rep(1:48, each=606006))

# rename matrix in list  

names_of_samples <- c(1:48)

for (n in (1:48)) {

	number_of_sample <- as.numeric(unlist(strsplit(as.matrix(chip_data_list[[n]][1,1]), "No"))[2])

	names_of_samples[n] <- paste ("sample", number_of_sample, sep="_")

	}

names(chip_data_list) <- names_of_samples

# delete unused info

for ( n in (1:48)) {

	chip_data_list [[n]] <- chip_data_list[[n]] [,c("SNP.Name", "Allele1...Top", "Allele2...Top", "Chr", "Position")] 
	colnames(chip_data_list[[n]]) <- c("SNP", "A1", "A2", "Chr", "Pos")
	
	}

# save data as .Rdata

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/Rdata")

save(chip_data_list, file="chip_data_pre-processed.Rdata")