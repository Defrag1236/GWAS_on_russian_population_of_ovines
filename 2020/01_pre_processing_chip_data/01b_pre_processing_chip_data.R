### pre_processing_chip_data ###

#load data

library(data.table)

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/raw_data")

chip_data <- fread ("raw_data.txt", head=T, stringsAsFactors=F)

# rename colnames 

colnames(chip_data) <- c("SNP", "Sample_ID", "Allele1_top", "Allele2_top", "Allele1_AB", "Allele2_AB", "GC_score", "GT_score", "Chr", "Position", "Sample_index")

# parse the full table for separate tables of samples

chip_data_list <- split(chip_data, rep(1:96, each=606006))

# rename matrix in list  

names_of_samples <- c(1:96)

for (n in (1:96)) {

	names_of_samples[n] <- chip_data_list[[n]][1,2]

	}

names(chip_data_list) <- names_of_samples

# delete unused info

for ( n in (1:96)) {

	chip_data_list [[n]] <- chip_data_list[[n]] [,c("SNP", "Allele1_AB", "Allele2_AB", "Chr", "Position")] 
	colnames(chip_data_list[[n]]) <- c("SNP", "A1", "A2", "Chr", "Pos")
	
	}


# save data as .Rdata

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

save(chip_data_list, file="chip_data_pre-processed_96_animals_AB_alleles.Rdata")