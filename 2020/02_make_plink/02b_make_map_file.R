### make map file ###

# load data

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

load("chip_data_pre-processed_96_animals_forward_alleles.Rdata")

# do map file

map <- data.frame(
  CHR = chip_data_list[[1]]$Chr,
  SNP = chip_data_list[[1]]$SNP,
  CM = 0,
  BP = chip_data_list[[1]]$Pos)

# save map file for each trait and each measurement

map_names_to_save_6d <- list()

for (n in (1:7)) {

	x <- paste("ped_6d_", n, "_index.map", sep="") 
	map_names_to_save_6d[[n]] <- x
	
	}

map_names_to_save_42d <- list()

for (n in (1:7)) {

	x <- paste("ped_42d_", n, "_index.map", sep="") 
	map_names_to_save_42d[[n]] <- x
	
	}

map_names_to_save_3m <- list()

for (n in (1:7)) {

	x <- paste("ped_3m_", n, "_index.map", sep="") 
	map_names_to_save_3m[[n]] <- x
	
	}


library (data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_6d")

for (n in (1:7)) {

	fwrite(map, map_names_to_save_6d[[n]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

	}

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_42d")

for (n in (1:7)) {

	fwrite(map, map_names_to_save_42d[[n]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

	}

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_3m")

for (n in (1:7)) {

	fwrite(map, map_names_to_save_3m[[n]], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

	}
