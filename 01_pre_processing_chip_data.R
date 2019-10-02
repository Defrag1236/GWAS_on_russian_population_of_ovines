### pre_processing_chip_data ###

# load data
library(data.table)

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data")

chip_data <- fread ("F1_OAM_RMNV.txt", head=T, stringsAsFactors=F)

# parse the full table for separate tables of samples

chip_data_list <- list()

for (n in (1:48)) {

	number <- paste("No", n, sep="")
	sample <- chip_data [grepl(pattern=number, x=chip_data[,1]),]
	chip_data_list[[n]] <- sample

	}

