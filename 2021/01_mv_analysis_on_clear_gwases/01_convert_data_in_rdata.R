### convert data in rdata ###


# load data

library(data.table)

metal_6d <- list()

for (n in (1:7)) {

	folder_name <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/index_", n, "_6d", sep="")
	setwd(folder_name)
	gwas_name <- paste("index_", n, "_6d_done.csv", sep="")
	metal_6d[[n]] <- read.csv(gwas_name, head=T, stringsAsFactors=F)

}

metal_42d <- list()

for (n in (1:7)) {

	folder_name <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/index_", n, "_42d", sep="")
	setwd(folder_name)
	gwas_name <- paste("index_", n, "_42d_done.csv", sep="")
	metal_42d[[n]] <- read.csv(gwas_name, head=T, stringsAsFactors=F)

}

metal_3m <- list()

for (n in (1:7)) {

	folder_name <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/index_", n, "_3m", sep="")
	setwd(folder_name)
	gwas_name <- paste("index_", n, "_3m_done.csv", sep="")
	metal_3m[[n]] <- read.csv(gwas_name, head=T, stringsAsFactors=F)

}

metal_mass <- list()

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/mass_6d")
metal_mass[[1]] <- read.csv("mass_6d_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/mass_42d")
metal_mass[[2]] <- read.csv("mass_42d_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/mass_3m")
metal_mass[[3]] <- read.csv("mass_3m_done.csv", head=T, stringsAsFactors=F)

# save results as .Rdata

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/Rdata")

save(metal_6d, metal_42d, metal_3m, metal_mass, file="metal_results_24_traits_from_database.Rdata")