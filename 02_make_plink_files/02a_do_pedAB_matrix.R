### do pedAB matrix ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/Rdata")

load("chip_data_pre-processed.Rdata")

# make AB column for each individual

library(tidyr)

for (n in (1:48)) {

	chip_data_list[[n]] <- unite(data=chip_data_list[[n]], c("A1","A2"), col=AB, sep=" ")

	}


# make AB matrix for ped files 

pedAB <- matrix(ncol=nrow(chip_data_list[[1]]), nrow=48)
rownames(pedAB) <- names(chip_data_list)

for (n in (c(1:48))) {

	pedAB_n <- t(as.matrix((unlist(chip_data_list[[n]]$AB))))

	pedAB[n,] <- pedAB_n
}

# save pedAB as .Rdata

save(pedAB, file="pedAB.Rdata")