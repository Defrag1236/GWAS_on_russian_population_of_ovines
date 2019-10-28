### extract snps from gwas for 2nd measurement of phenotypes ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/results/gemma_results/42d/output")

trait_list <- list()

for (n in (1:11)) {

	file_name <- paste ("42d_", n, "_trait.assoc.txt", sep="")

	trait_list[[n]] <- fread(file_name, head=T, stringsAsFactors=F)

	}

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d_key <- fread("pheno_42d_key.txt", head=F, stringsAsFactors=F)


# extract significant snps from each trait 

trait_list_snps <- list()

for (n in (1:11)) {

	trait_list_snps[[n]] <-  subset(trait_list[[n]], trait_list[[n]]$p_wald<0.05/(33*nrow(trait_list[[n]])), select=1:12)

	}


# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/results/gemma_results")

save(trait_list_snps, file="significant_snps_for_42d_measurement_of_phenotypes.Rdata")