### extract snps from gwas for 1st measurement of phenotypes ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/results/gemma_results/6d/output")

trait_list <- list()

for (n in (1:11)) {

	file_name <- paste ("6d_", n, "_trait.assoc.txt", sep="")

	trait_list[[n]] <- fread(file_name, head=T, stringsAsFactors=F)

	}

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d_key <- fread("pheno_6d_key.txt", head=F, stringsAsFactors=F)


# extract significant snps from each trait 

trait1_snps <- subset(trait_list[[1]], trait_list[[1]]$p_wald<)
