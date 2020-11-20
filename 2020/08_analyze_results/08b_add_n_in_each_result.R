### add n in each result of metal ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_21_traits_unified.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

pool_1_n <- read.table("sample_size_pool_1.txt", head=F, stringsAsFactors=F)
pool_2_n <- read.table("sample_size_pool_2.txt", head=F, stringsAsFactors=F)

# make vectors with n for each trait

n_6d <- pool_1_n[1:7,2] + pool_2_n[1:7,2]
n_42d <- pool_1_n[8:14,2] + pool_2_n[8:14,2] 
n_3m <- pool_1_n[15:21,2] + pool_2_n[15:21,2]  

# add n columns to each trait

for (n in (1:7)) {

	N <- rep (n_6d[n], nrow(metal_6d_clear[[n]]))

	metal_6d_clear[[n]] <- cbind(metal_6d_clear[[n]], N)

}

for (n in (1:7)) {

	N <- rep (n_42d[n], nrow(metal_42d_clear[[n]]))

	metal_42d_clear[[n]] <- cbind(metal_42d_clear[[n]], N)

}

for (n in (1:7)) {

	N <- rep (n_3m[n], nrow(metal_3m_clear[[n]]))

	metal_3m_clear[[n]] <- cbind(metal_3m_clear[[n]], N)

}

# save traits for multiabel

library(data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/for_multiabel")

for (n in (1:7)) {

	name_to_save <- paste("index_", n, "_6d.txt", sep="")

	fwrite(metal_6d_clear[[n]], name_to_save, row.names=F, col.names=T, quote=F)

}

for (n in (1:7)) {

	name_to_save <- paste("index_", n, "_42d.txt", sep="")

	fwrite(metal_42d_clear[[n]], name_to_save, row.names=F, col.names=T, quote=F)

}

for (n in (1:7)) {

	name_to_save <- paste("index_", n, "_3m.txt", sep="")

	fwrite(metal_3m_clear[[n]], name_to_save, row.names=F, col.names=T, quote=F)

}