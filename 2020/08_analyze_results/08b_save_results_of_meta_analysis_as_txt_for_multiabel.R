### add n in each result of metal ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

pool_1_n <- read.table("sample_size_pool_1.txt", head=F, stringsAsFactors=F)
pool_2_n <- read.table("sample_size_pool_2.txt", head=F, stringsAsFactors=F)

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


fwrite(metal_mass[[1]], "mass_6d.txt", row.names=F, col.names=T, quote=F)
fwrite(metal_mass[[2]], "mass_42d.txt", row.names=F, col.names=T, quote=F)
fwrite(metal_mass[[3]], "mass_3m.txt", row.names=F, col.names=T, quote=F)