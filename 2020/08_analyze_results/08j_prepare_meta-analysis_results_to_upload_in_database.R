### preare meta_analysis results to upload in database ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/mv_stephens/filtered")

index_1 <- read.table("p_value_index_4_multi_stephens_maf_0.01_filtered.txt", head=F, stringsAsFactors=F)

# make coord file to match

coord <- snp_info[match(index_1[,1], snp_info$SNP_name),c(2,4:6)]
coord <-  coord[!(grepl("X", coord$chromosome)),]
coord <-  coord[!(grepl("99", coord$chromosome)),]

# filter and save data

## 6 days

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/raw_data/for_database")

for (n in (1:7)) {

	metal_6d_clear[[n]] <- metal_6d_clear[[n]][match(coord$SNP_name, metal_6d_clear[[n]][,1]),]
	metal_6d_clear[[n]][,1] <- coord[match(coord$SNP_name, metal_6d_clear[[n]][,1]),"rs"]
	to_save <- metal_6d_clear[[n]][,c(1:4,6:8,10)]

	name_to_save <- paste("index_", n, "_6d_clear.txt", sep="")

	fwrite(to_save, name_to_save, row.names=F, col.names=T, quote=F)

}

## 42 days

for (n in (1:7)) {

	metal_42d_clear[[n]] <- metal_42d_clear[[n]][match(coord$SNP_name, metal_42d_clear[[n]][,1]),]
	metal_42d_clear[[n]][,1] <- coord[match(coord$SNP_name, metal_42d_clear[[n]][,1]),"rs"]
	to_save <- metal_42d_clear[[n]][,c(1:4,6:8,10)]

	name_to_save <- paste("index_", n, "_42d_clear.txt", sep="")

	fwrite(to_save, name_to_save, row.names=F, col.names=T, quote=F)

}

## 3 months

for (n in (1:7)) {

	metal_3m_clear[[n]] <- metal_3m_clear[[n]][match(coord$SNP_name, metal_3m_clear[[n]][,1]),]
	metal_3m_clear[[n]][,1] <- coord[match(coord$SNP_name, metal_3m_clear[[n]][,1]),"rs"]
	to_save <- metal_3m_clear[[n]][,c(1:4,6:8,10)]

	name_to_save <- paste("index_", n, "_3m_clear.txt", sep="")

	fwrite(to_save, name_to_save, row.names=F, col.names=T, quote=F)

}

## mass 

metal_mass[[1]] <- metal_mass[[1]][match(coord$SNP_name, metal_mass[[1]][,1]),]
metal_mass[[1]][,1] <- coord[match(coord$SNP_name, metal_mass[[1]][,1]),"rs"]

metal_mass[[2]] <- metal_mass[[2]][match(coord$SNP_name, metal_mass[[2]][,1]),]
metal_mass[[2]][,1] <- coord[match(coord$SNP_name, metal_mass[[2]][,1]),"rs"]

metal_mass[[3]] <- metal_mass[[3]][match(coord$SNP_name, metal_mass[[3]][,1]),]
metal_mass[[3]][,1] <- coord[match(coord$SNP_name, metal_mass[[3]][,1]),"rs"]

fwrite(metal_mass[[1]][,c(1:4,6:8,10)], "mass_6d_clear.txt", row.names=F, col.names=T, quote=F)
fwrite(metal_mass[[2]][,c(1:4,6:8,10)], "mass_42d_clear.txt", row.names=F, col.names=T, quote=F)
fwrite(metal_mass[[3]][,c(1:4,6:8,10)], "mass_3m_clear.txt", row.names=F, col.names=T, quote=F)