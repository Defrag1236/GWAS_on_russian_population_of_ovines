### extract additional info for mv snps ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/clumping")

snps <- read.table("clumping_for_all_mv_stephens.txt", head=T, stringsAsFactors=F, sep="\t")

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

argali_eaf <- read.table("mv_snps_argali_eaf.frq", head=T, stringsAsFactors=F)
romanovka_eaf <- read.table("mv_snps_romanovka_eaf.frq", head=T, stringsAsFactors=F)
katahdin_eaf <- read.table("mv_snps_katahdin_eaf.frq", head=T, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered")

index_1 <- read.table("p_value_index_1_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_2 <- read.table("p_value_index_2_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_3 <- read.table("p_value_index_3_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_4 <- read.table("p_value_index_4_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_5 <- read.table("p_value_index_5_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_6 <- read.table("p_value_index_6_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_7 <- read.table("p_value_index_7_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
mass <- read.table("p_value_mass_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)

# make snp table for plink 

for_plink <- snp_info[(snp_info$rs %in% snps[,1]),"SNP_name"]

# make big common table

snps <- as.data.frame(snps[order(snps$CHR),])

snps[,9] <- snp_info[match( snps[,1], snp_info$rs),"SNP_name"]

metal_grep_6d <- list()
metal_grep_42d <- list()
metal_grep_3m <- list()
metal_grep_mass <- list()


for (n in (1:7)) {

	grep_6d <- metal_6d_clear[[n]][match(snps[,9], metal_6d_clear[[n]]$MarkerName),]
	colnames(grep_6d) <- paste("index_", n, "_6d_", colnames(grep_6d), sep="")
	metal_grep_6d[[n]] <- grep_6d

	grep_42d <- metal_42d_clear[[n]][match(snps[,9], metal_42d_clear[[n]]$MarkerName),]
	colnames(grep_42d) <- paste("index_", n, "_42d_", colnames(grep_42d), sep="")
	metal_grep_42d[[n]] <- grep_42d

	grep_3m <- metal_3m_clear[[n]][match(snps[,9], metal_3m_clear[[n]]$MarkerName),]
	colnames(grep_3m) <- paste("index_", n, "_3m_", colnames(grep_3m), sep="")
	metal_grep_3m[[n]] <- grep_3m

	print(n)

}

metal_grep_mass[[1]] <-  metal_mass[[1]][match(snps[,9], metal_mass[[1]]$MarkerName),]
colnames(metal_grep_mass[[1]]) <- paste("mass_6d_", colnames(metal_grep_mass[[1]]), sep="")

metal_grep_mass[[2]] <-  metal_mass[[2]][match(snps[,9], metal_mass[[2]]$MarkerName),]
colnames(metal_grep_mass[[2]]) <- paste("mass_42d_", colnames(metal_grep_mass[[2]]), sep="")

metal_grep_mass[[3]] <-  metal_mass[[3]][match(snps[,9], metal_mass[[3]]$MarkerName),]
colnames(metal_grep_mass[[3]]) <- paste("mass_3m_", colnames(metal_grep_mass[[3]]), sep="")



metal_grep_6d <- do.call(cbind, metal_grep_6d)
metal_grep_42d <- do.call(cbind, metal_grep_42d)
metal_grep_3m <- do.call(cbind, metal_grep_3m)
metal_grep_mass <- do.call(cbind, metal_grep_mass)	

general_table <- cbind(snps, metal_grep_6d, metal_grep_42d, metal_grep_3m, metal_grep_mass)

colnames(argali_eaf) <- paste("argali_", colnames(argali_eaf), sep="")
colnames(romanovka_eaf) <- paste("romanovka_", colnames(romanovka_eaf), sep="")
colnames(katahdin_eaf) <- paste("katahdin_", colnames(katahdin_eaf), sep="")

argali_eaf <- rbind(argali_eaf[1:3,], NA, argali_eaf[4:8,])
romanovka_eaf <- rbind(romanovka_eaf[1:3,], NA, romanovka_eaf[4:8,])
katahdin_eaf <- rbind(katahdin_eaf[1:3,], NA, katahdin_eaf[4:8,])

general_table <- cbind(general_table, argali_eaf, romanovka_eaf, katahdin_eaf)


index_1_mv <- index_1[match(snps[,9], index_1[,1]), 2]
index_2_mv <- index_2[match(snps[,9], index_2[,1]), 2]
index_3_mv <- index_3[match(snps[,9], index_3[,1]), 2]
index_4_mv <- index_4[match(snps[,9], index_4[,1]), 2]
index_5_mv <- index_5[match(snps[,9], index_5[,1]), 2]
index_6_mv <- index_6[match(snps[,9], index_6[,1]), 2]
index_7_mv <- index_7[match(snps[,9], index_7[,1]), 2]
mass_mv <- mass[match(snps[,9], mass[,1]), 2]

general_table <- cbind(general_table, index_1_mv=index_1_mv, index_2_mv=index_2_mv, index_3_mv=index_3_mv, index_4_mv=index_4_mv, index_5_mv=index_5_mv,
						index_6_mv=index_6_mv, index_7_mv=index_7_mv, mass_mv=mass_mv)

# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

write.table(for_plink, "mv_snps_for_plink.txt", col.names=F, row.names=F, quote=F)

write.table(general_table, "mv_snps_with_additional_info.txt", col.names=T, row.names=F, quote=F)