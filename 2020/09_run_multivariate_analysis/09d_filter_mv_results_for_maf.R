### filter mv results for maf ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens")

index_1 <- read.table("p_value_index_1_multi_stephens.txt", head=F, stringsAsFactors=F)
index_2 <- read.table("p_value_index_2_multi_stephens.txt", head=F, stringsAsFactors=F)
index_3 <- read.table("p_value_index_3_multi_stephens.txt", head=F, stringsAsFactors=F)
index_4 <- read.table("p_value_index_4_multi_stephens.txt", head=F, stringsAsFactors=F)
index_5 <- read.table("p_value_index_5_multi_stephens.txt", head=F, stringsAsFactors=F)
index_6 <- read.table("p_value_index_6_multi_stephens.txt", head=F, stringsAsFactors=F)
index_7 <- read.table("p_value_index_7_multi_stephens.txt", head=F, stringsAsFactors=F)
mass <- read.table("p_value_mass_multi_stephens.txt", head=F, stringsAsFactors=F)

# count maf for each trait

maf_6d <- list()
maf_42d <- list()
maf_3m <- list()
maf_mass <- list()

for (n in (1:7)) {

	maf_6d[[n]] <- pmin(1-metal_6d_clear[[n]]$Freq1, metal_6d_clear[[n]]$Freq1)
	names(maf_6d[[n]]) <- metal_6d_clear[[n]]$MarkerName

}

for (n in (1:7)) {

	maf_42d[[n]] <- pmin(1-metal_42d_clear[[n]]$Freq1, metal_42d_clear[[1]]$Freq1)
	names(maf_42d[[n]]) <- metal_42d_clear[[n]]$MarkerName

}

for (n in (1:7)) {

	maf_3m[[n]] <- pmin(1-metal_3m_clear[[n]]$Freq1, metal_3m_clear[[1]]$Freq1)
	names(maf_3m[[n]]) <- metal_3m_clear[[n]]$MarkerName

}


for (n in (1:3)) {

	maf_mass[[n]] <- pmin(1-metal_mass[[n]]$Freq1, metal_mass[[1]]$Freq1)
	names(maf_mass[[n]]) <- metal_mass[[n]]$MarkerName

}

maf_6d <- do.call(cbind, maf_6d)
maf_42d <- do.call(cbind, maf_42d)
maf_3m <- do.call(cbind, maf_3m)
maf_mass <- do.call(cbind, maf_mass)

maf_all <- cbind(maf_6d, maf_42d, maf_3m, maf_mass)

# filter to maf >0.05

maf <- apply(maf_all, 1, min)

maf <- as.matrix(maf)

maf_f <- subset(maf, maf[,1]>0.05, select=1)

index_1_f <- index_1[match(rownames(maf_f), index_1[,1]),]
index_2_f <- index_2[match(rownames(maf_f), index_2[,1]),]
index_3_f <- index_3[match(rownames(maf_f), index_3[,1]),]
index_4_f <- index_4[match(rownames(maf_f), index_4[,1]),]
index_5_f <- index_5[match(rownames(maf_f), index_5[,1]),]
index_6_f <- index_6[match(rownames(maf_f), index_6[,1]),]
index_7_f <- index_7[match(rownames(maf_f), index_7[,1]),]
mass_f <- mass[match(rownames(maf_f), mass[,1]),]

colnames(index_1_f) <- c("SNP", "P")
colnames(index_2_f) <- c("SNP", "P")
colnames(index_3_f) <- c("SNP", "P")
colnames(index_4_f) <- c("SNP", "P")
colnames(index_5_f) <- c("SNP", "P")
colnames(index_6_f) <- c("SNP", "P")
colnames(index_7_f) <- c("SNP", "P")
colnames(mass_f) <- c("SNP", "P")

# save results 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis")

write.table(maf, "maf_for_each_snp.txt", col.names=F, row.names=T, quote=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered")

write.table(index_1_f, "p_value_index_1_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(index_2_f, "p_value_index_2_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(index_3_f, "p_value_index_3_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(index_4_f, "p_value_index_4_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(index_5_f, "p_value_index_5_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(index_6_f, "p_value_index_6_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(index_7_f, "p_value_index_7_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)
write.table(mass_f, "p_value_mass_multi_stephens_maf_0.05_filtered.txt", col.names=F, row.names=F, quote=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink")

write.table(index_1_f, "p_value_index_1_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(index_2_f, "p_value_index_2_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(index_3_f, "p_value_index_3_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(index_4_f, "p_value_index_4_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(index_5_f, "p_value_index_5_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(index_6_f, "p_value_index_6_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(index_7_f, "p_value_index_7_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
write.table(mass_f, "p_value_mass_multi_stephens_maf_0.05_filtered.txt", col.names=T, row.names=F, quote=F)
