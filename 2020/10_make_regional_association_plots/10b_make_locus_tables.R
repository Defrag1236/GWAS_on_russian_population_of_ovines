### make locus tables for regional association plot ### 

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/clumping")

new_snps <- read.table("clumping_for_all_mv_stephens.txt", head=T, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered")

index_4 <- read.table("p_value_index_4_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_5 <- read.table("p_value_index_5_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_7 <- read.table("p_value_index_7_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
mass <- read.table("p_value_mass_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/r2")

library(data.table)

ld <- fread("r2_pool_1.ld", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

# make chr_pos for locus table 

coord <- snp_info[match(index_4[,1], snp_info$SNP_name),c(2, 4:6)]
coord <-  coord[!(grepl("X", coord$chromosome)),]
coord <-  coord[!(grepl("99", coord$chromosome)),]

library(plyr)

coord <- arrange(coord, chromosome, position)

chr_pos <- coord[,2:3]
chr_pos<- apply(chr_pos, 2, as.numeric) 

# order p_value matrix

index_4 <- index_4[match(coord$SNP_name, index_4[,1]),]
index_5 <- index_5[match(coord$SNP_name, index_5[,1]),]
index_7 <- index_7[match(coord$SNP_name, index_7[,1]),]
mass <- mass[match(coord$SNP_name, mass[,1]),]

# replace snp_name for rs_id in p-value matrix

index_4[,1] <- coord[,1]
index_5[,1] <- coord[,1]
index_7[,1] <- coord[,1]
mass[,1] <- coord[,1]

# make locus table for each new snp associated with index_4

for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=index_4[,1],
						 "PVAL"=index_4[,2], "TYPE"=rep("imputed",length(index_4[,1])), stringsAsFactors=F)

top <- new_snps[7,1]
top_snp_name <- coord[grepl(top, coord$rs), 4]

top_chr_pos <- coord[grepl(top, coord$rs), 2:3]
top_chr_pos <- apply(top_chr_pos, 2, as.numeric)

snp <- subset(for_plot,  for_plot[,2]==top_chr_pos[1], select=1:7)
snp <- subset(snp, snp[,1]<(top_chr_pos[2]+500000), select=1:7)
snp <- subset(snp, snp[,1]>(top_chr_pos[2]-500000), select=1:7)



r2 <- ld[grepl(top_snp_name, ld$SNP_A),]

for (n in (1:nrow(r2))) {

	snp[which(r2$BP_B[n]==snp[,1]),"RSquared"] <- r2[n, "R2"]

}
	
write.table(snp, paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/for_regional_association_plots/locus_table_", top, ".txt", sep=""), row.names=F, quote=F)

# make locus table for each new snp associated with index_5


for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=index_5[,1],
						 "PVAL"=index_5[,2], "TYPE"=rep("imputed",length(index_5[,1])), stringsAsFactors=F)

top <- new_snps[9,1]
top_snp_name <- coord[grepl(top, coord$rs), 4]

top_chr_pos <- coord[grepl(top, coord$rs), 2:3]
top_chr_pos <- apply(top_chr_pos, 2, as.numeric)

snp <- subset(for_plot,  for_plot[,2]==top_chr_pos[1], select=1:7)
snp <- subset(snp, snp[,1]<(top_chr_pos[2]+500000), select=1:7)
snp <- subset(snp, snp[,1]>(top_chr_pos[2]-500000), select=1:7)



r2 <- ld[grepl(top_snp_name, ld$SNP_A),]

for (n in (1:nrow(r2))) {

	snp[which(r2$BP_B[n]==snp[,1]),"RSquared"] <- r2[n, "R2"]

}
	
write.table(snp, paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/for_regional_association_plots/locus_table_", top, ".txt", sep=""), row.names=F, quote=F)

# make locus table for each new snp associated with index_7

for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=index_7[,1],
						 "PVAL"=index_7[,2], "TYPE"=rep("imputed",length(index_7[,1])), stringsAsFactors=F)

for (i in (c(5,6))) {

	top <- new_snps[i,1]
	top_snp_name <- coord[grepl(top, coord$rs), 4]

	top_chr_pos <- coord[grepl(top, coord$rs), 2:3]
	top_chr_pos <- apply(top_chr_pos, 2, as.numeric)

	snp <- subset(for_plot,  for_plot[,2]==top_chr_pos[1], select=1:7)
	snp <- subset(snp, snp[,1]<(top_chr_pos[2]+500000), select=1:7)
	snp <- subset(snp, snp[,1]>(top_chr_pos[2]-500000), select=1:7)



	r2 <- ld[grepl(top_snp_name, ld$SNP_A),]

	for (n in (1:nrow(r2))) {

		snp[which(r2$BP_B[n]==snp[,1]),"RSquared"] <- r2[n, "R2"]

	}

write.table(snp, paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/for_regional_association_plots/locus_table_", top, ".txt", sep=""), row.names=F, quote=F)


}

# make locus table for each new snp associated with mass

for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=mass[,1],
						 "PVAL"=mass[,2], "TYPE"=rep("imputed",length(mass[,1])), stringsAsFactors=F)

for (i in (c(1:4,8))) {

	top <- new_snps[i,1]
	top_snp_name <- coord[grepl(top, coord$rs), 4]

	top_chr_pos <- coord[grepl(top, coord$rs), 2:3]
	top_chr_pos <- apply(top_chr_pos, 2, as.numeric)

	snp <- subset(for_plot,  for_plot[,2]==top_chr_pos[1], select=1:7)
	snp <- subset(snp, snp[,1]<(top_chr_pos[2]+500000), select=1:7)
	snp <- subset(snp, snp[,1]>(top_chr_pos[2]-500000), select=1:7)



	r2 <- ld[grepl(top_snp_name, ld$SNP_A),]

	for (n in (1:nrow(r2))) {

		snp[which(r2$BP_B[n]==snp[,1]),"RSquared"] <- r2[n, "R2"]

	}

write.table(snp, paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/for_regional_association_plots/locus_table_", top, ".txt", sep=""), row.names=F, quote=F)


}
