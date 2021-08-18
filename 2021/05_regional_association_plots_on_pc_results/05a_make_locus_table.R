### make locus table ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/clumping/pc")

new_snps <- read.table("pc_2.txt", head=T, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/r2")

library(data.table)

ld <- fread("r2_122_animals_with_pheno.ld", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/pc/pc_2")

pc <- fread("pc_2_done.csv", head=T, stringsAsFactors=F, data.table=F)

# make chr_pos for locus table 

coord <- snp_info[match(pc$rs_id, snp_info$rs),c(2, 4:6)]
coord <-  coord[!(grepl("X", coord$chromosome)),]
coord <-  coord[!(grepl("99", coord$chromosome)),]

library(plyr)

coord <- arrange(coord, chromosome, position)

chr_pos <- coord[,2:3]
chr_pos<- apply(chr_pos, 2, as.numeric) 

# match coord and pc 

pc <- pc[match(coord$rs, pc$rs_id),]

# make locus table for each new snp associated with pc

for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=pc$rs,
						 "PVAL"=pc$p, "TYPE"=rep("imputed",length(pc[,1])), stringsAsFactors=F)

for (i in 1:nrow(new_snps)) {

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

	write.table(snp, paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/for_regional_association_plot/locus_table_", top, ".txt", sep=""), row.names=F, quote=F)


}