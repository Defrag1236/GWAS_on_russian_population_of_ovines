### compare r^2 ### 

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/clumping/pc")

new_snps <- read.table("pc_2.txt", head=T, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/r2")

library(data.table)

ld <- fread("r2_122_animals_with_pheno.ld", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/pc/pc_2/")

pc_2 <- fread("pc_2_done.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/rdata")

load("supplement_from_Bolormaa.Rdata")

# make chr_pos for locus table 

coord <- snp_info[match(pc_2$rs_id, snp_info$rs),c(2, 4:6)]
coord <-  coord[!(grepl("X", coord$chromosome)),]
coord <-  coord[!(grepl("99", coord$chromosome)),]
rownames(coord) <- paste(coord$chromosome, coord$position, sep=":")

library(plyr)

coord <- arrange(coord, chromosome, position)

chr_pos <- coord[,2:3]
chr_pos<- apply(chr_pos, 2, as.numeric) 

# make gwas data in order with coord

pc_2 <- pc_2[match(coord$rs, pc_2$rs_id),]

# make r^2 matrix for top snps in window +- 100 kb

for_plot <- data.frame("Coordinate_HG18"=chr_pos[,2],"Chromosome"= chr_pos[,1],"RSquared"=0, "RecombinationRate"=0, "Proxy"=pc_2$rs,
						 "PVAL"=pc_2$p, "TYPE"=rep("imputed",length(pc_2[,1])), stringsAsFactors=F)

r2_pc <- list()

for (n in 1:nrow(new_snps)) {

	top <- new_snps[n,1]
	top_snp_name <- coord[grepl(top, coord$rs), 4]

	top_chr_pos <- coord[grepl(top, coord$rs), 2:3]
	top_chr_pos <- apply(top_chr_pos, 2, as.numeric)

	snp <- subset(for_plot,  for_plot[,2]==top_chr_pos[1], select=1:7)
	snp <- subset(snp, snp[,1]<(top_chr_pos[2]+500000), select=1:7)
	snp <- subset(snp, snp[,1]>(top_chr_pos[2]-500000), select=1:7)



	r2 <- ld[grepl(top_snp_name, ld$SNP_A),]

	for (i in (1:nrow(r2))) {

		snp[which(r2$BP_B[i]==snp[,1]),"RSquared"] <- r2[i, "R2"]

	}

	r2_pc[[n]] <- snp

	rownames(r2_pc[[n]]) <- paste(as.matrix(r2_pc[[n]][2]), as.matrix(r2_pc[[n]][1]), sep=":")
	
	x <- subset(r2_pc[[n]], r2_pc[[n]]$RSquared>0, select=1:7)

	r2_pc[[n]] <- rbind(x, r2_pc[[n]][grepl(top, r2_pc[[n]]$Proxy),])

}


# make r^2 matrix for top snps on Bolormaa data

r2_b <- list()

for (n in 1:nrow(new_snps)) {

	top <- paste(new_snps$CHR[n], new_snps$POS[n], sep=":")


	z_snp <- s_all[rownames(s_all) %in% rownames(r2_pc[[n]]),]
	z_top <- z_snp[which(rownames(z_snp) %in% top),]

	r_snp <- cor(t(z_top), t(z_snp))

	r2_snp <- r_snp^2 

	r2_b[[n]] <- r2_snp

}

