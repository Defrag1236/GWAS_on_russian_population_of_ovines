### make manhattan on uv results ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered")

index_1 <- read.table("p_value_index_1_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis")

maf <- read.table("maf_for_each_snp.txt", head=F, stringsAsFactors=F)


# make for_plot lists 

maf_f <- subset(maf, maf[,2]>0.05, select=1)

for_plot_6d <- list()
for_plot_42d <- list()
for_plot_3m <- list()
for_plot_mass <- list()

for (n in (1:7)) {

	for_plot_6d_tmp <- as.data.frame(metal_6d_clear[[n]][metal_6d_clear[[n]][,1] %in% index_1[,1],c(1,8)])
	for_plot_6d_tmp[,3:4] <- snp_info[match(for_plot_6d_tmp[,1], snp_info$SNP_name), 4:5]
	colnames(for_plot_6d_tmp) <- c("SNP", "P", "CHR", "BP")
	for_plot_6d_tmp <- for_plot_6d_tmp[!grepl("X", for_plot_6d_tmp$CHR),]
	for_plot_6d_tmp <- for_plot_6d_tmp[!grepl("99", for_plot_6d_tmp$CHR),]
	for_plot_6d_tmp$CHR <- as.numeric(for_plot_6d_tmp$CHR)
	for_plot_6d_tmp <- for_plot_6d_tmp[match(maf_f[,1], for_plot_6d_tmp[,1]),]
	for_plot_6d_tmp <- for_plot_6d_tmp[complete.cases(for_plot_6d_tmp),]
	for_plot_6d[[n]] <- for_plot_6d_tmp

	for_plot_42d_tmp <- as.data.frame(metal_42d_clear[[n]][metal_42d_clear[[n]][,1] %in% index_1[,1],c(1,8)])
	for_plot_42d_tmp[,3:4] <- snp_info[match(for_plot_42d_tmp[,1], snp_info$SNP_name), 4:5]
	colnames(for_plot_42d_tmp) <- c("SNP", "P", "CHR", "BP")
	for_plot_42d_tmp <- for_plot_42d_tmp[!grepl("X", for_plot_42d_tmp$CHR),]
	for_plot_42d_tmp <- for_plot_42d_tmp[!grepl("99", for_plot_42d_tmp$CHR),]
	for_plot_42d_tmp$CHR <- as.numeric(for_plot_42d_tmp$CHR)
	for_plot_42d_tmp <- for_plot_42d_tmp[match(maf_f[,1], for_plot_42d_tmp[,1]),]
	for_plot_42d_tmp <- for_plot_42d_tmp[complete.cases(for_plot_42d_tmp),]
	for_plot_42d[[n]] <- for_plot_42d_tmp


	for_plot_3m_tmp <- as.data.frame(metal_3m_clear[[n]][metal_3m_clear[[n]][,1] %in% index_1[,1],c(1,8)])
	for_plot_3m_tmp[,3:4] <- snp_info[match(for_plot_3m_tmp[,1], snp_info$SNP_name), 4:5]
	colnames(for_plot_3m_tmp) <- c("SNP", "P", "CHR", "BP")
	for_plot_3m_tmp <- for_plot_3m_tmp[!grepl("X", for_plot_3m_tmp$CHR),]
	for_plot_3m_tmp <- for_plot_3m_tmp[!grepl("99", for_plot_3m_tmp$CHR),]
	for_plot_3m_tmp$CHR <- as.numeric(for_plot_3m_tmp$CHR)
	for_plot_3m_tmp <- for_plot_3m_tmp[match(maf_f[,1], for_plot_3m_tmp[,1]),]
	for_plot_3m_tmp <- for_plot_3m_tmp[complete.cases(for_plot_3m_tmp),]
	for_plot_3m[[n]] <- for_plot_3m_tmp

	print(n)

}

for (n in (1:3)) {

	for_plot_mass_tmp <- as.data.frame(metal_mass[[n]][metal_mass[[n]][,1] %in% index_1[,1],c(1,8)])
	for_plot_mass_tmp[,3:4] <- snp_info[match(for_plot_mass_tmp[,1], snp_info$SNP_name), 4:5]
	colnames(for_plot_mass_tmp) <- c("SNP", "P", "CHR", "BP")
	for_plot_mass_tmp <- for_plot_mass_tmp[!grepl("X", for_plot_mass_tmp$CHR),]
	for_plot_mass_tmp <- for_plot_mass_tmp[!grepl("99", for_plot_mass_tmp$CHR),]
	for_plot_mass_tmp$CHR <- as.numeric(for_plot_mass_tmp$CHR)
	for_plot_mass_tmp <- for_plot_mass_tmp[match(maf_f[,1], for_plot_mass_tmp[,1]),]
	for_plot_mass_tmp <- for_plot_mass_tmp[complete.cases(for_plot_mass_tmp),]
	for_plot_mass[[n]] <- for_plot_mass_tmp


}

# make manhattan plots 

library(qqman)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots/manhattan")

names_mv <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")

pdf("manhattan_metal_results.pdf", width=16, height=10)

for (n in (1:7)) {

	manhattan(for_plot_6d[[n]], main=paste(names_mv[n], "_6d", sep=""), ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))

}

manhattan(for_plot_mass[[1]], main=paste(names_mv[8], "_6d", sep=""), ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))


for (n in (1:7)) {

	manhattan(for_plot_42d[[n]], main=paste(names_mv[n], "_42d", sep=""), ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))

}

manhattan(for_plot_mass[[2]], main=paste(names_mv[8], "_42d", sep=""), ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))


for (n in (1:7)) {

	manhattan(for_plot_3m[[n]], main=paste(names_mv[n], "_3m", sep=""), ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))

}

manhattan(for_plot_mass[[3]], main=paste(names_mv[8], "_3m", sep=""), ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))


dev.off()
