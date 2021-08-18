### make manhattan plots on multivariate results ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered")

index_1 <- read.table("p_value_index_1_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_2 <- read.table("p_value_index_2_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_3 <- read.table("p_value_index_3_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_4 <- read.table("p_value_index_4_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_5 <- read.table("p_value_index_5_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_6 <- read.table("p_value_index_6_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_7 <- read.table("p_value_index_7_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
mass <- read.table("p_value_mass_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis")

# make for_plot matrix

m37 <- snp_info[match(index_1[,1], snp_info$SNP_name),c(2, 4:5)]
colnames(m37) <- c ("SNP", "CHR", "BP")

mv_stephens_list_p_value <- list(index_1[,2], index_2[,2], index_3[,2], index_4[,2], index_5[,2], index_6[,2], index_7[,2], mass[,2])

for_plot <- list()

for (n in (1:8)) {

	for_plot[[n]] <- cbind(m37, P=mv_stephens_list_p_value[[n]])

}

## delete X and 99 from for_plot matrix

for (n in (1:8)) {

	for_plot[[n]] <- for_plot[[n]][!grepl("X", for_plot[[n]]$CHR),]
	for_plot[[n]] <- for_plot[[n]][!grepl("99", for_plot[[n]]$CHR),]
	for_plot[[n]]$CHR <- as.numeric(for_plot[[n]]$CHR)
}



# make manhattan plots

library(qqman)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots/manhattan")

names_mv <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")

pdf("manhattan_mv.pdf", width=16, height=10)

for (n in (1:8)) {

	manhattan(for_plot[[n]], main=names_mv[[n]], ylim=c(4,20), genomewideline=-log10(0.05/(nrow(index_1)*8)), suggestiveline=F, cex=1.2, 
    cex.axis=1.2, col=c("blue4", "orange3"))

}

dev.off()