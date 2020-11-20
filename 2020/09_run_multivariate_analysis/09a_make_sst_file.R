### make sst file ###

library(MultiABEL)

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_21_traits_unified.Rdata")

# make additional objectivies

indep.snps <- metal_6d_clear[[1]]$MarkerName

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/for_multiabel")

file.names <- list.files()

# make sst file

sst = load.summary(files=file.names,
	cor.pheno = NULL, indep.snps = as.character(indep.snps), est.var = TRUE,
	columnNames = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "N"), fixedN = NULL)

# do corrplot for sst$cor.pheno

for_order <- c(3, 6, 9, 12, 15, 18, 21, 2, 5, 8, 11, 14, 17, 20, 1, 4, 7, 10, 13, 16, 19)

for_plot <- sst$cor.pheno[for_order, for_order]

library(corrplot)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("corrplot_for_21_metal_traits_multiabel.pdf", width=16, height=10)

corrplot(for_plot, tl.cex=0.8)

dev.off()

# save sst file

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

save(sst, file="MV_sst_file_for_21_metal_traits.Rdata")