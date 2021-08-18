### make kinship and histogramms ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship")

emmax <- read.table("pool_1.aBN.kinf", head=F, stringsAsFactors=F)

# make histogramms for emmax

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/plots")

emmax_tri <- emmax[lower.tri(emmax, diag=F)]

emmax_diag <- emmax[col(emmax)==row(emmax)]

pdf("emmax_hist.pdf")

hist(emmax_tri, n=50)

dev.off()

pdf("emmax_hist_diag.pdf")

hist(emmax_diag, n=20)

dev.off()

# make kinship with popkin 

library("popkin")
library("BEDMatrix")

m_bed <- BEDMatrix("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1") 

kinship_popkin <- popkin(m_bed)

write.csv(kinship_popkin, "/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/kinship_popkin.csv")

# make histogramms for popkin

pop_tri <- kinship_popkin[lower.tri(kinship_popkin, diag=F)]

pop_diag <- emmax[col(kinship_popkin)==row(kinship_popkin)]

pdf("popkin_hist.pdf")

hist(pop_tri, n=50)

dev.off()

pdf("popkin_hist_diag.pdf")

hist(pop_diag, n=50)

dev.off()

# count correlation between emmax and popkin

kinship_cor <- cor(c(as.matrix(emmax)), c(as.matrix(kinship_popkin)))