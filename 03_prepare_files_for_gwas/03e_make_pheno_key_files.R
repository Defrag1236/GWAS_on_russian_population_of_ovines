### make pheno key files ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d <- read.csv("phenotypes_for_48_sheeps_6d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")
pheno_42d <- read.csv("phenotypes_for_48_sheeps_42d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")
pheno_3m <- read.csv("phenotypes_for_48_sheeps_3m.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# make pheno key files 

pheno_6d_key <- matrix(ncol=2, nrow=11)
pheno_6d_key[,1] <- seq(1:11)
pheno_6d_key[,2] <- colnames(pheno_6d)[3:13]

pheno_42d_key <- matrix(ncol=2, nrow=11)
pheno_42d_key[,1] <- seq(1:11)
pheno_42d_key[,2] <- colnames(pheno_42d)[3:13]

pheno_3m_key <- matrix(ncol=2, nrow=11)
pheno_3m_key[,1] <- seq(1:11)
pheno_3m_key[,2] <- colnames(pheno_3m)[3:13]

# save pheno key files

write.table(pheno_6d_key, "pheno_6d_key.txt", col.names=F, row.names=F, quote=F)
write.table(pheno_42d_key, "pheno_42d_key.txt", col.names=F, row.names=F, quote=F)
write.table(pheno_3m_key, "pheno_3m_key.txt", col.names=F, row.names=F, quote=F)