### make_pheno_files ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_6d <- read.csv("phenotypes_for_48_sheeps_6d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")
pheno_42d <- read.csv("phenotypes_for_48_sheeps_42d.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")
pheno_3m <- read.csv("phenotypes_for_48_sheeps_3m.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# re write pheno files without colnames and individual id

library(data.table)

fwrite (pheno_6d[3:13], "pheno_6d_for_gwas", col.names=F, row.names=F, quote=F)
fwrite (pheno_42d[3:13], "pheno_42d_for_gwas", col.names=F, row.names=F, quote=F)
fwrite (pheno_3m[3:13], "pheno_3m_for_gwas", col.names=F, row.names=F, quote=F)
