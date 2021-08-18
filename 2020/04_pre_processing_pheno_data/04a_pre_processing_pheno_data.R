### pre_processing_pheno_data ### 

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

library(data.table)

pheno <- fread("pheno_2020.csv", head=T, stringsAsFactors=F, sep=";", na.strings=c("","NA"), data.table=F)

pheno[,-c(1:4)] <-  apply(pheno[,-c(1:4)], 2, as.numeric)

rownames(pheno) <- pheno[,1]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

pool_1 <-read.table("pool_1.txt", stringsAsFactors=F, head=F)
pool_2 <-read.table("pool_2.txt", stringsAsFactors=F, head=F)
pool_3 <-read.table("pool_3.txt", stringsAsFactors=F, head=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pca")

pca <- read.table("pca_192_animals_filtered.eigenvec", head=F, stringsAsFactors=F)

# make hist for each trait 

for_plot <- pheno[,-c(1:4)]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf ("pheno_hist.pdf")

for (n in (1:ncol(for_plot))) {

	hist(for_plot[,n], main=colnames(for_plot)[n])

}

dev.off()


# make qc

## remove values which out +-3 iqr

pheno_qc_pool_1 <- pheno[(pheno[,1] %in% pool_1[,1]), -c(1:4)]
pheno_qc_pool_2 <- pheno[(pheno[,1] %in% pool_2[,1]), -c(1:4)]

for (n in (1:ncol(pheno_qc_pool_1))) {

	pheno_qc_pool_1[which(pheno_qc_pool_1[,n]>mean(pheno_qc_pool_1[,n], na.rm=T)+(IQR(pheno_qc_pool_1[,n], na.rm=T)*3)),n] <- NA
	pheno_qc_pool_1[which(pheno_qc_pool_1[,n]<mean(pheno_qc_pool_1[,n], na.rm=T)-(IQR(pheno_qc_pool_1[,n], na.rm=T)*3)),n] <- NA

	pheno_qc_pool_2[which(pheno_qc_pool_2[,n]>mean(pheno_qc_pool_2[,n], na.rm=T)+(IQR(pheno_qc_pool_2[,n], na.rm=T)*3)),n] <- NA
	pheno_qc_pool_2[which(pheno_qc_pool_2[,n]<mean(pheno_qc_pool_2[,n], na.rm=T)-(IQR(pheno_qc_pool_2[,n], na.rm=T)*3)),n] <- NA


}

## make plot after qc

for_plot_qc <- pheno_qc

pdf ("pheno_hist_qc.pdf")

for (n in (1:ncol(for_plot_qc))) {

	hist(for_plot_qc[,n], main=colnames(for_plot_qc)[n])

}

dev.off()


# make plots for each pool of animals 

pdf ("pheno_hist_qc_pool_1.pdf")

for (n in (1:ncol(pheno_qc_pool_1))) {

	hist(pheno_qc_pool_1[,n], main=colnames(pheno_qc_pool_1)[n], n=20)

}

dev.off()

pdf ("pheno_hist_qc_pool_2.pdf")

for (n in (1:ncol(pheno_qc_pool_2))) {

	hist(pheno_qc_pool_2[,n], main=colnames(pheno_qc_pool_2)[n], n=20)

}

dev.off()

# prepare pheno_qc for save

pheno_qc_pool_1 <- pheno_qc_pool_1[match(rownames(pheno_qc_pool_1), pool_1[,1]),]

pheno_qc_pool_2 <- pheno_qc_pool_2[match(rownames(pheno_qc_pool_2), pool_2[,1]),]

pheno_qc_pool_3 <- pheno_qc[(pheno[,1] %in% pool_3[,1]),]
pheno_qc_pool_3 <- pheno_qc_pool_3[match(rownames(pheno_qc_pool_3), pool_3[,1]),]

# calculate sample size for each trait

sample_size_pool_1 <- colSums(!(is.na(pheno_qc_pool_1)))
sample_size_pool_2 <- colSums(!(is.na(pheno_qc_pool_2)))
sample_size_pool_3 <- colSums(!(is.na(pheno_qc_pool_3)))

sample_size_pool_1_2 <- sample_size_pool_1 + sample_size_pool_2

# save qc pheno 


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data/pheno_qc")

for (n in (1:ncol(pheno_qc_pool_1))) {

	name_to_save <- paste(colnames(pheno_qc_pool_1)[n], "_pool_1.txt", sep="")

	write.table(pheno_qc_pool_1[,n], name_to_save, col.names=F, row.names=F, quote=F)

}

for (n in (1:ncol(pheno_qc_pool_2))) {

	name_to_save <- paste(colnames(pheno_qc_pool_2)[n], "_pool_2.txt", sep="")

	write.table(pheno_qc_pool_2[,n], name_to_save, col.names=F, row.names=F, quote=F)

}

for (n in (1:ncol(pheno_qc_pool_3))) {

	name_to_save <- paste(colnames(pheno_qc)[n], "_pool_3.txt", sep="")

	write.table(pheno_qc_pool_3[,n], name_to_save, col.names=F, row.names=F, quote=F)

}

# make covariates files 

#pca_pool_1 <- pca[pca[,2] %in% pool_1[,1],2:4]
#pca_pool_2 <- pca[pca[,2] %in% pool_2[,1],2:4]

pool_1_cov <- as.data.frame(pheno[(pheno[,1] %in% pool_1[,1]),1:3])
pool_1_cov <- pool_1_cov[match(rownames(pool_1_cov), pool_1[,1]),]
#pool_1_cov[,4:5] <- pca_pool_1[match(pca_pool_1[,1], pool_1[,1]),2:3]
pool_1_cov[1:46,4] <- 0
pool_1_cov[47:85,4] <- 1  
pool_1_cov[,1] <- 1

pool_1_cov_gcta <- cbind(0, rownames(pool_1_cov), pool_1_cov[,3:4])
pool_1_cov_gcta_q <- cbind(0, rownames(pool_1_cov), pool_1_cov[,2])


pool_2_cov <- as.data.frame(pheno[(pheno[,1] %in% pool_2[,1]),1:3])
pool_2_cov <- pool_2_cov[match(rownames(pool_2_cov), pool_2[,1]),]
#pool_2_cov[,4:5] <- pca_pool_2[match(pca_pool_2[,1], pool_2[,1]),2:3]
pool_2_cov[,1] <- 1

pool_3_cov <- pheno[(pheno[,1] %in% pool_3[,1]),1:3] 
pool_3_cov <- pool_3_cov[match(rownames(pool_3_cov), pool_3[,1]),]
pool_3_cov[,1] <- 1


# save .phen files for gcta (reml analysis)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data/pheno_for_gcta")

for (n in (1:ncol(pheno_qc_pool_1))) {

	name_to_save <- paste(colnames(pheno_qc_pool_1)[n], "_pool_1.phen", sep="")

	x <- cbind(0, rownames(pheno_qc_pool_1), pheno_qc_pool_1[,n])

	x[is.na(x[,3]),3] <- -9

	write.table(x, name_to_save, col.names=F, row.names=F, quote=F)

}


# save covariates files 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

write.table(pool_1_cov, "pool_1_cov.txt", col.names=F, row.names=F, quote=F)
write.table(pool_2_cov, "pool_2_cov.txt", col.names=F, row.names=F, quote=F)
write.table(pool_3_cov, "pool_3_cov.txt", col.names=F, row.names=F, quote=F)
write.table(pool_1_cov_gcta, "pool_1_cov_gcta.covar", col.names=F, row.names=F, quote=F)
write.table(pool_1_cov_gcta_q, "pool_1_cov_gcta_q.covar", col.names=F, row.names=F, quote=F)

# save sample size files

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

write.table(as.matrix(sample_size_pool_1), "sample_size_pool_1.txt", col.names=F, row.names=T, quote=F)
write.table(as.matrix(sample_size_pool_2), "sample_size_pool_2.txt", col.names=F, row.names=T, quote=F)
write.table(as.matrix(sample_size_pool_3), "sample_size_pool_3.txt", col.names=F, row.names=T, quote=F)
write.table(as.matrix(sample_size_pool_1_2), "sample_size_pool_1_2.txt", col.names=F, row.names=T, quote=F)