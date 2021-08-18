### prepare pheno data ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

library(data.table)

pheno <- fread("pheno_2020.csv", head=T, stringsAsFactors=F, sep=";", na.strings=c("","NA"), data.table=F)

pheno[,-c(1:4)] <-  apply(pheno[,-c(1:4)], 2, as.numeric)

rownames(pheno) <- pheno[,1]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/pca")
pca <- read.table("122_animals_with_pheno.eigenvec", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

breed_info <- read.csv("pheno_index_6d_2020.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files")

fam <- read.table("122_animals_with_pheno.fam", head=F, stringsAsFactors=F)

# make qc

## remove values which out +-3 iqr

pheno_qc <- pheno[(pheno[,1] %in% fam[,2]), -c(1:4)]

	for (n in (1:ncol(pheno_qc))) {

		pheno_qc[which(pheno_qc[,n]>mean(pheno_qc[,n], na.rm=T)+(IQR(pheno_qc[,n], na.rm=T)*3)),n] <- NA
		pheno_qc[which(pheno_qc[,n]<mean(pheno_qc[,n], na.rm=T)-(IQR(pheno_qc[,n], na.rm=T)*3)),n] <- NA

	}

# prepare pheno_qc for save

pheno_qc <- pheno_qc[match(fam[,2], rownames(pheno_qc)),]

# make covariate file 

cov <- as.data.frame(pheno[(pheno[,1] %in% fam[,2]),1:3])
cov <- cov[match(rownames(pheno_qc), rownames(cov)),] 
cov[,4] <- NA
cov[,5] <- NA
cov[,6] <- NA
cov[,7] <- NA
cov[,8] <- NA
cov[,9] <- NA
cov[,10] <- NA

## add 5 vectors with code for each breed rom_arg=10000, rom_arg_kat=01000, rom_mouf_kat=00100, rom_arg_kat_mouf=00010, katahdin=00001

rom_arg <- breed_info[grepl("^rom_arg$", breed_info$Breed),1]
rom_arg_kat <-breed_info[grepl("^rom_arg_kat$", breed_info$Breed),1]
rom_mouf_kat <-breed_info[grepl("^rom_mouf_kat$", breed_info$Breed),1]
rom_arg_kat_mouf <-breed_info[grepl("^rom_arg_kat_mouf$", breed_info$Breed),1]
katahdin <-breed_info[grepl("^katahdin$", breed_info$Breed),1]

### rom_arg=10000

cov[1:47,4] <- 1
cov[(cov[,1] %in% rom_arg),4] <- 1
cov[1:47,5] <- 0
cov[(cov[,1] %in% rom_arg),5] <- 0 
cov[1:47,6] <- 0
cov[(cov[,1] %in% rom_arg),6] <- 0
cov[1:47,7] <- 0
cov[(cov[,1] %in% rom_arg),7] <- 0
cov[1:47,8] <- 0
cov[(cov[,1] %in% rom_arg),8] <- 0

### rom_arg_kat=01000

cov[(cov[,1] %in% rom_arg_kat),4] <- 0
cov[(cov[,1] %in% rom_arg_kat),5] <- 1 
cov[(cov[,1] %in% rom_arg_kat),6] <- 0
cov[(cov[,1] %in% rom_arg_kat),7] <- 0
cov[(cov[,1] %in% rom_arg_kat),8] <- 0

### rom_mouf_kat=00100

cov[(cov[,1] %in% rom_mouf_kat),4] <- 0
cov[(cov[,1] %in% rom_mouf_kat),5] <- 0 
cov[(cov[,1] %in% rom_mouf_kat),6] <- 1
cov[(cov[,1] %in% rom_mouf_kat),7] <- 0
cov[(cov[,1] %in% rom_mouf_kat),8] <- 0

### rom_arg_kat_mouf=00010

cov[(cov[,1] %in% rom_arg_kat_mouf),4] <- 0
cov[(cov[,1] %in% rom_arg_kat_mouf),5] <- 0 
cov[(cov[,1] %in% rom_arg_kat_mouf),6] <- 0
cov[(cov[,1] %in% rom_arg_kat_mouf),7] <- 1
cov[(cov[,1] %in% rom_arg_kat_mouf),8] <- 0

### katahdin=00001

cov[(cov[,1] %in% katahdin),4] <- 0
cov[(cov[,1] %in% katahdin),5] <- 0 
cov[(cov[,1] %in% katahdin),6] <- 0
cov[(cov[,1] %in% katahdin),7] <- 0
cov[(cov[,1] %in% katahdin),8] <- 1

## add 2 first pca values

cov[,9] <- pca[,3]
cov[,10] <- pca[,4]

## add intercept as 1st column

cov[,1] <- 1

# adjust phenotypes on covariates 

data <- cbind(cov[,-1], pheno_qc[,c(1:21,29:31)])

colnames(data)[1:9] <- c("litter_size", "sex", "rom_arg", "rom_arg_kat", "rom_mouf_kat", "rom_arg_kat_mouf", "katahdin", "pc1", "pc2")

pheno_adj <- matrix(ncol=ncol(data), nrow=nrow(data))
colnames(pheno_adj) <- colnames(data)
rownames(pheno_adj) <- rownames(data)

for (n in 10:33) {

	l <- lm(data[,n]~litter_size+sex+rom_arg+rom_arg_kat+rom_mouf_kat+rom_arg_kat_mouf+katahdin+pc1+pc2, data=data)

	res <- data[,n]-predict(l,data)

 	pheno_adj [,n] <- as.matrix(res)

}

pheno_adj <- pheno_adj[,-c(1:9)]

# normal transformation to adjust phenotypes animals

pheno_adj_norm <- matrix(ncol=ncol(pheno_adj), nrow=nrow(pheno_adj))
colnames(pheno_adj_norm) <- colnames(pheno_adj)
rownames(pheno_adj_norm) <- rownames(pheno_adj)

for (n in 1:ncol(pheno_adj)) {

	pheno_adj_norm[,n] <- qnorm((rank(pheno_adj[,n],na.last="keep")-0.5)/sum(!is.na(pheno_adj[,n])))


}


# make hist for pheno_adjust

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/plots")

pdf("pheno_adj_122_animals.pdf")

for (n in 1:ncol(pheno_adj)) {

	hist(pheno_adj[,n], main=colnames(pheno_adj)[n], n=20)

}

dev.off()

# count sample size

sample_size <- apply(pheno_adj_norm, 2, function(x) table(!is.na(x))[2])
sample_size[is.na(sample_size)] <- 122
sample_size <- as.matrix(sample_size)

# save results 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pheno_data/pheno_qc_122_animals")

for (n in (1:ncol(pheno_adj_norm))) {

	name_to_save <- paste(colnames(pheno_adj_norm)[n], ".txt", sep="")

	write.table(pheno_adj_norm[,n], name_to_save, col.names=F, row.names=F, quote=F)

}


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pheno_data")

write.table(sample_size, "sample_size_122_animals.txt", col.names=F, row.names=T, quote=F)

write.table(pheno_qc[,c(1:21, 29:31)], "pheno_qc_122_animals.txt", col.names=T, row.names=T, quote=F)