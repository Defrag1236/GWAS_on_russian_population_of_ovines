### make pheno and covariates files for 107 animals ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

library(data.table)

pheno <- fread("pheno_2020.csv", head=T, stringsAsFactors=F, sep=";", na.strings=c("","NA"), data.table=F)

pheno[,-c(1:4)] <-  apply(pheno[,-c(1:4)], 2, as.numeric)

rownames(pheno) <- pheno[,1]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data")

pool <- read.table("pool_1_pool_2_id_without_katahdin_for_plink.txt", stringsAsFactors=F, head=F)[,2]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/pca")

pca <- read.table("pool_1_and_pool_2_without_katahdin.eigenvec", head=F, stringsAsFactors=F)
pca_122 <- read.table("122_animals_with_pheno.eigenvec", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

breed_info <- read.csv("pheno_index_6d_2020.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files")

fam <- read.table("pool_1_and_pool_2_without_katahdin.fam", head=F, stringsAsFactors=F)
fam_122 <- read.table("122_animals_with_pheno.fam", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

pool_1 <-read.table("pool_1.txt", stringsAsFactors=F, head=F)
pool_2 <-read.table("pool_2.txt", stringsAsFactors=F, head=F)

# make qc

## remove values which out +-3 iqr

pheno_qc_pool <- pheno[(pheno[,1] %in% pool), -c(1:4)]

	for (n in (1:ncol(pheno_qc_pool))) {

		pheno_qc_pool[which(pheno_qc_pool[,n]>mean(pheno_qc_pool[,n], na.rm=T)+(IQR(pheno_qc_pool[,n], na.rm=T)*3)),n] <- NA
		pheno_qc_pool[which(pheno_qc_pool[,n]<mean(pheno_qc_pool[,n], na.rm=T)-(IQR(pheno_qc_pool[,n], na.rm=T)*3)),n] <- NA

	}

pheno_qc_122 <- pheno[(pheno[,1] %in% fam_122[,2]), -c(1:4)]

	for (n in (1:ncol(pheno_qc_122))) {

		pheno_qc_122[which(pheno_qc_122[,n]>mean(pheno_qc_122[,n], na.rm=T)+(IQR(pheno_qc_122[,n], na.rm=T)*3)),n] <- NA
		pheno_qc_122[which(pheno_qc_122[,n]<mean(pheno_qc_122[,n], na.rm=T)-(IQR(pheno_qc_122[,n], na.rm=T)*3)),n] <- NA

	}

# prepare pheno_qc for save

pheno_qc_pool <- pheno_qc_pool[match(fam[,2], rownames(pheno_qc_pool)),]
pheno_qc_122 <- pheno_qc_122[match(fam_122[,2], rownames(pheno_qc_122)),]

index_4_6d <- pheno_qc_pool[,"index_4_6d"]
index_4_6d_122 <- pheno_qc_122[,"index_4_6d"]

# make covariate file for 107 animals

cov <- as.data.frame(pheno[(pheno[,1] %in% pool),1:3])
cov <- cov[match(rownames(pheno_qc_pool), rownames(cov)),] 
cov[,4] <- NA
cov[,5] <- NA
cov[,6] <- NA
cov[,7] <- NA
cov[,8] <- NA
cov[,9] <- NA

## add 4 vectors with code for each breed rom_arg=1000, rom_arg_kat=0100, rom_mouf_kat=0010, rom_arg_kat_mouf=0001

rom_arg <- breed_info[grepl("^rom_arg$", breed_info$Breed),1]
rom_arg_kat <-breed_info[grepl("^rom_arg_kat$", breed_info$Breed),1]
rom_mouf_kat <-breed_info[grepl("^rom_mouf_kat$", breed_info$Breed),1]
rom_arg_kat_mouf <-breed_info[grepl("^rom_arg_kat_mouf$", breed_info$Breed),1]

### rom_arg=1000

cov[1:47,4] <- 1
cov[(cov[,1] %in% rom_arg),4] <- 1
cov[1:47,5] <- 0
cov[(cov[,1] %in% rom_arg),5] <- 0 
cov[1:47,6] <- 0
cov[(cov[,1] %in% rom_arg),6] <- 0
cov[1:47,7] <- 0
cov[(cov[,1] %in% rom_arg),7] <- 0

### rom_arg_kat=0100

cov[(cov[,1] %in% rom_arg_kat),4] <- 0
cov[(cov[,1] %in% rom_arg_kat),5] <- 1 
cov[(cov[,1] %in% rom_arg_kat),6] <- 0
cov[(cov[,1] %in% rom_arg_kat),7] <- 0

### rom_mouf_kat=0010

cov[(cov[,1] %in% rom_mouf_kat),4] <- 0
cov[(cov[,1] %in% rom_mouf_kat),5] <- 0 
cov[(cov[,1] %in% rom_mouf_kat),6] <- 1
cov[(cov[,1] %in% rom_mouf_kat),7] <- 0

### rom_arg_kat_mouf=0001

cov[(cov[,1] %in% rom_arg_kat_mouf),4] <- 0
cov[(cov[,1] %in% rom_arg_kat_mouf),5] <- 0 
cov[(cov[,1] %in% rom_arg_kat_mouf),6] <- 0
cov[(cov[,1] %in% rom_arg_kat_mouf),7] <- 1

## add 2 first pca values

cov[,8] <- pca[,3]
cov[,9] <- pca[,4]

## add intercept as 1st column

cov[,1] <- 1

# make cov_ariate file for 122 animals

cov_122 <- as.data.frame(pheno[(pheno[,1] %in% fam_122[,2]),1:3])
cov_122 <- cov_122[match(rownames(pheno_qc_122), rownames(cov_122)),] 
cov_122[,4] <- NA
cov_122[,5] <- NA
cov_122[,6] <- NA
cov_122[,7] <- NA
cov_122[,8] <- NA
cov_122[,9] <- NA
cov_122[,10] <- NA

## add 4 vectors with code for each breed rom_arg=10000, rom_arg_kat=01000, rom_mouf_kat=00100, rom_arg_kat_mouf=00010, katahdin=00001

rom_arg <- breed_info[grepl("^rom_arg$", breed_info$Breed),1]
rom_arg_kat <-breed_info[grepl("^rom_arg_kat$", breed_info$Breed),1]
rom_mouf_kat <-breed_info[grepl("^rom_mouf_kat$", breed_info$Breed),1]
rom_arg_kat_mouf <-breed_info[grepl("^rom_arg_kat_mouf$", breed_info$Breed),1]
katahdin <-breed_info[grepl("^katahdin$", breed_info$Breed),1]

### rom_arg=10000

cov_122[1:47,4] <- 1
cov_122[(cov_122[,1] %in% rom_arg),4] <- 1
cov_122[1:47,5] <- 0
cov_122[(cov_122[,1] %in% rom_arg),5] <- 0 
cov_122[1:47,6] <- 0
cov_122[(cov_122[,1] %in% rom_arg),6] <- 0
cov_122[1:47,7] <- 0
cov_122[(cov_122[,1] %in% rom_arg),7] <- 0
cov_122[1:47,8] <- 0
cov_122[(cov_122[,1] %in% rom_arg),8] <- 0

### rom_arg_kat=01000

cov_122[(cov_122[,1] %in% rom_arg_kat),4] <- 0
cov_122[(cov_122[,1] %in% rom_arg_kat),5] <- 1 
cov_122[(cov_122[,1] %in% rom_arg_kat),6] <- 0
cov_122[(cov_122[,1] %in% rom_arg_kat),7] <- 0
cov_122[(cov_122[,1] %in% rom_arg_kat),8] <- 0

### rom_mouf_kat=00100

cov_122[(cov_122[,1] %in% rom_mouf_kat),4] <- 0
cov_122[(cov_122[,1] %in% rom_mouf_kat),5] <- 0 
cov_122[(cov_122[,1] %in% rom_mouf_kat),6] <- 1
cov_122[(cov_122[,1] %in% rom_mouf_kat),7] <- 0
cov_122[(cov_122[,1] %in% rom_mouf_kat),8] <- 0

### rom_arg_kat_mouf=00010

cov_122[(cov_122[,1] %in% rom_arg_kat_mouf),4] <- 0
cov_122[(cov_122[,1] %in% rom_arg_kat_mouf),5] <- 0 
cov_122[(cov_122[,1] %in% rom_arg_kat_mouf),6] <- 0
cov_122[(cov_122[,1] %in% rom_arg_kat_mouf),7] <- 1
cov_122[(cov_122[,1] %in% rom_arg_kat_mouf),8] <- 0

### katahdin=00001

cov_122[(cov_122[,1] %in% katahdin),4] <- 0
cov_122[(cov_122[,1] %in% katahdin),5] <- 0 
cov_122[(cov_122[,1] %in% katahdin),6] <- 0
cov_122[(cov_122[,1] %in% katahdin),7] <- 0
cov_122[(cov_122[,1] %in% katahdin),8] <- 1

## add 2 first pca values

cov_122[,9] <- pca_122[,3]
cov_122[,10] <- pca_122[,4]

## add intercept as 1st column

cov_122[,1] <- 1

# adjust phenotypes on covariates 107 animals

data <- cbind(cov[,-1], index_4_6d)

colnames(data) <- c("litter_size", "sex", "rom_arg", "rom_arg_kat", "rom_mouf_kat", "rom_arg_kat_mouf", "pc1", "pc2", "index_4_6d")

l <- lm(index_4_6d~litter_size+sex+rom_arg+rom_arg_kat+rom_mouf_kat+rom_arg_kat_mouf+pc1+pc2, data=data)

res <- data$index_4_6d-predict(l,data)

index_4_6d_res <- as.matrix(res)

# normal transformation to adjust phenotypes 107 animals

index_4_6d_res_norm <- qnorm((rank(index_4_6d_res,na.last="keep")-0.5)/sum(!is.na(index_4_6d_res)))

names(index_4_6d_res_norm) <- rownames(pheno_qc_pool)

# make hist of original and adjust pheno

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/plots")

pdf("index_4_6d_hist.pdf")

hist(index_4_6d, main="Original_index_4_6d", n=20)
hist(index_4_6d_res, main="Adjust_index_4_6d", n=20)
hist(index_4_6d_res_norm, main="Adjust_norm_index_4_6d", n=20)

dev.off()

# adjust phenotypes on covariates 122 animals

data_122 <- cbind(cov_122[,-1], index_4_6d_122)

colnames(data_122) <- c("litter_size", "sex", "rom_arg", "rom_arg_kat", "rom_mouf_kat", "rom_arg_kat_mouf", "katahdin", "pc1", "pc2", "index_4_6d")

l_122 <- lm(index_4_6d~litter_size+sex+rom_arg+rom_arg_kat+rom_mouf_kat+rom_arg_kat_mouf+katahdin+pc1+pc2, data=data_122)

res <- data_122$index_4_6d-predict(l_122,data_122)

index_4_6d_res_122 <- as.matrix(res)

# normal transformation to adjust phenotypes 122 animals

index_4_6d_res_norm_122 <- qnorm((rank(index_4_6d_res_122,na.last="keep")-0.5)/sum(!is.na(index_4_6d_res_122)))

names(index_4_6d_res_norm_122) <- rownames(pheno_qc_122)

# make hist of original and adjust pheno for 122 animals

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/plots")

pdf("index_4_6d_hist_122.pdf")

hist(pheno[(pheno[,1] %in% fam_122[,2]), "index_4_6d"], main="Original_index_4_6d_122", n=20)
hist(index_4_6d_122, main="Filtered_original_index_4_6d_122", n=20)
hist(index_4_6d_res_122, main="Adjust_index_4_6d_122", n=20)
hist(index_4_6d_res_norm_122, main="Adjust_norm_index_4_6d", n=20)

dev.off()

# make pheno files for each pool

pool_1_phe <- c()
pool_2_phe <- c()

pool_1_phe[which(noquote(names(index_4_6d_res_norm)) %in% pool_1[,1])] <- index_4_6d_res_norm[which(noquote(names(index_4_6d_res_norm)) %in% pool_1[,1])]
pool_2_phe[which(noquote(names(index_4_6d_res_norm)) %in% pool_2[,1])] <- index_4_6d_res_norm[which(noquote(names(index_4_6d_res_norm)) %in% pool_2[,1])]
pool_2_phe[69:107] <- NA

# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas")

# 107 animals

write.table(index_4_6d, "index_4_6d.txt", col.names=F, row.names=F, quote=F) 
write.table(index_4_6d_res, "index_4_6d_adjust.txt", col.names=F, row.names=F, quote=F)
write.table(index_4_6d_res_norm, "index_4_6d_adjust_norm.txt", col.names=F, row.names=F, quote=F)
write.table(pool_1_phe, "index_4_6d_adjust_norm_pool_1.txt", col.names=F, row.names=F, quote=F) 
write.table(pool_2_phe, "index_4_6d_adjust_norm_pool_2.txt", col.names=F, row.names=F, quote=F)   
write.table(cov, "cov_all.txt", col.names=F, row.names=F, quote=F) 
write.table(cov[,c(1,8:9)], "cov.txt", col.names=F, row.names=F, quote=F)

# 122 animals

write.table(index_4_6d_res_norm_122, "index_4_6d_adjust_norm_122.txt", col.names=F, row.names=F, quote=F)
write.table(cov_122[,c(1,9:10)], "cov_122.txt", col.names=F, row.names=F, quote=F)
