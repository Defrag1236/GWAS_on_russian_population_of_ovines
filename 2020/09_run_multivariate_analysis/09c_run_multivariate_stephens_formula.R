### run multivariate stephens formula ###


# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("MV_sst_file_for_21_metal_traits.Rdata")

# calculate z_score 

z_score_6d <- list()

for (n in c(5, 11, 17, 23, 29, 35, 41)) {

	z_score_6d[[n]] <- sst$gwa[,n]/sst$gwa[,n+1]

}

z_score_42d <- list()

for (n in c(3, 9, 15, 21, 27, 33, 39)) {

	z_score_42d[[n]] <- sst$gwa[,n]/sst$gwa[,n+1]

}

z_score_3m <- list()

for (n in c(1, 7, 13, 19, 25, 31, 37)) {

	z_score_3m[[n]] <- sst$gwa[,n]/sst$gwa[,n+1]

}



z_6d <- do.call(cbind, z_score_6d)
z_42d <- do.call(cbind, z_score_42d)
z_3m <- do.call(cbind, z_score_3m)


colnames(z_6d) <- c("index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d")
colnames(z_42d) <- c("index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d")
colnames(z_3m) <- c("index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m")

z_all <- cbind(z_6d, z_42d, z_3m)

cor_all <- cor(z_all)

cor_6d <- sst$cor.pheno[c(3, 6, 9, 12, 15, 18, 21), c(3, 6, 9, 12, 15, 18, 21)] 
cor_42d <- sst$cor.pheno[c(2, 5, 8, 11, 14, 17, 20), c(2, 5, 8, 11, 14, 17, 20)]
cor_3m <- sst$cor.pheno[c(1, 4, 7, 10, 13, 16, 19), c(1, 4, 7, 10, 13, 16, 19)]

# calculate lamda from sst statistics

p_sst_6d <- list()
p_sst_42d <- list()
p_sst_3m <- list()


for (n in (1:7)) {

	p_sst_6d[[n]] <- pchisq((z_6d[,n])^2, df=1, lower.tail=F)

}

for (n in (1:7)) {

	p_sst_42d[[n]] <- pchisq((z_42d[,n])^2, df=1, lower.tail=F)

}

for (n in (1:7)) {

	p_sst_3m[[n]] <- pchisq((z_3m[,n])^2, df=1, lower.tail=F)

}



p_sst_6d <- do.call(cbind, p_sst_6d)
p_sst_42d <- do.call(cbind, p_sst_42d)
p_sst_3m <- do.call(cbind, p_sst_3m)

lambda_sst_6d <- matrix(ncol=7, nrow=1)
lambda_sst_42d <- matrix(ncol=7, nrow=1)
lambda_sst_3m <- matrix(ncol=7, nrow=1)


for (n in (1:7)) {

	lambda_sst_6d[1,n] <- median (qchisq(p_sst_6d[,n], df=1, lower.tail=F)/qchisq(0.5,df=1))

}

for (n in (1:7)) {

	lambda_sst_42d[1,n] <- median (qchisq(p_sst_42d[,n], df=1, lower.tail=F)/qchisq(0.5,df=1))

}

for (n in (1:7)) {

	lambda_sst_3m[1,n] <- median (qchisq(p_sst_3m[,n], df=1, lower.tail=F)/qchisq(0.5,df=1))

}


lambda_sst <- rbind(lambda_sst_6d, lambda_sst_3m, lambda_sst_42d)

colnames(lambda_sst) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7")
rownames(lambda_sst) <- c("6d", "42d", "3m")

# calculate mv with stephens formula

## 1st trait (6d)

chi_6d_multi <- matrix(ncol=1, nrow=nrow(z_6d))
rownames(chi_6d_multi) <- rownames(z_6d)

for (n  in (1:nrow(z_6d))) {

	Z <- z_6d[n,]
	M <- cor_6d

	chi_6d_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_6d_multi <- matrix(ncol=1, nrow=nrow(z_6d))
p_value_6d_multi[,1] <- pchisq(chi_6d_multi, df=7, low=F)
rownames(p_value_6d_multi) <- rownames(z_6d)

## 2nd trait (42d)

chi_42d_multi <- matrix(ncol=1, nrow=nrow(z_42d))
rownames(chi_42d_multi) <- rownames(z_42d)

for (n  in (1:nrow(z_42d))) {

	Z <- z_42d[n,]
	M <- cor_42d

	chi_42d_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_42d_multi <- matrix(ncol=1, nrow=nrow(z_42d))
p_value_42d_multi[,1] <- pchisq(chi_42d_multi, df=7, low=F)
rownames(p_value_42d_multi) <- rownames(z_42d)

## 3rd trait (3m)

chi_3m_multi <- matrix(ncol=1, nrow=nrow(z_3m))
rownames(chi_3m_multi) <- rownames(z_3m)

for (n  in (1:nrow(z_3m))) {

	Z <- z_3m[n,]
	M <- cor_3m

	chi_3m_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_3m_multi <- matrix(ncol=1, nrow=nrow(z_3m))
p_value_3m_multi[,1] <- pchisq(chi_3m_multi, df=7, low=F)
rownames(p_value_3m_multi) <- rownames(z_3m)

## 4th trait (all)

chi_all_multi <- matrix(ncol=1, nrow=nrow(z_all))
rownames(chi_all_multi) <- rownames(z_all)

for (n  in (1:nrow(z_all))) {

	Z <- z_all[n,]
	M <- sst$cor.pheno

	chi_all_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_all_multi <- matrix(ncol=1, nrow=nrow(z_all))
p_value_all_multi[,1] <- pchisq(chi_all_multi, df=21, low=F)
rownames(p_value_all_multi) <- rownames(z_all)


# calculate lambda

lambda_6d <- c()
lambda_42d <- c()
lambda_3ь <- c()
lambda_all <- c()

lambda_6d <- median (qchisq(p_value_6d_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_42d <- median (qchisq(p_value_42d_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_3m <- median (qchisq(p_value_3m_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_all <- median (qchisq(p_value_all_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))

lambda <- rbind(lambda_6d, lambda_42d, lambda_3m, lambda_all)

# save results 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens")

write.table(p_value_6d_multi, "p_value_6d_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_42d_multi, "p_value_42d_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_3m_multi, "p_value_3m_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_all_multi, "p_value_all_multi_stephens.txt", row.names=T, col.names=F, quote=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

write.table(lambda_sst, "lambda_from_sst_21_trait.txt", row.names=T, col.names=T, quote=F)