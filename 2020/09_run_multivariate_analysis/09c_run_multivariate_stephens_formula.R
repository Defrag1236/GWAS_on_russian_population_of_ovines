### run multivariate stephens formula ###


# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("MV_sst_file_for_24_metal_traits.Rdata")

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

z_score_mass <- list()

for (n in c(47, 45, 43)) {

	z_score_mass[[n]] <- sst$gwa[,n]/sst$gwa[,n+1]

}

z_6d <- do.call(cbind, z_score_6d)
z_42d <- do.call(cbind, z_score_42d)
z_3m <- do.call(cbind, z_score_3m)
z_mass <- do.call(cbind, z_score_mass)

colnames(z_6d) <- c("index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d")
colnames(z_42d) <- c("index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d")
colnames(z_3m) <- c("index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m")
colnames(z_mass) <- c("mass_6d", "mass_42d", "mass_3m")
 
z_all <- cbind(z_6d, z_42d, z_3m, z_mass)

cor_all <- cor(z_all)

cor_6d <- sst$cor.pheno[c(3, 6, 9, 12, 15, 18, 21), c(3, 6, 9, 12, 15, 18, 21)] 
cor_42d <- sst$cor.pheno[c(2, 5, 8, 11, 14, 17, 20), c(2, 5, 8, 11, 14, 17, 20)]
cor_3m <- sst$cor.pheno[c(1, 4, 7, 10, 13, 16, 19), c(1, 4, 7, 10, 13, 16, 19)]


z_index_1 <- cbind(index_1_6d=z_6d[,1], index_1_42d=z_42d[,1], index_1_3m=z_3m[,1])
z_index_2 <- cbind(index_2_6d=z_6d[,2], index_2_42d=z_42d[,2], index_2_3m=z_3m[,2])
z_index_3 <- cbind(index_3_6d=z_6d[,3], index_3_42d=z_42d[,3], index_3_3m=z_3m[,3])
z_index_4 <- cbind(index_4_6d=z_6d[,4], index_4_42d=z_42d[,4], index_4_3m=z_3m[,4])
z_index_5 <- cbind(index_5_6d=z_6d[,5], index_5_42d=z_42d[,5], index_5_3m=z_3m[,5])
z_index_6 <- cbind(index_6_6d=z_6d[,6], index_6_42d=z_42d[,6], index_6_3m=z_3m[,6])
z_index_7 <- cbind(index_7_6d=z_6d[,7], index_7_42d=z_42d[,7], index_7_3m=z_3m[,7])

cor_index_1 <- sst$cor.pheno[c(3, 2, 1), c(3, 2, 1)] 
cor_index_2 <- sst$cor.pheno[c(6, 5, 4), c(6, 5, 4)] 
cor_index_3 <- sst$cor.pheno[c(9, 8, 7), c(9, 8, 7)] 
cor_index_4 <- sst$cor.pheno[c(12, 11, 10), c(12, 11, 10)] 
cor_index_5 <- sst$cor.pheno[c(15, 14, 13), c(15, 14, 13)] 
cor_index_6 <- sst$cor.pheno[c(18, 17, 16), c(18, 17, 16)] 
cor_index_7 <- sst$cor.pheno[c(21, 20, 19), c(21, 20, 19)] 
cor_mass <- sst$cor.pheno[c(24, 23, 22), c(24, 23, 22)] 

# calculate lambda from sst statistics

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

# calculate mv stephens formula for each index

## index_1

chi_index_1_multi <- matrix(ncol=1, nrow=nrow(z_index_1))
rownames(chi_index_1_multi) <- rownames(z_index_1)

for (n  in (1:nrow(z_index_1))) {

	Z <- z_index_1[n,]
	M <- cor_index_1

	chi_index_1_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_1_multi <- matrix(ncol=1, nrow=nrow(z_index_1))
p_value_index_1_multi[,1] <- pchisq(chi_index_1_multi, df=3, low=F)
rownames(p_value_index_1_multi) <- rownames(z_index_1)

## index_2

chi_index_2_multi <- matrix(ncol=1, nrow=nrow(z_index_2))
rownames(chi_index_2_multi) <- rownames(z_index_2)

for (n  in (1:nrow(z_index_2))) {

	Z <- z_index_2[n,]
	M <- cor_index_2

	chi_index_2_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_2_multi <- matrix(ncol=1, nrow=nrow(z_index_2))
p_value_index_2_multi[,1] <- pchisq(chi_index_2_multi, df=3, low=F)
rownames(p_value_index_2_multi) <- rownames(z_index_2)

## index_3

chi_index_3_multi <- matrix(ncol=1, nrow=nrow(z_index_3))
rownames(chi_index_3_multi) <- rownames(z_index_3)

for (n  in (1:nrow(z_index_3))) {

	Z <- z_index_3[n,]
	M <- cor_index_3

	chi_index_3_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_3_multi <- matrix(ncol=1, nrow=nrow(z_index_3))
p_value_index_3_multi[,1] <- pchisq(chi_index_3_multi, df=3, low=F)
rownames(p_value_index_3_multi) <- rownames(z_index_3)

## index_4

chi_index_4_multi <- matrix(ncol=1, nrow=nrow(z_index_4))
rownames(chi_index_4_multi) <- rownames(z_index_4)

for (n  in (1:nrow(z_index_4))) {

	Z <- z_index_4[n,]
	M <- cor_index_4

	chi_index_4_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_4_multi <- matrix(ncol=1, nrow=nrow(z_index_4))
p_value_index_4_multi[,1] <- pchisq(chi_index_4_multi, df=3, low=F)
rownames(p_value_index_4_multi) <- rownames(z_index_4)

## index_5

chi_index_5_multi <- matrix(ncol=1, nrow=nrow(z_index_5))
rownames(chi_index_5_multi) <- rownames(z_index_5)

for (n  in (1:nrow(z_index_5))) {

	Z <- z_index_5[n,]
	M <- cor_index_5

	chi_index_5_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_5_multi <- matrix(ncol=1, nrow=nrow(z_index_5))
p_value_index_5_multi[,1] <- pchisq(chi_index_5_multi, df=3, low=F)
rownames(p_value_index_5_multi) <- rownames(z_index_5)

## index_6

chi_index_6_multi <- matrix(ncol=1, nrow=nrow(z_index_6))
rownames(chi_index_6_multi) <- rownames(z_index_6)

for (n  in (1:nrow(z_index_6))) {

	Z <- z_index_6[n,]
	M <- cor_index_6

	chi_index_6_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_6_multi <- matrix(ncol=1, nrow=nrow(z_index_6))
p_value_index_6_multi[,1] <- pchisq(chi_index_6_multi, df=3, low=F)
rownames(p_value_index_6_multi) <- rownames(z_index_6)

## index_7

chi_index_7_multi <- matrix(ncol=1, nrow=nrow(z_index_7))
rownames(chi_index_7_multi) <- rownames(z_index_7)

for (n  in (1:nrow(z_index_7))) {

	Z <- z_index_7[n,]
	M <- cor_index_7

	chi_index_7_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_index_7_multi <- matrix(ncol=1, nrow=nrow(z_index_7))
p_value_index_7_multi[,1] <- pchisq(chi_index_7_multi, df=3, low=F)
rownames(p_value_index_7_multi) <- rownames(z_index_7)

## mass 

chi_mass_multi <- matrix(ncol=1, nrow=nrow(z_mass))
rownames(chi_mass_multi) <- rownames(z_mass)

for (n  in (1:nrow(z_mass))) {

	Z <- z_mass[n,]
	M <- cor_mass

	chi_mass_multi[n] <- t(Z) %*% solve(M) %*% Z

	}

p_value_mass_multi <- matrix(ncol=1, nrow=nrow(z_mass))
p_value_mass_multi[,1] <- pchisq(chi_mass_multi, df=3, low=F)
rownames(p_value_mass_multi) <- rownames(z_mass)

# calculate lambda for index traits

lambda_index <- matrix(ncol=1, nrow=8)

colnames(lambda_index) <- "lambda"
rownames(lambda_index) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")

lambda_index[1,] <- median (qchisq(p_value_index_1_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[2,] <- median (qchisq(p_value_index_2_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[3,] <- median (qchisq(p_value_index_3_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[4,] <- median (qchisq(p_value_index_4_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[5,] <- median (qchisq(p_value_index_5_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[6,] <- median (qchisq(p_value_index_6_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[7,] <- median (qchisq(p_value_index_7_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_index[8,] <- median (qchisq(p_value_mass_multi[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))


# save results 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens")

write.table(p_value_index_1_multi, "p_value_index_1_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_2_multi, "p_value_index_2_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_3_multi, "p_value_index_3_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_4_multi, "p_value_index_4_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_5_multi, "p_value_index_5_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_6_multi, "p_value_index_6_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_7_multi, "p_value_index_7_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_mass_multi, "p_value_mass_multi_stephens.txt", row.names=T, col.names=F, quote=F)

write.table(p_value_6d_multi, "p_value_6d_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_42d_multi, "p_value_42d_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_3m_multi, "p_value_3m_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_all_multi, "p_value_all_multi_stephens.txt", row.names=T, col.names=F, quote=F)

write.table(lambda_index, "lambda_index.txt", row.names=T, col.names=T, quote=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

write.table(lambda_sst, "lambda_from_sst_21_trait.txt", row.names=T, col.names=T, quote=F)
