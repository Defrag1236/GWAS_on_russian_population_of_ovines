### run mv analysis stephens ###


# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/Rdata")

load("metal_results_24_traits_from_database.Rdata")

# make z_score lists

z_score_6d <- list()

for (n in 1:7) {

	z_score_6d[[n]] <- metal_6d[[n]]$z

}

z_score_42d <- list()

for (n in 1:7) {

	z_score_42d[[n]] <- metal_42d[[n]]$z

}

z_score_3m <- list()

for (n in 1:7) {

	z_score_3m[[n]] <- metal_3m[[n]]$z

}

z_score_mass <- list()

for (n in 1:3) {

	z_score_mass[[n]] <- metal_mass[[n]]$z

}

z_6d <- do.call(cbind, z_score_6d)
z_42d <- do.call(cbind, z_score_42d)
z_3m <- do.call(cbind, z_score_3m)
z_mass <- do.call(cbind, z_score_mass)

colnames(z_6d) <- c("index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d")
colnames(z_42d) <- c("index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d")
colnames(z_3m) <- c("index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m")
colnames(z_mass) <- c("mass_6d", "mass_42d", "mass_3m")

rownames(z_6d) <- metal_6d[[1]]$rs_id
rownames(z_42d) <- metal_6d[[1]]$rs_id
rownames(z_3m) <- metal_6d[[1]]$rs_id
rownames(z_mass) <- metal_6d[[1]]$rs_id

z_all <- cbind(z_6d, z_42d, z_3m, z_mass)

cor_all <- cor(z_all)

z_index_1 <- cbind(index_1_6d=z_6d[,1], index_1_42d=z_42d[,1], index_1_3m=z_3m[,1])
z_index_2 <- cbind(index_2_6d=z_6d[,2], index_2_42d=z_42d[,2], index_2_3m=z_3m[,2])
z_index_3 <- cbind(index_3_6d=z_6d[,3], index_3_42d=z_42d[,3], index_3_3m=z_3m[,3])
z_index_4 <- cbind(index_4_6d=z_6d[,4], index_4_42d=z_42d[,4], index_4_3m=z_3m[,4])
z_index_5 <- cbind(index_5_6d=z_6d[,5], index_5_42d=z_42d[,5], index_5_3m=z_3m[,5])
z_index_6 <- cbind(index_6_6d=z_6d[,6], index_6_42d=z_42d[,6], index_6_3m=z_3m[,6])
z_index_7 <- cbind(index_7_6d=z_6d[,7], index_7_42d=z_42d[,7], index_7_3m=z_3m[,7])

cor_index_1 <- cor_all[c("index_1_6d", "index_1_42d", "index_1_3m"), c("index_1_6d", "index_1_42d", "index_1_3m")]
cor_index_2 <- cor_all[c("index_2_6d", "index_2_42d", "index_2_3m"), c("index_2_6d", "index_2_42d", "index_2_3m")]
cor_index_3 <- cor_all[c("index_3_6d", "index_3_42d", "index_3_3m"), c("index_3_6d", "index_3_42d", "index_3_3m")]
cor_index_4 <- cor_all[c("index_4_6d", "index_4_42d", "index_4_3m"), c("index_4_6d", "index_4_42d", "index_4_3m")]
cor_index_5 <- cor_all[c("index_5_6d", "index_5_42d", "index_5_3m"), c("index_5_6d", "index_5_42d", "index_5_3m")]
cor_index_6 <- cor_all[c("index_6_6d", "index_6_42d", "index_6_3m"), c("index_6_6d", "index_6_42d", "index_6_3m")]
cor_index_7 <- cor_all[c("index_7_6d", "index_7_42d", "index_7_3m"), c("index_7_6d", "index_7_42d", "index_7_3m")]
cor_mass <- cor_all[c("mass_6d", "mass_42d", "mass_3m"), c("mass_6d", "mass_42d", "mass_3m")]

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

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/mv_stephens/on_database_gwas")

write.table(p_value_index_1_multi, "p_value_index_1_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_2_multi, "p_value_index_2_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_3_multi, "p_value_index_3_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_4_multi, "p_value_index_4_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_5_multi, "p_value_index_5_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_6_multi, "p_value_index_6_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_index_7_multi, "p_value_index_7_multi_stephens.txt", row.names=T, col.names=F, quote=F)
write.table(p_value_mass_multi, "p_value_mass_multi_stephens.txt", row.names=T, col.names=F, quote=F)

write.table(lambda_index, "lambda_index.txt", row.names=T, col.names=T, quote=F)
