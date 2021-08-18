### add n to gwas_results ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pheno_data")

sample_size <- read.table("sample_size_122_animals.txt", head=F, stringsAsFactors=F)

library(data.table)

gwas_6d <- list()
gwas_42d <- list()
gwas_3m <- list()


for (n in (1:7)) {

	name_to_read_6d <- paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/6d/output/index_", n, "_6d_adjust_norm.assoc.txt", sep="")
	name_to_list_6d <- paste("index_", n, "_6d", sep="")
	gwas_6d[[n]] <- fread(name_to_read_6d, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_6d)[[n]] <- name_to_list_6d

	name_to_read_42d <- paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/42d/output/index_", n, "_42d_adjust_norm.assoc.txt", sep="")
	name_to_list_42d <- paste("index_", n, "_42d", sep="")
	gwas_42d[[n]] <- fread(name_to_read_42d, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_42d)[[n]] <- name_to_list_42d

	name_to_read_3m <- paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/3m/output/index_", n, "_3m_adjust_norm.assoc.txt", sep="")
	name_to_list_3m <- paste("index_", n, "_3m", sep="")
	gwas_3m[[n]] <- fread(name_to_read_3m, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_3m)[[n]] <- name_to_list_3m
 
}

gwas_m <- list()

gwas_m[[1]] <- fread("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/6d/output/mass_6d_adjust_norm.assoc.txt", head=T, stringsAsFactors=F, data.table=F)
gwas_m[[2]] <- fread("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/42d/output/mass_42d_adjust_norm.assoc.txt", head=T, stringsAsFactors=F, data.table=F)
gwas_m[[3]] <- fread("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/3m/output/mass_3m_adjust_norm.assoc.txt", head=T, stringsAsFactors=F, data.table=F)


names(gwas_m) <- c("mass_6d", "mass_42d", "mass_3m")


# make n file dor each time

n_6d <- sample_size[1:7,]
n_42d <- sample_size[8:14,]
n_3m <- sample_size[15:21,]
n_mass <- sample_size[22:24,]

# make n_obs column for each trait and filter by CR>=0.95

gwas_6d_with_n_obs <- list()
gwas_42d_with_n_obs <- list()
gwas_3m_with_n_obs <- list()

CR <- 0.95

for (n in (1:7)) {

	n_obs_6d <- seq(n_6d[n,2], n_6d[n,2], nrow(gwas_6d[[n]])) - gwas_6d[[n]]$n_miss
	gwas_6d_with_n_obs[[n]] <- cbind(gwas_6d[[n]], n_obs=n_obs_6d)
	gwas_6d_with_n_obs[[n]] <- subset(gwas_6d_with_n_obs[[n]], gwas_6d_with_n_obs[[n]]$n_obs>=(max(gwas_6d_with_n_obs[[n]]$n_obs)*CR), select=1:13)

	n_obs_42d <- seq(n_42d[n,2], n_42d[n,2], nrow(gwas_42d[[n]])) - gwas_42d[[n]]$n_miss
	gwas_42d_with_n_obs[[n]] <- cbind(gwas_42d[[n]], n_obs=n_obs_42d)
	gwas_42d_with_n_obs[[n]] <- subset(gwas_42d_with_n_obs[[n]], gwas_42d_with_n_obs[[n]]$n_obs>=(max(gwas_42d_with_n_obs[[n]]$n_obs)*CR), select=1:13)

	n_obs_3m <- seq(n_3m[n,2], n_3m[n,2], nrow(gwas_3m[[n]])) - gwas_3m[[n]]$n_miss
	gwas_3m_with_n_obs[[n]] <- cbind(gwas_3m[[n]], n_obs=n_obs_3m)
	gwas_3m_with_n_obs[[n]] <- subset(gwas_3m_with_n_obs[[n]], gwas_3m_with_n_obs[[n]]$n_obs>=(max(gwas_3m_with_n_obs[[n]]$n_obs)*CR), select=1:13)

	print(n)

}

gwas_mass_with_n_obs <- list()

gwas_mass_with_n_obs[[1]] <- cbind(gwas_m[[1]], n_obs=(seq(n_mass[1,2], n_mass[1,2], nrow(gwas_m[[1]])) - gwas_m[[1]]$n_miss))
gwas_mass_with_n_obs[[2]] <- cbind(gwas_m[[2]], n_obs=(seq(n_mass[2,2], n_mass[2,2], nrow(gwas_m[[2]])) - gwas_m[[2]]$n_miss))
gwas_mass_with_n_obs[[3]] <- cbind(gwas_m[[3]], n_obs=(seq(n_mass[3,2], n_mass[3,2], nrow(gwas_m[[3]])) - gwas_m[[3]]$n_miss))

gwas_mass_with_n_obs[[1]] <- subset(gwas_mass_with_n_obs[[1]], gwas_mass_with_n_obs[[1]]$n_obs>=(max(gwas_mass_with_n_obs[[1]]$n_obs)*CR), select=1:13)
gwas_mass_with_n_obs[[2]] <- subset(gwas_mass_with_n_obs[[2]], gwas_mass_with_n_obs[[2]]$n_obs>=(max(gwas_mass_with_n_obs[[2]]$n_obs)*CR), select=1:13)
gwas_mass_with_n_obs[[3]] <- subset(gwas_mass_with_n_obs[[3]], gwas_mass_with_n_obs[[3]]$n_obs>=(max(gwas_mass_with_n_obs[[3]]$n_obs)*CR), select=1:13)


# count lambda 

lambda_6d <- matrix(ncol=1, nrow=7)
lambda_42d <- matrix(ncol=1, nrow=7)
lambda_3m <- matrix(ncol=1, nrow=7)
lambda_mass <- matrix(ncol=1, nrow=3)

for (n in (1:7)) {

	lambda_6d[n,1] <- median (qchisq(gwas_6d_with_n_obs[[n]]$p_wald, df=1, lower.tail=F)/qchisq(0.5,df=1))
	lambda_42d[n,1] <- median (qchisq(gwas_42d_with_n_obs[[n]]$p_wald, df=1, lower.tail=F)/qchisq(0.5,df=1))
	lambda_3m[n,1] <- median (qchisq(gwas_3m_with_n_obs[[n]]$p_wald, df=1, lower.tail=F)/qchisq(0.5,df=1))

	print(n)

}

for (n in (1:3)) {

	lambda_mass[n,1] <-  median (qchisq(gwas_mass_with_n_obs[[n]]$p_wald, df=1, lower.tail=F)/qchisq(0.5,df=1))

}

lambda <- cbind(lambda_6d, lambda_42d, lambda_3m)
lambda <- rbind(lambda, t(lambda_mass))

rownames(lambda) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")
colnames(lambda) <- c("6d", "42d", "3m")


# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/with_n")

for (n in 1:7) {

	name_to_save_6d <- paste("index_", n, "_6d_adj_norm_with_n.assoc.txt", sep="")
	fwrite(gwas_6d_with_n_obs[[n]], name_to_save_6d, col.names=T, row.names=F, quote=F, sep=",")

	name_to_save_42d <- paste("index_", n, "_42d_adj_norm_with_n.assoc.txt", sep="")
	fwrite(gwas_42d_with_n_obs[[n]], name_to_save_42d, col.names=T, row.names=F, quote=F, sep=",")

	name_to_save_3m <- paste("index_", n, "_3m_adj_norm_with_n.assoc.txt", sep="")
	fwrite(gwas_3m_with_n_obs[[n]], name_to_save_3m, col.names=T, row.names=F, quote=F, sep=",")

}

fwrite(gwas_mass_with_n_obs[[1]], "mass_6d_adj_norm_with_n.assoc.txt", col.names=T, row.names=F, quote=F, sep=",")
fwrite(gwas_mass_with_n_obs[[2]], "mass_42d_adj_norm_with_n.assoc.txt", col.names=T, row.names=F, quote=F, sep=",")
fwrite(gwas_mass_with_n_obs[[3]], "mass_3m_adj_norm_with_n.assoc.txt", col.names=T, row.names=F, quote=F, sep=",")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results")

write.table(lambda, "lambda_122_animals_before_db.txt", col.names=T, row.names=T, quote=F)