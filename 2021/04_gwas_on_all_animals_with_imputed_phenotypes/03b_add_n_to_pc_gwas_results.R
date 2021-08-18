### add n to pc gwas results ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/pc_on_imputed_pheno/output")

pc_gwas <- list()

for (n in 1:24) {

	name_to_read <- paste("pc_", n, ".assoc.txt", sep="")
	name_to_list <- paste("pc_", n, sep="")
	pc_gwas[[n]] <- fread(name_to_read, head=T, stringsAsFactors=F, data.table=F)
	names(pc_gwas)[[n]] <- name_to_list

}

# add n and filter for CR>=0.95

N=122
CR=0.95

pc_gwas_with_n_obs <- list()

for (n in 1:24) {

	n_obs <- seq(N, N, nrow(pc_gwas[[n]])) - pc_gwas[[n]]$n_miss
	pc_gwas_with_n_obs[[n]] <- cbind(pc_gwas[[n]], n_obs=n_obs)
	pc_gwas_with_n_obs[[n]] <- subset(pc_gwas_with_n_obs[[n]], pc_gwas_with_n_obs[[n]]$n_obs>=(max(pc_gwas_with_n_obs[[n]]$n_obs)*CR), select=1:13)

}

# count lambda

lambda <- matrix(ncol=1, nrow=24)

for(n in 1:24) {

	lambda[n,1] <- median (qchisq(pc_gwas_with_n_obs[[n]]$p_wald, df=1, lower.tail=F)/qchisq(0.5,df=1))

}

# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/animals_with_pheno/pc_on_imputed_pheno/with_n")

for (n in 1:24) {

	name_to_save <- paste("pc_", n, "_with_n.assoc.txt", sep="")
	fwrite(pc_gwas_with_n_obs[[n]], name_to_save, col.names=T, row.names=F, quote=F, sep=",")

}

