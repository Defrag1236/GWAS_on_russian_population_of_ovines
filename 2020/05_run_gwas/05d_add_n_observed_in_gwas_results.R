### add n observed in gwas results ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/pool_1/output")

library(data.table)

pool_1_6d <- list()
pool_1_42d <- list()
pool_1_3m <- list()


for (n in (1:7)) {

	name_to_read_6d <- paste("index_", n, "_6d_pool_1.assoc.txt", sep="")
	name_to_list_6d <- paste("index_", n, "_6d_pool_1", sep="")
	pool_1_6d[[n]] <- fread(name_to_read_6d, head=T, stringsAsFactors=F, data.table=F)
	names(pool_1_6d)[[n]] <- name_to_list_6d

	name_to_read_42d <- paste("index_", n, "_42d_pool_1.assoc.txt", sep="")
	name_to_list_42d <- paste("index_", n, "_42d_pool_1", sep="")
	pool_1_42d[[n]] <- fread(name_to_read_42d, head=T, stringsAsFactors=F, data.table=F)
	names(pool_1_42d)[[n]] <- name_to_list_42d

	name_to_read_3m <- paste("index_", n, "_3m_pool_1.assoc.txt", sep="")
	name_to_list_3m <- paste("index_", n, "_3m_pool_1", sep="")
	pool_1_3m[[n]] <- fread(name_to_read_3m, head=T, stringsAsFactors=F, data.table=F)
	names(pool_1_3m)[[n]] <- name_to_list_3m
 
}

pool_1_m <- list()

pool_1_m[[1]] <- fread("mass_6d_pool_1.assoc.txt", head=T, stringsAsFactors=F, data.table=F)
pool_1_m[[2]] <- fread("mass_42d_pool_1.assoc.txt", head=T, stringsAsFactors=F, data.table=F)
pool_1_m[[3]] <- fread("mass_3m_pool_1.assoc.txt", head=T, stringsAsFactors=F, data.table=F)

names(pool_1_m) <- c("mass_6d", "mass_42d", "mass_3m")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/pool_2/output")

pool_2_6d <- list()
pool_2_42d <- list()
pool_2_3m <- list()


for (n in (1:7)) {

	name_to_read_6d <- paste("index_", n, "_6d_pool_2.assoc.txt", sep="")
	name_to_list_6d <- paste("index_", n, "_6d_pool_2", sep="")
	pool_2_6d[[n]] <- fread(name_to_read_6d, head=T, stringsAsFactors=F, data.table=F)
	names(pool_2_6d)[[n]] <- name_to_list_6d

	name_to_read_42d <- paste("index_", n, "_42d_pool_2.assoc.txt", sep="")
	name_to_list_42d <- paste("index_", n, "_42d_pool_2", sep="")
	pool_2_42d[[n]] <- fread(name_to_read_42d, head=T, stringsAsFactors=F, data.table=F)
	names(pool_2_42d)[[n]] <- name_to_list_42d

	name_to_read_3m <- paste("index_", n, "_3m_pool_2.assoc.txt", sep="")
	name_to_list_3m <- paste("index_", n, "_3m_pool_2", sep="")
	pool_2_3m[[n]] <- fread(name_to_read_3m, head=T, stringsAsFactors=F, data.table=F)
	names(pool_2_3m)[[n]] <- name_to_list_3m
 
}

pool_2_m <- list()

pool_2_m[[1]] <- fread("mass_6d_pool_2.assoc.txt", head=T, stringsAsFactors=F, data.table=F)
pool_2_m[[2]] <- fread("mass_42d_pool_2.assoc.txt", head=T, stringsAsFactors=F, data.table=F)
pool_2_m[[3]] <- fread("mass_3m_pool_2.assoc.txt", head=T, stringsAsFactors=F, data.table=F)

names(pool_2_m) <- c("mass_6d", "mass_42d", "mass_3m")


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pheno_data")

pool_1_n <- read.table("sample_size_pool_1.txt", head=F, stringsAsFactors=F)
pool_2_n <- read.table("sample_size_pool_2.txt", head=F, stringsAsFactors=F)

# make vectors with n for each trait

## pool 1

n_6d_pool_1 <- pool_1_n[1:7,2] 
n_42d_pool_1 <- pool_1_n[8:14,2]  
n_3m_pool_1 <- pool_1_n[15:21,2] 
n_6d_m_pool_1 <- pool_1_n[29,2]
n_42d_m_pool_1 <- pool_1_n[30,2]  
n_3m_m_pool_1 <- pool_1_n[31,2]   

## pool 2

n_6d_pool_2 <- pool_2_n[1:7,2] 
n_42d_pool_2 <- pool_2_n[8:14,2]  
n_3m_pool_2 <- pool_2_n[15:21,2] 
n_6d_m_pool_2 <- pool_2_n[29,2]
n_42d_m_pool_2 <- pool_2_n[30,2]  
n_3m_m_pool_2 <- pool_2_n[31,2] 

# make n_obs column for each trait

## pool 1

pool_1_6d_with_n_obs <- list()
pool_1_42d_with_n_obs <- list()
pool_1_3m_with_n_obs <- list()

for (n in (1:7)) {

	n_obs_6d <- seq(n_6d_pool_1[n], n_6d_pool_1[n], nrow(pool_1_6d[[n]])) - pool_1_6d[[n]]$n_miss
	pool_1_6d_with_n_obs[[n]] <- cbind(pool_1_6d[[n]], n_obs=n_obs_6d)

	n_obs_42d <- seq(n_42d_pool_1[n], n_42d_pool_1[n], nrow(pool_1_42d[[n]])) - pool_1_42d[[n]]$n_miss
	pool_1_42d_with_n_obs[[n]] <- cbind(pool_1_42d[[n]], n_obs=n_obs_42d)

	n_obs_3m <- seq(n_3m_pool_1[n], n_3m_pool_1[n], nrow(pool_1_3m[[n]])) - pool_1_3m[[n]]$n_miss
	pool_1_3m_with_n_obs[[n]] <- cbind(pool_1_3m[[n]], n_obs=n_obs_3m)

	print(n)

}

pool_1_mass_with_n_obs <- list()

pool_1_mass_with_n_obs[[1]] <- cbind(pool_1_m[[1]], n_obs=(seq(n_6d_m_pool_1, n_6d_m_pool_1, nrow(pool_1_m[[1]])) - pool_1_m[[1]]$n_miss))
pool_1_mass_with_n_obs[[2]] <- cbind(pool_1_m[[2]], n_obs=(seq(n_42d_m_pool_1, n_42d_m_pool_1, nrow(pool_1_m[[2]])) - pool_1_m[[2]]$n_miss))
pool_1_mass_with_n_obs[[3]] <- cbind(pool_1_m[[3]], n_obs=(seq(n_3m_m_pool_1, n_3m_m_pool_1, nrow(pool_1_m[[3]])) - pool_1_m[[3]]$n_miss))

## pool 2

pool_2_6d_with_n_obs <- list()
pool_2_42d_with_n_obs <- list()
pool_2_3m_with_n_obs <- list()

for (n in (1:7)) {

	n_obs_6d <- seq(n_6d_pool_2[n], n_6d_pool_2[n], nrow(pool_2_6d[[n]])) - pool_2_6d[[n]]$n_miss
	pool_2_6d_with_n_obs[[n]] <- cbind(pool_2_6d[[n]], n_obs=n_obs_6d)

	n_obs_42d <- seq(n_42d_pool_2[n], n_42d_pool_2[n], nrow(pool_2_42d[[n]])) - pool_2_42d[[n]]$n_miss
	pool_2_42d_with_n_obs[[n]] <- cbind(pool_2_42d[[n]], n_obs=n_obs_42d)

	n_obs_3m <- seq(n_3m_pool_2[n], n_3m_pool_2[n], nrow(pool_2_3m[[n]])) - pool_2_3m[[n]]$n_miss
	pool_2_3m_with_n_obs[[n]] <- cbind(pool_2_3m[[n]], n_obs=n_obs_3m)

	print(n)

}

pool_2_mass_with_n_obs <- list()

pool_2_mass_with_n_obs[[1]] <- cbind(pool_2_m[[1]], n_obs=(seq(n_6d_m_pool_2, n_6d_m_pool_2, nrow(pool_2_m[[1]])) - pool_2_m[[1]]$n_miss))
pool_2_mass_with_n_obs[[2]] <- cbind(pool_2_m[[2]], n_obs=(seq(n_42d_m_pool_2, n_42d_m_pool_2, nrow(pool_2_m[[2]])) - pool_2_m[[2]]$n_miss))
pool_2_mass_with_n_obs[[3]] <- cbind(pool_2_m[[3]], n_obs=(seq(n_3m_m_pool_2, n_3m_m_pool_2, nrow(pool_2_m[[3]])) - pool_2_m[[3]]$n_miss))

# save results 



## pool 1

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/pool_1/pool_1_with_n_obs")

for (n in (1:7)) {

	name_to_save_6d <- paste("index_", n, "_6d_pool_1.assoc.txt", sep="")
	fwrite(pool_1_6d_with_n_obs[[n]], name_to_save_6d, col.names=T, row.names=F, quote=F, sep="\t")

	name_to_save_42d <- paste("index_", n, "_42d_pool_1.assoc.txt", sep="")
	fwrite(pool_1_42d_with_n_obs[[n]], name_to_save_42d, col.names=T, row.names=F, quote=F, sep="\t")  

	name_to_save_3m <- paste("index_", n, "_3m_pool_1.assoc.txt", sep="")
	fwrite(pool_1_3m_with_n_obs[[n]], name_to_save_3m, col.names=T, row.names=F, quote=F, sep="\t")  

}

fwrite(pool_1_mass_with_n_obs[[1]], "mass_6d_pool_1.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")
fwrite(pool_1_mass_with_n_obs[[2]], "mass_42d_pool_1.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")
fwrite(pool_1_mass_with_n_obs[[3]], "mass_3m_pool_1.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")

## pool 2

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas/pool_2/pool_2_with_n_obs")

for (n in (1:7)) {

	name_to_save_6d <- paste("index_", n, "_6d_pool_2.assoc.txt", sep="")
	fwrite(pool_2_6d_with_n_obs[[n]], name_to_save_6d, col.names=T, row.names=F, quote=F, sep="\t")

	name_to_save_42d <- paste("index_", n, "_42d_pool_2.assoc.txt", sep="")
	fwrite(pool_2_42d_with_n_obs[[n]], name_to_save_42d, col.names=T, row.names=F, quote=F, sep="\t")  

	name_to_save_3m <- paste("index_", n, "_3m_pool_2.assoc.txt", sep="")
	fwrite(pool_2_3m_with_n_obs[[n]], name_to_save_3m, col.names=T, row.names=F, quote=F, sep="\t")  

}

fwrite(pool_2_mass_with_n_obs[[1]], "mass_6d_pool_2.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")
fwrite(pool_2_mass_with_n_obs[[2]], "mass_42d_pool_2.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")
fwrite(pool_2_mass_with_n_obs[[3]], "mass_3m_pool_2.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")
