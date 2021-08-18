### add n observed in gwas results ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/gwas_results/output")

pool_1 <- read.table("index_4_6d_adjust_norm_pool_1.assoc.txt", head=T, stringsAsFactors=F)
pool_2 <- read.table("index_4_6d_adjust_norm_pool_2.assoc.txt", head=T, stringsAsFactors=F)

# add n 

pool_1_n <- 86
pool_2_n <- 20

n_obs_pool_1 <-  pool_1_n-pool_1$n_miss
n_obs_pool_2 <-  pool_2_n-pool_2$n_miss

pool_1 <- cbind(pool_1, n_obs=n_obs_pool_1)
pool_2 <- cbind(pool_2, n_obs=n_obs_pool_2)

# save results 

write.table(pool_1, "index_4_6d_adjust_norm_pool_1_with_n.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(pool_2, "index_4_6d_adjust_norm_pool_2_with_n.assoc.txt", col.names=T, row.names=F, quote=F, sep="\t")