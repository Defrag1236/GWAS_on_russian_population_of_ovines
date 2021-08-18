### count lambda after uplaod in db ###

# load data

library(data.table)

gwas_6d <- list()
gwas_42d <- list()
gwas_3m <- list()

for (n in 1:7) {

	name_to_read_6d <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/adj_norm/index_", n, "_6d_adj_norm/index_", n, "_6d_adj_norm_done.csv", sep="")
	name_to_list_6d <- paste("index_", n, "_6d", sep="")
	gwas_6d[[n]] <- fread(name_to_read_6d, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_6d)[[n]] <- name_to_list_6d

	name_to_read_42d <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/adj_norm/index_", n, "_42d_adj_norm/index_", n, "_42d_adj_norm_done.csv", sep="")
	name_to_list_42d <- paste("index_", n, "_42d", sep="")
	gwas_42d[[n]] <- fread(name_to_read_42d, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_42d)[[n]] <- name_to_list_42d

	name_to_read_3m <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/adj_norm/index_", n, "_3m_adj_norm/index_", n, "_3m_adj_norm_done.csv", sep="")
	name_to_list_3m <- paste("index_", n, "_3m", sep="")
	gwas_3m[[n]] <- fread(name_to_read_3m, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_3m)[[n]] <- name_to_list_3m

}

gwas_m <- list()

gwas_m[[1]] <- fread("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/adj_norm/mass_6d_adj_norm/mass_6d_adj_norm_done.csv", head=T, stringsAsFactors=F, data.table=F)
gwas_m[[2]] <- fread("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/adj_norm/mass_42d_adj_norm/mass_42d_adj_norm_done.csv", head=T, stringsAsFactors=F, data.table=F)
gwas_m[[3]] <- fread("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/adj_norm/mass_3m_adj_norm/mass_3m_adj_norm_done.csv", head=T, stringsAsFactors=F, data.table=F)

# count lambda 

lambda_6d <- matrix(ncol=1, nrow=7)
lambda_42d <- matrix(ncol=1, nrow=7)
lambda_3m <- matrix(ncol=1, nrow=7)
lambda_mass <- matrix(ncol=1, nrow=3)

for (n in (1:7)) {

	lambda_6d[n,1] <- median (qchisq(gwas_6d[[n]]$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
	lambda_42d[n,1] <- median (qchisq(gwas_42d[[n]]$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
	lambda_3m[n,1] <- median (qchisq(gwas_3m[[n]]$p, df=1, lower.tail=F)/qchisq(0.5,df=1))

	print(n)

}

for (n in (1:3)) {

	lambda_mass[n,1] <-  median (qchisq(gwas_m[[n]]$p, df=1, lower.tail=F)/qchisq(0.5,df=1))

}

lambda <- cbind(lambda_6d, lambda_42d, lambda_3m)
lambda <- rbind(lambda, t(lambda_mass))

rownames(lambda) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")
colnames(lambda) <- c("6d", "42d", "3m")
 

# save result

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results")

write.table(lambda, "lambda_122_animals_after_db.txt", col.names=T, row.names=T, quote=F)
