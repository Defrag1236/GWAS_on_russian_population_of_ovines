### count lambda for 2nd measurement of phenotypes ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/")

pheno_3m_key <- fread("pheno_3m_key.txt", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/results/gemma_results/3m/output")

trait_list <- list()

for (n in (1:11)) {

	file_name <- paste ("3m_", n, "_trait.assoc.txt", sep="")

	trait_list[[n]] <- fread(file_name, head=T, stringsAsFactors=F)

	}

# make lambda table

lambda_table <- matrix(ncol=2, nrow=11)
lambda_table[,1] <- as.matrix(pheno_3m_key[,2])

for (n in (1:11)) {

	lambda_table[n,2] <- median (qchisq(trait_list[[n]]$p_wald, df=1, lower.tail=F)/qchisq(0.5,df=1))

	}


# save lambda_table

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/results")

fwrite(lambda_table, "lambda_table_3m.txt", col.names=F, row.names=F, quote=F, sep="\t")