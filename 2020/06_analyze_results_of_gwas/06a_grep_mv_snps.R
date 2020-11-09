### grep mv snps ###


# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/gwas/pool_1/output")

library(data.table)

pool_1_6d <- list()
pool_1_42d <- list()
pool_1_3m <- list()
pool_1_6m <- list()

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

for (n in (c(1:2,7))) {
	
	name_to_read_6m <- paste("index_", n, "_6m_pool_1.assoc.txt", sep="")
	name_to_list_6m <- paste("index_", n, "_6m_pool_1", sep="")
	pool_1_6m[[n]] <- fread(name_to_read_6m, head=T, stringsAsFactors=F, data.table=F)
	names(pool_1_6m)[[n]] <- name_to_list_6m

}

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)

# grep new mv_snps

grep_table <- matrix(ncol=31, nrow=8)

grep_table[,1] <- mv_snps[-c(8:11),1]
grep_table[,2] <- mv_snps[-c(8:11),2]
grep_table[,3] <- mv_snps[-c(8:11),3]

p_index_6d <- matrix(ncol=7, nrow=8)
p_index_42d <- matrix(ncol=7, nrow=8)
p_index_3m <- matrix(ncol=7, nrow=8)
p_index_6m <- matrix(ncol=7, nrow=8)


for (n in (c(1:7))) {

	for (i in (c(1:8))) {

		if (identical(pool_1_6d[[n]][grepl(grep_table[i,1], pool_1_6d[[n]]$rs), "p_wald"], numeric(0))==T) {

			p_index_6d[i,n] <- NA

		} 

		else {

		p_index_6d[i,n] <- pool_1_6d[[n]][grepl(grep_table[i,1], pool_1_6d[[n]]$rs), "p_wald"]

		}


	}

}

for (n in (c(1:7))) {

	for (i in (c(1:8))) {

		if (identical(pool_1_42d[[n]][grepl(grep_table[i,1], pool_1_42d[[n]]$rs), "p_wald"], numeric(0))==T) {

			p_index_42d[i,n] <- NA

		} 

		else {

		p_index_42d[i,n] <- pool_1_42d[[n]][grepl(grep_table[i,1], pool_1_42d[[n]]$rs), "p_wald"]

		}


	}

}


for (n in (c(1:6))) {

	for (i in (c(1:8))) {

		if (identical(pool_1_3m[[n]][grepl(grep_table[i,1], pool_1_3m[[n]]$rs), "p_wald"], numeric(0))==T) {

			p_index_3m[i,n] <- NA

		} 

		else {

		p_index_3m[i,n] <- pool_1_3m[[n]][grepl(grep_table[i,1], pool_1_3m[[n]]$rs), "p_wald"]

		}


	}

}


for (n in (c(1:2, 7))) {

	for (i in (c(1:8))) {

		if (identical(pool_1_6m[[n]][grepl(grep_table[i,1], pool_1_6m[[n]]$rs), "p_wald"], numeric(0))==T) {

			p_index_6m[i,n] <- NA

		} 

		else {

		p_index_6m[i,n] <- pool_1_6m[[n]][grepl(grep_table[i,1], pool_1_6m[[n]]$rs), "p_wald"]

		}


	}

}

grep_table[,4:10] <- p_index_6d
grep_table[,11:17] <- p_index_42d
grep_table[,18:24] <- p_index_3m
grep_table[,25:31] <- p_index_6m

colnames(grep_table) <- c("SNP", "Chr", "Pos", "index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d",
							"index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d",
							"index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m",
							"index_1_6m", "index_2_6m", "index_3_6m", "index_4_6m", "index_5_6m", "index_6_6m", "index_7_6m")


# save result

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

write.table(grep_table, "grep_mv_snps_on_gwas_data.txt", col.names=T, row.names=F, quote=F)
