### grep mv snps from pc results ###

# load data 

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)

pc <- list()

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/pc/")

for (n in 1:24) {

	name_to_read <- paste("pc_", n, "/pc_", n, "_done.csv", sep="")
	name_to_list <- paste("pc_", n, sep="")
	pc[[n]] <- fread(name_to_read, head=T, stringsAsFactors=F, data.table=F)
	names(pc)[[n]] <- name_to_list

}

# grep mv snps

grep_mv_snps <- list()

for (n in (1:24)) {

	grep_mv_snps[[n]] <- matrix(ncol=1, nrow=12)

	for (i in (c(1:2, 4:12))) {

		grep_mv_snps[[n]][i,1] <- pc[[n]]$p[grepl(mv_snps$rs[i], pc[[n]]$rs_id)]

	}


	print(n)

}

# make matrix with rs and minimal p_value

final_table  <- matrix(ncol=5, nrow=24)

colnames(final_table) <- c("trait", "rs", "chr", "pos", "p")


for (n in (1:24)) {

	final_table[n,1] <- paste("pc_", n, sep="")
	final_table[n,2] <- mv_snps[which.min(grep_mv_snps[[n]]),"rs"]
	final_table[n,3] <- mv_snps[which.min(grep_mv_snps[[n]]),"CHR"]
	final_table[n,4] <- mv_snps[which.min(grep_mv_snps[[n]]),"POS"]
	final_table[n,5] <- min(grep_mv_snps[[n]], na.rm=T)
}

# save result

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results")

write.table(final_table, "mv_snps_grep_on_pc_results.txt", row.names=F, col.names=T, quote=F)