### grep mv snps from mv results ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("multi_6d.Rdata")
load("multi_42d.Rdata")
load("multi_3m.Rdata")
load("multi_all.Rdata")

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens")
 
p_6d <- fread("p_value_6d_multi_stephens.txt", head=F, stringsAsFactors=F, data.table=F)
p_42d <- fread("p_value_42d_multi_stephens.txt", head=F, stringsAsFactors=F, data.table=F)
p_3m <- fread("p_value_3m_multi_stephens.txt", head=F, stringsAsFactors=F, data.table=F)
p_all <- fread("p_value_all_multi_stephens.txt", head=F, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered")

index_1 <- read.table("p_value_index_1_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_2 <- read.table("p_value_index_2_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_3 <- read.table("p_value_index_3_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_4 <- read.table("p_value_index_4_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_5 <- read.table("p_value_index_5_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_6 <- read.table("p_value_index_6_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
index_7 <- read.table("p_value_index_7_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)
mass <- read.table("p_value_mass_multi_stephens_maf_0.05_filtered.txt", head=F, stringsAsFactors=F)

# count lambda for mv multiabel

median (qchisq(multi_6d$scan$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
median (qchisq(multi_42d$scan$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
median (qchisq(multi_3m$scan$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
median (qchisq(multi_all$scan$p, df=1, lower.tail=F)/qchisq(0.5,df=1))

# grep mv snps from multiabel

mv_6d <- list()
mv_42d <- list()
mv_3m <- list()
mv_all <- list()

for (n in(1:nrow(mv_snps))) {

	mv_6d[[n]] <- multi_6d$scan[grepl(mv_snps[n,1], multi_6d$scan$marker), 1:6]
	mv_42d[[n]] <- multi_42d$scan[grepl(mv_snps[n,1], multi_42d$scan$marker), 1:6] 
	mv_3m[[n]] <- multi_3m$scan[grepl(mv_snps[n,1], multi_3m$scan$marker), 1:6] 
	mv_all[[n]] <- multi_all$scan[grepl(mv_snps[n,1], multi_all$scan$marker), 1:6]  

}

mv_6d <- do.call(rbind, mv_6d)
mv_42d <- do.call(rbind, mv_3m)
mv_3m <- do.call(rbind, mv_3m)
mv_all <- do.call(rbind, mv_all)

# grep mv snps from stephens

## for rime traits

mv_6d_s <- list()
mv_42d_s <- list()
mv_3m_s <- list()
mv_all_s <- list()

for (n in(1:nrow(mv_snps))) {

	mv_6d_s[[n]] <- p_6d[grepl(mv_snps[n,1], p_6d[,1]),]
	mv_42d_s[[n]] <- p_42d[grepl(mv_snps[n,1], p_42d[,1]),]
	mv_3m_s[[n]] <- p_3m[grepl(mv_snps[n,1], p_3m[,1]),]
	mv_all_s[[n]] <- p_all[grepl(mv_snps[n,1], p_all[,1]),] 

}

mv_6d_s <- do.call(rbind, mv_6d_s)
mv_42d_s <- do.call(rbind, mv_42d_s)
mv_3m_s <- do.call(rbind, mv_3m_s)
mv_all_s <- do.call(rbind, mv_all_s)

## for index traits 

mv_index_1 <- list()
mv_index_2 <- list()
mv_index_3 <- list()
mv_index_4 <- list()
mv_index_5 <- list()
mv_index_6 <- list()
mv_index_7 <- list()
mv_mass <- list()

for (n in(1:nrow(mv_snps))) {

	mv_index_1[[n]] <- index_1[grepl(mv_snps[n,1], index_1[,1]),]
	mv_index_2[[n]] <- index_2[grepl(mv_snps[n,1], index_2[,1]),]
	mv_index_3[[n]] <- index_3[grepl(mv_snps[n,1], index_3[,1]),]
	mv_index_4[[n]] <- index_4[grepl(mv_snps[n,1], index_4[,1]),]
	mv_index_5[[n]] <- index_5[grepl(mv_snps[n,1], index_5[,1]),]
	mv_index_6[[n]] <- index_6[grepl(mv_snps[n,1], index_6[,1]),]
	mv_index_7[[n]] <- index_7[grepl(mv_snps[n,1], index_7[,1]),]
	mv_mass[[n]] <- mass[grepl(mv_snps[n,1], mass[,1]),]

	print (n)

}



mv_index_1 <- do.call(rbind, mv_index_1)
mv_index_2 <- do.call(rbind, mv_index_2)
mv_index_3 <- do.call(rbind, mv_index_3)
mv_index_4 <- do.call(rbind, mv_index_4)
mv_index_5 <- do.call(rbind, mv_index_5)
mv_index_6 <- do.call(rbind, mv_index_6)
mv_index_7 <- do.call(rbind, mv_index_7)
mv_mass <- do.call(rbind, mv_mass)

mv_index <- cbind(SNP=mv_index_1[,1], index_1=mv_index_1[,2], index_2=mv_index_2[,2], index_3=mv_index_3[,2], index_4=mv_index_4[,2], index_5=mv_index_5[,2], index_6=mv_index_6[,2], index_7=mv_index_7[,2], mass=mv_mass[,2])

# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_results")

fwrite(mv_6d, "mv_6d_mv_snps_grep.txt", col.names=T, row.names=F, quote=F)
fwrite(mv_42d, "mv_42d_mv_snps_grep.txt", col.names=T, row.names=F, quote=F)
fwrite(mv_3m, "mv_3m_mv_snps_grep.txt", col.names=T, row.names=F, quote=F)
fwrite(mv_all, "mv_all_mv_snps_grep.txt", col.names=T, row.names=F, quote=F)


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_results")

fwrite(mv_6d_s, "mv_6d_stephens_mv_snps_grep_p.txt", col.names=F, row.names=F, quote=F)
fwrite(mv_42d_s, "mv_42d_stephens_mv_snps_grep_p.txt", col.names=F, row.names=F, quote=F)
fwrite(mv_3m_s, "mv_3m_stephens_mv_snps_grep_p.txt", col.names=F, row.names=F, quote=F)
fwrite(mv_all_s, "mv_all_stephens_mv_snps_grep_p.txt", col.names=F, row.names=F, quote=F)

fwrite(mv_index, "mv_index_mv_snps_grep_p.txt", col.names=T, row.names=F, quote=F)