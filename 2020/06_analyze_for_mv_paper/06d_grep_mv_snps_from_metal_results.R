### grep mv snps from metal results ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_21_traits.Rdata")

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)

# grep mv snps

## 6d

mv_snps_id  <- mv_snps[c(1:7, 12),1]


mv_snps_6d_index_1 <- list()
mv_snps_6d_index_2 <- list()
mv_snps_6d_index_3 <- list()
mv_snps_6d_index_4 <- list()
mv_snps_6d_index_5 <- list()
mv_snps_6d_index_6 <- list()
mv_snps_6d_index_7 <- list()

	
for (n in (1:8)) {

	mv_snps_6d_index_1[[n]] <- metal_6d[[1]][grepl(mv_snps_id[n], metal_6d[[1]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_6d_index_2[[n]] <- metal_6d[[2]][grepl(mv_snps_id[n], metal_6d[[2]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_6d_index_3[[n]] <- metal_6d[[3]][grepl(mv_snps_id[n], metal_6d[[3]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_6d_index_4[[n]] <- metal_6d[[4]][grepl(mv_snps_id[n], metal_6d[[4]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_6d_index_5[[n]] <- metal_6d[[5]][grepl(mv_snps_id[n], metal_6d[[5]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_6d_index_6[[n]] <- metal_6d[[6]][grepl(mv_snps_id[n], metal_6d[[6]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_6d_index_7[[n]] <- metal_6d[[7]][grepl(mv_snps_id[n], metal_6d[[7]]$MarkerName),]

}

mv_snps_6d_index_1 <- do.call(rbind, mv_snps_6d_index_1)
mv_snps_6d_index_2 <- do.call(rbind, mv_snps_6d_index_2)
mv_snps_6d_index_3 <- do.call(rbind, mv_snps_6d_index_3)
mv_snps_6d_index_4 <- do.call(rbind, mv_snps_6d_index_4)
mv_snps_6d_index_5 <- do.call(rbind, mv_snps_6d_index_5)
mv_snps_6d_index_6 <- do.call(rbind, mv_snps_6d_index_6)
mv_snps_6d_index_7 <- do.call(rbind, mv_snps_6d_index_7)

## 42d 

mv_snps_42d_index_1 <- list()
mv_snps_42d_index_2 <- list()
mv_snps_42d_index_3 <- list()
mv_snps_42d_index_4 <- list()
mv_snps_42d_index_5 <- list()
mv_snps_42d_index_6 <- list()
mv_snps_42d_index_7 <- list()

	
for (n in (1:8)) {

	mv_snps_42d_index_1[[n]] <- metal_42d[[1]][grepl(mv_snps_id[n], metal_42d[[1]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_42d_index_2[[n]] <- metal_42d[[2]][grepl(mv_snps_id[n], metal_42d[[2]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_42d_index_3[[n]] <- metal_42d[[3]][grepl(mv_snps_id[n], metal_42d[[3]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_42d_index_4[[n]] <- metal_42d[[4]][grepl(mv_snps_id[n], metal_42d[[4]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_42d_index_5[[n]] <- metal_42d[[5]][grepl(mv_snps_id[n], metal_42d[[5]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_42d_index_6[[n]] <- metal_42d[[6]][grepl(mv_snps_id[n], metal_42d[[6]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_42d_index_7[[n]] <- metal_42d[[7]][grepl(mv_snps_id[n], metal_42d[[7]]$MarkerName),]

}

mv_snps_42d_index_1 <- do.call(rbind, mv_snps_42d_index_1)
mv_snps_42d_index_2 <- do.call(rbind, mv_snps_42d_index_2)
mv_snps_42d_index_3 <- do.call(rbind, mv_snps_42d_index_3)
mv_snps_42d_index_4 <- do.call(rbind, mv_snps_42d_index_4)
mv_snps_42d_index_5 <- do.call(rbind, mv_snps_42d_index_5)
mv_snps_42d_index_6 <- do.call(rbind, mv_snps_42d_index_6)
mv_snps_42d_index_7 <- do.call(rbind, mv_snps_42d_index_7)


## 3m

mv_snps_id  <- mv_snps[c(1:7, 12),1]


mv_snps_3m_index_1 <- list()
mv_snps_3m_index_2 <- list()
mv_snps_3m_index_3 <- list()
mv_snps_3m_index_4 <- list()
mv_snps_3m_index_5 <- list()
mv_snps_3m_index_6 <- list()
mv_snps_3m_index_7 <- list()

	
for (n in (1:8)) {

	mv_snps_3m_index_1[[n]] <- metal_3m[[1]][grepl(mv_snps_id[n], metal_3m[[1]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_3m_index_2[[n]] <- metal_3m[[2]][grepl(mv_snps_id[n], metal_3m[[2]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_3m_index_3[[n]] <- metal_3m[[3]][grepl(mv_snps_id[n], metal_3m[[3]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_3m_index_4[[n]] <- metal_3m[[4]][grepl(mv_snps_id[n], metal_3m[[4]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_3m_index_5[[n]] <- metal_3m[[5]][grepl(mv_snps_id[n], metal_3m[[5]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_3m_index_6[[n]] <- metal_3m[[6]][grepl(mv_snps_id[n], metal_3m[[6]]$MarkerName),]

}

for (n in (1:8)) {

	mv_snps_3m_index_7[[n]] <- metal_3m[[7]][grepl(mv_snps_id[n], metal_3m[[7]]$MarkerName),]

}

mv_snps_3m_index_1 <- do.call(rbind, mv_snps_3m_index_1)
mv_snps_3m_index_2 <- do.call(rbind, mv_snps_3m_index_2)
mv_snps_3m_index_3 <- do.call(rbind, mv_snps_3m_index_3)
mv_snps_3m_index_4 <- do.call(rbind, mv_snps_3m_index_4)
mv_snps_3m_index_5 <- do.call(rbind, mv_snps_3m_index_5)
mv_snps_3m_index_6 <- do.call(rbind, mv_snps_3m_index_6)
mv_snps_3m_index_7 <- do.call(rbind, mv_snps_3m_index_7)

# make z_score matrix for mv snps for mv_analysis




# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/metal_results/6d")

write.table(mv_snps_6d_index_1, "mv_snps_6d_index_1.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_6d_index_2, "mv_snps_6d_index_2.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_6d_index_3, "mv_snps_6d_index_3.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_6d_index_4, "mv_snps_6d_index_4.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_6d_index_5, "mv_snps_6d_index_5.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_6d_index_6, "mv_snps_6d_index_6.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_6d_index_7, "mv_snps_6d_index_7.txt", col.names=T, row.names=F, quote=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/metal_results/42d")

write.table(mv_snps_42d_index_1, "mv_snps_42d_index_1.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_42d_index_2, "mv_snps_42d_index_2.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_42d_index_3, "mv_snps_42d_index_3.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_42d_index_4, "mv_snps_42d_index_4.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_42d_index_5, "mv_snps_42d_index_5.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_42d_index_6, "mv_snps_42d_index_6.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_42d_index_7, "mv_snps_42d_index_7.txt", col.names=T, row.names=F, quote=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/metal_results/3m")

write.table(mv_snps_3m_index_1, "mv_snps_3m_index_1.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_3m_index_2, "mv_snps_3m_index_2.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_3m_index_3, "mv_snps_3m_index_3.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_3m_index_4, "mv_snps_3m_index_4.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_3m_index_5, "mv_snps_3m_index_5.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_3m_index_6, "mv_snps_3m_index_6.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_3m_index_7, "mv_snps_3m_index_7.txt", col.names=T, row.names=F, quote=F)

