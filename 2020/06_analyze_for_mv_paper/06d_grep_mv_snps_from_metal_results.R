### grep mv snps from metal results ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")


library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)

# grep mv snps (need add _clear to all names (metal_6d_clear etc.))

## 6d

mv_snps_id  <- mv_snps[,1]


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

# grep only p from metal results

mv_snps_6d_p <- matrix(ncol=8, nrow=length(mv_snps_id))
mv_snps_6d_p[,1] <- mv_snps_id

for (n in (1:length(mv_snps_id))) {

	for (i in (1:7)) {

		if (identical(metal_6d_clear[[i]][which(metal_6d_clear[[i]]$MarkerName==mv_snps_id[n]),8], numeric(0))==T)

		{

			mv_snps_6d_p[n, i+1] <- NA


		} else {


		mv_snps_6d_p[n, i+1] <- metal_6d_clear[[i]][grepl(mv_snps_id[n], metal_6d_clear[[i]]$MarkerName),8]
			

		}	

	}
}

mv_snps_42d_p <- matrix(ncol=8, nrow=length(mv_snps_id))
mv_snps_42d_p[,1] <- mv_snps_id

for (n in (1:length(mv_snps_id))) {

	for (i in (1:7)) {

		if (identical(metal_42d_clear[[i]][which(metal_42d_clear[[i]]$MarkerName==mv_snps_id[n]),8], numeric(0))==T)

		{

			mv_snps_42d_p[n, i+1] <- NA


		} else {


		mv_snps_42d_p[n, i+1] <- metal_42d_clear[[i]][grepl(mv_snps_id[n], metal_42d_clear[[i]]$MarkerName),8]
			

		}	

	}
}


mv_snps_3m_p <- matrix(ncol=8, nrow=length(mv_snps_id))
mv_snps_3m_p[,1] <- mv_snps_id

for (n in (1:length(mv_snps_id))) {

	for (i in (1:7)) {

		if (identical(metal_3m_clear[[i]][which(metal_3m_clear[[i]]$MarkerName==mv_snps_id[n]),8], numeric(0))==T)

		{

			mv_snps_3m_p[n, i+1] <- NA


		} else {


		mv_snps_3m_p[n, i+1] <- metal_3m_clear[[i]][grepl(mv_snps_id[n], metal_3m_clear[[i]]$MarkerName),8]
			

		}	

	}
}


mv_snps_p <- cbind(mv_snps_6d_p, mv_snps_42d_p[,-1], mv_snps_3m_p[,-1])
colnames(mv_snps_p) <- c("SNP", "index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d",
						"index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d",
						"index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m")


mv_snps_mass_p <-  matrix(ncol=4, nrow=length(mv_snps_id))
mv_snps_mass_p[,1]  <- mv_snps_id

for (n in (1:length(mv_snps_id))) {

	for (i in (1:3)) {

		if (identical(metal_mass[[i]][which(metal_mass[[i]]$MarkerName==mv_snps_id[n]),8], numeric(0))==T)

		{

			mv_snps_mass_p[n, i+1] <- NA


		} else {


		mv_snps_mass_p[n, i+1] <- metal_mass[[i]][grepl(mv_snps_id[n], metal_mass[[i]]$MarkerName),8]
			

		}	

	}
}

colnames(mv_snps_mass_p) <- c("SNP", "mass_6d", "mass_42d", "mass_3m")

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

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/metal_results")

write.table(mv_snps_p, "mv_snps_p.txt", col.names=T, row.names=F, quote=F)
write.table(mv_snps_mass_p, "mv_snps_mass_p.txt", col.names=T, row.names=F, quote=F)