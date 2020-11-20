### count lambda and do corr plot for metal results ###

# load data 

library(data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/6d")

metal_6d <- list()

for (n in (1:7)) {

	name <- paste("index_", n, "_6d.TBL", sep="")
	metal_6d[[n]] <- fread(name, head=T, stringsAsFactors=F, data.table=F)

}

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d")

metal_42d <- list()

for (n in (1:7)) {

	name <- paste("index_", n, "_42d.TBL", sep="")
	metal_42d[[n]] <- fread(name, head=T, stringsAsFactors=F, data.table=F)

}

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m")

metal_3m <- list()

for (n in (1:7)) {

	name <- paste("index_", n, "_3m.TBL", sep="")
	metal_3m[[n]] <- fread(name, head=T, stringsAsFactors=F, data.table=F)

}


# count lambda for each trait

lambda_6d <- matrix(ncol=1, nrow=7)
lambda_42d <- matrix(ncol=1, nrow=7)
lambda_3m <- matrix(ncol=1, nrow=7)

for (n in (1:7)) {

	lambda_6d[n,1] <- median (qchisq(metal_6d[[n]][,8], df=1, lower.tail=F)/qchisq(0.5,df=1))
	lambda_42d[n,1] <- median (qchisq(metal_42d[[n]][,8], df=1, lower.tail=F)/qchisq(0.5,df=1))
	lambda_3m[n,1] <- median (qchisq(as.numeric(metal_3m[[n]][,8]), df=1, lower.tail=F)/qchisq(0.5,df=1))

	print(n)

}

lambda_metal <- cbind(lambda_6d, lambda_42d, lambda_3m)

rownames(lambda_metal) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7")
colnames(lambda_metal) <- c("6d", "42d", "3m")


# make common matrix with z_score

## calculate z-score

for (n in (1:7)) {

	metal_6d[[n]]$z_score <- metal_6d[[n]]$Effect/metal_6d[[n]]$StdErr
	metal_42d[[n]]$z_score <- metal_42d[[n]]$Effect/metal_42d[[n]]$StdErr
	metal_3m[[n]]$z_score <- metal_3m[[n]]$Effect/metal_3m[[n]]$StdErr

}

## found interesct snps for all traits

snps_6d <- list()
snps_42d <- list()
snps_3m <- list()

for (n in(1:7)) {


	snps_6d[[n]] <- metal_6d[[n]]$MarkerName
	snps_42d[[n]] <- metal_42d[[n]]$MarkerName
	snps_3m[[n]] <- metal_3m[[n]]$MarkerName

	print(n)

}


snps_clear_6d <- Reduce(intersect, snps_6d) 
snps_clear_42d <- Reduce(intersect, snps_42d) 
snps_clear_3m <- Reduce(intersect, snps_3m) 

snps_clear <- Reduce(intersect, list(snps_clear_6d, snps_clear_42d, snps_clear_3m))

## make matrix with z-score for 21 traits

z_6d <- list()
z_42d <- list()
z_3m <- list()

for (n in (1:7)) {

	z_6d[[n]] <- metal_6d[[n]][(metal_6d[[n]]$MarkerName %in% snps_clear),10]
	z_42d[[n]] <- metal_42d[[n]][(metal_42d[[n]]$MarkerName %in% snps_clear),10]
	z_3m[[n]] <- metal_3m[[n]][(metal_3m[[n]]$MarkerName %in% snps_clear),10]

	print(n)

}

z_6d <- do.call(cbind, z_6d)
z_42d <- do.call(cbind, z_42d)
z_3m <- do.call(cbind, z_3m)

colnames(z_6d) <- c("index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d")
colnames(z_42d) <- c("index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d")
colnames(z_3m) <- c("index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m")

z_all <- cbind(z_6d, z_42d, z_3m)
z_all[abs(z_all)>=2] <- NA
z_all <- na.omit(z_all)

# make corr plot for 21 traits

z_corr <- cor(z_all)

library(corrplot)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("corrplot_for_21_traits_after_metal.pdf", width=16, height=10)

corrplot(z_corr, tl.cex=0.8)

dev.off()

# make results unified by snps 

metal_6d_clear <- list()
metal_42d_clear <- list()
metal_3m_clear <- list()


for (n in (1:7)) {

	metal_6d_clear[[n]] <- metal_6d[[n]][(metal_6d[[n]]$MarkerName %in% snps_clear),]
	metal_42d_clear[[n]] <- metal_42d[[n]][(metal_42d[[n]]$MarkerName %in% snps_clear),]
	metal_3m_clear[[n]] <- metal_3m[[n]][(metal_3m[[n]]$MarkerName %in% snps_clear),]

	print(n)

}

# save results

## save lambda file

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

write.table(lambda_metal, "lambda_meta_analysis_21_trait.txt", row.names=T, col.names=T, quote=F)

## save metal results ad .Rdata

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

save(metal_6d, metal_42d, metal_3m, file="metal_results_21_traits.Rdata")

## save metal results unified by snps as .Rdata

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

save(metal_6d_clear, metal_42d_clear, metal_3m_clear, file="metal_results_21_traits_unified.Rdata")