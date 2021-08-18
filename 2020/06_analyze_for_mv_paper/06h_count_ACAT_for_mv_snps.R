### count ACAT for mv snps ###

# load data

library("sumFREGAT")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/metal_results")

index_p <- read.table("mv_snps_p.txt", head=T, stringsAsFactors=F)
mass_p <- read.table("mv_snps_mass_p.txt", head=T, stringsAsFactors=F)

index_p <- index_p[-c(4:5),]
mass_p <- mass_p[-c(4:5),]

# count ACAT 

acat_p <- matrix(nrow=10, ncol=8)
colnames(acat_p) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")
rownames(acat_p) <- index_p$SNP


for (n in (1:7)) {

	for (i in (1:10)) {

		x <- index_p[i,c(n+1, n+8, n+15)]

		x <- as.numeric(x)

		y <- ACATO(x)

		acat_p[i,n] <- y

		print(i)
	}

print(n)

}


for (n in (1:10)) {

	x <- mass_p[n,2:4]

	x <- as.numeric(x)

	y <- ACATO(x)

	acat_p[n,8] <- y


}


# save result

write.table(acat_p, "acat_p.txt", col.names=T, row.names=T, quote=F)