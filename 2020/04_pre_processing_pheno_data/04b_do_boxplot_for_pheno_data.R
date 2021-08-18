### do boxplot for pheno data ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

library(data.table)

pheno <- fread("pheno_2020.csv", head=T, stringsAsFactors=F, sep=";", na.strings=c("","NA"), data.table=F)

pheno[,-c(1:4)] <-  apply(pheno[,-c(1:4)], 2, as.numeric)

rownames(pheno) <- pheno[,1]

pheno <- pheno[,c(1:25)]

# count na for each animal

trait_count <- matrix(ncol=1, nrow=nrow(pheno))

for (n in (1:nrow(pheno))) {

	trait_count[n,1] <- table(is.na(pheno[n, -c(1:4)]))[1]

}

## remove values which out +-5 standart deviations

pheno_qc <- pheno[-c(1:4)]

for (n in (1:ncol(pheno_qc))) {

	pheno_qc[which(pheno_qc[,n]>mean(pheno_qc[,n], na.rm=T)+(sd(pheno_qc[,n], na.rm=T)*5)),n] <- NA
	pheno_qc[which(pheno_qc[,n]<mean(pheno_qc[,n], na.rm=T)-(sd(pheno_qc[,n], na.rm=T)*5)),n] <- NA

}


# make linear plot for all animals

index_1 <- pheno_qc[,c(1,8,15)]
index_1 <- index_1[complete.cases(index_1),]

index_2 <- pheno_qc[,c(2,9,16)]
index_2 <- index_2[complete.cases(index_2),]

index_3 <- pheno_qc[,c(3,10,17)]
index_3 <- index_3[complete.cases(index_3),]

index_4 <- pheno_qc[,c(4,11,18)]
index_4 <- index_4[complete.cases(index_4),]

index_5 <- pheno_qc[,c(5,12,19)]
index_5 <- index_5[complete.cases(index_5),]

index_6 <- pheno_qc[,c(6,13,20)]
index_6 <- index_6[complete.cases(index_6),]

index_7 <- pheno_qc[,c(7,14,2)]
index_7 <- index_7[complete.cases(index_7),]

for_plot_linear <- list(index_1, index_2, index_3, index_4, index_5, index_6, index_7)

pdf("linear_plot_for_all_animals.pdf", width=16, height=10)

for (n in (1:7)) {

	plot(x=c(6, 90),y=c(min(for_plot_linear[[n]]), max(for_plot_linear[[n]])))

	 for (i in (1:nrow(for_plot_linear[[n]]))) {

	 	points(x=c(6,42,90), y=for_plot_linear[[n]][i,], type="l", col="red")

	 }

}

dev.off()


# make boxplot

index_1_long <- t(index_1)
t_index_1 <- c(seq(6, 6, length.out=(ncol(index_1_long))),  seq(42, 42, length.out=(ncol(index_1_long))), seq(90, 90, length.out=(ncol(index_1_long))))
index_1_long <- c(index_1_long[1,], index_1_long[2,], index_1_long[3,])

index_2_long <- t(index_2)
t_index_2 <- c(seq(6, 6, length.out=(ncol(index_2_long))),  seq(42, 42, length.out=(ncol(index_2_long))), seq(90, 90, length.out=(ncol(index_2_long))))
index_2_long <- c(index_2_long[1,], index_2_long[2,], index_2_long[3,])

index_3_long <- t(index_3)
t_index_3 <- c(seq(6, 6, length.out=(ncol(index_3_long))),  seq(42, 42, length.out=(ncol(index_3_long))), seq(90, 90, length.out=(ncol(index_3_long))))
index_3_long <- c(index_3_long[1,], index_3_long[2,], index_3_long[3,])

index_4_long <- t(index_4)
t_index_4 <- c(seq(6, 6, length.out=(ncol(index_4_long))),  seq(42, 42, length.out=(ncol(index_4_long))), seq(90, 90, length.out=(ncol(index_4_long))))
index_4_long <- c(index_4_long[1,], index_4_long[2,], index_4_long[3,])

index_5_long <- t(index_5)
t_index_5 <- c(seq(6, 6, length.out=(ncol(index_5_long))),  seq(42, 42, length.out=(ncol(index_5_long))), seq(90, 90, length.out=(ncol(index_5_long))))
index_5_long <- c(index_5_long[1,], index_5_long[2,], index_5_long[3,])

index_6_long <- t(index_6)
t_index_6 <- c(seq(6, 6, length.out=(ncol(index_6_long))),  seq(42, 42, length.out=(ncol(index_6_long))), seq(90, 90, length.out=(ncol(index_6_long))))
index_6_long <- c(index_6_long[1,], index_6_long[2,], index_6_long[3,])

index_7_long <- t(index_7)
t_index_7 <- c(seq(6, 6, length.out=(ncol(index_7_long))),  seq(42, 42, length.out=(ncol(index_7_long))), seq(90, 90, length.out=(ncol(index_7_long))))
index_7_long <- c(index_7_long[1,], index_7_long[2,], index_7_long[3,])

for_boxplot <- list(index_1_long, index_2_long, index_3_long, index_4_long, index_5_long, index_6_long, index_7_long)
for_boxplot_t <- list(t_index_1, t_index_2, t_index_3, t_index_4, t_index_5, t_index_6, t_index_7)


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("boxplot_for_all_animals.pdf", width=16, height=10)

for (n in (1:7)) {

	boxplot(for_boxplot[[n]]~for_boxplot_t[[n]])

}

dev.off()

# make regression plot

pdf("regression_plot.pdf", width=16, height=10)

for (n in (1:7)) {

	plot(for_boxplot_t[[n]], for_boxplot[[n]])

	abline(lm(for_boxplot[[n]]~for_boxplot_t[[n]]))

}

dev.off()


# count sd and mean and cor

mean_and_sd <- matrix(ncol=2, nrow=7)

colnames(mean_and_sd) <- c("mean", "sd")
rownames(mean_and_sd) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7")

for (n in (1:7)) {

	mean_and_sd[n,1] <- mean(for_boxplot[[n]])
	mean_and_sd[n,2] <- sd(for_boxplot[[n]])

}


cor_index <- matrix(ncol=1, nrow=7)

rownames(cor_index) <-  c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7")
colnames(cor_index) <- c("correlation")

for (n in (1:7)) {

	cor_index[n,1] <- cor(for_boxplot_t[[n]], for_boxplot[[n]])

}

# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

fwrite(mean_and_sd, "mean_and_sd_for_animals_with_all_measurements_for_1_index.txt", col.names=T, row.names=T, quote=F, sep="\t")
fwrite(as.data.frame(cor_index), "correlation_between_t_and_index.txt", col.names=T, row.names=T, quote=F, sep="\t")