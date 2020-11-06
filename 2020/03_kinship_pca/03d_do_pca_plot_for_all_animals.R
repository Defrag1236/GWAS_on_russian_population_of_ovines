### do pca plot for all animals ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pca")

pca <- read.table("pca_192_animals_filtered.eigenvec", head=F, stringsAsFactors=F)
pca_96 <- read.table("pca_96_animals_filtered.eigenvec", head=F, stringsAsFactors=F)
pca_qc <- read.table("pca_rom_rom_plus_arg_kat_arg_filtered.eigenvec", head=F, stringsAsFactors=F)
pca_with_pheno <- read.table("pca_animals_with_pheno_filtered.eigenvec", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

pheno_6d <- read.csv("pheno_index_6d_2020.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

setwd("/home/common/projects/ovine_selection/ovines_gwas_map/Data")

breed_info <- read.csv("sex_for_1st_48_sheeps.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# make for_plot matrix

for_plot <- pca[,2:4]

breed <- c(1:192)

for_plot <- cbind(for_plot, breed)

for_plot$breed <- NA

for_plot[(for_plot[,1] %in% breed_info[,1]),4] <- breed_info[match(for_plot[(for_plot[,1] %in% breed_info[,1]),1],breed_info[,1]),3]
for_plot[(for_plot[,1] %in% pheno_6d[,1]),4] <- pheno_6d[match(for_plot[(for_plot[,1] %in% pheno_6d[,1]),1],pheno_6d[,1]),4]


for_plot[1:48,4] <- "rom_arg"

for_plot[is.na(for_plot[,4]),4] <- "Unknown"

for_plot[grepl("^Katahdin$", for_plot$breed),4] <- "katahdin" 

# make plot 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("pca_plot.pdf", width=16, height=10)

plot(for_plot[,2], for_plot[,3], cex=1.8)

points(for_plot[grepl("Argali", for_plot$breed),2], for_plot[grepl("Argali", for_plot$breed),3], col="red", pch=19)
points(for_plot[grepl("katahdin", for_plot$breed),2], for_plot[grepl("katahdin", for_plot$breed),3], col="orange", pch=19)
points(for_plot[grepl("F1- KAT/ROM", for_plot$breed),2], for_plot[grepl("F1- KAT/ROM", for_plot$breed),3], col="yellow", pch=19)
points(for_plot[grepl("rom_arg", for_plot$breed),2], for_plot[grepl("rom_arg", for_plot$breed),3], col="green", pch=19)
points(for_plot[grepl("rom_arg_kat", for_plot$breed),2], for_plot[grepl("rom_arg_kat", for_plot$breed),3], col="blue", pch=19)
points(for_plot[grepl("rom_arg_kat_mouf", for_plot$breed),2], for_plot[grepl("rom_arg_kat_mouf", for_plot$breed),3], col="pink", pch=19)
points(for_plot[grepl("rom_mouf_kat", for_plot$breed),2], for_plot[grepl("rom_mouf_kat", for_plot$breed),3], col="brown", pch=19)
points(for_plot[grepl("Romanov", for_plot$breed),2], for_plot[grepl("Romanov", for_plot$breed),3], col="grey", pch=19)
points(for_plot[grepl("Unknown", for_plot$breed),2], for_plot[grepl("Unknown", for_plot$breed),3], col="black", pch=19)


legend("topleft", legend=c("Argali", "katahdin", "F1- KAT/ROM", "rom_arg", "rom_arg_kat", "rom_arg_kat_mouf", "rom_mouf_kat", "Romanov", "Unknown"), 
	col=c("red", "orange", "yellow", "green", "blue", "pink", "brown", "grey", "black"), cex=1.8, pch=19)

dev.off()


# do pca plot for new 96 (new) animals

## make for_plot_96 matrix

for_plot_96 <- pca_96[,2:4]

breed <- c(1:96)

for_plot_96 <- cbind(for_plot_96, breed)

for_plot_96$breed <- NA

for_plot_96[(for_plot_96[,1] %in% pheno_6d[,1]),4] <- pheno_6d[match(for_plot_96[(for_plot_96[,1] %in% pheno_6d[,1]),1],pheno_6d[,1]),4]

for_plot_96[is.na(for_plot_96[,4]),4] <- "Unknown"

for_plot_96[grepl("^Katahdin$", for_plot_96$breed),4] <- "katahdin" 

# make plot fo 96 new (2020) animals

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("pca_plot_96_animals_2020.pdf", width=16, height=10)

plot(for_plot_96[,2], for_plot_96[,3], cex=1.8)

points(for_plot_96[grepl("katahdin", for_plot_96$breed),2], for_plot_96[grepl("katahdin", for_plot_96$breed),3], col="orange", pch=19)
points(for_plot_96[grepl("rom_arg", for_plot_96$breed),2], for_plot_96[grepl("rom_arg", for_plot_96$breed),3], col="green", pch=19)
points(for_plot_96[grepl("rom_arg_kat", for_plot_96$breed),2], for_plot_96[grepl("rom_arg_kat", for_plot_96$breed),3], col="blue", pch=19)
points(for_plot_96[grepl("rom_arg_kat_mouf", for_plot_96$breed),2], for_plot_96[grepl("rom_arg_kat_mouf", for_plot_96$breed),3], col="pink", pch=19)
points(for_plot_96[grepl("rom_mouf_kat", for_plot_96$breed),2], for_plot_96[grepl("rom_mouf_kat", for_plot_96$breed),3], col="brown", pch=19)
points(for_plot_96[grepl("Unknown", for_plot_96$breed),2], for_plot_96[grepl("Unknown", for_plot_96$breed),3], col="black", pch=19)


legend("topleft", legend=c("katahdin", "rom_arg", "rom_arg_kat", "rom_arg_kat_mouf", "rom_mouf_kat", "Unknown"), 
	col=c("orange", "green", "blue", "pink", "brown", "black"), cex=1.8, pch=19)

dev.off()


# pca plot for rom, rom+arg, kat and argali

## make for_plot_qc matrix


for_plot_qc <- pca_qc[,2:4]

breed <- c(1:nrow(for_plot_qc))

for_plot_qc <- cbind(for_plot_qc, breed)


for_plot_qc[(for_plot_qc[,1] %in% breed_info[,1]),4] <- breed_info[match(for_plot_qc[(for_plot_qc[,1] %in% breed_info[,1]),1],breed_info[,1]),3]
for_plot_qc[(for_plot_qc[,1] %in% pheno_6d[,1]),4] <- pheno_6d[match(for_plot_qc[(for_plot_qc[,1] %in% pheno_6d[,1]),1],pheno_6d[,1]),4]


for_plot_qc[1:48,4] <- "rom_arg"

for_plot_qc[grepl("^Katahdin$", for_plot_qc$breed),4] <- "katahdin" 


new_samples <- for_plot_qc[(grepl("sample", for_plot_qc[,1])),]

## make plot rom, rom+arg, kat and argali

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("pca_plot_rom_rom_plus_arg_kat_arg.pdf", width=16, height=10)

plot(for_plot_qc[,2], for_plot_qc[,3], cex=1.8)

points(for_plot_qc[grepl("Argali", for_plot_qc$breed),2], for_plot_qc[grepl("Argali", for_plot_qc$breed),3], col="red", pch=19)
points(for_plot_qc[(grepl("katahdin", for_plot_qc$breed) & grepl("sample", for_plot_qc[,1])),2], for_plot_qc[(grepl("katahdin", for_plot_qc$breed) & grepl("sample", for_plot_qc[,1])),3], col="orange", pch=17)
points(for_plot_qc[(grepl("katahdin", for_plot_qc$breed) & !grepl("sample", for_plot_qc[,1])),2], for_plot_qc[(grepl("katahdin", for_plot_qc$breed) & !grepl("sample", for_plot_qc[,1])),3], col="orange", pch=19)
points(for_plot_qc[(grepl("rom_arg", for_plot_qc$breed) & grepl("sample", for_plot_qc[,1])),2], for_plot_qc[(grepl("rom_arg", for_plot_qc$breed) & grepl("sample", for_plot_qc[,1])),3], col="green", pch=17)
points(for_plot_qc[(grepl("rom_arg", for_plot_qc$breed) & !grepl("sample", for_plot_qc[,1])),2], for_plot_qc[(grepl("rom_arg", for_plot_qc$breed) & !grepl("sample", for_plot_qc[,1])),3], col="green", pch=19)
points(for_plot_qc[grepl("Romanov", for_plot_qc$breed),2], for_plot_qc[grepl("Romanov", for_plot_qc$breed),3], col="grey", pch=19)


legend("topleft", legend=c("Argali", "katahdin_new", "katahdin", "rom_arg_new", "rom_arg", "Romanov"), 
  col=c("red", "orange", "orange", "green", "green", "grey"), cex=1.8, pch=c(19,17,19,17,19,19))

dev.off()

# pca plot for all animals with pheno

for_plot_pheno <- pca_with_pheno[,2:4]

breed <- c(1:nrow(for_plot_pheno))

for_plot_pheno <- cbind(for_plot_pheno, breed)


for_plot_pheno[(for_plot_pheno[,1] %in% breed_info[,1]),4] <- breed_info[match(for_plot_pheno[(for_plot_pheno[,1] %in% breed_info[,1]),1],breed_info[,1]),3]
for_plot_pheno[(for_plot_pheno[,1] %in% pheno_6d[,1]),4] <- pheno_6d[match(for_plot_pheno[(for_plot_pheno[,1] %in% pheno_6d[,1]),1],pheno_6d[,1]),4]


for_plot_pheno[1:46,4] <- "rom_arg"

for_plot_pheno[grepl("^Katahdin$", for_plot_pheno$breed),4] <- "katahdin" 


## make plot pca for animals with pheno

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots")

pdf("pca_plot_with_pheno.pdf", width=16, height=10)

plot(for_plot_pheno[,2], for_plot_pheno[,3], cex=1.8)

points(for_plot_pheno[grepl("katahdin", for_plot_pheno$breed),2], for_plot_pheno[grepl("katahdin", for_plot_pheno$breed),3], col="orange", pch=19)
points(for_plot_pheno[grepl("^rom_arg$", for_plot_pheno$breed),2], for_plot_pheno[grepl("^rom_arg$", for_plot_pheno$breed),3], col="green", pch=19)
points(for_plot_pheno[grepl("rom_arg_kat", for_plot_pheno$breed),2], for_plot_pheno[grepl("rom_arg_kat", for_plot_pheno$breed),3], col="blue", pch=19)
points(for_plot_pheno[grepl("rom_arg_kat_mouf", for_plot_pheno$breed),2], for_plot_pheno[grepl("rom_arg_kat_mouf", for_plot_pheno$breed),3], col="pink", pch=19)
points(for_plot_pheno[grepl("rom_mouf_kat", for_plot_pheno$breed),2], for_plot_pheno[grepl("rom_mouf_kat", for_plot_pheno$breed),3], col="brown", pch=19)

legend("topright", legend=c( "katahdin",  "rom_arg", "rom_arg_kat", "rom_arg_kat_mouf", "rom_mouf_kat"), 
	col=c("orange", "green", "blue", "pink", "brown"), cex=1.8, pch=19)

dev.off()


# make pools of animals based on for_plot_pheno pool of animals | sample_2, sample_35, sample_45, sample_47, sample_48 are outliers

outliers <- c("sample_2", "sample_35",  "sample_45", "sample_47", "sample_48")

## make pool_1, rom_arg, right cluster on plot | also delete sample_44, because he is in another pool on pca_plot_with_pheno

pool_1 <- for_plot_pheno[grepl("^rom_arg$", for_plot_pheno$breed),1]

pool_1 <- pool_1[!(pool_1 %in% outliers)]
pool_1 <- pool_1[-47]

## make pool_2, leftdown cluster on pca_plot_with_pheno and pool_3 leftup cluster 

pool_2 <- for_plot_pheno[!(for_plot_pheno[,1] %in% pool_1),1]
pool_2 <- pool_2[!(pool_2 %in% outliers)]

pool_3 <- c("sample_4", "sample_5", "sample_7", "sample_28", "sample_30", "sample_42", "sample_38", "sample_57")

pool_2 <- pool_2[!(pool_2 %in% pool_3)]


## make pool files for plink 

pool_1_for_plink <- matrix(ncol=2, nrow=length(pool_1)) 
pool_1_for_plink[,1] <- 0
pool_1_for_plink[,2] <- pool_1

pool_2_for_plink <- matrix(ncol=2, nrow=length(pool_2)) 
pool_2_for_plink[,1] <- 0
pool_2_for_plink[,2] <- pool_2

pool_3_for_plink <- matrix(ncol=2, nrow=length(pool_3)) 
pool_3_for_plink[,1] <- 0
pool_3_for_plink[,2] <- pool_3

## save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

write.table(pool_1, "pool_1.txt", col.names=F, row.names=F, quote=F)
write.table(pool_2, "pool_2.txt", col.names=F, row.names=F, quote=F)
write.table(pool_3, "pool_3.txt", col.names=F, row.names=F, quote=F)

write.table(pool_1_for_plink, "pool_1_for_plink.txt", col.names=F, row.names=F, quote=F)
write.table(pool_2_for_plink, "pool_2_for_plink.txt", col.names=F, row.names=F, quote=F)
write.table(pool_3_for_plink, "pool_3_for_plink.txt", col.names=F, row.names=F, quote=F)