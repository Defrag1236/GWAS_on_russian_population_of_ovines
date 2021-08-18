### replicate db snps on pc results ###

# load data

db_snps <- read.table("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/snps_from_database.txt", head=T, stringsAsFactors=F, sep="\t")

library(data.table)

pc <- list()

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/pc/")

for (n in 1:6) {

	name_to_read <- paste("pc_", n, "/pc_", n, "_done.csv", sep="")
	name_to_list <- paste("pc_", n, sep="")
	pc[[n]] <- fread(name_to_read, head=T, stringsAsFactors=F, data.table=F)
	names(pc)[[n]] <- name_to_list

}

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

# makee corrd matrix to match 

coord <- snp_info[match(pc[[1]]$rs_id, snp_info$rs),c(2, 4:6)]
coord <-  coord[!(grepl("X", coord$chromosome)),]
coord <-  coord[!(grepl("99", coord$chromosome)),]
rownames(coord) <- paste(coord$chromosome, coord$position, sep=":")

# grep db snps in chip data

db_snps_clear <- list()

for (n in 1:nrow(db_snps)) {

	db_snps_clear[[n]] <- coord[grepl(db_snps[n,1], coord$SNP_name),1]	

	print(n)

}

db_snps_clear <- do.call(rbind, db_snps_clear)

# grep clear db snps on 6 pc gwas

replicated_snps_pc_1 <- list()
replicated_snps_pc_2 <- list()
replicated_snps_pc_3 <- list()
replicated_snps_pc_4 <- list()
replicated_snps_pc_5 <- list()
replicated_snps_pc_6 <- list()

replicated <- list()

replicated <- list(replicated_snps_pc_1, replicated_snps_pc_2, replicated_snps_pc_3, replicated_snps_pc_4, replicated_snps_pc_5, replicated_snps_pc_6)

for (n in 1:nrow(db_snps_clear)) {

	for (i in 1:6) {

		replicated[[i]][[n]] <- pc[[i]][grepl(db_snps_clear[n,1], pc[[i]]$rs_id), c(3, 5:6, 13)]	

		

		}

	print (n)

}

pc_1_replicated <- do.call(rbind, replicated[[1]])
pc_2_replicated <- do.call(rbind, replicated[[2]])
pc_3_replicated <- do.call(rbind, replicated[[3]])
pc_4_replicated <- do.call(rbind, replicated[[4]])
pc_5_replicated <- do.call(rbind, replicated[[5]])
pc_6_replicated <- do.call(rbind, replicated[[6]])

# folter result s for p<0.05/(6*nrow(db_snps_clear))

pc_1_replicated_c <- subset(pc_1_replicated, pc_1_replicated$p<=0.05/(6*nrow(db_snps_clear)), select=1:4)
pc_2_replicated_c <- subset(pc_2_replicated, pc_2_replicated$p<=0.05/(6*nrow(db_snps_clear)), select=1:4)
pc_3_replicated_c <- subset(pc_3_replicated, pc_3_replicated$p<=0.05/(6*nrow(db_snps_clear)), select=1:4)
pc_4_replicated_c <- subset(pc_4_replicated, pc_4_replicated$p<=0.05/(6*nrow(db_snps_clear)), select=1:4)
pc_5_replicated_c <- subset(pc_5_replicated, pc_5_replicated$p<=0.05/(6*nrow(db_snps_clear)), select=1:4)
pc_6_replicated_c <- subset(pc_6_replicated, pc_6_replicated$p<=0.05/(6*nrow(db_snps_clear)), select=1:4)

# count lambda for grep snps 

lambda <- matrix(ncol=1, nrow=6)

lambda[1,1] <- median (qchisq(pc_1_replicated$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda[2,1] <- median (qchisq(pc_2_replicated$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda[3,1] <- median (qchisq(pc_3_replicated$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda[4,1] <- median (qchisq(pc_4_replicated$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda[5,1] <- median (qchisq(pc_5_replicated$p, df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda[6,1] <- median (qchisq(pc_6_replicated$p, df=1, lower.tail=F)/qchisq(0.5,df=1))