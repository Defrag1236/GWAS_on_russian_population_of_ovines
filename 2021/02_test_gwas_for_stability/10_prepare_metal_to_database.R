### prepare metal results to save in db ###


# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/gwas_results")

metal <- read.table("metal_with_n_weights_n.TBL", head=T, stringsAsFactors=F)

# make coord file to match

coord <- snp_info[match(metal[,1], snp_info$SNP_name),c(2,4:6)]
coord <-  coord[!(grepl("X", coord$chromosome)),]
coord <-  coord[!(grepl("99", coord$chromosome)),]

# prepare metal results

metal_to_db <- data.frame()

metal_to_db <- metal[match(coord$SNP_name, metal[,1]),]
metal_to_db[,1] <- coord[match(coord$SNP_name, metal_to_db[,1]),"rs"]
metal_to_db <- metal_to_db[,c(1:4,6:8,10)]

# save results 

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/gwas_stability_test")

write.table(metal_to_db, "metal_to_db_with_n_weights_n.txt", row.names=F, col.names=T, quote=F, sep=",")