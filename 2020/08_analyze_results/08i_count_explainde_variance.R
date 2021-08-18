### discovery how many traits explained 90% of variance for 8 traits (6d) ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

# add all z_score 6d in 1 table

all_z_score <- cbind(metal_6d_clear[[1]]$z_score, metal_6d_clear[[2]]$z_score, metal_6d_clear[[3]]$z_score, metal_6d_clear[[4]]$z_score, metal_6d_clear[[5]]$z_score, metal_6d_clear[[6]]$z_score, metal_6d_clear[[7]]$z_score, metal_mass[[1]]$z_score )
rownames(all_z_score) <- metal_6d_clear[[1]]$markerName
colnames(all_z_score) <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")

# delete snps with abs(z_score) > 2 

all_z_score <- all_z_score[!rowSums(abs(all_z_score)>2),]

# count correlation and eigen()

all_cor <- cor(all_z_score)

eigen(all_cor)