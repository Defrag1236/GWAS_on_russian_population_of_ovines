### run multivariate ###

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("MV_sst_file_for_21_metal_traits.Rdata")

# run multivariate and save results

library(MultiABEL)

multi_6d <- MultiSummary(sst, index=c(3, 6, 9, 12, 15, 18, 21), type = "outbred")
save(file="/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata/multi_6d.Rdata", list=c("multi_6d"))

multi_42d <- MultiSummary(sst, index=c(2, 5, 8, 11, 14, 17, 20), type = "outbred")
save(file="/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata/multi_42d.Rdata", list=c("multi_42d"))

multi_3m <- MultiSummary(sst, index=c(1, 4, 7, 10, 13, 16, 19), type = "outbred")
save(file="/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata/multi_3m.Rdata", list=c("multi_3m"))

multi_all <- MultiSummary(sst, index=c(1:21), type = "outbred")
save(file="/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata/multi_all.Rdata", list=c("multi_all"))

