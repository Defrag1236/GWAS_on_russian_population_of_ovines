### grep mv snps and compare beta ###

# load data

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)


setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/Rdata")

load("metal_results_24_traits_unified.Rdata")

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/data/supplement_from_Bolormaa_et_al_2016")

s4 <- fread("12864_2016_2538_MOESM4_ESM", head=T, stringsAsFactors=F, data.table=F)
s5 <- fread("12864_2016_2538_MOESM5_ESM", head=T, stringsAsFactors=F, data.table=F)
s6 <- fread("12864_2016_2538_MOESM6_ESM", head=T, stringsAsFactors=F, data.table=F)
s7 <- fread("12864_2016_2538_MOESM7_ESM", head=T, stringsAsFactors=F, data.table=F)

# rbind data in 1 data.frame

s_all <- rbind(s4, s5, s6, s7)
rownames(s_all) <- s_all [,1]
s_all <- s_all[,-c(1)]

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/metal_results")
metal_p <- fread("mv_snps_p.txt", head=T, stringsAsFactors=F, data.table=F)

rownames(metal_p) <- metal_p[,1]
metal_p <- metal_p[8:11,-c(1)]

top_uv <- fread("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/top_uv_traits_for_general_table_new_threshold.txt", head=F, stringsAsFactors=F, data.table=F)

top_uv <- top_uv[8:11,1] 

# extract beta for 4 known mv snps

snps <- mv_snps[8:11,1]
snps_chr_pos <- paste(mv_snps[8:11,2], mv_snps[8:11,3], sep=":")

metal_p_min_names <- colnames(metal_p)[apply(metal_p, 1, which.min)]

## grep from traits in metal_p_names

beta_our <-list()

beta_our[[1]] <- metal_42d_clear[[3]][grepl(snps[1], metal_42d_clear[[3]]$MarkerName),c(1:3,6,8)]
beta_our[[2]] <- metal_3m_clear[[2]][grepl(snps[2], metal_3m_clear[[2]]$MarkerName),c(1:3,6,8)]
beta_our[[3]] <- metal_6d_clear[[3]][grepl(snps[3], metal_6d_clear[[3]]$MarkerName),c(1:3,6,8)]
beta_our[[4]] <- metal_6d_clear[[1]][grepl(snps[4], metal_6d_clear[[1]]$MarkerName),c(1:3,6,8)]

beta_our <- do.call(rbind, beta_our)

## grep from bolorrma original z-score

beta_bolormaa <-list()

beta_bolormaa[[1]] <- s_all[grepl(snps_chr_pos[1], rownames(s_all)),top_uv[1]]
beta_bolormaa[[2]] <- s_all[grepl(snps_chr_pos[2], rownames(s_all)),top_uv[2]]
beta_bolormaa[[3]] <- s_all[grepl(snps_chr_pos[3], rownames(s_all)),top_uv[3]]
beta_bolormaa[[4]] <- s_all[grepl(snps_chr_pos[4], rownames(s_all)),top_uv[4]]

beta_bolormaa <- do.call(rbind, beta_bolormaa)

## grep alleles for bolormaa

alleles_ref <- list()

for (n in (1:4)) {

	alleles_ref[[n]] <- snp_info[grepl(snps[n], snp_info$SNP_name),c(3,6)]

}

alleles_ref <- do.call(rbind, alleles_ref)

# make final table

beta_check <- cbind(beta_our, beta_bolormaa, alleles_ref[,1])

colnames(beta_check) <- c("SNP", "A1_metal", "A2_metal", "b_metal", "p-metal", "z_score_bolormaa", "alleles_AB_bolormaa")

# save results

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper")

fwrite(beta_check, "beta_check.txt", col.names=T, row.names=F, quote=F)