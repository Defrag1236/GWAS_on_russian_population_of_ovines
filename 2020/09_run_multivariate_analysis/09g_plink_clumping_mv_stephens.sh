

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/clumping/plink_clumping 

plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_1_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_1
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_2_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_2
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_3_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_3
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_4_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_4
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_5_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_5
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_6_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_6
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_index_7_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out index_7
plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/filtered/filtered_for_plink/p_value_mass_multi_stephens_maf_0.05_filtered.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out mass

R

index_1 <- read.table("index_1.clumped", head=T, stringsAsFactors=F) 
index_2 <- read.table("index_2.clumped", head=T, stringsAsFactors=F) #empty
index_3 <- read.table("index_3.clumped", head=T, stringsAsFactors=F) 
index_4 <- read.table("index_4.clumped", head=T, stringsAsFactors=F) 
index_5 <- read.table("index_5.clumped", head=T, stringsAsFactors=F) 
index_6 <- read.table("index_6.clumped", head=T, stringsAsFactors=F) 
index_7 <- read.table("index_7.clumped", head=T, stringsAsFactors=F) 
mass <- read.table("mass.clumped", head=T, stringsAsFactors=F) 

mv_all <- rbind(index_1[,c(3,5)], index_3[,c(3,5)], index_4[,c(3,5)], index_5[,c(3,5)], index_6[,c(3,5)], index_7[,c(3,5)], mass[,c(3,5)])

write.table(mv_all, "mv_all_clumped.txt", col.names=T, row.names=F, quote=F)

q()


plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --clump mv_all_clumped.txt  --autosome-num 26 --clump-p1 1.67e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 1.67e-8 --out mv_all
