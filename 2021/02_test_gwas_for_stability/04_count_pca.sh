#!/usr/bin/bash 

# make pca for 107 animals pool_1+pool_2 without katahdin

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files

plink --bfile pool_1_and_pool_2_without_katahdin --autosome-num 26 --pca --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/pca/pool_1_and_pool_2_without_katahdin

# make pca for 122 animals with pheno data

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files

plink --bfile 122_animals_with_pheno --autosome-num 26 --pca --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/pca/122_animals_with_pheno
