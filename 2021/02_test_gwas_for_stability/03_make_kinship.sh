#!/usr/bin/bash 

# make kinship for 107 animals pool_1+pool_2 without katahdin

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files

plink --bfile pool_1_and_pool_2_without_katahdin --recode12 --output-missing-genotype 0 --transpose --autosome-num 26 --out pool_1_and_pool_2_without_katahdin

/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin

mv /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin.aBN.kinf /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship

# make kinship for 122 animals with pheno data

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files

plink --bfile 122_animals_with_pheno --recode12 --output-missing-genotype 0 --transpose --autosome-num 26 --out 122_animals_with_pheno

/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/122_animals_with_pheno

mv /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/122_animals_with_pheno.aBN.kinf /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship
