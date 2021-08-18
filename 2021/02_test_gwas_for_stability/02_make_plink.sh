#!/usr/bin/bash 

# make plink for 107 animals pool_1+pool_2 without katahdin

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 2nd_and_3rd_batch_raw --make-bed --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pool_1_pool_2_id_without_katahdin_for_plink.txt --autosome-num 26 --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin

# make plink for 122 animals with pheno data

plink --bfile 2nd_and_3rd_batch_raw --make-bed --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/122_animals_with_pheno.txt --autosome-num 26 --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/122_animals_with_pheno
	