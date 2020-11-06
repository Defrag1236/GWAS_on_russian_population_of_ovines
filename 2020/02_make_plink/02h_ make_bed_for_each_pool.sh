#!/usr/bin/bash 


cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 192_animals --make-bed --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pool_1_for_plink.txt --autosome-num 26 --out pool_1
plink --bfile 192_animals --make-bed --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pool_2_for_plink.txt --autosome-num 26 --out pool_2
plink --bfile 192_animals --make-bed --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pool_3_for_plink.txt --autosome-num 26 --out pool_3