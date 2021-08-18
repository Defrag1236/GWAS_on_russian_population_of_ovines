#!/usr/bin/bas

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/r2


plink -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 --r2 --autosome-num 26 --out r2_pool_1