#!/usr/bin/bas

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/r2


plink -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/122_animals_with_pheno --r2 --autosome-num 26 --ld-window-r2 0 --out r2_122_animals_with_pheno