#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --file 2nd_batch_48_animals_2020_raw --make-bed --autosome-num 26 --out 2nd_batch_48_animals_2020_raw
