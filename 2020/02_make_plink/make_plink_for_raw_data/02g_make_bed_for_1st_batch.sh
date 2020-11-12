#!/usr/bin/bash  

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --file /home/common/projects/ovine_selection/ovines_gwas_map/Data/plink_files/plink_for_1st_48_sheeps --make-bed --autosome-num 26 --out 1st_batch_48_animals_2020_raw
