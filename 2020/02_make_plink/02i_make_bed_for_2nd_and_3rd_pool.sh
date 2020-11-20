#!/usr/bin/bash

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 96_animals_2020_raw  --exclude 2nd_and_3rd_batch_raw-merge.missnp --make-bed --autosome-num 26 --out 96_animals_2020_raw_excluded_v1

plink --bfile 96_animals_2020_raw_excluded_v1 --bmerge 2nd_batch_48_animals_2020_raw --make-bed --autosome-num 26 --out  2nd_and_3rd_batch_raw
