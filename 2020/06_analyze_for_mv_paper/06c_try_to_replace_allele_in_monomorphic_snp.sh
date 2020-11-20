#!/usr/bin/bash

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 2nd_batch_48_animals_2020_raw --reference-allele /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snp_to_replace_allele.txt  --autosome-num 26 --make-bed --out 2nd_batch_48_animals_2020_raw
