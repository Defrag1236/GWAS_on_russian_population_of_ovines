#!/usr/bin/bas

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile pool_1 --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_freq_callrate

