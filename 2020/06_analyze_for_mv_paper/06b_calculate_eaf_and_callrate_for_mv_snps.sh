#!/usr/bin/bash

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

# eaf and callrate on 1st pool

plink --bfile pool_1 --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_freq_callrate

# eaf and callrate on 1st batch

plink --bfile 1st_batch_48_animals_2020_raw --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_snps_freq_callrate_batch_1

# eaf and callrate on 2nd batch

plink --bfile 2nd_batch_48_animals_2020_raw --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_snps_freq_callrate_batch_2

# eaf and callrate on 3rd batch

plink --bfile 96_animals_2020_raw --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_snps_freq_callrate_batch_3

# eaf and callrate on 192_animals_new

plink --bfile 192_animals_raw --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_snps_freq_callrate_192_animals_raw

# eaf and callrate on 2nd_and_3rd_batch_raw

plink --bfile 2nd_and_3rd_batch_raw --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/mv_snps_to_replicate.txt --autosome-num 26 --freq --missing --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/for_mv_paper/mv_snps_freq_callrate_2nd_and_3rd_batch_raw
