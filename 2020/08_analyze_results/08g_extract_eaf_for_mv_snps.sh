#!/usr/bin/bash

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

## replace 1 snp allele

#plink --bfile 1st_batch_48_animals_2020_raw --reference-allele /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/new_mv_snp_to_replace_allele.txt  --autosome-num 26 --make-bed --out 1st_batch_48_animals_2020_raw

plink --bfile 192_animals_raw --keep /home/common/projects/ovine_selection/ovines_gwas_map/Data/argali_id.txt --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_for_plink.txt --autosome-num 26 --freq --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_argali_eaf

plink --bfile 192_animals_raw --keep /home/common/projects/ovine_selection/ovines_gwas_map/Data/romanovka_id.txt --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_for_plink.txt --autosome-num 26 --freq  --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_romanovka_eaf

plink --bfile 192_animals_raw --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/katahdin_id.txt --extract /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_for_plink.txt --autosome-num 26 --freq  --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_snps_katahdin_eaf