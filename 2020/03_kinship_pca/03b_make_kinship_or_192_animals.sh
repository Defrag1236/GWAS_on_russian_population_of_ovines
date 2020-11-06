#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

# make tped/tfam files 

plink --bfile 192_animals --recode12 --output-missing-genotype 0 --transpose --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/for_all_samples_192_filtered --maf 0.1 --geno 0.05 --autosome-num 26

# make kinship IBS matrix 

/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/for_all_samples_192_filtered


mv /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/for_all_samples_192_filtered.aBN.kinf /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/kinship