#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files

# make tped/tfam files 

plink --bfile for_all_samples --recode12 --output-missing-genotype 0 --transpose --out for_all_samples_filtered --maf 0.1 --geno 0.05

# make kinship IBS matrix 

/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files/for_all_samples_filtered