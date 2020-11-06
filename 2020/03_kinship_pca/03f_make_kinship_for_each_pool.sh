#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

# make tped/tfam files 

plink --bfile pool_1 --recode12 --output-missing-genotype 0 --transpose --autosome-num 26 --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1 
plink --bfile pool_2 --recode12 --output-missing-genotype 0 --transpose --autosome-num 26 --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_2 
plink --bfile pool_3 --recode12 --output-missing-genotype 0 --transpose --autosome-num 26 --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_3 


# make kinship IBS matrix 

/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1
/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_2
/home/common/projects/ovine_selection/Soft/emmax/emmax-kin-intel64 -v -d 10 /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_3



mv /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_1.aBN.kinf /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/kinship
mv /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_2.aBN.kinf /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/kinship
mv /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/pool_3.aBN.kinf /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/kinship