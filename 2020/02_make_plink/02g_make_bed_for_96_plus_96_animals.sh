#!/usr/bin/bash 

# make plink with filtered data

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_6d/bed

plink --bfile ped_6d_1_index_filtered --exclude /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/192_animals-merge.missnp --make-bed --autosome-num 26 --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/96_animals_2020_excluded

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 96_animals_2020_excluded --bmerge /home/common/projects/ovine_selection/ovines_gwas_map/Data/plink_files/plink_96_sheeps_reference_excluded --make-bed --autosome-num 26 --out  192_animals

# make plink with raw_data

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 96_animals_2020_raw  --exclude 192_animals_new-merge.missnp --make-bed --autosome-num 26 --out 96_animals_2020_raw_excluded


plink --bfile 96_animals_2020_raw_excluded --bmerge /home/common/projects/ovine_selection/ovines_gwas_map/Data/plink_files/plink_96_sheeps_reference_excluded --make-bed --autosome-num 26 --out  192_animals_new


