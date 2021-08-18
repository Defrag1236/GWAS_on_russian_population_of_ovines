#!/usr/bin/bas

# make bed

plink2 --pfile /mnt/polyomica/projects/RuDDS/tmp/sheep_13.05.2021/quality_control/merged_top_v2_qc_passed --make-bed --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/122_animals_with_pheno_new_qc.txt  --chr-override --sheep --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/new_qc_data/animals_with_pheno

# make tped for emmax 

plink --bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/new_qc_data/animals_with_pheno --autosome-num 26 --recode12 --output-missing-genotype 0 --transpose --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/new_qc_data/animals_with_pheno

# make kinship 