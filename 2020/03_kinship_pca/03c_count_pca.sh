#!/usr/bin/bash 

# pca for 192 animals

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

plink --bfile 192_animals_raw --maf 0.1 --geno 0.05 --autosome-num 26 --pca --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pca/pca_192_animals_filtered

# pca for new (2020) animals

plink --bfile 192_animals --maf 0.1 --geno 0.05 --autosome-num 26 --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/96_animals_id_2020.txt --pca --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pca/pca_96_animals_filtered

# pca for rom_arg, arg and kat 

plink --bfile 192_animals --maf 0.1 --geno 0.05 --autosome-num 26 --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/id_rom_rom_plus_arg_kat_arg.txt --pca --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pca/pca_rom_rom_plus_arg_kat_arg_filtered

# pca for animals with pheno

plink --bfile 192_animals_raw  --maf 0.1 --geno 0.05 --keep /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/id_animals_with_pheno.txt  --autosome-num 26 --pca --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/pca/pca_animals_with_pheno_filtered


