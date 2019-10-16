#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files/ped_42d

plink --file ped_42d_1_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_1_trait_filtered
plink --file ped_42d_2_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_2_trait_filtered
plink --file ped_42d_3_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_3_trait_filtered
plink --file ped_42d_4_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_4_trait_filtered
plink --file ped_42d_5_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_5_trait_filtered
plink --file ped_42d_6_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_6_trait_filtered
plink --file ped_42d_7_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_7_trait_filtered
plink --file ped_42d_8_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_8_trait_filtered
plink --file ped_42d_9_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_9_trait_filtered
plink --file ped_42d_10_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_10_trait_filtered
plink --file ped_42d_11_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_42d_11_trait_filtered