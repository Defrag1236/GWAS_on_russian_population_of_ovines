#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files/ped_6d

plink --file ped_6d_1_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_1_trait_filtered
plink --file ped_6d_2_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_2_trait_filtered
plink --file ped_6d_3_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_3_trait_filtered
plink --file ped_6d_4_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_4_trait_filtered
plink --file ped_6d_5_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_5_trait_filtered
plink --file ped_6d_6_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_6_trait_filtered
plink --file ped_6d_7_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_7_trait_filtered
plink --file ped_6d_8_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_8_trait_filtered
plink --file ped_6d_9_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_9_trait_filtered
plink --file ped_6d_10_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_10_trait_filtered
plink --file ped_6d_11_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_6d_11_trait_filtered