#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/data/plink_files/ped_3m

plink --file ped_3m_1_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_1_trait_filtered
plink --file ped_3m_2_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_2_trait_filtered
plink --file ped_3m_3_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_3_trait_filtered
plink --file ped_3m_4_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_4_trait_filtered
plink --file ped_3m_5_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_5_trait_filtered
plink --file ped_3m_6_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_6_trait_filtered
plink --file ped_3m_7_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_7_trait_filtered
plink --file ped_3m_8_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_8_trait_filtered
plink --file ped_3m_9_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_9_trait_filtered
plink --file ped_3m_10_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_10_trait_filtered
plink --file ped_3m_11_trait --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --out ped_3m_11_trait_filtered