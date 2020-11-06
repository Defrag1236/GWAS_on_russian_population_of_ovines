#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_6d

plink --file ped_6d_1_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_1_index_filtered
plink --file ped_6d_2_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_2_index_filtered
plink --file ped_6d_3_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_3_index_filtered
plink --file ped_6d_4_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_4_index_filtered
plink --file ped_6d_5_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_5_index_filtered
plink --file ped_6d_6_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_6_index_filtered
plink --file ped_6d_7_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_6d_7_index_filtered

mv *.bed *.bim *.fam *.hh *.nosex /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_6d/bed/
rm *.log 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_42d

plink --file ped_42d_1_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_1_index_filtered
plink --file ped_42d_2_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_2_index_filtered
plink --file ped_42d_3_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_3_index_filtered
plink --file ped_42d_4_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_4_index_filtered
plink --file ped_42d_5_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_5_index_filtered
plink --file ped_42d_6_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_6_index_filtered
plink --file ped_42d_7_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_42d_7_index_filtered

mv *.bed *.bim *.fam *.hh *.nosex /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_42d/bed/
rm *.log 


cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_3m

plink --file ped_3m_1_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_1_index_filtered
plink --file ped_3m_2_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_2_index_filtered
plink --file ped_3m_3_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_3_index_filtered
plink --file ped_3m_4_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_4_index_filtered
plink --file ped_3m_5_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_5_index_filtered
plink --file ped_3m_6_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_6_index_filtered
plink --file ped_3m_7_index --make-bed --maf 0.1 --geno 0.05 --output-missing-genotype 0 --autosome-num 26 --out ped_3m_7_index_filtered

mv *.bed *.bim *.fam *.hh *.nosex /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/ped_3m/bed/
rm *.log 
