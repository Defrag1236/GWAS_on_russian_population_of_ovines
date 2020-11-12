#!/usr/bin/bash 

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files

# exclude triallelic snps 

plink --bfile 96_animals_2020_raw  --exclude 192_animals_raw-merge.missnp --make-bed --autosome-num 26 --out 96_animals_2020_raw_excluded
plink --bfile 1st_batch_48_animals_2020_raw  --exclude 192_animals_raw-merge.missnp --make-bed --autosome-num 26 --out 1st_batch_48_animals_2020_raw_excluded
plink --bfile 2nd_batch_48_animals_2020_raw  --exclude 192_animals_raw-merge.missnp --make-bed --autosome-num 26 --out 2nd_batch_48_animals_2020_raw_excluded

# exclude again 

plink --bfile 96_animals_2020_raw_excluded  --exclude 192_animals_raw-merge.missnp --make-bed --autosome-num 26 --out 96_animals_2020_raw_excluded
plink --bfile 1st_batch_48_animals_2020_raw_excluded  --exclude 192_animals_raw-merge.missnp --make-bed --autosome-num 26 --out 1st_batch_48_animals_2020_raw_excluded
plink --bfile 2nd_batch_48_animals_2020_raw_excluded  --exclude 192_animals_raw-merge.missnp --make-bed --autosome-num 26 --out 2nd_batch_48_animals_2020_raw_excluded

plink --bfile 96_animals_2020_raw_excluded --merge-list /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/merge_list_bed_files.txt --make-bed --autosome-num 26 --out 192_animals_raw