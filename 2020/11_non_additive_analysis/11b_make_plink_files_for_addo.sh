### make 2 plink files for addo: 1st: 1-23 chromosomes 2nd:24-26 chromosmes ###

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files

plink --bfile pool_1 --chr 1-23 --autosome-num 26 --make-bed --out /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/for_addo/pool_1_1_23