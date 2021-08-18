#!/usr/bin/bas

cd /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/gwas_results

# adjust

/home/common/projects/ovine_selection/Soft/gemma-0.98.1-linux-static -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin -k /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship/pool_1_and_pool_2_without_katahdin.aBN.kinf -c /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/cov.txt -p /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/index_4_6d_adjust.txt -maf 0 -lmm 1 -o index_4_6d_adjust

# adjust norm

/home/common/projects/ovine_selection/Soft/gemma-0.98.1-linux-static -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin -k /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship/pool_1_and_pool_2_without_katahdin.aBN.kinf -c /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/cov.txt -p /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/index_4_6d_adjust_norm.txt -maf 0 -lmm 1 -o index_4_6d_adjust_norm

# run for each pool 

/home/common/projects/ovine_selection/Soft/gemma-0.98.1-linux-static -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin -k /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship/pool_1_and_pool_2_without_katahdin.aBN.kinf -c /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/cov.txt -p /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/index_4_6d_adjust_norm_pool_1.txt -maf 0 -lmm 1 -o index_4_6d_adjust_norm_pool_1
/home/common/projects/ovine_selection/Soft/gemma-0.98.1-linux-static -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/pool_1_and_pool_2_without_katahdin -k /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship/pool_1_and_pool_2_without_katahdin.aBN.kinf -c /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/cov.txt -p /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/index_4_6d_adjust_norm_pool_2.txt -maf 0 -lmm 1 -o index_4_6d_adjust_norm_pool_2

# run for 122 animals adjust norm

/home/common/projects/ovine_selection/Soft/gemma-0.98.1-linux-static -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/122_animals_with_pheno -k /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship/122_animals_with_pheno.aBN.kinf -c /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/cov_122.txt -p /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/index_4_6d_adjust_norm_122.txt -maf 0 -lmm 4 -o index_4_6d_adjust_norm_122
/home/common/projects/ovine_selection/Soft/gemma-0.98.1-linux-static -bfile /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/122_animals_with_pheno -k /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/kinship/122_animals_with_pheno.aBN.kinf -p /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/index_4_6d_adjust_norm_122.txt -maf 0 -lmm 4 -o index_4_6d_adjust_norm_122_without_cov
