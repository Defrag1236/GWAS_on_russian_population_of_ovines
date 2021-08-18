# run metal

cd /home/common/projects/ovine_selection/Soft/metal/generic-metal


./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/data_for_gwas/config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/gwas_results/metal_with_n_weights_n.TBL
rm -r METAANALYSIS1.TBL.info
