#!/usr/bin/bash

cd /home/common/projects/ovine_selection/Soft/metal/generic-metal

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_1_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_1_42d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_2_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_2_42d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_3_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_3_42d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_4_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_4_42d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_5_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_5_42d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_6_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_6_42d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/42d/index_7_42d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/42d/index_7_42d.TBL
rm -r METAANALYSIS1.TBL.info
