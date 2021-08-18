#!/usr/bin/bash

cd /home/common/projects/ovine_selection/Soft/metal/generic-metal

# indexes

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_1_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_1_6d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_2_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_2_6d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_3_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_3_6d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_4_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_4_6d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_5_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_5_6d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_6_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_6_6d.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/index_7_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/index_7_6d.TBL
rm -r METAANALYSIS1.TBL.info

# mass

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/6d/mass_6d_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/with_n/6d/mass_6d.TBL
rm -r METAANALYSIS1.TBL.info
