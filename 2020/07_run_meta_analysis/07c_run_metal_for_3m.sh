#!/usr/bin/bash

cd /home/common/projects/ovine_selection/Soft/metal/generic-metal

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_1_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_1_3m.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_2_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_2_3m.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_3_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_3_3m.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_4_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_4_3m.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_5_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_5_3m.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_6_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_6_3m.TBL
rm -r METAANALYSIS1.TBL.info

./metal /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/07_run_meta_analysis/config_files/3m/index_7_3m_config.txt

mv METAANALYSIS1.TBL /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/meta_analysis/3m/index_7_3m.TBL
rm -r METAANALYSIS1.TBL.info
