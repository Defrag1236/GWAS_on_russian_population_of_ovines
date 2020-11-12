# GWAS_on_russian_population_of_ovines

Изначально, есть 3 бэтча данных: 
1) 48 животных, в том числе архары, катадины, романовки и гибриды первого поления. Адекватных прмоеров нет, в гвасы не идут, но нужны для PCA и прочих вещей. 
2) 48 гибридов бэккроссов романовки и архара с промерами. По факту в гвасы берем 46, так как для двух образцов непонятны фенотипы. 
3) 96 животных, бэккросов романовки с архаром, гибридов романовки с катадином и фулоно разной кровности, гибридов с архаром разной кровности. Для 75 есть фенотипы.

Изначально, при переводе в плинк был сделан QC и в некоторых случаях фильтрация. Скрипты это отражают, однако, это не верный подход. Также имеются скрипты и/или строчки в уже существующих скриптах для перевода в плинк без фильтраций.

Полные нефильтроаванные файлы плинков для каждого бэтча:  
1) /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/1st_batch_48_animals_2020_raw.*
2) /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/2nd_batch_48_animals_2020_raw.*
3) /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/96_animals_2020_raw.*

Все 3 батча были объеднены в один плинковский файл -> /home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files/192_animals_raw.*

Начиная со скрипта 02h в папке 2020 для расчетов используется указанный выше общий плинковский файл.
