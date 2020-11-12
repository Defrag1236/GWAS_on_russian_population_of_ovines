### make map for raw data ###

# load_data

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

load ("chip_data_pre-processed_96_animals_forward_alleles_all_snps.Rdata")

# do map file

map <- data.frame(
  CHR = chip_data_list[[1]]$Chr,
  SNP = chip_data_list[[1]]$SNP,
  CM = 0,
  BP = chip_data_list[[1]]$Pos)

# save map file

library (data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files")

fwrite (map, "96_animals_2020_raw.map", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")