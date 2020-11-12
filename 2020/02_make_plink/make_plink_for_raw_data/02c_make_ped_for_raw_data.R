### make ped for raw data ###

# load_data

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

load ("chip_data_pre-processed_96_animals_forward_alleles_all_snps.Rdata")

load ("pedAB_96_animals_raw_data.Rdata")

# load pheno data (it's need for IID and sex info)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/pheno_data")

pheno_6d <- read.csv("pheno_index_6d_2020.csv", head=T, stringsAsFactors=F, sep=";", fileEncoding="latin1")

# add unknown samples to pheno data

unknown_samples <- matrix(nrow=21, ncol=ncol(pheno_6d))

unknown_samples[,1] <- names(chip_data_list)[!(names(chip_data_list) %in% pheno_6d[,1])]

unknown_samples[,2:4] <- 0

unknown_samples[,5:11] <- NA 

colnames(unknown_samples) <- colnames(pheno_6d)

pheno_6d <- rbind(pheno_6d, unknown_samples)


# make first 6 columns ped file

ped6 <- data.frame(
  		FID = 0,
  		IID = pheno_6d[,1],
  		FatherID = 0,
  		MotherID = 0,
  		Sex = pheno_6d$Sex,
  		Phenotype = -9)

to_match <- match(rownames(pedAB), ped6$IID)	

# make ped file

ped_full <- cbind(ped6[to_match,], pedAB)

# save result

library (data.table)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/plink_files")

fwrite(ped_full, "96_animals_2020_raw.ped", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

	

