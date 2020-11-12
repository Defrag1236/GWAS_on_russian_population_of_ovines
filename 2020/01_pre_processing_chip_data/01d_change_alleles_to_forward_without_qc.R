### change alleles to forward without qc ###



setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

load("chip_data_pre-processed_96_animals_AB_alleles.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/reference_data")

reference <- read.table("SNPchimp_result_400473265.tsv", head=T)

# change AB alleles to FORWARD alleles

## prepare reference and chip tables

reference_forward_alleles <- strsplit(reference$Alleles_A_B_FORWARD, "/")

reference_forward_alleles <- do.call(rbind, reference_forward_alleles)

reference_forward_alleles <- cbind(reference[,c(2,11:13)], reference_forward_alleles)
colnames(reference_forward_alleles) <- c("rs", "Chr", "Pos", "SNP_name", "A1", "A2")

reference_forward_alleles <- reference_forward_alleles[order(reference_forward_alleles$SNP_name, chip_data_list[[1]]$SNP), ]

for (n in (1:96)) {

	chip_data_list[[n]]$A1[(which(chip_data_list[[n]]$A1=="A"))] <- reference_forward_alleles[which(chip_data_list[[n]]$A1=="A"), "A1"]
	chip_data_list[[n]]$A1[(which(chip_data_list[[n]]$A1=="B"))] <- reference_forward_alleles[which(chip_data_list[[n]]$A1=="B"), "A2"]
	chip_data_list[[n]]$A2[(which(chip_data_list[[n]]$A2=="A"))] <- reference_forward_alleles[which(chip_data_list[[n]]$A2=="A"), "A1"]
	chip_data_list[[n]]$A2[(which(chip_data_list[[n]]$A2=="B"))] <- reference_forward_alleles[which(chip_data_list[[n]]$A2=="B"), "A2"]
	
	print(n)

}

# save results 

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

save(chip_data_list, file="chip_data_pre-processed_96_animals_forward_alleles_all_snps.Rdata")