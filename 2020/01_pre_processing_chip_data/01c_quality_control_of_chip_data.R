### quality control of chip data ###


# load data

setwd ("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/Rdata")

load("chip_data_pre-processed_96_animals_AB_alleles.Rdata")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/reference_data")

reference <- read.table("SNPchimp_result_400473265.tsv", head=T)


# delete 0 from chr column and missing values from allele columns

for (n in(1:96)) {

	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^0$" ,x=chip_data_list[[n]]$Chr)),]
	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^-$" ,x=chip_data_list[[n]]$A1)),]
	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^-$" ,x=chip_data_list[[n]]$A2)),]
	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^I$" ,x=chip_data_list[[n]]$A1)),]
	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^I$" ,x=chip_data_list[[n]]$A2)),]
	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^D$" ,x=chip_data_list[[n]]$A1)),]
	chip_data_list[[n]] <- chip_data_list[[n]][!(grepl(pattern="^D$" ,x=chip_data_list[[n]]$A2)),]


	print(n)

}


# change AB alleles to FORWARD alleles

## prepare reference and chip tables

reference_forward_alleles <- strsplit(reference$Alleles_A_B_FORWARD, "/")

reference_forward_alleles <- do.call(rbind, reference_forward_alleles)

reference_forward_alleles <- cbind(reference[,c(2,11:13)], reference_forward_alleles)
colnames(reference_forward_alleles) <- c("rs", "Chr", "Pos", "SNP_name", "A1", "A2")

reference_forward_alleles <- reference_forward_alleles[!(grepl(pattern="^99$", x=reference_forward_alleles$Chr)),]
reference_forward_alleles <- reference_forward_alleles[!(grepl(pattern="IN", x=reference_forward_alleles$A1)),]
reference_forward_alleles <- reference_forward_alleles[!(grepl(pattern="DE", x=reference_forward_alleles$A2)),]

reference_forward_alleles <- reference_forward_alleles[(reference_forward_alleles$SNP_name %in% chip_data_list[[1]]$SNP),]

## unify dimension for all samples 

snps <- list()

for (n in(1:96)) {


	snps[[n]] <- chip_data_list[[n]]$SNP

	print(n)

}

snps[[97]] <- reference_forward_alleles$SNP_name

snps_clear <- Reduce(intersect, snps)

for (n in(1:96)) {


	chip_data_list[[n]] <- chip_data_list[[n]][match(snps_clear, chip_data_list[[n]]$SNP),]

	print(n)

}



reference_forward_alleles <- reference_forward_alleles[match(snps_clear, reference_forward_alleles$SNP_name),]
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

save(chip_data_list, file="chip_data_pre-processed_96_animals_forward_alleles.Rdata")