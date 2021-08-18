### run addo heterotic model ###

# load data

library('ADDO')

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results")

id <- read.table("pool_1.txt", head=F, stringsAsFactors=F)

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pheno_data")

phe_cov <- list()

phe_cov[[1]] <- id
phe_cov[[2]] <- read.table("pool_1_cov.txt", head=F, stringsAsFactors=F)[,1]
phe_cov[[3]] <- read.table("pool_1_cov.txt", head=F, stringsAsFactors=F)[,2]
phe_cov[[4]] <- read.table("pool_1_cov.txt", head=F, stringsAsFactors=F)[,3]
phe_cov[[5]] <- read.table("pool_1_cov.txt", head=F, stringsAsFactors=F)[,4]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/pheno_data/pheno_qc")

phe_cov[[6]] <- read.table("index_1_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[7]] <- read.table("index_2_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[8]] <- read.table("index_3_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[9]] <- read.table("index_4_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[10]] <- read.table("index_5_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[11]] <- read.table("index_6_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[12]] <- read.table("index_7_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[13]] <- read.table("mass_6d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[14]] <- read.table("index_1_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[15]] <- read.table("index_2_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[16]] <- read.table("index_3_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[17]] <- read.table("index_4_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[18]] <- read.table("index_5_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[19]] <- read.table("index_6_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[20]] <- read.table("index_7_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[21]] <- read.table("mass_42d_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[22]] <- read.table("index_1_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[23]] <- read.table("index_2_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[24]] <- read.table("index_3_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[25]] <- read.table("index_4_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[26]] <- read.table("index_5_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[27]] <- read.table("index_6_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[28]] <- read.table("index_7_3m_pool_1.txt", head=F, stringsAsFactors=F)
phe_cov[[29]] <- read.table("mass_3m_pool_1.txt", head=F, stringsAsFactors=F)


# correct sex fom 0 (female) and 1 (male) 

phe_cov[[4]][which(phe_cov[[4]]==2)] <- 0

# save .phe and .cov files

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/for_addo")


phe_cov_to_save <- do.call(cbind, phe_cov)
 
colnames(phe_cov_to_save) <- c("id", "intercept", "fetus", "sex", "batch", "index_1_6d", "index_2_6d", "index_3_6d", "index_4_6d", "index_5_6d", "index_6_6d", "index_7_6d", "mass_6d",
					"index_1_42d", "index_2_42d", "index_3_42d", "index_4_42d", "index_5_42d", "index_6_42d", "index_7_42d", "mass_42d",
					"index_1_3m", "index_2_3m", "index_3_3m", "index_4_3m", "index_5_3m", "index_6_3m", "index_7_3m", "mass_3m")


write.table(phe_cov_to_save[,1:6], "pool_1_1_23.phe", col.names=T, row.names=F, quote=F)

# Run 1st step (QC)

covariates_types <- c("n", "n", "n", "f")
names(covariates_types) = c("intercept", "fetus", "sex", "batch")

ADDO_Heterotic1_QC(indir = "/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/data/plink_files/for_addo" , 
	outdir = "/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/ADDO_results",
	Input_name = "pool_1_1_23", Input_type = "PLINK",
	Kinship_type = "EMMAX", covariates_types = covariates_types, GT_maf = 0.01, GT_missing = 0.05, num_nodes = 10)
