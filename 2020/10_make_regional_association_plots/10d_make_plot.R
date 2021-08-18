### make plot ###

# path to function
source("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/zlobin_src/GWAS_on_russian_population_of_ovines/2020/10_make_regional_association_plots/10c_plot_function.R")

# load data 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/mv_stephens/clumping")
new_snps <- read.table("clumping_for_all_mv_stephens.txt", head=T, stringsAsFactors=F)

# make plots for new snps associated with index_4

  # interesting (best) SNP
  snp       <- new_snps[7,1]
  # title of plot
  locusname <- paste("SNP", snp,"for index_4", sep=" ")
  # locus file (e.g. output of "building_locus_file.r")
  locus     <- read.table(paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/data/for_regional_association_plots/locus_table_", snp, ".txt", sep=""), header=T, stringsAsFactors=F)
  # known genes path
	known_genes_path <- "na"	
  # path for plot
  pfad_ausgabe <- paste("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2020/results/plots/regional_association_plots/assocplot_LDregions_", snp,".pdf", sep="")
  # range of y-axis (P-values) 
  minimaler_pwert <- min(locus[,6])
  range_graph     <- round(-(log10(minimaler_pwert))+1)
  # pvalue of interesting (best) SNP
  pval_snp    <- locus[locus[,5]==snp,6]
  # Build dependend on used HapMap version for SNAP informations
  build_ncbi       <- "b36"
  # plattform used 
  plattform   <- "linux"

  # open plot
  pdf(pfad_ausgabe, width=6, height=6)
  # call function
  make.fancy.locus.plot(snp, locusname, locus, range_graph, pval_snp, build_ncbi, plattform, known_genes_path)
  # close plot
  dev.off()




