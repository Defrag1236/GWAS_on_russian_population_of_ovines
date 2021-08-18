### clumping on imputed pheno ###

# load data

library(data.table)

gwas_6d <- list()
gwas_42d <- list()
gwas_3m <- list()

for (n in 1:7) {

	name_to_read_6d <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/imp_adj_norm/index_", n, "_6d_imp_adj_norm/index_", n, "_6d_imp_adj_norm_done.csv", sep="")
	name_to_list_6d <- paste("index_", n, "_6d", sep="")
	gwas_6d[[n]] <- fread(name_to_read_6d, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_6d)[[n]] <- name_to_list_6d

	name_to_read_42d <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/imp_adj_norm/index_", n, "_42d_imp_adj_norm/index_", n, "_42d_imp_adj_norm_done.csv", sep="")
	name_to_list_42d <- paste("index_", n, "_42d", sep="")
	gwas_42d[[n]] <- fread(name_to_read_42d, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_42d)[[n]] <- name_to_list_42d

	name_to_read_3m <- paste("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/imp_adj_norm/index_", n, "_3m_imp_adj_norm/index_", n, "_3m_imp_adj_norm_done.csv", sep="")
	name_to_list_3m <- paste("index_", n, "_3m", sep="")
	gwas_3m[[n]] <- fread(name_to_read_3m, head=T, stringsAsFactors=F, data.table=F)
	names(gwas_3m)[[n]] <- name_to_list_3m

}

gwas_mass <- list()

gwas_mass[[1]] <- fread("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/imp_adj_norm/mass_6d_imp_adj_norm/mass_6d_imp_adj_norm_done.csv", head=T, stringsAsFactors=F, data.table=F)
gwas_mass[[2]] <- fread("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/imp_adj_norm/mass_42d_imp_adj_norm/mass_42d_imp_adj_norm_done.csv", head=T, stringsAsFactors=F, data.table=F)
gwas_mass[[3]] <- fread("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/imp_adj_norm/mass_3m_imp_adj_norm/mass_3m_imp_adj_norm_done.csv", head=T, stringsAsFactors=F, data.table=F)

# clumping function

function_for_shlop_28_12_2017=function(locus_table,p_value="P",pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=5e-8,trait=NULL){
 #locus_table=bt
    locus_table[,p_value]=as.numeric(locus_table[,p_value])
    
    if (!is.null(trait)){
        traits="traits"
        locus_table=cbind(locus_table,traits=locus_table[,trait])
        locus_table[,traits]=as.character(locus_table[,traits])
    }
    
    out=locus_table[0,]
    
    locus_table=locus_table[locus_table[,p_value]<=thr,]
    
    i=1
    if (nrow(locus_table)>0){
        locus_table[,pos]=as.numeric(locus_table[,pos])
        locus_table[,p_value]=as.numeric(locus_table[,p_value])
        Zx <-locus_table
        Zx=Zx[order(Zx[,p_value]),]
        #n_traits=1
        #Zx=cbind(Zx,n_traits)
        i=1
        while (nrow(Zx)>0){       
            ind=which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))
                      
            if (!is.null(trait)){
                Zx[i,traits]=paste(unique(Zx[ind,trait]),collapse = ";")
            }
            
            out=rbind(out,Zx[i,])
            Zx=Zx[-ind,]            
        }
        rownames(out)=as.character(out[,snp])
    }
    
    if (!is.null(trait)){
        j=1
        out=cbind(out,Ntraits=1)
        out[,"Ntraits"]=as.numeric(out[,"Ntraits"])
        for (j in 1:nrow(out)){
            trs=unique(unlist(strsplit(out[j,traits],split = ";")))
            out[j,traits]=paste(trs,collapse = ";")
            out[j,"Ntraits"]=length(trs)
        }
    }
    
 return(out)
}

# make clumping

clumped_6d <- list()
clumped_42d <- list()
clumped_3m <- list()
clumped_mass <- list()

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/clumping/imputed_pheno")

for (n in 1:7) {

	for_clumping_6d <- gwas_6d[[n]][,c("rs_id", "chr", "bp", "p")]
	for_clumping_6d[,5] <- paste("index_", n, "_6d", sep="") 
	colnames(for_clumping_6d)  <- c("SNP", "CHR", "POS", "P", "trait")

	if (min(for_clumping_6d$P)<(0.05/nrow(for_clumping_6d))) {

    clumped <- function_for_shlop_28_12_2017(locus_table=for_clumping_6d,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=(0.05/nrow(for_clumping_6d)), trait="trait")

 	clumped_6d[[n]] <- clumped

    colnames(clumped_6d[[n]]) <- c("SNP", "CHR", "POS", "P", "trait", "traits", "Ntrait")

    fwrite(file=paste("index_", n, "_6d.txt", sep=""),x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") }

}

for (n in 1:7) {

	for_clumping_42d <- gwas_42d[[n]][,c("rs_id", "chr", "bp", "p")]
	for_clumping_42d[,5] <- paste("index_", n, "_42d", sep="") 
	colnames(for_clumping_42d)  <- c("SNP", "CHR", "POS", "P", "trait")

	if (min(for_clumping_42d$P)<(0.05/nrow(for_clumping_42d))) {

    clumped <-
    function_for_shlop_28_12_2017(locus_table=for_clumping_42d,p_value="P",pos="POS",snp="SNP",delta=5e5,
    chr="CHR",thr=(0.05/nrow(for_clumping_42d)), trait="trait")

 	clumped_42d[[n]] <- clumped

    colnames(clumped_42d[[n]]) <- c("SNP", "CHR", "POS", "P", "trait",
    "traits", "Ntrait")

    fwrite(file=paste("index_", n, "_42d.txt",
    sep=""),x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") }

}

for (n in 1:7) {

	for_clumping_3m <- gwas_3m[[n]][,c("rs_id", "chr", "bp", "p")]
	for_clumping_3m[,5] <- paste("index_", n, "_3m", sep="") 
	colnames(for_clumping_3m)  <- c("SNP", "CHR", "POS", "P", "trait")

	if (min(for_clumping_3m$P)<(0.05/nrow(for_clumping_3m))) {

    clumped <- function_for_shlop_28_12_2017(locus_table=for_clumping_3m,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=(0.05/nrow(for_clumping_3m)), trait="trait")

 	clumped_3m[[n]] <- clumped

    colnames(clumped_3m[[n]]) <- c("SNP", "CHR", "POS", "P", "trait", "traits", "Ntrait")

    fwrite(file=paste("index_", n, "_3m.txt", sep=""),x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") }

}

for (n in 1:3) {

	for_clumping_mass <- gwas_mass[[n]][,c("rs_id", "chr", "bp", "p")]
	for_clumping_mass[,5] <- paste(n, "_mass", sep="") 
	colnames(for_clumping_mass)  <- c("SNP", "CHR", "POS", "P", "trait")

	if (min(for_clumping_mass$P)<(0.05/nrow(for_clumping_mass))) {

    clumped <- function_for_shlop_28_12_2017(locus_table=for_clumping_mass,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=(0.05/nrow(for_clumping_mass)), trait="trait")

 	clumped_mass[[n]] <- clumped

    colnames(clumped_mass[[n]]) <- c("SNP", "CHR", "POS", "P", "trait", "traits", "Ntrait")

    fwrite(file=paste(n, "_mass.txt", sep=""),x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") }

}

# clumping for all traits

clumped_6d <- do.call(rbind, clumped_6d)
clumped_42d <- do.call(rbind, clumped_42d)
clumped_3m <- do.call(rbind, clumped_3m)
clumped_mass <- do.call(rbind, clumped_mass)

## rename 2_mass to mass_42d

clumped_mass$trait <- "mass_42d"

## do clumping

clumped_all <- rbind(clumped_6d[,1:5], clumped_42d[,1:5], clumped_3m[,1:5], clumped_mass[,1:5])

clumped_all_final <- function_for_shlop_28_12_2017(locus_table=clumped_all,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=(0.05/nrow(for_clumping_mass)), trait="trait")

fwrite("clumped_all.txt",x=clumped_all_final,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")