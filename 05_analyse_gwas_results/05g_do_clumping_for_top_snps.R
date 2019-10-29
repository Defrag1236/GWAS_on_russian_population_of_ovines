### do clumping for top snps ###

library(data.table)

# load data and unlist 1st measurement of phenotypes

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/results/gemma_results")

load("6d_top5_snps.Rdata")

snps_6d <- rbindlist(trait_list_top_snps)

# load data and unlist 1st measurement of phenotypes

load("42d_top5_snps.Rdata")

snps_42d <- rbindlist(trait_list_top_snps)

# load data and unlist 1st measurement of phenotypes

load("3m_top5_snps.Rdata")

snps_3m <- rbindlist(trait_list_top_snps)

# extract and rename columns

snps_6d <- snps_6d[,c(1,3,2,12:13)]
colnames(snps_6d) <- c("CHR", "POS", "SNP", "P", "trait")

snps_42d <- snps_42d[,c(1,3,2,12:13)]
colnames(snps_42d) <- c("CHR", "POS", "SNP", "P", "trait")

snps_3m <- snps_3m[,c(1,3,2,12:13)]
colnames(snps_3m) <- c("CHR", "POS", "SNP", "P", "trait")

# do clumping

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

snps_all <- rbind(snps_6d, snps_42d, snps_3m)

snps_all <- as.data.frame(snps_all)

clumping_result <- function_for_shlop_28_12_2017(snps_all, p_value="P", pos="POS",snp="SNP",
                                       delta=2.5e5,chr="CHR",thr=0.05,trait="trait")

# save results 

write.table(clumping_result, "clumping_result.csv", col.names=T, row.names=F, quote=F, sep="\t")
