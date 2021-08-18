### clumping mv results ###

# load data

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/mv_stephens/on_database_gwas")

index_1 <- read.table("p_value_index_1_multi_stephens.txt", head=F, stringsAsFactors=F)
index_2 <- read.table("p_value_index_2_multi_stephens.txt", head=F, stringsAsFactors=F)
index_3 <- read.table("p_value_index_3_multi_stephens.txt", head=F, stringsAsFactors=F)
index_4 <- read.table("p_value_index_4_multi_stephens.txt", head=F, stringsAsFactors=F)
index_5 <- read.table("p_value_index_5_multi_stephens.txt", head=F, stringsAsFactors=F)
index_6 <- read.table("p_value_index_6_multi_stephens.txt", head=F, stringsAsFactors=F)
index_7 <- read.table("p_value_index_7_multi_stephens.txt", head=F, stringsAsFactors=F)
mass <- read.table("p_value_mass_multi_stephens.txt", head=F, stringsAsFactors=F)

lambda <- read.table("lambda_index.txt", head=T, stringsAsFactors=F)

library(data.table)

setwd("/home/common/projects/ovine_selection/Data/ilumina_600K_forward_forward")

snp_info <- fread("SNPchimp_result_3250871638.csv", head=T, stringsAsFactors=F, data.table=F)

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


# do clumping for mv results 

mv_stephens_list_p_value <- list(index_1[,2], index_2[,2], index_3[,2], index_4[,2], index_5[,2], index_6[,2], index_7[,2], mass[,2])

mv_stephens_names_for_save <- c("index_1_stephens.txt", "index_2_stephens.txt", "index_3_stephens.txt", "index_4_stephens.txt", "index_5_stephens.txt",
									"index_6_stephens.txt", "index_7_stephens.txt", "mass_stephens.txt")

names_mv <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6", "index_7", "mass")

shlop_list_stephens_mv <- list()

m37 <- snp_info[match(index_1[,1], snp_info$rs),c(2, 4:5)]
colnames(m37) <- c ("SNP", "Chr", "Pos")

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/mv_stephens/on_database_gwas/clumping")

for (n in (1:8)) {

	y <- mv_stephens_list_p_value[[n]]
	x <- names_mv[n]

    cycle_chr_pos <- cbind(m37, y, x)
    colnames(cycle_chr_pos) <- c("SNP", "CHR", "POS", "P", "trait")
    
    if (min(mv_stephens_list_p_value[[n]])<(0.05/(nrow(index_1)*8))) {

    clumped <- function_for_shlop_28_12_2017(locus_table=cycle_chr_pos,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=(0.05/(nrow(index_1)*8)), trait="trait")

  

    shlop_list_stephens_mv[[n]] <- clumped

    colnames(shlop_list_stephens_mv[[n]]) <- c("SNP", "CHR", "POS", "P", "trait", "traits", "Ntrait")

    fwrite(file=mv_stephens_names_for_save[n],x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") }

	}

# do clumping for all mv stephens

shlop_all_stephens_mv <- do.call (rbind, shlop_list_stephens_mv)

shlop_result_all_stephens_mv <- function_for_shlop_28_12_2017(shlop_all_stephens_mv[,1:5], p_value="P", pos="POS",snp="SNP",
                                       delta=5e5,chr="CHR",thr=(0.05/(nrow(index_1)*8)),trait="trait")

fwrite("clumping_for_all_mv_stephens.txt", x=shlop_result_all_stephens_mv,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# make genomic control for clumping results 

chi <- qchisq(shlop_result_all_stephens_mv$P, df=1, low=F)

corrected_chi <- c(1:14)

corrected_chi[c(10,12)] <- chi[c(10,12)]/lambda[1,] 
corrected_chi[c(1:8,13)] <- chi[c(1:8,13)]/lambda[4,] 
corrected_chi[c(11,14)] <- chi[c(11,14)]/lambda[5,]
corrected_chi[9] <- chi[9]/lambda[8,] 

p_corrected <- pchisq(corrected_chi, df=1, low=F) 

gc_results <- shlop_result_all_stephens_mv

gc_results$P <- p_corrected

gc_results <- subset(gc_results, gc_results$P<=(0.05/(nrow(index_1)*8)), select=1:7)

fwrite("clumping_for_all_mv_stephens_gc.txt", x=gc_results,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")