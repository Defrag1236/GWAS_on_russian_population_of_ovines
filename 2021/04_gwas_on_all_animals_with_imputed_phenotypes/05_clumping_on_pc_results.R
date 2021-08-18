# load data

library(data.table)

pc <- list()

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/pc/")

for (n in 1:24) {

	name_to_read <- paste("pc_", n, "/pc_", n, "_done.csv", sep="")
	name_to_list <- paste("pc_", n, sep="")
	pc[[n]] <- fread(name_to_read, head=T, stringsAsFactors=F, data.table=F)
	names(pc)[[n]] <- name_to_list

}

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

pc_clumped <- list()

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/clumping/pc")

for (n in 1:24) {

	for_clumping <- pc[[n]][,c("rs_id", "chr", "bp", "p")]
	for_clumping[,5] <- paste("pc_", n, sep="") 
	colnames(for_clumping)  <- c("SNP", "CHR", "POS", "P", "trait")

	if (min(for_clumping$P)<(0.05/nrow(for_clumping))) {

    clumped <- function_for_shlop_28_12_2017(locus_table=for_clumping,p_value="P",pos="POS",snp="SNP",delta=5e5, chr="CHR",thr=(0.05/nrow(for_clumping)), trait="trait")

 	pc_clumped[[n]] <- clumped

    colnames(pc_clumped[[n]]) <- c("SNP", "CHR", "POS", "P", "trait", "traits", "Ntrait")

    fwrite(file=paste("pc_", n, ".txt", sep=""),x=clumped,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t") }

}