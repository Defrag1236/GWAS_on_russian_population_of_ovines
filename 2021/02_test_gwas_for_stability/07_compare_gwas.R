### compare 2 gwases ###

# load data

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/unified_data/index_4_6d")

orig <- read.csv("index_4_6d_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/gwas_stability_test/index_4_6d_adjust")

adj <- read.csv("index_4_6d_adjust_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/gwas_stability_test/index_4_6d_adjust_norm")

adj_norm <- read.csv("index_4_6d_adjust_norm_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/gwas_stability_test/index_4_6d_metal")

metal <- read.csv("index_4_6d_metal_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/gwas_stability_test/index_4_6d_adjust_norm_122")

adj_norm_122 <- read.csv("index_4_6d_adjust_norm_122_done.csv", head=T, stringsAsFactors=F)

setwd("/mnt/polyomica/ovis/gwas_storage/MARH/gwas_stability_test/index_4_6d_adjust_norm_122_without_cov")

adj_norm_122_wc <- read.csv("index_4_6d_adjust_norm_122_without_cov_done.csv", head=T, stringsAsFactors=F)

# make z-z plot

common <- intersect(adj$rs_id, orig$rs_id)

orig <- orig[orig$rs_id %in% common,]
adj <- adj[adj$rs_id %in% common,]
adj_norm <- adj_norm[adj_norm$rs_id %in% common,]
metal <- metal[metal$rs_id %in% common,]
adj_norm_122 <- adj_norm_122[adj_norm_122$rs_id %in% common,]
adj_norm_122_wc <- adj_norm_122_wc[adj_norm_122_wc$rs_id %in% common,]


orig <- orig[match(adj$rs_id, orig$rs_i),]
metal <- metal[match(adj$rs_id, metal$rs_i),]

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test/plots")

pdf("z_z_plot_orig_vs_adj.pdf")

plot(orig$z, adj$z, main="orig_vs_adj")

dev.off()

pdf("z_z_plot_orig_vs_adj_norm.pdf")

plot(orig$z, adj_norm$z, main="orig_vs_adj_norm")

dev.off()

pdf("z_z_plot_adj_vs_adj_norm.pdf")

plot(adj$z, adj_norm$z, main="adj_vs_adj_normm")

dev.off()

pdf("z_z_plot_adj_norm_vs_metal.pdf")

plot(adj_norm$z, metal$z, main="adj_norm_vs_metal")

dev.off()

pdf("z_z_plot_orig_vs_metal.pdf")

plot(adj_norm$z, metal$z, main="orig_vs_metal")

dev.off()

jpeg("z_z_plot_122_vs_122_without_cov.jpg", height=1500, width=1500)

plot(adj_norm_122$z, adj_norm_122_wc$z)

dev.off()

cor(orig$z, adj$z)
cor(orig$z, adj_norm$z)
cor(adj$z, adj_norm$z)
cor(adj_norm$z, metal$z)
cor(orig$z, metal$z)

## make af-af plot for adj_norm and metal

pdf("af_af_plot_adj_norm_vs_metal.pdf")

plot(adj_norm$eaf, metal$eaf, main="adj_norm_vs_metal")

dev.off()


# make clumping 

## clumping function

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


orig_for_clumping <- as.data.frame(orig[,c(3,5:6,13)])
orig_for_clumping[,5] <- "index_4_6d"
colnames(orig_for_clumping) <- c("SNP", "CHR", "POS", "P", "trait")

adj_for_clumping <- as.data.frame(adj[,c(3,5:6,13)])
adj_for_clumping[,5] <- "index_4_6d_adj"
colnames(adj_for_clumping) <- c("SNP", "CHR", "POS", "P", "trait")

adj_norm_for_clumping <- as.data.frame(adj_norm[,c(3,5:6,13)])
adj_norm_for_clumping[,5] <- "index_4_6d_adj_norm"
colnames(adj_norm_for_clumping) <- c("SNP", "CHR", "POS", "P", "trait")

metal_for_clumping <- as.data.frame(metal[,c(3,5:6,13)])
metal_for_clumping[,5] <- "index_4_6d_metal"
colnames(metal_for_clumping) <- c("SNP", "CHR", "POS", "P", "trait")

orig_clumping <- function_for_shlop_28_12_2017(locus_table=orig_for_clumping, p_value="P", pos="POS", snp="SNP", delta=5e5, chr="CHR", thr=(0.05/nrow(orig_for_clumping)), trait="trait")
adj_clumping <- function_for_shlop_28_12_2017(locus_table=adj_for_clumping, p_value="P", pos="POS", snp="SNP", delta=5e5, chr="CHR", thr=(0.05/nrow(adj_for_clumping)), trait="trait")
adj_norm_clumping <- function_for_shlop_28_12_2017(locus_table=adj_norm_for_clumping, p_value="P", pos="POS", snp="SNP", delta=5e5, chr="CHR", thr=(0.05/nrow(adj_norm_for_clumping)), trait="trait")
metal_clumping <- function_for_shlop_28_12_2017(locus_table=metal_for_clumping, p_value="P", pos="POS", snp="SNP", delta=5e5, chr="CHR", thr=(0.05/nrow(metal_for_clumping)), trait="trait")


common_for_clumping <- rbind(orig_clumping, adj_clumping)
common_for_clumping <- common_for_clumping[,-c(6:7)]

common_clumping <- function_for_shlop_28_12_2017(locus_table=common_for_clumping, p_value="P", pos="POS", snp="SNP", delta=5e5, chr="CHR", thr=(0.05/nrow(adj_for_clumping)), trait="trait")		

# save results 

setwd("/home/common/projects/ovine_selection/GWAS_on_russian_population_of_ovines/2021/results/gwas_staibility_test")

write.table(common_clumping, "common_clumping.txt", row.names=F, col.names=T, quote=F)