### count new p-value for mv snps ###

# load data 

library(data.table)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

mv_snps <- fread("general_table_for_paper_new_threshold.txt", head=T, stringsAsFactors=F, data.table=F)

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results/mv_rdata/mv_stephens_formula")

load("multi_weight.Rdata")
load("multi_fat.Rdata")
load("multi_muscle.Rdata")

# count lambda

lambda_w <- median (qchisq(p_value_weight_multivariate[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_f <- median (qchisq(p_value_fat_multivariate[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))
lambda_m <- median (qchisq(p_value_muscle_multivariate[,1], df=1, lower.tail=F)/qchisq(0.5,df=1))

# count new mv p 

new_mv_p <- mv_snps[,1:4]

new_mv_p[,4] <- qchisq(new_mv_p[,4], df=1, low=F)
new_mv_p[c(1:3, 5, 8:12),4] <- pchisq(new_mv_p[c(1:3, 5, 8:12),4]/lambda_w, df=1, low=F)
new_mv_p[c(4, 6:7),4] <- pchisq(new_mv_p[c(4, 6:7),4]/lambda_f, df=1, low=F)

# save result 

setwd("/home/common/projects/ovine_selection/ovines_multivariate_analysis/results")

fwrite(new_mv_p, "mv_p_new.txt", col.names=T, row.names=F, quote=F, sep="\t")