library(survival)
library(pROC)
library(stringr)
library(data.table)
library(dplyr)
source("utility.R")

# loading scores and phenotype dataset
eur_case <- fread('~/scratch60/prs_comparison/case_list/eur_case_for_prs_comparison.tsv',header=T)
eur_cv1_id <- fread('~/scratch60/prs_comparison/id_list/eur_cv1_id.tsv', header=F)
eur_cv2_id <- fread('~/scratch60/prs_comparison/id_list/eur_cv2_id.tsv', header=F)

# helper functions
# preprocess
data_process <- function(data,disease){
    data$prs <- scale(data$prs)
    data$prs_top0.1 <- 0
    data$prs_top0.1[which(data$prs > quantile(data$prs,0.9,na.rm=T))] <- 1
    data$prs_top0.2 <- 0
    data$prs_top0.2[which(data$prs > quantile(data$prs,0.8,na.rm=T))] <- 1
    data$prs_top0.05 <- 0
    data$prs_top0.05[which(data$prs > quantile(data$prs,0.95,na.rm=T))] <- 1
    data$prs_top0.01 <- 0
    data$prs_top0.01[which(data$prs > quantile(data$prs,0.99,na.rm=T))] <- 1
    
    cencus <- paste0(disease,'_age')
    data$cencus <- apply(cbind(data[,..cencus],data$age_end),1,my.min)
    data$status <- data[,..disease] + 1
    return(data)
}

# min
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)

# roc
cal_roc <- function(dataset) {
    res.cox.prs <- coxph(Surv(cencus, status) ~ prs , data =  dataset)
    prediction.prs <- predict(res.cox.prs,newdata=dataset,type='risk')
    roc.prs <- roc(dataset$status,prediction.prs,col='1',lwd=1,main='ROC Curves',plot = F)
    roc.prs$auc
    round(ci.auc(roc.prs)[1:3],4)
    return(roc.prs)
}


# input files
# loading scores and phenotype dataset
eur_case <- fread('~/scratch60/prs_comparison/case_list/eur_case_for_prs_comparison.tsv',header=T)
eur_cv1_id <- fread("~/project/TiPred/data//cad_2000.txt")
score_path <- "~/scratch60/prs_comparison/shared_files/LDpred/LDpred_qc1/CADscore/p1.0000e+00score.txt"
score_file <- fread(score_path)
disease <- "CAD"

# preprocess (the first column of score file is ID)
score_file <- score_file[order(score_file[,1]),]
score_file_matched <- score_file[which(score_file$sample %in% eur_cv1_id$V1),]
score_file_ini <- score_file

# score
str_score <- "score"
eur_cv1 <- eur_case[which(eur_case$eid %in% eur_cv1_id$V1),]
eur_cv1$prs <- score_file_matched[, ..str_score]
eur_cv1 <- data_process(eur_cv1, disease)

# cal
tmp_roc <- cal_roc(eur_cv1)




