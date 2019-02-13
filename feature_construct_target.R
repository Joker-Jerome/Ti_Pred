# library
library(dplyr)
library(ggplot2)
library(ROCR)
library(data.table)
library(doParallel)
library(svMisc)
library(stringr)
library(Matrix)
library(RSQLite)
library(DBI)
library(IRdisplay)

# helper function snp
extract_snp <- function(str) {
    snp <- strsplit(str, "_")[[1]][1]
    return(snp)
}

# select gene
select_gene <- function(r2_df, cut_off) {
    gene_id <- r2_df %>%
        filter(rsq >= cut_off) %>%
        select(gene)
    gene_id <- as.vector(sapply(gene_id$gene, function(x) strsplit(x, "\\.")[[1]][1]))
    return(gene_id)
}

# filter match gene
filter_gene <- function(gene_vec, ref_table) {
    output_vec <- c()
    for (i in 1:length(gene_vec)) {
        if (gene_vec[i] %in% ref_table$'Gene stable ID') {
            output_vec <- c(output_vec, gene_vec[i])
        }
    }
    return(output_vec)
}   

# add ensembl
add_ensembl_id <- function(r2_df) {
    r2_df$ensembl_id <- as.vector(sapply(r2_df$gene, function(x) strsplit(x, "\\.")[[1]][1]))
    return(r2_df)
}


impute_exp <- function(gene, dosage_list_sig, weight_list, ref_table, i) {
    chr <- ref_table %>%
        filter(`Gene stable ID` == !!gene) %>%
        select(`Chromosome/scaffold name`)
    weight_df <- weight_list[[i]]
    tar_snp <- weight_df %>%
        filter(ensembl_id == !!gene) %>%
        select(rsid, weight)
    # target
    snp_vec <- tar_snp$rsid
    weight <- tar_snp$weight
    print(snp_vec)
    print(weight)
    
    # dosage
    chr <- as.numeric(chr[1,1])
    cur_snp_vec <- dosage_list_sig[[i]][[chr]][["snp_list"]]
    cur_dosage <- dosage_list_sig[[i]][[chr]][["dosage"]]
    cur_dosage[is.na(cur_dosage)] <- 0
    
    
    # id match
    weight_vec <- rep(0, length(cur_snp_vec))
    midx <- match(snp_vec, cur_snp_vec)
    nnmidx <- !is.na(midx)
    nnwidx <- midx[nnmidx]
    weight_vec[nnwidx] <- weight[nnmidx] 
    
    # feature
    #display_markdown(paste0(dim(cur_dosage)[1]))
    #display_markdown(paste0(dim(cur_dosage)[2]))
    #display_markdown(paste0(length(weight_vec)))
    #display_markdown(paste0(weight_vec))
    #cat(weight_vec)
    feature_vec <- as.matrix(cur_dosage) %*% weight_vec
    return(feature_vec)
}



# reading dosage
dosage_dir <- "~/scratch60/UKB/ibd/"
setwd(dosage_dir)

snp_list <- list()
dosage_list <- list()

for (i in 1:22) {
    idx <- i
    print(paste0("INFO: loading chr", i))
    dosage_file_prefix <-"ukb_qc1_ibd_2000_chr"
    dosage_file <- paste0(dosage_file_prefix, idx, ".RData")
    load(dosage_file)
    core_snp_vec <- colnames(dose_df)[2:ncol(dose_df)]
    #colnames(dose_df) <- c("ID", core_snp_vec)
    dosage_list[[i]] <- dose_df
    snp_list[[i]] <- core_snp_vec
}

ref_table <- fread("~/project/TiPred/mart_export-2.txt")
sig_table <- fread("~/project/TiPred/ibd_signature_JCC_v1.txt", header = T)
sig_table_c <- sig_table[which(sig_table$`Chromosome/scaffold name` != "X"),]
sig_vec_ens <- sig_table_c$`Gene stable ID`

# weights
setwd("~/project/TiPred/weight/")

weight_si = fread("Small_Intestine_Terminal_Ileum.elnt_allTissues_weight.txt")
weight_ct = fread("Colon_Transverse.elnt_allTissues_weight.txt")
weight_cs = fread("Colon_Sigmoid.elnt_allTissues_weight.txt")
weight_bl = fread("Whole_Blood.elnt_allTissues_weight.txt")

weight_list <- list()

weight_list[[1]] <- add_ensembl_id(weight_si)
weight_list[[2]] <- add_ensembl_id(weight_ct)
weight_list[[3]] <- add_ensembl_id(weight_cs)
weight_list[[4]] <- add_ensembl_id(weight_bl)

# sebset of dosage
gene_pool <- sig_vec_ens
weight_list_sig <- list()
dosage_list_sig <- list()
dosage_df_list <- list()
for (i in 1:22) {
    print(paste0("INFO: data frame chr",i))
    dosage_df_list[[i]] <- as.data.frame(dosage_list[[i]]) 
    dosage_list[[i]] <- NA
}

for (i in 1:4) {
    weight_df <- weight_list[[i]] %>% 
        filter(ensembl_id %in% gene_pool)
    weight_list_sig[[i]] <- weight_df
    dosage_list_sig[[i]] <- list()
    for (j in 1:22) {
        print(paste0("INFO: collect info chr",j))
        snp_set <- unique(weight_df$rsid)
        cur_dose_df <- dosage_df_list[[j]]
        midx <- match(snp_set, c("ID", snp_list[[j]]))
        nnmidx <- midx[!is.na(midx)]
        nnsnp <- snp_set[!is.na(midx)]
        dosage_list_sig[[i]][[j]] <- list()
        dosage_list_sig[[i]][[j]][["snp_list"]] <- nnsnp
        dosage_list_sig[[i]][[j]][["dosage"]] <- cur_dose_df[, nnmidx]
        
    }
}

feature_mtx <- c(rep(1, 2000), rep(0, 2000))

for (i in 1:length(gene_pool)) {
    display_markdown(paste0("gene ", i))
    for (j in 1:4) {
        feature_vec <- impute_exp(gene_pool[i], dosage_list_sig, weight_list_sig, ref_table, j)   
        feature_mtx <- cbind(feature_mtx, feature_vec)
    }
}

save(feature_mtx, file = "~/project/TiPred/data/ibd_2000_jccsig_en.RData")







