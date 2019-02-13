#  library
library(dplyr)
library(ggplot2)
library(ROCR)
library(data.table)
library(doParallel)
library(stringr)
library(Matrix)
library(RSQLite)
library(DBI)

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

extract_hq_gene <- function(r2_df, cut_off) {
    tmp_df <- r2_df %>%
        filter(rsq >= cut_off) %>%
        select(ensembl_id)
    return(tmp_df$ensembl_id)
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

# reading dosage
dosage_dir <- "~/scratch60/UKB/cad/"
setwd(dosage_dir)

snp_list <- list()
dosage_list <- list()

for (i in 1:22) {
    idx <- i
    #print(i)
    display_markdown(paste0("chr", i))
    dosage_file_prefix <-"ukb_qc1_cad_2000_chr"
    dosage_file <- paste0(dosage_file_prefix, idx, ".RData")
    #dosage_file_c <- paste0(dosage_file_prefix, idx, "_c.RData")
    start_time <- Sys.time()
    display_markdown("loading")
    load(dosage_file)
    end_time <- Sys.time()
    print(end_time - start_time)
    start_time <- Sys.time()
    display_markdown("extracting")
    core_snp_vec <- colnames(dose_df)[2:ncol(dose_df)]
    #dose_df[is.na(dose_df)] <- 0
    end_time <- Sys.time()
    print(end_time - start_time)
    start_time <- Sys.time()
    end_time <- Sys.time()
    print(end_time - start_time)
    dosage_list[[i]] <- dose_df
    snp_list[[i]] <- core_snp_vec
}


# ref table
ref_table <- fread("~/project/TiPred/mart_export-2.txt")

# r2_path
sig_ens_list <- list()
r2_df <- list()
r2_path <- "~/project/TiPred/model_r2/selected_models_elnt_"
tissue_vec <- c("Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Whole_Blood", "Liver")
file_vec <- as.character(sapply(tissue_vec, function(x) paste0(r2_path, x, ".txt")))

# reading r2
for (i in 1:length(file_vec)) {
    r2_df[[i]] <- fread(file_vec[i])
    r2_df[[i]] <- add_ensembl_id(r2_df[[i]])
    sig_ens_list[[i]] <- extract_hq_gene(r2_df[[i]], 0.8)
    sig_ens_list[[i]] <- filter_gene(sig_ens_list[[i]], ref_table)
}

# sebset of dosage
#gene_pool <- sig_vec_ens
weight_list_sig <- list()
dosage_list_sig <- list()
dosage_df_list <- list()
for (i in 1:22) {
    display_markdown(paste0("chr",i))
    dosage_df_list[[i]] <- as.data.frame(dosage_list[[i]]) 
}

# weight
weight_path <- "~/project/TiPred/weight/"
tissue_vec <- c("Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Whole_Blood", "Liver")
weight_file_vec <- as.character(sapply(tissue_vec, function(x) paste0(weight_path, x, ".elnt_allTissues_weight.txt")))

# reading weights
weight_df <- list()
for (i in 1:length(weight_file_vec)) {
    weight_df[[i]] <- fread(weight_file_vec[i])
    weight_df[[i]] <- add_ensembl_id(weight_df[[i]])
    
}

for (i in 1:4) {
    gene_pool <- sig_ens_list[[i]]
    weight_list_sig[[i]] <- weight_df[[i]] %>% 
        filter(ensembl_id %in% gene_pool)
    dosage_list_sig[[i]] <- list()
    for (j in 1:22) {
        display_markdown(paste0("chr",j))
        snp_set <- unique(weight_list_sig[[i]]$rsid)
        cur_dose_df <- dosage_df_list[[j]]
        midx <- match(snp_set, c("ID", snp_list[[j]]))
        nnmidx <- midx[!is.na(midx)]
        nnsnp <- snp_set[!is.na(midx)]
        dosage_list_sig[[i]][[j]] <- list()
        dosage_list_sig[[i]][[j]][["snp_list"]] <- nnsnp
        dosage_list_sig[[i]][[j]][["dosage"]] <- cur_dose_df[, nnmidx]
        
    }
}

# merge all the ens id
sig_ens_vec <- c()
for (i in 1:length(sig_ens_list)) {
    sig_ens_vec <- union(sig_ens_vec, sig_ens_list[[i]])
}

# sig ref_table
ref_table_sig <- ref_table %>%
    filter(`Gene stable ID` %in% sig_ens_vec)

# weight matrix function
build_weight_mtx_hq <- function(weight_list_sig, sig_ens_list, dosage_list, ref_table, n_tissue) {
    #feature_tmx
    feature_mx <- c(rep(1,2000), rep(1, 2000))
    for (i in 1:n_tissue) {
        cur_weight_df <- weight_list_sig[[i]]
        for (j in 1:length(sig_ens_list[[i]])) {
            display_markdown(as.character(i))
            gene_id <- sig_ens_list[[i]][j]
            
            # match the gene
            chr_info <- ref_table %>%
                filter(`Gene stable ID` == gene_id)
            #print(gene_id) 
            chr <- as.numeric(chr_info$`Chromosome/scaffold name`[1])
            #print(chr)
            
            
            cur_snp_vec <- dosage_list[[i]][[chr]][["snp_list"]]
            
            cur_dosage <- dosage_list[[i]][[chr]][["dosage"]]
            start_time <- Sys.time()
            cur_weight <- cur_weight_df %>% filter(ensembl_id == gene_id)
            end_time <- Sys.time()
            print(end_time - start_time)
            start_time <- Sys.time()
            map_idx <- match(cur_weight$rsid, cur_snp_vec)
            end_time <- Sys.time()
            print(end_time - start_time)
            vec_idx <- !is.na(map_idx)
            map_idx_nn <- map_idx[vec_idx]
            print(paste0("match: ", length(map_idx_nn)))
            # start_time <- Sys.time()
            # 
            # dose_mtx <- cur_dosage[, 2:ncol(cur_dosage)]
            # end_time <- Sys.time()
            print(end_time - start_time)
            weight_vec <- rep(0, length(cur_snp_vec))
            start_time <- Sys.time()
            
            weight_vec[map_idx_nn] <- cur_weight$weight[vec_idx]
            end_time <- Sys.time()
            print(end_time - start_time)
            print(length(weight_vec))
            #print(dim(dose_mtx))
            print(length(weight_vec))
            start_time <- Sys.time()
            feature_mtx <- as.matrix(cur_dosage) %*% (weight_vec)
            end_time <- Sys.time()
            print(end_time - start_time)
            #feature_mtx <- as.matrix(dose_mtx) %*% (weight_vec)
            print(dim(feature_mx))
            print(dim(feature_mtx))
            feature_mx <- cbind(feature_mx, feature_mtx)
            
            
        }
    }
    
    return(feature_mx)
}

tmp <- build_weight_mtx_hq(weight_list_sig, sig_ens_list, dosage_list_sig, ref_table, 4)


