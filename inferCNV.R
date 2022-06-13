rm(list=ls())

Sys.setenv(lang="en")
library(rjags)
library(infercnv)
library(biomaRt)
library(dplyr)
library(Seurat)
library(tidyverse)
library(data.table)

#set working directory


#InferCNV is an algortihm to infer if a given cell is tumor or normal by using single cell RNAseq data. 
#The decsription and documentation can be found at: https://github.com/broadinstitute/inferCNV/wiki
#The algorithm requires: 
# - a count matrix that this script extracts from a Seurat objects and provides it formatted without writing to a file
# - a metadata file containing CellID and cell type or another label. Cells should be a known normal used as a refence (cand also be firborblast or endothelial cells)
#and the suspected tumor cells. 
# - a file containing the genes and their chromosomes locations ordered by chromosomes generated from ensemble
#Both files are tab separated files with no column titles

#It is recommended where possible that the algortihm is run by patient in order to mitigate inter-patient variability and long excution times.  


#----------------------------------------------------------------------------
# #Create Infercnv Object # -----------------------------------------------
#----------------------------------------------------------------------------


# 1. raw_counts matrix 
# extract from Seurat object 
# Read rds dataset

setwd("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/ToTransfer")
source("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/SingleCellSeuratBased/scRNAseq-unimodal-Seurat-based-functions.R")
source("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/Rprojects/CITEseq-workflows/CITEseq-functions-seurat-based.R")

DotPlot(SeuratObj, features=c("TPBG", "EPCAM", "SLC39A7"), split.by = "PatientNumber",cols=c("Red", "Blue", "Green", "Gray", "Orange", "Yellow","Cyan"))


SeuratObj <-readRDS("H:/data/10x datasets/Seurat objects old/Lambrechts_CRC-CL-original-cell-type.RDS")
ref_group_names=c("B cell","CD4 T cell", "CD8 T cell", "endothelial cell", "fibroblast",  "monocyte/macrophage", "natural killer cell")
tumor_group_names=c("epithelial cell", "other")

SeuratObj <- subset(SeuratObj, cell.type %in% c(ref_group_names, tumor_group_names))

raw_counts_matrix <- as.matrix(SeuratObj@assays[["RNA"]]@counts)

# 2. cell annotation files

MetaData <- SeuratObj@meta.data

#------------------------------------------------------
# 3. gene order file (order on Chromosomes)
#------------------------------------------------------
all.genes <- data.frame(hgnc_symbol=rownames(raw_counts_matrix))

# retrieve chromosomes positions given a list of genes 

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results_id <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'chromosome_name',
                     'start_position', 'end_position'),
      filters = "hgnc_symbol", 
      values = all.genes$hgnc_symbol, 
      mart = ensembl)

chromo_list <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19",
                 "20", "21", "22", "X", "Y")

results <- results_id %>% filter(chromosome_name %in% chromo_list)
results <- results %>% dplyr::select(hgnc_symbol, chromosome_name, start_position, end_position)

#check if any duplicates in gene position
rep_gene <- data.frame(table(results$hgnc_symbol))
results[results$hgnc_symbol %in% rep_gene[rep_gene$Freq>1, ]$Var1, ]

#clear replicates
results_unique <- results[!duplicated(results$hgnc_symbol), ]

# write table of gene notations
#write.table(results_unique, paste(getwd(), "data/output/gene_chromopos.txt", sep = "/"),
#            col.names = FALSE, row.names = FALSE,  sep = "\t", quote = FALSE)

# filter the counts matrix according to results of chromosome positions
counts_matrix <- raw_counts_matrix[c(results$hgnc_symbol), ]
# write.table(counts_matrix, file = paste(getwd(), "data/output/cnt_matrix", sep = "/")
#              , sep = "\t", col.names= FALSE, row.names = FALSE)

rm(raw_counts_matrix)
rm(SeuratObj)
out_dir="C:/data/inferCNV/GSE131907-batch1/"
#-------------------------------------------------------------------------------
# Create InferCNV object and run ------------------------------------------
#-------------------------------------------------------------------------------




#MetaData$Sample_type_patient <- paste(MetaData$Sample_type, MetaData$orig.ident, sep="-")
MetaData$Cell_ID <- rownames(MetaData)
#for (i in 1:nrow(MetaData)) {if (MetaData$cell.type[i]=="common lymphoid progenitor" & MetaData$mutation_annotation[i]=="mutation detected") MetaData$cell.type[i]<-"common lymphoid progenitor tumor"}
for (i in 1:nrow(MetaData)) {if (MetaData$cell.type[i] %in% ref_group_names ) MetaData$cell.type[i]<-"reference"}
for (i in 1:nrow(MetaData)) {if (MetaData$cell.type[i] %in% tumor_group_names ) MetaData$cell.type[i]<-"observation"}
#write.table(MetaData[,c(which(colnames(MetaData) %in% c("Cell_ID")), which(colnames(MetaData) %in% c("cell.type")))],"C:/data/inferCNV/metadata-for-infer-cnv.txt", row.names = FALSE, sep="\t", quote=FALSE)

results_unique$chromosome_name <- paste("chr", results_unique$chromosome_name, sep="")
#write.table(results_unique,"C:/data/inferCNV/gene-annotation-for-inferCNV.txt", row.names = FALSE, sep="\t", quote=FALSE)
#fwrite(counts_matrix,"C:/data/10x datasets/Seurat objects/count-for-inferCNV.txt", row.names = TRUE, sep="\t", quote=FALSE)



closeAllConnections()



#out_dir <- paste(getwd(), "data/output/InferCNV/", sep = "/")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="C:/data/inferCNV/metadata-for-infer-cnv.txt",
                                    gene_order_file= "C:/data/inferCNV/gene-annotation-for-inferCNV.txt", 
                                    ref_group_names=c("reference")) 

rm(counts_matrix)
gc()

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             analysis_mode = "samples",
                             num_threads = 3,
                             cluster_by_groups=FALSE, #if true clusters by cell types as called previously rather than by CNv 
                             denoise=TRUE,
                             HMM=TRUE)


