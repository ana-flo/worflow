library(Seurat)
library(tidyverse)
library(dplyr)
library(SingleR)
library(SingleCellExperiment)
library(SeuratDisk)
library(scran)
library(DropletUtils)


source("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/RLibraries/ClusteringFunctionLibraryRNAseq.R")



#========================================================
#use DropletUtils to remove empty dropplets and return the count matrix only containing the "full" droplets
#==================================================
RemoveEmptyDrops.Barcodes <- function(count.matrix){

  #as this is only for uni modal data, count matrix should be a matrix, not an element of a list of matrices
  
  set.seed(100)
  e.out <- emptyDrops(count.matrix)
  is.cell <- e.out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  cell.barcodes <- colnames(count.matrix)[which(is.cell==TRUE)]
  
  
  plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
       xlab="Total UMI count", ylab="-Log Probability")
  
  count.matrix <- count.matrix[,which(colnames(count.matrix) %in% cell.barcodes)]  
  return(count.matrix)
}



#==============================================
#call SingleR on clusters to assign cell types (data is SCTransformed)
#===================================================

Seurat.Singler.SCT <- function(SeuratObj, ref) {
  
  #requires reference file input for SingleR so the same function can be used for human and mice 
  #also requires seurat clusters to be present
  
  tt <-SeuratObj@assays[["SCT"]]@data
  singler1<-SingleR(tt, ref=df.encode, labels = df.encode$cell.type,  method = c("cluster"),
                    clusters=SeuratObj@meta.data$seurat_cluster,genes = "de", quantile = 0.8, fine.tune = TRUE,
                    tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                    check.missing = TRUE)
 
  ClusterCellTypes <- data.frame(singler1@listData[["pruned.labels"]])
  Clusters <<- data.frame(seurat_cluster=as.numeric(rownames(ClusterCellTypes))-1, cell.type=ClusterCellTypes[1])

  colnames(Clusters)[2]<-"cell.type"
  rm(tt)
  
  
  
  MetaDataM <- data.frame(CellID=rownames(SeuratObj@meta.data),SeuratObj@meta.data)
  #colnames(MetaDataM)[29]<-"cluster_subtype"
  MetaDataM$cluster.type.singler <- "unclassified"
  for (i in 1: nrow(MetaDataM)){
    for (k in 1: nrow(Clusters)){
      if (is.na(Clusters$cell.type[k])) Clusters$cell.type[k]<-"unclassified"
      if (MetaDataM$seurat_cluster[i]==Clusters$seurat_cluster[k]) MetaDataM$cluster.type.singler[i]<-Clusters$cell.type[k]
      
      
    } 
  }
  SeuratObj <- AddMetaData(SeuratObj, MetaDataM)

  
  return(SeuratObj)
  
}


#==============================================
#call SingleR on  to assign cell types (data is SCTransformed)
#===================================================

Seurat.Singler.SCT.cells <- function(SeuratObj, ref) {
  
  #requires reference file input for SingleR so the same function can be used for human and mice 
  #also requires seurat clusters to be present
  
  tt <-SeuratObj@assays[["SCT"]]@data
  singler1<-SingleR(tt, ref=df.encode, labels = df.encode$cell.type,  method = c("single"),
                    clusters=SeuratObj@meta.data$seurat_cluster,genes = "de", quantile = 0.8, fine.tune = TRUE,
                    tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                    check.missing = TRUE)
  
  ClusterCellTypes <- data.frame(singler1@listData[["pruned.labels"]])
  Clusters <<- data.frame(seurat_cluster=as.numeric(rownames(ClusterCellTypes))-1, cell.type=ClusterCellTypes[1])
  
  colnames(Clusters)[2]<-"cell.type"
  rm(tt)
  
  
  
  MetaDataM <- data.frame(CellID=rownames(SeuratObj@meta.data),SeuratObj@meta.data)
  #colnames(MetaDataM)[29]<-"cluster_subtype"
  MetaDataM$cluster.type.singler <- "unclassified"
  for (i in 1: nrow(MetaDataM)){
    for (k in 1: nrow(Clusters)){
      if (is.na(Clusters$cell.type[k])) Clusters$cell.type[k]<-"unclassified"
      if (MetaDataM$seurat_cluster[i]==Clusters$seurat_cluster[k]) MetaDataM$cluster.type.singler[i]<-Clusters$cell.type[k]
      
      
    } 
  }
  SeuratObj <- AddMetaData(SeuratObj, MetaDataM)
  
  
  return(SeuratObj)
  
}



#==============================================
#call SingleR on clusters to assign cell types (data is NOT SCTransformed)
#===================================================

Seurat.Singler <- function(SeuratObj, ref) {
  
  #requires reference file input for SingleR so the same function can be used for human and mice 
  #also requires seurat clusters to be present
  
  
  tt <-SeuratObj@assays[["RNA"]]@data
  singler1<-SingleR(tt, ref=df.encode, labels = df.encode$cell.type,  method = c("cluster"),
                    clusters=SeuratObj@meta.data$seurat_cluster,genes = "de", quantile = 0.8, fine.tune = TRUE,
                    tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                    check.missing = TRUE)
  
  
  ClusterCellTypes <- data.frame(singler1@listData[["pruned.labels"]])
  Clusters <- data.frame(seurat_cluster=as.numeric(rownames(ClusterCellTypes))-1, cell.type=ClusterCellTypes[1])
  
  colnames(Clusters)[2]<-"cell.type"
  rm(tt)
  
  
  
  MetaDataM <- data.frame(CellID=rownames(SeuratObj@meta.data),SeuratObj@meta.data)
  #colnames(MetaDataM)[29]<-"cluster_subtype"
  MetaDataM$cluster.type.singler <- "unclassified"
  for (i in 1: nrow(MetaDataM)){
    for (k in 1: nrow(Clusters)){
      if (MetaDataM$seurat_cluster[i]==Clusters$seurat_cluster[k]) MetaDataM$cluster.type.singler[i]<-Clusters$cell.type[k]
      
      
    } 
  }
  SeuratObj <- AddMetaData(SeuratObj, MetaDataM)
  
  
  return(SeuratObj)
  
}


#==================================================
# call Azimuth to identify cell types
#===============================================
#================================================
Seurat.Azimuth.celltypes <- function(SeuratObj){
  
  
  reference <- LoadH5Seurat("H:/data/10x datasets/Seurat objects/multi.h5seurat")

  DefaultAssay(object = SeuratObj) <- "SCT"
    
  anchors <- FindTransferAnchors(
    reference = reference,
    query = SeuratObj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:50
  )
  
  SeuratObj <- MapQuery(
    anchorset = anchors,
    query = SeuratObj,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  
  #rename major cell types to align to cell ontology

  SeuratObj@meta.data$predicted.celltype.l1 <- recode(SeuratObj@meta.data$predicted.celltype.l1, "mono"="monocyte", "CD4 T" = "CD4 T cell", "B" = "B cell", "other"="unclassified", "DC"="myeloid dendritic cell", "other T"= "T cell other", "CD8 T"="CD8 T cell", "NK"="natural killer cell")
  return(SeuratObj) 
  
}


#==========================================
# do clustering and resampling on a pre.normalized seurat object
#==========================================================

Seurat.clustering.resampling <- function(SeuratObj, cluster.res, Nresamp, filename){

 
  
  reference1 <- data.frame(cell.name=rownames(SeuratObj@meta.data), cluster=SeuratObj@meta.data$seurat_clusters)
  
  tt0<-c(0:(length(unique(reference1$cluster))-1))
  #clusterRes <-0.2
  
  for (iresamp in 1:Nresamp){
    
    set.seed(iresamp)
    gc()
    #dfrs <-SampleDataFrameFraction(transposeBigData(DataTumor), 0.2)
    cellList <-  colnames(SeuratObj)
    frac.resampled <- 0.8*length(cellList)
    dfrs <-  sample(cellList,frac.resampled, replace = FALSE)
    
    cells.use <- colnames(SeuratObj)[ colnames(SeuratObj) %in% dfrs]
    SeuratObj_subset <- subset(SeuratObj, cells = cells.use)
    
    SeuratObj <- SCTransform(SeuratObj,  vars.to.regress = "percent.mt", verbose = FALSE)
    #SeuratObj_subset <- NormalizeData(SeuratObj_subset)
    SeuratObj_subset <- FindVariableFeatures(SeuratObj_subset, selection.method = "vst", nfeatures = 2000)
    SeuratObj_subset <- ScaleData(object = SeuratObj_subset, use.umi = TRUE)
    SeuratObj_subset <- RunPCA(object = SeuratObj_subset, pc.genes = SeuratObj_subset@var.genes, do.print = TRUE, pcs.print = 1:5, 
                               genes.print = 5)
    SeuratObj_subset<- RunTSNE(object = SeuratObj_subset, dims.use = 1:15, do.fast = TRUE,check_duplicates = FALSE)
    SeuratObj_subset <- FindNeighbors(SeuratObj_subset, dims = 1:15) 
    SeuratObj_subset <- FindClusters(object = SeuratObj_subset, resolution=clusterRes)
    
    print("interation no:")
    print(iresamp)
    reference2 <-data.frame(cell.name=rownames(SeuratObj_subset@meta.data), cluster=SeuratObj_subset@meta.data$seurat_clusters)
    
    tt <- Do_Jaccard_Distribution_Clustered(reference1, reference2)
    
    attach(tt)
    tt <- tt[order(tt$Cluster),]
    detach(tt)
    tt0 <-rbind(tt0,tt$MaxJ)
    #this is just a dump for every iteration -  no need for it in a working pipeline
    #write.table(tt0[-1,],"test_clusters_Jaccards_RNA_res0.1.txt", sep="\t")
    
  }
  
  tt0 <- tt0[-1,]
  dfJ <- apply(tt0,2,median)
  
  ########
  #this is the final output
  dfJ1 <- data.frame(Cluster=tt$Cluster, Jaccard=dfJ )
  #======================================
 # filename <- paste("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Projects/PTCL/analysis/jaccardIndex_all-t-and-nk_res", clusterRes,".txt", sep="")
  
  write.table(dfJ1, filename, sep="\t", row.names = FALSE)
  
  return(dfJ1)
  
}
