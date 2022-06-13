library(Seurat)
library(tidyverse)
library(dplyr)
library(SingleR)
library(SingleCellExperiment)
library(SeuratDisk)
library(scran)
library(DropletUtils)
library(Azimuth)


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
# call Azimuth to identify cell types -
#function has l1 and l2 resolution cell types renamed according to the cell ontology: output will be cell.type and cell.type.detailed
#mapping to cell ontology is the one  on the Azimuth webiste: https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC
#markers of each cell type can be found at the same link and the original azimuth cell types
#Here we will return cell types according to cell ontology and MP ontology, with the cell ontology carrying CL at the end. 
#als the path to the reference file is included here
#===============================================
#================================================
Seurat.Azimuth.celltypes <- function(SeuratObj){
  
  
  reference <- LoadH5Seurat("H:/data/10x datasets/Seurat objects/multi.h5seurat")

  DefaultAssay(object = SeuratObj) <- "SCT"
    
  anchors <- FindTransferAnchors(
    reference = reference,
    query = SeuratObj,
    recompute.residuals = FALSE,
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
  SeuratObj@meta.data$cell.type <-  SeuratObj@meta.data$predicted.celltype.l1
  SeuratObj@meta.data$cell.type.detailed <-  SeuratObj@meta.data$predicted.celltype.l2
  
  

  SeuratObj@meta.data$cell.type <- recode(SeuratObj@meta.data$cell.type, "mono"="monocyte", "CD4 T" = "CD4 T cell", "B" = "B cell", "other"="other", "DC"="myeloid dendritic cell", "other T"= "other T cell", "CD8 T"="CD8 T cell", "NK"="natural killer cell")
  SeuratObj@meta.data$cell.type.detailed <- recode(SeuratObj@meta.data$cell.type.detailed, "ASDC"= "plasmacytoid dendritic cell", "B intermediate"="mature B cell", 
                                                      "B memory"="memory B cell", "B naive"="naive B cell", "CD14 Mono"="CD14-positive, CD16-negative classical monocyte", "CD16 Mono"=
                                                        "CD14-low, CD16-positive monocyte", "CD4 CTL"="cytotoxic CD4 T cell", "CD4 Naive"="naive CD4 T cell", "CD4 Proliferating"="proliferating CD4 T cell", 
                                                   "central memory CD4 T cell"= "CD4-positive, alpha-beta memory T cell", "CD4 TEM"="effector memory CD4 T cell", "CD8 Naive"="naive CD8 T cell", "CD8 Proliferating"="proliferating CD8 T cell", 
                                                   "CD8 TCM"="central memory CD8 T cell", "CD8 TEM"="effector memory CD8 T cell", "cDC1"="CD141-positive myeloid dendritic cell", "cDC2"="CD1c-positive myeloid dendritic cell","dnT"="double negative T cell", 
                                                   "Doublet"="doublet", "Eryth"="erythroid lineage cell", "gdT"="gamma-delta T cell", "HSPC"="hematopoietic precursor cell", "ILC"="innate lymphoid cell", "MAIT"="mucosal invariant T cell", "NK"="natural killer cell", "NK Proliferating"="proliferating natural killer cell", 
                                                   "NK_CD56bright"="natural killer cell", "pDC"="plasmacytoid dendritic cell", 
                                                   "Plasmablast"="plasmablast", "Platelet"="platelet", "Treg"="regulatory T cell")
  
  
   return(SeuratObj) 
  
}


#==================================================
# call Azimuth to identify cell types for bone marrow, using bmcite as reference
#function has l1 and l2 resolution cell types renamed according to the cell ontology: output will be cell.type and cell.type.detailed
#mapping to cell ontology is the one  on the Azimuth webiste: https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC
#markers of each cell type can be found at the same link and the original azimuth cell types
#Here we will return cell types according to cell ontology and MP ontology, with the cell ontology carrying CL at the end. 
#als the path to the reference file is included here
#===============================================
#================================================
Seurat.Azimuth.bonemarrow <- function(query){
  
  #query <- subset(query,cell_names %in% MetaData$cell_names[150001:200000] )
  
 
  
  
  
  reference <-  LoadReference(path = "H:/data/10x datasets/reference-data/human_bonemarrow/")
  
  DefaultAssay(object = query) <- "SCT"
  
  # Find anchors between query and reference
  anchors <- FindTransferAnchors(
    reference = reference$map,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "SCT",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
    dims = 1:50,
    n.trees = 20,
    mapping.score.k = 100
  )
  
  # Transfer cell type labels and impute protein expression
  #
  # Transferred labels are in metadata columns named "predicted.*"
  # The maximum prediction score is in a metadata column named "predicted.*.score"
  # The prediction scores for each class are in an assay named "prediction.score.*"
  # The imputed assay is named "impADT" if computed
  
  refdata <- lapply(X = "celltype.l2", function(x) {
    reference$map[[x, drop = TRUE]]
  })
  names(x = refdata) <- "celltype.l2"
  if (FALSE) {
    refdata[["impADT"]] <- GetAssayData(
      object = reference$map[['ADT']],
      slot = 'data'
    )
  }
  query <- TransferData(
    reference = reference$map,
    query = query,
    dims = 1:50,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE
  )
  
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference$map,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE
  )
  
  # Calculate the query neighbors in the reference
  # with respect to the integrated embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference$map[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = reference$map[[]]
  )
  
  query@meta.data$cell.type <- recode(query@meta.data$predicted.celltype.l2, "ASDC"= "plasmacytoid dendritic cell", "B intermediate"="B cell", 
                                      "Memory B"="B cell", "Naive B"="B cell", "CD14 Mono"="monocyte", "CD16 Mono"=
                                        "monocyte", "CD4 Effector"="cytotoxic CD4 T cell", "CD4 Naive"="naive CD4 T cell", "CD4 Proliferating"="proliferating CD4 T cell", 
                                      "CD4 Memory"="effector memory CD4 T cell", "CD4 TEM"="effector memory CD4 T cell", "CD8 Naive"="naive CD8 T cell", "T Proliferating"="proliferating CD8 T cell", 
                                      "CD8 TCM"="central memory CD8 T cell", "CD8 Memory"="effector memory CD8 T cell", "cDC1"="myeloid dendritic cell", "cDC2"="myeloid dendritic cell","dnT"="double negative T cell", 
                                      "Doublet"="doublet", "Eryth"="erythroid lineage cell", "gdT"="gamma-delta T cell", "HSPC"="hematopoietic precursor cell", "ILC"="innate lymphoid cell", "MAIT"="mucosal invariant T cell", "NK"="natural killer cell", "NK Proliferating"="proliferating natural killer cell", 
                                      "NK CD56+"="natural killer cell", "pDC"="plasmacytoid dendritic cell", 
                                      "Plasmablast"="plasmablast", "Platelet"="platelet", "Treg"="regulatory T cell",  "BaEoMa"= "basophil mast progenitor cell", "CD8 Effector_1"="effector  CD8 T cell","CD8 Effector_2"="effector  CD8 T cell", "CLP"="common lymphoid progenitor", "Early Eryth"="erythroid lineage cell","EMP"="megakaryocyte-erythroid progenitor cell", 
                                      "GMP"="granulocyte monocyte progenitor", "HSC"="hematopoietic stem cell", "ILC"="innate lmyphoid cell","Late Eryth"="erythroid lineage cell",
                                      "LMPP"="multipotent progenitor cell", "Macrophage"="macrophage", "Plasma"="plasma cell", "pre-mDC"="common dendritic progenitor", 
                                      "pre-pDC"="immature plasmacytoid dendritic cell", "pre B"="pro-B cell", "pro B"="pro-B cell", "Prog Mk"="megakaryocyte progenitor cell",
                                      "transitional B"="transitional stage B cell", "Stromal"="stromal cell")
  
  query@meta.data$cell.type <- recode(query@meta.data$cell.type, "CD14-low, CD16-positive monocyte"="monocyte/macrophage", "CD14-positive, CD16-negative classical monocyte"=
                                        "monocyte/macrophage", "CD141-positive myeloid dendritic cell"="myeloid dendritic cell", "CD1c-positive myeloid dendritic cell"="myeloid dendritic cell", 
                                      "cytotoxic CD4 T cell"="CD4 T cell", "effector  CD8 T cell"="CD8 T cell", "effector memory CD4 T cell"="CD4 T cell", "effector memory CD8 T cell"="CD8 T cell", 
                                      "macrophage"="monocyte/macrophage","CD8 Effector_3"="CD8 T cell" , "memory B cell"="B cell", "naive B cell"="B cell", "naive CD4 T cel"="CD4 T cell", "naive CD8 T cell"="CD8 T cell", "naive CD4 T cell"="CD4 T cell")
  
  
  
  query@meta.data$cell.type.coarse <- recode(query@meta.data$cell.type, "basophil mast progenitor cell"="myeloid progenitor", " common dendritic progenitor"="myeloid progenitor", 
                                                    "granulocyte monocyte progenitor"="myeloid progenitor", "immature plasmacytoid dendritic cell"="myeloid progenitor", 
                                                    "innate lymphoid cell"="natural killer cell", "megakaryocyte-erythroid progenitor cell"="myeloid progenitor", "megakaryocyte progenitor cell"="myeloid progenitor", 
                                                    "mucosal invariant T cell"="T cell", "multipotent progenitor cell"="hematopoietic stem cell","common lymphoid progenitor"="lymphoid progenitor", 
                                                    "pro-B cell"="lymphoid progenitor", "proliferating CD8 T cell"="T cell", "CD4 T cell"="T cell", "CD8 T cell"="T cell", "common dendritic progenitor"="myeloid progenitor", 
                                                    "monocyte/macrophage"="monocyte", "proliferating natural killer cell"="natural killer cell", "transitional stage B cell"="B cell")
  
  
  return(query)
}




#==========================================
# do clustering and resampling on a pre.normalized seurat object
#==========================================================

Seurat.clustering.resampling <- function(SeuratObj, clusterRes, Nresamp, filename){

 
  
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
    
    #SeuratObj <- SCTransform(SeuratObj,  vars.to.regress = "percent.mt", verbose = FALSE)
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


#=========================================
  # do clustering and resampling on a pre.normalized sand MNN aligned eurat object
  #==========================================================

Seurat.clustering.resampling.mnn <- function(SeuratObj, clusterRes, Nresamp, filename){
  
  
  
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
    
    #SeuratObj <- SCTransform(SeuratObj,  vars.to.regress = "percent.mt", verbose = FALSE)
    #SeuratObj_subset <- NormalizeData(SeuratObj_subset)
    #SeuratObj_subset <- FindVariableFeatures(SeuratObj_subset, selection.method = "vst", nfeatures = 2000)
    #SeuratObj_subset <- ScaleData(object = SeuratObj_subset, use.umi = TRUE)
    #SeuratObj_subset <- RunPCA(object = SeuratObj_subset, pc.genes = SeuratObj_subset@var.genes, do.print = TRUE, pcs.print = 1:5, 
    #                           genes.print = 5)
    #SeuratObj_subset<- RunTSNE(object = SeuratObj_subset, dims.use = 1:15, do.fast = TRUE,check_duplicates = FALSE)
    SeuratObj_subset <- FindNeighbors(SeuratObj_subset,reduction = "mnn",, dims = 1:30) 
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
