
# sce subset is a large sce object with all the samples combined if two_cohorts == TRUE

cluster_protein_data<- function(sce_subset, 
                                 protein_normalization_method = "DSB_norm"){ 

  set.seed(1333)
  
  # create the sce object with the data as whatever protein normalization method. 
  # Normalized data stored at: (sce_subset@int_colData@listData[["altExps"]]@listData[["Protein"]]@se@assays@data@listData[["DSB_norm"]])
  # we use normalized data as main data in seurat.
  s <- Seurat::as.Seurat(x = altExp(sce_subset),counts="Protein",data=protein_normalization_method) # DSB_norm or CLR_norm
  s <- SeuratObject::RenameAssays(object = s, originalexp = 'Protein')
  s <- Seurat::ScaleData(s, assay = "Protein")
  s <- Seurat::FindVariableFeatures(s,assay = "Protein")
  num_variable_features<- length(VariableFeatures(object = s))
  cat(silver("Using", num_variable_features, "variable features \n"))
  
  if (protein_normalization_method == "CLR_norm"){
    cat(blue("Using CLR normalized protein data in dimension reduction and clustering... \n"))
    num_principle_components<-6
    s <- Seurat::RunPCA(s, features = VariableFeatures(object = s), npcs = num_principle_components, approx = FALSE)
    cat(silver("Using", num_principle_components-1, "principle components \n"))
    elbow<- ElbowPlot(s, ndims = 50); print(elbow)
    s <- Seurat::FindNeighbors(s, dims = 1:6) # dims specifies which principal components (PCs) or dimensions of your reduced dimensionality data should be used to calculate the nearest neighbors between cells;
    s <- Seurat::FindClusters(s, resolution = 1.5, algorithm =4)
    s <- Seurat::RunUMAP(s, dims = 1:6, n.neighbors = 70, min.dist = 0.2)
    p <- FeaturePlot(s, features = "CD11b"); print(p)
    DimHeatmap(s, dims = 1:6, cells = 500, balanced = TRUE)
    DimPlot(s,reduction = "umap",label = TRUE, label.size = 6)
  } 
  if (protein_normalization_method == "DSB_norm"){
    cat(blue("Using DSB normalized protein data in dimension reduction and clustering... \n"))
    s <- Seurat::FindVariableFeatures(s,assay = "Protein")
    # cluster and run umap (based directly on dsb normalized values without isotype controls)
    s <- FindNeighbors(object = s, dims = NULL, k.param = 30)
    s <- FindClusters(object = s, # direct graph clustering
                      assay = "Protein",
                      resolution = 0.5,
                      algorithm = 3,
                      verbose = TRUE)
    s <- RunUMAP(object = s,
                 seed.use = 1333,
                 min.dist = 0.4, # if increase, Points are more spread out, emphasizing the broader topology of the data and preserving global relationships.
                 n.neighbors = 30, # if increase Captures more global structures, potentially merging smaller clusters and revealing overarching patterns.
                 features=s@assays$Protein@var.features,
                 verbose = TRUE)
    cat(green("----Finished Clutering----"))
  } else {cat(red("Didn't perform that normalization."))}
  return(s)
}
