
#----------------------------------------------------------------------------------------------------
#' Combines the sces in the sce list 
#' Helper function to dsb_normalize_protein_data. We store the normalized protein data in sce objects
#' so had to create them
#----------------------------------------------------------------------------------------------------
make_group_sce <- function(sce_list){
  convert_assay_to_matrix <- function(sce, assay_name) { # to convert specified assay to matrix class bc scMerge needs it that way.
    assay_data <- assay(sce, assay_name)
    assay(sce, assay_name) <- as.matrix(assay_data)
    return(sce)
  }
  
  
  sce_list <- lapply(sce_list, convert_assay_to_matrix, assay_name = "NGT")
  # Combine the main experiments using sce_cbind
  group_sce <- sce_cbind(sce_list = sce_list, method = "union", exprs = "NGT")
  protein_altExps_list <- list() #  a list to store 'Protein' alternative experiments
  for (i in seq_along(sce_list)) {# Extract 'Protein' alternative experiments from each SingleCellExperiment object
    sce <- sce_list[[i]]
    if ("Protein" %in% altExpNames(sce)) {
      protein_altExps_list[[i]] <- altExp(sce, "Protein")
    } else {
      stop(paste("The 'Protein' altExp is missing in SCE object at index", i))
    }
  }
  protein_features <- lapply(protein_altExps_list, rownames)  # Check if all 'Protein' alternative experiments have the same features
  if (!all(sapply(protein_features, identical, protein_features[[1]]))) {
    print("Not all 'Protein' altExps have the same features.")
  }
  combined_protein_altExp <- sce_cbind(  # Combine the 'Protein' alternative experiments
    sce_list = protein_altExps_list,method = "union", exprs = "Protein" 
  )
  altExp(group_sce, "Protein") <- combined_protein_altExp  # Add the combined 'Protein' alternative experiment to the main SingleCellExperiment object
  return(group_sce)
}

#----------------------------------------------------------------------------------------------------
# metadata is droplet metadata. Needed for DSB normalization
# Denoised and Scaled by Background
# CLR AND DSB
#----------------------------------------------------------------------------------------------------
dsb_normalize_protein_data<-function(sce_list, group_metadata, group_background){ 

  group_sce<- make_group_sce(sce_list)# constructed to store the resulting normalized protein data.
  group_protein_sce <- SingleCellExperiment::altExp(group_sce,"Protein")#   in order to store the resulting normalized protein data.
  group_protein_matrix<- list()
  group_empty_drops_matrix_input<- list()

  # gets the protein matrix and identifies background droplets for each sample sce. 
  for(i in seq_along(sce_list)){
    sample_name<- sce_list[[i]]@metadata[["sample_name"]]
    protein_sce <- SingleCellExperiment::altExp(sce_list[[i]],"Protein")
    protein_mat <- protein_sce@assays@data$Protein # cells are the columns and the rows are the proteins like CD11b, CD11c, etc. 
    cat(blue("Getting protein matrix for",  sce_list[[i]]@metadata[["sample_name"]], "which has", length(colnames(sce_list[[i]])),"cells: \n"))
    print(protein_mat[1:3, 1:3])
    group_protein_matrix[[i]] <- protein_mat
    #=========================identify and extract background (empty droplets)===================
    file_name<-sce_list[[i]]@metadata$file; file_name
    cells_of_interest<-colnames(protein_mat); cells_of_interest
    all_protein_droplets <- sce_list[[i]]@metadata[["protein_droplet_subset"]] # rhdf5::h5read(file=file, name="/all_barcodes/protein_read_counts/layers/read_counts")
    
    all_protein_droplets<- as.data.frame(all_protein_droplets)
    
    empty_drops_matrix_input <- data.frame(all_protein_droplets) %>% 
      dplyr::select(any_of(group_background)) # select columns corresponding to background_droplets
    rownames(empty_drops_matrix_input)<- rownames(protein_mat)
    group_empty_drops_matrix_input[[i]] <- empty_drops_matrix_input
 
    #=================================================================================================
  }
  # combine the samples' protein data.
  #final_group_protein_matrix <- bind_cols(lapply(group_protein_matrix, as.data.frame)) # can get directly from roup_protein_sce@assays@data$Protein? why update? 
  final_group_protein_matrix<- do.call(cbind, group_protein_matrix)
  rownames(final_group_protein_matrix) <- rownames(group_protein_matrix[[1]])
  # update the group_protein_sce
  group_protein_sce@assays@data$Protein<- final_group_protein_matrix # i think can get directly. no need to update. check
  final_group_empty_drops_matrix_input <- do.call(cbind, group_empty_drops_matrix_input)

  rownames(final_group_empty_drops_matrix_input) <- rownames(group_empty_drops_matrix_input[[1]])
  
  isotype <- grep("IgG",rownames(final_group_protein_matrix),value=TRUE) # vector of isotype control names
  cat(blue("Performing DSB Normalization..."))
  adt_norm <- dsb::DSBNormalizeProtein( # remove ambient protien noise reflected in counts from empty droplets
    cell_protein_matrix = final_group_protein_matrix, # cell-containing droplet raw protein count matrix
    empty_drop_matrix = final_group_empty_drops_matrix_input, # empty/background droplet raw protein count
    denoise.counts = TRUE, # (default = TRUE); run step II
    use.isotype.control = TRUE, # (default = TRUE): use isotype controls to define technical components.
    isotype.control.name.vec = isotype # isotype controls (IgG markers) to adjust for non-specific antibody binding
  )
  # add DSB norm data to group_protein_sce which is then copied to group_sce in the Protein assay of altExp
  SummarizedExperiment::assay(group_protein_sce, "DSB_norm")<-adt_norm
  
  cat(blue("Performing CLR normalization..."))
  group_protein_matrix<- group_protein_sce@assays@data$Protein 
  s <- Seurat::CreateAssayObject(group_protein_matrix,assay="Protein")
  s <- Seurat::CreateSeuratObject(s,assay="Protein")
  s <- Seurat::NormalizeData(s, normalization.method = "CLR")
  SummarizedExperiment::assay(group_protein_sce, "CLR_norm") <- s@assays$Protein@data # add the CLR data to the protein sce. the protein sce is then added to the main group_sce later. 
  
  group_protein_sce@metadata<-group_metadata
  group_protein_sce@colData<-S4Vectors::DataFrame(group_metadata%>% dplyr::filter(Cell%in%colnames(final_group_protein_matrix)))
  SingleCellExperiment::altExp(group_sce, "Protein") <- group_protein_sce
  cat(green("✅ Normalization Complete."))
  return(group_sce)
}
