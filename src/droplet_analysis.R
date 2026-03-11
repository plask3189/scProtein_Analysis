#----------------------------------------------------------------------------------------------------
# Helper function to group_droplete_data
# only thing that is different is that it adds the sample name to the barcode.
#----------------------------------------------------------------------------------------------------
extract_droplet_size_kp<- function(sce){
  sample_name<-sce@metadata[["sample_name"]]
  file<-sce@metadata$file
  all_protein_droplets <- rhdf5::h5read(file=file,name="/all_barcodes/protein_read_counts/layers/read_counts")
  all_dna_droplets <- rhdf5::h5read(file=file,name="/all_barcodes/dna_read_counts/layers/read_counts")
  colnames(all_dna_droplets) <-rhdf5::h5read(file=file,name="/all_barcodes/dna_read_counts/ra/barcode")
  colnames(all_dna_droplets) <- paste(sample_name, colnames(all_dna_droplets), sep = "_")
  colnames(all_protein_droplets) <-rhdf5::h5read(file=file,name="/all_barcodes/protein_read_counts/ra/barcode")
  colnames(all_protein_droplets) <- paste(sample_name, colnames(all_protein_droplets), sep = "_")
  # create a metadata dataframe of simple qc stats for each droplet
  dna_size <- data.frame("Cell"=colnames(all_dna_droplets),
                         "dna_size"=log10(Matrix::colSums(all_dna_droplets)),
                         "amplicons"=Matrix::colSums(all_dna_droplets > 0))
  protein_size <- data.frame("Cell"=colnames(all_protein_droplets),
                             "proteins"=Matrix::colSums(all_protein_droplets > 0),
                             "protein_size"=log10(Matrix::colSums(all_protein_droplets)))
  md <- dplyr::inner_join(dna_size,protein_size)%>% dplyr::mutate(Cell=gsub("-1","",Cell))
  md<- md%>% mutate(Droplet_type=case_when( Cell%in%colnames(sce)~"Cell", TRUE~"Empty" ))%>% full_join(sce@metadata$NGT,by="Cell")
  return(md)
}#----------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# merge each sce object's droplet metadata. 
#----------------------------------------------------------------------------------------------------
group_droplete_data<- function(sce_list){
  droplet_metadata_list <- list()
  merged_sample_names<- list()
  h5closeAll()
  for (i in seq_along(sce_list)) {
    group <- sce_list[[i]]@metadata[["group"]]
    cat(silver("   Analyzing", sce_list[[i]]@metadata[["sample_name"]], "metadata in Group", group,  "... \n"))
    droplet_metadata <- extract_droplet_size_kp(sce_list[[i]]) %>%
      dplyr::select(c("Cell", "dna_size", "amplicons", "proteins", "protein_size", "Droplet_type", "Group", "Clone"))
    cat(silver("   Preview of droplet metadata for Sample",sce_list[[i]]@metadata[["sample_name"]] , "\n"))
    print(head(droplet_metadata))
    droplet_metadata$Group <- group
    merged_sample_names<- c(merged_sample_names, sce_list[[i]]@metadata[["sample_name"]])
    droplet_metadata_list[[i]] <- droplet_metadata
    rm(droplet_metadata)
    gc()  # Garbage collection
  }
  cat(blue("✅ Merged droplet metadata from ", paste(merged_sample_names, collapse = ", "), "\n"))
  group_droplet_metadata <- dplyr::bind_rows(droplet_metadata_list)
  cat(blue("Classifying if each droplet is 'Empty' or is a 'Cell' \n"))
  group_droplet_metadata<-GetMembershipFaster(group_droplet_metadata) # Classifies each droplet into “Cell” or “Empty” based on DNA size and protein size.
  cat(blue("Computing stain index \n"))
  group_droplet_metadata<-GetStainIndex(group_droplet_metadata) # computes the likelihood that it belongs to the cell distribution or empty droplet distribution using the multivariate normal probability density function (PDF)
  print(head(group_droplet_metadata))
  cat(green("----Finished making the group droplet metadata---- \n"))
  return(group_droplet_metadata)
}#----------------------------------------------------------------------------------------------------
