rename_cells_with_group_and_sample_identifier<- function(sample_name, group_identifier, sce){
  # like from this "ACAACCTAGCTTACTAG" to this "A0386_2_AACAACCTAGCTTACTAG"
  sce@metadata[["sample_name"]]<-sample_name 
  sce@metadata[["group"]]<-group_identifier 
  new_names <- paste(sample_name, colnames(sce), sep = "_")
  for (assay_name in assayNames(sce)) { # rename the cells (columns) of the assays with the samplename added to the barcode : 
    colnames(assay(sce, assay_name, withDimnames=FALSE)) <- new_names
  }
  for (alt_exp_name in altExpNames(sce)) { # rename the cells (columns) of the altExpNames ("Protein" "CNV" ) with the samplename added to the barcode : 
    alt_exp <- altExp(sce, alt_exp_name)
    for (assay_name in assayNames(alt_exp)) {
      colnames(assay(alt_exp, assay_name, withDimnames=FALSE)) <- new_names
    }
    altExp(sce, alt_exp_name) <- alt_exp
  }
  rownames(colData(sce)) <- new_names
  return(sce)
}

identify_missing_markers<- function(sce_list){
  min_length_marker_list<-rownames(altExp(sce_list[[1]], "Protein")) # INIT
  max_length_marker_list<-rownames(altExp(sce_list[[1]], "Protein")) # INIT
  for(i in seq_along(sce_list)){
    marker_list<-rownames(altExp(sce_list[[i]], "Protein"));marker_list
    if(length(marker_list) < length(min_length_marker_list)){
      cat(blue("Found sample with min length marker list:", sce_list[[i]]@metadata[["sample_name"]],"\n"))
      
      min_length_marker_list <- marker_list
    } 
    if(length(marker_list) > length(max_length_marker_list)){
      cat(blue("Found max length marker list \n"))
      max_length_marker_list <- marker_list
    } 
  }
  if (length(max_length_marker_list) != length(min_length_marker_list)){
    missing_markers<- setdiff(max_length_marker_list, min_length_marker_list)
  } else {
    missing_markers <- NA
  }
  
  return(missing_markers)
}



handle_missing_marker_data<- function(sce_list, missing_markers){
  
  if (length(missing_markers)==0){
    cat(green("No markers to remove \n "))
  } 
  for(i in seq_along(sce_list)){
      row_nums_of_missing_markers <- which(rownames(altExp(sce_list[[i]], "Protein")) %in% missing_markers)
      # handle the missing marker in protein data first
      markers_to_keep <- setdiff(rownames(altExp(sce_list[[i]], "Protein")), missing_markers)
      altExp(sce_list[[i]], "Protein") <- altExp(sce_list[[i]], "Protein")[markers_to_keep, ]
      # now handle the droplet data in the h5s
      h5_file_path<- sce_list[[i]]@metadata[["file"]]
      protein_droplets <- rhdf5::h5read(file=h5_file_path, name="/all_barcodes/protein_read_counts/layers/read_counts")
      colnames(protein_droplets) <-rhdf5::h5read(file=h5_file_path, name="/all_barcodes/protein_read_counts/ra/barcode")
      colnames(protein_droplets) <- gsub("-1","",colnames(protein_droplets))
      sample_name<- sce_list[[i]]@metadata[["sample_name"]]
      colnames(protein_droplets) <- paste(sample_name, colnames(protein_droplets), sep = "_")
      # remove marker_to_remove
      beware_sample_names<- c("LAM_1962", "Sample17020", "LAM_3307")
      #if(!(sce_list[[i]]@metadata[["sample_name"]] %in% beware_sample_names)){
      
      if(length(row_nums_of_missing_markers)>0){
        protein_droplets_subset <-protein_droplets[-row_nums_of_missing_markers, , drop=FALSE]
      } else{
        cat(blue("There are no missing markers so we're using all", nrow(protein_droplets), "of them \n"))
        protein_droplets_subset <-protein_droplets
      }
      #} 
      sce_list[[i]]@metadata[["protein_droplet_subset"]]<-protein_droplets_subset
    }
  cat(red("Removed marker:", missing_markers, "\n"))
  markers_using <- rownames(altExp(sce_list[[1]], "Protein"))
  cat(blue("Using markers:", paste(markers_using, collapse = ", "), "\n"))
  h5closeAll()
  return(sce_list)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
process_sces<- function(sample_files, group_identifier = 0, path_to_save_variant_data, base_variant_output_list = NA){
  blacklist <- c( # exclude these variants. not pathogenic.
    "TET2.I1762V", "TET2.961*", "TET2.Y1579*", "TET2.A1283T", "TET2.L1721W",
    "DNMT3A.F772C", "NRAS.L56P", "NRAS.T58A", "NRAS.L56Q", "NRAS.D57N",
    "DNMT3A.I310S", "NRAS.T58I", "DNMT3A.I292S", "DNMT3A.L888Q", "DNMT3A.L888P",
    "PTPN11.L525R", "DNMT3A.N489T", "DNMT3A.K429T", "DNMT3A.N757T", "NRAS.T58P",
    "NRAS.D57Y", "NRAS.Q61P", "TET2.Q618H", "TET2.L1819F", "DNMT3A.Y481S",
    "FLT3.N847T", "TET2.A584T", "TET2.A584P", "DNMT3A.F290L", "TP53.P72R", 
    "BRAF.L79Q"
    # "ASXL1.L815P","TP53.N288T","TP53.K291T", "RUNX1.1193L", "TP53.P72R",
    # "TET2.585L","ASXL1.1292S", "TET2.S585L", "EZH2.K500T", "ASXL1.G1292S", 
    # "DNMT3A.M801L", "ASXL1.Y974H", "ASXL1.Q733R", "EZH2.K500N", "RUNX1.R250H" ,# benign
    # "ASXL1.L1043M",  "ASXL1.G1292S",  "TP53.150S" , "DNMT3A.M801L",  "TET2.S585L",  "BRAF.S36A" , 
    # "TET2.Y819*", 
  )
  sce_list <- vector("list", length(sample_files)) # Pre-allocate memory
  
  for(sample_file in sample_files){  #sample_file<- "BRAF/A0634braf.dna+protein.h5"
    sample_name <- sub("\\..*", "", basename(sample_file)) 
    cat(blue("------ Creating sce object for", sample_name, "------\n"))
    if (length(base_variant_output_list)<1){
      cat(red("Computing Base Variants..."))
      suppressWarnings(variant_output<-variant_ID(file=sample_file,panel="Myeloid",GT_cutoff=90,
                                                  VAF_cutoff=2))
    } else {
      cat("using given base variants thx.\n")
      variant_output<- base_variant_output_list[[sample_file]]
    }

    braf_variants_to_select <- c("BRAF.V600E", "BRAF.G469", "BRAF.D594")
    braf_pattern <- paste(braf_variants_to_select, collapse = "|")
    braf_variants_of_interest <- as.data.frame(variant_output)%>%
      dplyr::filter(VAF >= 0.05)%>%
      dplyr::filter(!AA_change %in% blacklist) %>%
      dplyr::filter(!CONSEQUENCE == "synonymous" ) %>%
      dplyr::filter(str_detect(AA_change, braf_pattern)) %>%
      dplyr::group_by(AA_change) %>%
      dplyr::slice_max(VAF, n = 1) %>%  # Keep only the row with the highest VAF per AA_change
      dplyr::ungroup()  %>%
      arrange(desc(VAF)) %>%
      dplyr::slice(1:3)
    cat(blue("BRAF vars:", paste0(braf_variants_of_interest$AA_change, collapse = ","), "\n" ))
    # for just braf cohort look for these ras vars. 
    #ras_vars<-  c("KRAS.T58", "KRAS.Q61", "NRAS.G12", "NRAS.G13","PTPN11.A72V")

    ras_vars<-  c("KRAS.T58", "KRAS.Q61", "NRAS.G12", "NRAS.G13","PTPN11.A72V",
                  "PTPN11.T507K", "NRAS.Q61R", "PTPN11.G503A")
    non_ras_non_braf<- c("TET2.H1778R",  "TET2.P1723S",  "TP53.R342*" , "TP53.Y205*",  "ASXL1.F1238C", # 7456
                         "TET2.Y1724S", # 1912
                         "RUNX1.G165R", # 0292
                         "NPM1.287", # 2459
                         "EZH2.D185H",#5330
                         "FLT3.N676K", "TET2.Q886*", #0634 
                         "ASXL1.L1393M",
                         "TET2.961", "FLT3.N676S" ,# 0290
                         "TET2.R1269I", "TET2.I1873T", # 3307
                         "TET2.R1516*", #1962
                         "TET2.C1221S", #17020
                         "TET2.H1219Q", "FLT3.594", "FLT3.F594*I*SQMV",
                         "TET2.P1723S", "TET2.Y1608*" #A1107
                         )
    
    all_non_braf_variants<- c(ras_vars, non_ras_non_braf)
    non_braf_vars_pattern <- paste(all_non_braf_variants, collapse = "|")
    
    
    non_braf_variants_of_interest <- variant_output %>%
      dplyr::filter(VAF >= 1)%>%
      dplyr::filter(!CONSEQUENCE == "synonymous" ) %>%
      dplyr::filter(str_detect(AA_change, non_braf_vars_pattern)) %>%
      dplyr::group_by(AA_change) %>%
      dplyr::slice_max(VAF, n = 1) %>%  # Keep only the row with the highest VAF per AA_change
      dplyr::ungroup()  %>%
      arrange(desc(VAF))# %>%# dplyr::slice(1:4)
    
    cat(blue("Other vars:", paste0(non_braf_variants_of_interest$AA_change, collapse = ","), "\n" ))
    
    variants_of_interest<- rbind(braf_variants_of_interest, non_braf_variants_of_interest) %>% arrange(desc(VAF))
    save_variant_data_pdf_used_in_protein_processing(variants_of_interest, sample_name, path_to_save_variant_data)
    
    sce<-tapestri_h5_to_sce(file=sample_file,variant_set = variants_of_interest,GT_cutoff=90, 
                            VAF_cutoff=0.01, DP_cutoff=10, AF_cutoff=20)
    cat("Rows:", paste0(rownames(sce), collapse = ", "), "\n")

    sce<- rename_cells_with_group_and_sample_identifier(sample_name, group_identifier, sce)
    gc()
    invisible(capture.output(sce <- enumerate_clones_2(sce)))
    
    cat("New NGT colnames: \n")
    print(head(colnames(sce@assays@data@listData[["NGT"]])))
    sce_list[[length(sce_list) + 1]] <- sce
    rm(variant_output, variants_of_interest, sce)
    gc() 
  }
  cat(blue("✅ We now have a list of sce objects that represent a grouping.\n"))
  head(sce_list)
  sce_list <- Filter(Negate(is.null), sce_list)
  return(sce_list)
}#----------------------------------------------------------------------------------------------------
# sample_files <-c("BRAF/A0634braf.dna+protein.h5", "./data/A0290.dna+protein.h5")
# for (obj in sce_list){
#   cat(obj@metadata[["sample_name"]], "\n")
#   cat(unique(rownames(obj)), "\n")
#   cat(length(unique(rownames(obj))), "\n")
# }


save_variant_data_pdf_used_in_protein_processing <- function(variants_of_interest, sample_name, path_to_save_variant_data) {
  variants_of_interest_to_save<- variants_of_interest
  variants_of_interest_to_save$Sample<- sample_name
  
  variants_data_to_save <- variants_of_interest_to_save %>%
    dplyr::select(Sample, id, AA_change, CDS_change, CONSEQUENCE, WT, Het, Hom, Missing, VAF, genotyping_rate)
  
  data_save_name <- file.path(path_to_save_variant_data, paste0(sample_name, "_used_variants.pdf"))
  pdf(data_save_name, width = 15, height = 7) 
  grid.table(variants_data_to_save)
  dev.off()
  
  cat(crayon::green("Saved variant data to", data_save_name, "\n"))
}





# get rid of cells that are marked as wt (don't have any of our special variants) but  have other pathogenic variants. 
remove_cells_not_true_wt <- function(sce, variant_output){
  alleged_wt_cells<- get_alleged_wt_cells(sce)
  pathogenic_cells<- get_cells_with_pathogenic_mutations(sce@metadata[["file"]], variant_output)

  # which of my alleged wt cells actually have pathogenic mutations
  cells_to_get_rid_of<- intersect(alleged_wt_cells, pathogenic_cells)
  cat(silver("getting rid of", length(cells_to_get_rid_of), "cells \n"))
  
  sce <- sce[, !colnames(sce) %in% cells_to_get_rid_of]
  return(sce)
}



get_alleged_wt_cells<- function(sce){
  ngt <- assay(sce, "NGT")
  ngt[ngt == 2] <- 1 # replace 2s with 1s b/c there or not. 
  ngt[ngt == 3] <- 0 # ignore 3s. they don't count as mutations
  alleged_wt_cells <- colnames(ngt)[colSums(ngt) == 0]
  cat(cyan(length(alleged_wt_cells),"cells out of",length(colnames(sce)) ,"don't have", paste((variants_of_interest$AA_change), collapse = ", ")))
  return(alleged_wt_cells)
}

get_cells_with_pathogenic_mutations<- function(sample_file, variant_output){
  blacklist <- c( # exclude these variants. not pathogenic.
    "TET2.I1762V", "TET2.961*", "TET2.Y1579*", "TET2.A1283T", "TET2.L1721W",
    "DNMT3A.F772C", "NRAS.L56P", "NRAS.T58A", "NRAS.L56Q", "NRAS.D57N",
    "DNMT3A.I310S", "NRAS.T58I", "DNMT3A.I292S", "DNMT3A.L888Q", "DNMT3A.L888P",
    "PTPN11.L525R", "DNMT3A.N489T", "DNMT3A.K429T", "DNMT3A.N757T", "NRAS.T58P",
    "NRAS.D57Y", "NRAS.Q61P", "TET2.Q618H", "TET2.L1819F", "DNMT3A.Y481S",
    "FLT3.N847T", "TET2.A584T", "TET2.A584P", "DNMT3A.F290L", "TP53.P72R",
    "ASXL1.L815P",
    "TP53.N288T",
    "TP53.K291T",
    "RUNX1.1193L"
  )
  potentially_pathogenic_genes <- c("IDH2", "NPM1", "TET2", "FLT3","IDH1", "DNMT3A", "SRSF2", "BCOR", "RUNX1", "U2AF1", "EZH2", "ASXL1", "FLT3", "TP53", "KMT2A")
  
  # identify variants that would make a cell pathogenic. 
  other_pathogenic_variants <- variant_output %>%
    dplyr::filter(!CONSEQUENCE == "synonymous" ) %>%
    dplyr::filter(SYMBOL %in% potentially_pathogenic_genes) %>%
    dplyr::filter(!AA_change %in% blacklist) %>%
    arrange(desc(VAF)) %>%
    dplyr::slice(1:4)
  other_pathogenic_variants$AA_change
  
  # get cell names that have these mutations. 
  sce_for_checking_wt<-tapestri_h5_to_sce(file=sample_file,variant_set = other_pathogenic_variants,GT_cutoff=90, VAF_cutoff=0.01,DP_cutoff=10,GQ_cutoff=20, AF_cutoff=20)
  ngt <- assay(sce_for_checking_wt, "NGT")
  ngt[ngt == 2] <- 1 # replace 2s with 1s b/c there or not. 
  ngt[ngt == 3] <- 0 # ignore 3s.
  cells_with_other_mutations <- ngt[, colSums(ngt == 1) > 0] # select columns that have a 1 for any of their values
  pathogenic_cell_names<- colnames(cells_with_other_mutations)

  cat(cyan(length(pathogenic_cell_names),"cells out of",length(colnames(sce_for_checking_wt)) ,"have other pathogenic mutations: ", paste((other_pathogenic_variants$AA_change), collapse = ", ")))
  return(pathogenic_cell_names)
}


