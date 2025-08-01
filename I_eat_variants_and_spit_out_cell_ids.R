# JUST BRAF
library(crayon)
library(scDNA)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
source("R/protein_kp/braf_protein_helper_functions.R")
source("R/variant_ID.R")
source("R/enumerate_clones_2.R")
source("R/protein_kp/cluster.R")
library(BiocGenerics)
showMethods("match")



want <- c(
  "TET2.H1778R", "TET2.P1723S", "TP53.R342*", "TP53.Y205*", "ASXL1.F1238C",
  "BRAF.G469A", "BRAF.V600E", "TET2.Y1724S", "RUNX1.G165R", "KRAS.T58I",
  "NPM1.287", "ASXL1.L1393M", "EZH2.D185H", "NRAS.G12R", "BRAF.D594N",
  "TET2.Q886*", "FLT3.N676K", "PTPN11.A72V", "BRAF.D594G", "PTPN11.T507K",
  "NRAS.Q61R", "TET2.961", "PTPN11.G503A", "FLT3.N676S", "TET2.Y1608*",
  "NRAS.G13V", "TET2.H1219Q", "NRAS.G12V", "FLT3.594", "TET2.C1221S",
  "KRAS.Q61H", "TET2.R1516*", "TET2.R1269I", "TET2.I1873T", "NRAS.G13D"
)

g1_var_list<- readRDS("two_cohorts/two_cohorts_base_variant_output_list_g1.rds")
g2_var_list<-readRDS("two_cohorts/two_cohorts_base_variant_output_list_g2.rds")

big_var_list<- c(g1_var_list,g2_var_list)

list_of_variant_outputs<- init(generate_base_variant_list= FALSE, big_var_list) # loads base variants. 
file_path = "two_cohorts" # to save var data
non_braf_non_ras_variants <-c("TP53", "NPM1", "TET2", "FLT3","IDH1","IDH2", "DNMT3A", "RUNX1", "ASXL1", "FLT3", "EZH2", "KMT2A")

sample_files_group1 <-c("BRAF/M7456braf.dna+protein.h5",
                        "BRAF/M1912braf.dna+protein.h5", 
                        "BRAF/M0292braf.dna+protein.h5",
                        "BRAF/2459_Braf.dna+protein.h5", #"BRAF/M0626braf.dna+protein.h5", "BRAF/6232_braf.dna+protein.h5",
                        "BRAF/4629_braf.dna+protein.h5",
                        "BRAF/A5330braf.dna+protein.h5",
                        "BRAF/A0634braf.dna+protein.h5")

sample_files_group2 = c("./data/A0290.dna+protein.h5", # TET2, NPM1, NRAS
                        "./data/A1107.dna+protein.h5", # TET2, NPM1, NRAS
                        "./data/A0259.dna+protein.h5",# TET2, NRAS, FLT3
                        "./data/Sample17020.dna+protein.h5", # TET2, NPM1, NRAS
                        "./data/LAM_1962.dna+protein.h5", # TET2, NPM1, NRAS
                        "./data/LAM_3307.dna+protein.h5")# TET2, NPM1, NRAS
both_cohorts_files<- c(sample_files_group1, sample_files_group2)

sce_subset <- readRDS("two_cohorts/two_cohorts_sce.rds")

setdiff(rownames(sce_subset))

new_variants<- c("KRAS.T58I",  "PTPN11.G503A", "FLT3.N676S")
sce_list_extra_variants<- which_cells_have_this_variant(new_variants, 
                                                        both_cohorts_files, 
                                                        list_of_variant_outputs); rownames(sce_list_extra_variants[[1]]) 

sce_subset_added_variants <- add_new_variant_to_sce(sce_list_extra_variants, new_variants, sce_subset)

saveRDS(sce_subset_added_variants, file.path("two_cohorts", "two_cohorts_sce_subset_v4.rds"))




add_new_variant_to_sce<- function(sce_list_extra_variants, new_variants, sce_subset){
  row_list <- list() 
  for (i in seq_along(sce_list_extra_variants)) {
    row_i <- sce_list_extra_variants[[i]]@assays@data@listData[["NGT"]][1, , drop = FALSE]
    row_list[[i]] <- row_i
  }
  merged_row <- do.call(cbind, row_list)# new ngt. 
  rownames(merged_row)<- new_variants; merged_row[1:5]
  
  little_new_sce<-SingleCellExperiment::SingleCellExperiment(list(
    NGT=as.data.frame(merged_row)))
  
  new_sce_list <- list(sce_subset,little_new_sce)
  convert_assay_to_matrix <- function(sce, assay_name) { # to convert specified assay to matrix class bc scMerge needs it that way.
    assay_data <- assay(sce, assay_name)
    assay(sce, assay_name) <- as.matrix(assay_data)
    return(sce)
  }
  sce_list <- lapply(new_sce_list, convert_assay_to_matrix, assay_name = "NGT")
  #group_sce <- sce_cbind(sce_list = sce_list, method = "union", exprs = "NGT") # Combine the main experiments 
  
  common_cells <- Reduce(intersect, lapply(sce_list, function(x) colnames(assay(x, "NGT"))))
  
  #   align each object's assay
  merged_counts <- do.call(rbind, lapply(sce_list, function(x) {
    assay(x, "NGT")[, common_cells, drop = FALSE]
  }))
  
  #  Create merged SingleCellExperiment object
  merged_sce <- SingleCellExperiment(
    assays = list(NGT = merged_counts)
  )
  
  altExp(merged_sce, "Protein") <- altExp(sce_subset, "Protein")  

cat(blue("Old dimensions:", paste0(dim(sce_subset), collapse = "x"), "\n"))
cat(blue("New dimensions:", paste0(dim(merged_sce), collapse = "x"), "\n"))
return(merged_sce)
}







# relabeled_umap<- relab_cells_for_variant(original_labeled_seurat_umap, cells_with_new_variant, variant_name= new_variants,
#                                          new_colors = c("yellow"),
#                                          dont_overwrite= c("BRAF")); relabeled_umap
# 
# new_cell_labels_s<- relab_seurat(s, cells_positive_for_variant, variant_name= new_variants, dont_overwrite= c("BRAF"))
# as.data.frame(new_cell_labels_s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
# #umap_plot <- DimPlot(s, reduction = "umap", group.by = "mutation_status", alpha=0.4, pt.size =0.05)
# color_palette<- set_color_pallette()
# original_violin_plot <- suppressWarnings(
#   VlnPlot(new_cell_labels_s, features = "CD34",
#           pt.size = 0.01, alpha = 0.4, 
#           group.by = "mutation_status", cols = color_palette)) + theme(axis.title.x = element_blank());original_violin_plot


init<- function(generate_base_variant_list= FALSE, big_variant_output_list = NA){
  if(generate_base_variant_list){
    list_of_variant_outputs<- make_list_of_variant_outputs(sample_files_group1)
    saveRDS(list_of_variant_outputs, paste0(path, "list_of_var_outputs.rds")) #base_variants<- readRDS("list_of_var_outputs.rds")
  } else{
    list_of_variant_outputs<-big_variant_output_list # readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copylist_of_var_outputs.rds")
  }
  return(list_of_variant_outputs)
}



which_cells_have_this_variant<- function(variant_list, sample_files, list_of_variant_outputs){
  blacklist <- c(
    "TET2.I1762V", "TET2.961*", "TET2.Y1579*", "TET2.A1283T", "TET2.L1721W",
    "DNMT3A.F772C", "NRAS.L56P", "NRAS.T58A", "NRAS.L56Q", "NRAS.D57N",
    "DNMT3A.I310S", "NRAS.T58I", "DNMT3A.I292S", "DNMT3A.L888Q", "DNMT3A.L888P",
    "PTPN11.L525R", "DNMT3A.N489T", "DNMT3A.K429T", "DNMT3A.N757T", "NRAS.T58P",
    "NRAS.D57Y", "NRAS.Q61P", "TET2.Q618H", "TET2.L1819F", "DNMT3A.Y481S",
    "FLT3.N847T", "TET2.A584T", "TET2.A584P", "DNMT3A.F290L", "RUNX1.R250H","EZH2.D185H", "ASXL1.V7511")
  
  sce_list <- vector("list", length(sample_files)) # Pre-allocate memory
  for(sample_file in sample_files){  #sample_file<- "BRAF/A0634braf.dna+protein.h5"
    sample_name <- sub("\\..*", "", basename(sample_file)) 
    cat(blue("------ Creating sce object for", sample_name, "------\n"))
    index<- which(sample_files == sample_file)
    variant_output<- list_of_variant_outputs[[index]]
    
    pattern <- paste(variant_list, collapse = "|")
    variants_of_interest <- variant_output %>%
      dplyr::filter(!AA_change %in% blacklist) %>%
      dplyr::filter(!CONSEQUENCE == "synonymous" ) %>%
      dplyr::filter(str_detect(AA_change, pattern)) %>%
      dplyr::group_by(AA_change) %>%
      dplyr::slice_max(VAF, n = 1) %>%  # Keep only the row with the highest VAF per AA_change
      dplyr::ungroup()  %>% 
      arrange(desc(VAF))
    
    #save_variant_data_pdf(variants_of_interest, sample_file, file_path)
    
    sce<-tapestri_h5_to_sce(file=sample_file,variant_set = variants_of_interest,GT_cutoff=90, VAF_cutoff=0.01,DP_cutoff=10,GQ_cutoff=20, AF_cutoff=20)

    sce<- rename_cells_with_group_and_sample_identifier(sample_name, group_identifier=0, sce)
    
    apply(variants_of_interest[, c("id", "AA_change", "VAF", "genotyping_rate")], 1, function(row) {
      cat(magenta("--- Found ---\n"))
      cat(green(paste0("ID: ", row["id"],
                        " | AA_change: ", row["AA_change"],
                        " | VAF: ", row["VAF"],
                        " | Genotyping rate: ", row["genotyping_rate"], "\n")))
    })
    cat("New NGT colnames: \n")
    print(head(colnames(sce@assays@data@listData[["NGT"]])))
    cat("New NGT rownames: \n")
    print(head(rownames(sce@assays@data@listData[["NGT"]])))
    sce_list[[length(sce_list) + 1]] <- sce
    rm(variant_output, variants_of_interest, sce)
    gc() 
  }
  sce_list <- Filter(Negate(is.null), sce_list) 
  saveRDS(sce_list, "sce_list_2.rds")
  return(sce_list)
}#----------------------------------------------------------------------------------------------------


get_list_of_cells<- function(sce_list, which_variant){
  cat(blue("---------Tallying cells for", which_variant, "---------\n"))
  which_samples_have_this_variant <- list()
  all_cells_with_this_variant<- list()
  for (sce in sce_list){
    cat(magenta("      -----Processing", sce@metadata[["sample_name"]], which_variant, "-----\n"))
    ngt <- assay(sce, "NGT")
    cat("      Looking for", which_variant ,"in available variants:", rownames(ngt), "\n")
    ngt[ngt == 2] <- 1 # replace 2s with 1s b/c there or not. 
    ngt[ngt == 3] <- 0 # ignore 3s. they don't count as mutations
    if (any(grepl(which_variant, rownames(ngt)))) {
      which_samples_have_this_variant<- c(which_samples_have_this_variant, sce@metadata[["sample_name"]])
      matching_rows <- rownames(ngt)[grepl(which_variant, rownames(ngt))];matching_rows
      variant_ngt <- ngt[matching_rows, , drop = FALSE]
      cells_with_variant <- unique(colnames(variant_ngt)[colSums(variant_ngt == 1, na.rm = TRUE) > 0])
      cat(blue("      Found",length(cells_with_variant) ,"with", paste0(rownames(variant_ngt), collapse = ", "),"cells in", sce@metadata[["sample_name"]], "with", which_variant, "\n"))
      #cat(cells_with_variant[1:2], "\n")
      all_cells_with_this_variant<-c(all_cells_with_this_variant, cells_with_variant)
    } else {
      cat(red("      ",sce@metadata[["sample_name"]], "doesn't have", which_variant, "\n"))
    }
    
  }
  cat(green("Finished tallying cells. There are", length(all_cells_with_this_variant), "cells with", which_variant, "in", paste0(which_samples_have_this_variant, collapse = ", ") ,"\n"))
  return(all_cells_with_this_variant)
}


save_variant_data_pdf <- function(variants_data_to_save, sample_file, path) {
  name<- paste0(sub("\\..*", "", basename(sample_file)))
  print(colnames(variants_data_to_save))
  variants_data_to_save <- variants_data_to_save %>%
    dplyr::select(
      id, AA_change,
      dplyr::matches("var_found_in"), 
      dplyr::matches("sample_count"), 
      CDS_change, CONSEQUENCE, WT, Het, Hom, Missing, VAF, genotyping_rate
    )
  data_save_name <- file.path(path, paste0(sub("\\..*", "", basename(sample_file)), "_variants.pdf"))
  save_to_pdf(data_save_name, variants_data_to_save)
}


save_to_pdf<- function(path, data_to_save){
  pdf(path, width = 17, height = 8)  
  grid.table(data_to_save)
  dev.off()
  cat(crayon::green("Saved variant data to", path, "\n"))
}


make_list_of_variant_outputs<- function(sample_files){
  list_of_variant_outputs<- list()
  for(sample_file1 in sample_files){
    cat(blue("--------Generating Variant Output for", sub("\\..*", "", basename(sample_file1)), " -----------------\n"))
    variant_output<-variant_ID(file=sample_file1,
                               panel="Myeloid", 
                               GT_cutoff=90,  # mimimum percent of cells where a successful genotyping call was made
                               VAF_cutoff=0.01) 
    list_of_variant_outputs[[length(list_of_variant_outputs)+1]]<- variant_output
  }
  return(list_of_variant_outputs)
}

relab_seurat<- function(s, cells_positive_for_variant, variant_name, dont_overwrite= c("BRAF")){

  metadata1 <- as.data.frame(s@meta.data)
  
  metadata2<- metadata1 %>%
    dplyr::mutate(
      mutation_status = ifelse(
        rownames(metadata1) %in% cells_positive_for_variant &
          !sapply(metadata1$mutation_status, function(x) any(grepl(paste0(dont_overwrite, collapse = "|"), x))),
        variant_name,
        mutation_status
      )
    )
  cat(paste0(unique(metadata2$mutation_status), collapse = ", "))
  s@meta.data$mutation_status <- metadata2$mutation_status
  return(s)
}

