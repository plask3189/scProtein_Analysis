
#RESTART R
#rstudioapi::restartSession()
# FLT3-ITD, NRAS, KRAS, PTPN11, KIT
# non signaling:  DNMT3A, TET2, IDH1/2, NPM1, SRSF2, ASXL1

library(crayon)
library(scDNA)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(gridExtra)
source("src/braf_protein_helper_functions.R")
source("src/variant_ID.R")
source("src/enumerate_clones_2.R")
source("src/protein_kp/cluster.R")
library(BiocGenerics)#showMethods("match")

#⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐
# ================================ Main Run  ================================ 
#⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐
comparing_cohorts<- FALSE

if(comparing_cohorts){ # if two different groups (like BRAF samples vs other AML samples)
  path_to_save_variant_data<- "two_cohorts"
  run_name<- "two_cohorts"
  sce_save_name <- paste0(run_name, "sce_subset")
  base_variant_output_list_save_name <- paste0(run_name, "_base_variant_output_list")
  cat(magenta("Comparing two cohorts: BRAF samples vs other AML samples \n"))
  sample_files_group1 <- Sys.glob("BRAF/*.h5")
  
  #sample_files_group2 = c()# TET2, NPM1, NRAS
  base_variant_output_list_g1 <- list()
  for (sample_file in sample_files_group1){
    base_variant_output_list_g1[[sample_file]]<- variant_ID(file=sample_file,panel="Myeloid", GT_cutoff=90,VAF_cutoff=1)
  }; saveRDS(base_variant_output_list_g1, paste0( file.path(path_to_save_variant_data, base_variant_output_list_save_name), "_g1.rds"))
  
  base_variant_output_list_g2<- list()
  for (sample_file in sample_files_group2){
    base_variant_output_list_g2[[sample_file]]<- variant_ID(file=sample_file,panel="MSK_RL", GT_cutoff=90,VAF_cutoff=1)
  }; saveRDS(base_variant_output_list_g2, paste0( file.path(path_to_save_variant_data,base_variant_output_list_save_name), "_g2.rds"))
  
  group1_sce_list<- process_sces(sample_files_group1, group_identifier = 1,path_to_save_variant_data, 
                                 base_variant_output_list_g1) # make the sce objects and add the group identifier to the metadata.

  group2_sce_list<- process_sces(sample_files_group2, group_identifier = 2, path_to_save_variant_data,
                                 base_variant_output_list_g2)
  
  sce_list<- c(group1_sce_list, group2_sce_list)
  sce_list <- Filter(Negate(is.null), sce_list) # take out rows with NAs
  missing_markers<- identify_missing_markers(sce_list)
  sce_list<- handle_missing_marker_data(sce_list, missing_markers)
  cat(green("✅ Merged the two cohorts sce lists into one."))
  saveRDS(sce_list, file.path(path_to_save_variant_data,"both_cohorts_sce_list.rds"))
} else{ # JUST ONE COHORT (LIKE JUST BRAF SAMPLES) ----------------------------------
  path_to_save_variant_data<- "just_braf"
  run_name<- "just_braf"
  save_name <- "just_braf_cohort_sce_subset_n7_v3"
  sce_save_name <- paste0(run_name, "sce_subset")
  base_variant_output_list_save_name <- paste0(run_name, "base_variant_output_list")
  

  sample_files_group1 <-  Sys.glob("BRAF/*.h5")
  for (sample_file in sample_files_group1){
    base_variant_output_list<-variant_ID(file=sample_file,panel="Myeloid",
                                         GT_cutoff=85,VAF_cutoff=2)
  }; saveRDS(base_variant_output_list, paste0(base_variant_output_list_save_name, ".rds"))

  sce_list<- process_sces(sample_files_group1, group_identifier = 1, 
                          path_to_save_variant_data = path_to_save_variant_data, 
                          base_variant_output_list = base_variant_output_list) 
  sce_list <- Filter(Negate(is.null), sce_list) # take out rows with NAs
  missing_markers<- identify_missing_markers(sce_list); missing_markers
  sce_list<- handle_missing_marker_data(sce_list, missing_markers)
  saveRDS(sce_list, paste0("just_braf_sce_list_v2.rds"))
}
#sce_list1<- readRDS("two_cohorts/both_cohorts_sce_list.rds")
all_unique_vars<- list()
for (obj in sce_list){
  cat(obj@metadata[["sample_name"]], "\n")
  variants_used<- unique(rownames(obj))
  cat(variants_used, "\n")
  all_unique_vars<- c(all_unique_vars, variants_used)
  cat(length(unique(rownames(obj))), "\n")
}
all_unique_vars<- unique(unlist(all_unique_vars)); all_unique_vars
source("src/load_libraries.R")
load_libraries_2()

group_droplet_metadata<- group_droplete_data(sce_list)

group_background_droplets<-group_droplet_metadata%>% # group background cells/droplets
  dplyr::filter(Droplet_type=="Empty")%>% 
  dplyr::filter(dna_size<1.5&dna_size>0.15)%>%
  pull(Cell)

# does DSB and CLR normalization and stores in the respective assay
time_taken <- system.time({
  sce_subset<-dsb_normalize_protein_data(sce_list,
                                         group_metadata=group_droplet_metadata,
                                         group_background=group_background_droplets)
}); cat("⏱️ Time to run grouping and normalization took", round(time_taken["elapsed"] / 60, 2), "minutes\n")

print(unique(rownames(sce_subset@assays@data@listData[["NGT"]])))
saveRDS(sce_subset, file.path(path_to_save_variant_data, paste0(run_name, "_sce.rds")))

seurat_obj<- cluster_protein_data(sce_subset,protein_normalization_method = "DSB_norm")

saveRDS(seurat_obj, file.path(path_to_save_variant_data, paste0(run_name, "_SEURAT.rds")))








