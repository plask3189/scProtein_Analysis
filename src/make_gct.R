library(data.table)
library(writexl)
source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/R/protein_kp/visualization_helper_functions.R", echo=TRUE)
library(stringr)
library(rlang)

just_BRAF_gct<- TRUE

if (just_BRAF_gct == TRUE){
  seurat <- readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/just_braf_cohort_sce_subset_n7_v2_SEURAT.rds")
  sce <- readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/just_braf_cohort_sce_subset_n7_v2.rds")
} else {
  seurat<- readRDS("two_cohorts/two_cohorts_SEURAT.rds")
  sce <- readRDS("two_cohorts/two_cohorts_sce_subset_v4.rds")
}

basic_seurat <- initialize_seurat_object_for_gct(seurat = seurat, sce =sce, muts = c("BRAF"))
table(Idents(basic_seurat))# view the current idents

if  (just_BRAF_gct == TRUE){ 
  #find markers
  basic_seurat@meta.data$mutation_status <- ifelse(
    grepl("D594", basic_seurat@meta.data$mutation_status), "BRAF.D594",
    ifelse(
      grepl("G469", basic_seurat@meta.data$mutation_status), "BRAF.G469",
      basic_seurat@meta.data$mutation_status
    )
  )
  Idents(basic_seurat) <-basic_seurat@meta.data[["mutation_status"]]; table(Idents(basic_seurat))# # set new indents
  #markers <- FindMarkers(basic_seurat, ident.1 = "BRAF.G469", ident.2 = "BRAF.D594") %>%arrange(p_val_adj)
  markers <- FindAllMarkers(basic_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%arrange(p_val_adj)
  markers <- markers %>% filter(cluster != "WT") %>% group_by(cluster) %>% top_n(10, avg_log2FC) 
  top_cds<- markers$gene[1:30]
  top_cds <- top_cds[!grepl("Ig", top_cds)]
  top_cds<- c(top_cds, "CD14", "CD123", "CD11b", "CD34")
  DoHeatmap(basic_seurat, features = top_cds) + NoLegend()
} else {
  Idents(basic_seurat) <-basic_seurat@meta.data[["Group"]] # set new indents
  markers <- FindMarkers(basic_seurat, ident.1 = "BRAF Cohort", ident.2 = "RAS Cohort") %>%arrange(p_val_adj)
  top_cds<- markers$gene[1:30]
  top_cds <- top_cds[!grepl("Ig", top_cds)]
  top_cds<- c(top_cds, "CD14", "CD123", "CD11b", "CD34")
}

basic_metadata<- as.data.frame(basic_seurat@meta.data) %>% dplyr::select(-c("mutation_status")) # will add different statuses later

mutations_to_label<- c("TET2", "TP53", "ASXL1", "BRAF", "RUNX1", "NPM1", "NRAS", "FLT3", "PTPN11", "KRAS")
new_mutation_status_columns<- make_many_mutation_status_columns(mutations_to_label, seurat, sce)

# new_mutation_status_columns$BRAF_mutation_status <- ifelse(
#   grepl("D594", new_mutation_status_columns$BRAF_mutation_status), "BRAF.D594",
#   ifelse(
#     grepl("G469", new_mutation_status_columns$BRAF_mutation_status), "BRAF.G469",
#     new_mutation_status_columns$BRAF_mutation_status
#   )
# )
two_cohorts_metadata<- left_join(basic_metadata, new_mutation_status_columns, by= "Cell")

# # remove braf wt cells.
# two_cohorts_metadata <- two_cohorts_metadata %>% filter(BRAF_mutation_status != "WT")
braf_cohort_data <- two_cohorts_metadata %>% filter(Group== "BRAF Cohort")
ras_cohort_data <- two_cohorts_metadata %>% filter(Group == "RAS Cohort")

data<- braf_cohort_data # 🎆🎆🎆🎆🎆 data to use 


make_gct(data)

make_gct<- function(data){
  data$sample_name <- str_sub(data$Cell, end = -20)
  data <- data[, names(data) %in% top_cds | names(data) %in% c("Cell", "Group", "Clone", "sample_name") |
                 grepl("mutation_status", names(data)) ]
  
  sample_metadata_header <- as.data.frame(t(data[, names(data) %in% c("Cell", "Group", "Clone", "sample_name") | grepl("mutation_status", names(data))]))
  print(sample_metadata_header[,1:2])
  marker_data <- as.data.frame(t(data[, !colnames(data) %in% rownames(sample_metadata_header)])); print(marker_data[1:3,1:2])

  selected_cells_all<- select_cells(sample_metadata_header)
  marker_data <- marker_data[, selected_cells_all, drop = FALSE]
  sample_metadata_header <- sample_metadata_header[, selected_cells_all, drop = FALSE]
  version<- "#1.3"
  gene_metadata_header<- as.data.frame(rownames(marker_data)); colnames(gene_metadata_header) <- "marker_name"; head(gene_metadata_header)
  num_markers<- nrow(marker_data); num_markers
  num_samples<- ncol(sample_metadata_header);num_samples
  metadata_rows<- nrow(sample_metadata_header)
  dimension_vals<- c(num_markers, num_samples, ncol(gene_metadata_header),  metadata_rows); dimension_vals
  #------------examine data-------------
  head(gene_metadata_header); head(marker_data[1:5, 1:3]); head(sample_metadata_header[,1:3]); print(dimension_vals)
  #-------------section construction----------------
  rbind(version, dimension_vals)
  version_row <- c(version, rep("", length(dimension_vals) - 1))
  file_header <- rbind(version_row, dimension_vals)
  #----------- make marker and id columns ------------
  sample_metadata_header_formatted<- sample_metadata_header
  new_col_name <- colnames(gene_metadata_header); new_col_name
  sample_metadata_header_formatted[[new_col_name]] <- NA
  sample_metadata_header_formatted <- sample_metadata_header_formatted[, c(new_col_name, setdiff(names(sample_metadata_header_formatted), new_col_name))]
  sample_metadata_header_formatted[["id"]] <- rownames(sample_metadata_header)
  sample_metadata_header_formatted <- sample_metadata_header_formatted[, c("id", setdiff(names(sample_metadata_header_formatted), "id"))]
  rownames(sample_metadata_header_formatted) <- NULL
  print(sample_metadata_header_formatted[,1:4])
  
  marker_data_formatted<- as.data.frame(marker_data)
  new_col_name <- colnames(gene_metadata_header); new_col_name
  marker_data_formatted[[new_col_name]] <- rownames(marker_data_formatted)
  marker_data_formatted <- marker_data_formatted[, c(new_col_name, setdiff(names(marker_data_formatted), new_col_name))]
  marker_data_formatted[["id"]] <- rownames(marker_data_formatted)
  marker_data_formatted <- marker_data_formatted[, c("id", setdiff(names(marker_data_formatted), "id"))]
  rownames(marker_data_formatted) <- NULL
  print(marker_data_formatted[1:4,1:4])
  
  #---------------------combine sections --------- 
  main_section<- rbind(sample_metadata_header_formatted, marker_data_formatted)
  file_header_df <- as.data.frame(file_header, stringsAsFactors = FALSE)
  main_section_df <- as.data.frame(main_section, stringsAsFactors = FALSE)
  
  write.table(file_header_df, "for_morpheus.gct",
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(main_section_df, "for_morpheus.gct",
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
}

# 🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄
make_many_mutation_status_columns<- function(mutations_to_label, seurat, sce){
  set_of_new_mutation_status_columns<- data.frame()
  for (mut in mutations_to_label){ #make mutation_status_column for each. 
    seurat_with_mut_label <- initialize_seurat_object_for_gct(seurat = seurat, sce =sce, muts = c(mut))
    meta_seurat_with_mut_label<- seurat_with_mut_label@meta.data
    new_mut_column_name<- paste0(mut, "_mutation_status")
    meta_seurat_with_mut_label <- meta_seurat_with_mut_label  %>%
      dplyr::rename_with(~new_mut_column_name, "mutation_status")
    new_mut_column_data <- as.data.frame(meta_seurat_with_mut_label) %>%
      dplyr::select("Cell", matches("mutation_status"))
    cat("Num cells", nrow(new_mut_column_data), "\n")
    if (length(set_of_new_mutation_status_columns)<=0) {# if is frst column, initialize df as new_mut_column_data
      set_of_new_mutation_status_columns <- new_mut_column_data
    } else {
      set_of_new_mutation_status_columns <- left_join(set_of_new_mutation_status_columns, new_mut_column_data, by = "Cell")
    }
  }
  return(set_of_new_mutation_status_columns)
}# 🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄🦄


select_cells<- function(sample_metadata_header){
  cell_sample_ids <- as.character(sample_metadata_header["sample_name", ])
  unique_samples <- unique(cell_sample_ids)
  selected_cells_all <- c()
  set.seed(1333)  
  for (sample_id in unique_samples) {
    cell_names <- colnames(marker_data)[cell_sample_ids == sample_id]
    if (length(cell_names) < 1000) {
      warning(paste("Sample", sample_id, "has fewer than 500 cells, keeping all"))
      selected_cells <- cell_names
    } else {
      # Randomly sample 500 cell names
      selected_cells <- sample(cell_names, 1000)
    }
    
    selected_cells_all <- c(selected_cells_all, selected_cells)
  }
  num_cells<- length(selected_cells_all)
  cat("Using",num_cells, "cells", "\n" )
  return(selected_cells_all)
}



initialize_seurat_object_for_gct<- function(seurat, sce , muts = c("BRAF")){ 
  cat(blue("Using given seurat object and sce subset \n"))
  if(! "RAS Cohort" %in% unique(seurat@meta.data[["Group"]]) ){
    seurat@meta.data[["Group"]][seurat@meta.data[["Group"]] == 1] <- "BRAF Cohort"
    seurat@meta.data[["Group"]][seurat@meta.data[["Group"]] == 2] <- "RAS Cohort"
  }
  for(cluster_of_differentiation in rownames(seurat@assays$Protein@data)){ # for all the markers: 
    seurat@meta.data[[cluster_of_differentiation]] <- seurat@assays$Protein@data[cluster_of_differentiation, ] # makes the marker a column in the seurat object
    seurat@meta.data[[cluster_of_differentiation]][seurat@meta.data[[cluster_of_differentiation]] < 0] <- 0 # https://www.nature.com/articles/s41467-022-29356-8
  }
  #variants_in_sce_clean <- unique(sub("\\..*", "", rownames(sce)))
  seurat<- all_mutations_plus_no_overwrite(sce, seurat, original_desired_variants = muts,
                                           include_individual_braf_muts = TRUE)
  Idents(seurat) <-seurat@meta.data[["Group"]]
  return(seurat)
}

two_cohorts_metadata %>%
  dplyr::count(seurat_clusters, BRAF_mutation_status) %>%
  ggplot(aes(x = factor(seurat_clusters), y = n, fill = BRAF_mutation_status)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Seurat Cluster", y = "Proportion", fill = "BRAF Status") +
  theme_minimal()
RidgePlot(basic_seurat, features = c("CD34", "CD33", "CD117"), group.by = "BRAF_mutation_status")
