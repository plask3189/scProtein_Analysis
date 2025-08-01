library(ComplexHeatmap)
library(scCustomize)
make_dot_plot<- function(s, cds_of_interest, path){
  # theme_colors<- set_color_pallette()
  # merica<- custom_palette <- colorRampPalette(c(theme_colors[["theme1"]], theme_colors[["gray_theme"]], theme_colors[["theme2"]]))(100)
  all_cds <- rownames(sce_subset@int_colData@listData[["altExps"]]@listData[["Protein"]]@se@assays@data@listData[["DSB_norm"]]);all_cds 
  to_remove <- c("Fc&#949;RI&#945;", "IgG1", "IgG2a", "IgG2b")
  all_cds <- setdiff(all_cds, to_remove)
  
  og_s<-  relabel_seurat_for_two_cohort_run(sce, s)
  s<- og_s
  print(unique( s@meta.data[["mutation_status"]]))
  Idents(s)<- s@meta.data[["mutation_status"]]
 #desired_order <- c("BRAF.D594G", "BRAF.G469A", "BRAF.V600E", "RAS", "Non-RAS Mutant WT", "WT")
  desired_order <- c("WT", "Non-RAS Mutant WT", "RAS", "BRAF.V600E", "BRAF.G469A", "BRAF.D594G")
  s@active.ident <- factor(Idents(s), levels = desired_order)
 
  summary(range(s@assays$Protein@data))  # full global range
  quantile(s@assays$Protein@data, probs = c(0.01, 0.25, 0.5, 0.75, 0.99))
  hist(as.numeric(s@assays$Protein@data[all_cds, ]), breaks = 100)
  abline(v = c(-2.5, 0, 2.5), col = c("blue", "gray", "red"), lty = 2)
  
  # dotplot <- Clustered_DotPlot(s, features = all_cds,# group.by = "mutation_status",
  #                              group.by = "mutation_status", #"Group",
  #                              split.by = 'Group',
  #                              flip = TRUE,
  #                              plot_km_elbow = FALSE, # otherwise outputs the elbow plot too w/ sse info
  #                              cluster_feature = TRUE,
  #                              cluster_ident = FALSE,
  #                              row_names_side = "left",
  #                              seed = 1333,
  # 
  #                              colors_use_exp = c("#4a84c2", "gray", "#9F352D") 
  #   
  #                              );dotplot; png(file.path(path, "dotplot3.png"), width = 10, height = 2.5, units = "in", res = 300);draw(dotplot);dev.off()
  # 
  n_mutations <- length(unique(og_s@meta.data$mutation_status))

  
  dotplot <- DotPlot(
    s,
    features = all_cds,
    group.by = "mutation_status", #"Group",
    split.by = 'Group',
    cols = colorRampPalette(c("#e3e3e3", "#9F352D"))(n_mutations) 
  )
  dotplot
  dot_data<- dotplot[["data"]]
 
custom_dotplot <- ggplot(dot_data, aes(x = features.plot, y = id, size = pct.exp,
                                       color = avg.exp.scaled)) +
  geom_point() +
  scale_size(range = c(1, 6), name = "Percent Exp") +
  RotatedAxis() +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10)) +
  theme(axis.text.y = element_text(size =11)) +
  theme(axis.text.y = element_text(size =11))+
        # legend.text = element_text(size = 15),
        # legend.title = element_text(size = 15),
        # legend.spacing.y = unit(1, "cm")) +
  scale_color_gradientn(
    colors = c("#4a84c2", "#e3e3e3", "#9F352D"),# c("#461248", "#e3e3e3", "#1A4320"),
    # values = scales::rescale(c(-2, 0, 2)),  # Rescale data range
    # limits = c(-2, 2),
    name = "Avg. Exp \n"
  ); custom_dotplot
saveRDS(custom_dotplot, file.path(path,paste0("dotplot.rds")))

ggsave(file.path(path, "dotplot2.png"), plot = custom_dotplot, width = 20, height = 5, dpi = umap_resolution)
  return(custom_dotplot)
}

# does the + situation 
# HANDLE BOTH
#original_desired_variants = c("NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594")
label_seurat_with_mutations<-function(sce_subset, s, original_desired_variants = NA){ #original_desired_variants = c("BRAF", "NRAS", "KRAS", "PTPN11")
  ngt <- assay(sce_subset, "NGT")
  ngt_sparse <- as(ngt, "dgCMatrix")
  ngt_sparse@x[ngt_sparse@x == 2] <- 1# replace all 2s with 1s bc like just BRAF present or not.
  ngt <- as.matrix(ngt_sparse)
  ngt <- as.data.frame(ngt) 
  mutations<- rownames(ngt)
  label_key<- data.frame(og_desired_variants = original_desired_variants)
  label_key$columns_to_choose <- list(NA)
  if ("RAS" %in% label_key$og_desired_variants) {
    #ras_mutations <- c("PTPN11", "NRAS", "KRAS")
    ras_mutations <- rownames(ngt)[grepl("NRAS|PTPN11|KRAS", rownames(ngt))]
    ras_matches <- unlist(lapply(ras_mutations, function(gene) {
      mutations[grepl(gene, mutations)]
    }))
    ras_vars2 <- unique(ras_matches); print(ras_vars2)
    label_key$columns_to_choose[label_key$og_desired_variants == "RAS"] <- list(ras_vars2)
  }
  
  for (var in label_key$og_desired_variants){
    if(var != "RAS"){
      for (i in mutations){ 
        matches <- mutations[grepl(var, mutations)]
      }
      col3<- unique((matches))
      cat("var:", var, "\n", "mut", col3, "\n")
      label_key$columns_to_choose[label_key$og_desired_variants == var]<- list(col3)
    }
  }
  label_key
  
  cat(blue("Identified cells with these variant mutations:", paste(mutations, collapse = ", ")), "\n")
  cat(yellow("But we will identify", paste(label_key$og_desired_variants,collapse = ", "), "\n"))
  
  cells_with_mutations <- ngt[, colSums(ngt == 1) > 0] # select columns that have a 1 for any of their values
  cells_with_muts<-as.data.frame(t(cells_with_mutations))
  
  s@meta.data$mutation_status <- 'WT' # default is WT. 
  
  for(label_index in seq_along(label_key$og_desired_variants)){
    columns<- label_key$columns_to_choose[[label_index]]; columns
    cat("labeling",paste(columns, collapse=","), "as", label_key$og_desired_variants[[label_index]] , "\n")
    for(mutation_name in columns){
      cells_with_mut <- rownames(cells_with_muts[cells_with_muts[[mutation_name]] == 1, ])
      #s@meta.data$mutation_status[rownames(s@meta.data) %in% cells_with_mut] <- label_key$og_desired_variants[[label_index]]
      # combo: 
      # s@meta.data[rownames(s@meta.data) %in% cells_with_mut, "mutation_status"] <- 
      #   paste0( s@meta.data[rownames(s@meta.data) %in% cells_with_mut, "mutation_status"], 
      #     "+", label_key$og_desired_variants[[label_index]])
      target_rows <- rownames(s@meta.data) %in% cells_with_mut #  Identify which rows in s@meta.data correspond to cells with the mutation
      current_status <- s@meta.data[target_rows, "mutation_status"];current_status # Get the current mutation status for those rows
      new_label <- label_key$og_desired_variants[[label_index]]; new_label# Append the new mutation label from label_key
      updated_status <- paste0(current_status, "+", new_label); updated_status #Create the updated status string
      #  Assign the updated status back to those rows in s@meta.data
      s@meta.data[target_rows, "mutation_status"] <- updated_status; s@meta.data[target_rows, "mutation_status"] 
    }
  }
  # REmove the WT+
  s@meta.data$mutation_status <- ifelse(s@meta.data$mutation_status == "WT",
    "WT",gsub("^WT\\+", "", s@meta.data$mutation_status))
  
  as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
  cat(magenta("Overwritting mutation profiles with BRAF in them as just BRAF \n"))
  # make any braf win, relabel as braf if it contains one.
  for (idx in seq_along(s@meta.data$mutation_status)) {
    i <- s@meta.data$mutation_status[idx]
    if (grepl("BRAF", i)) {
      #cat(blue(i), "\n")
      braf_mut <- regmatches(i, regexpr("BRAF\\.[^+]+", i))
     # cat(green("----", braf_mut, "\n"))
      
      if (length(braf_mut) == 0) { # if there's nothing after BRAF so liek just BRAF instead of "BRAF.G469E"
        braf_mut <- "BRAF"
      }
      s@meta.data$mutation_status[idx] <- braf_mut
    } #else {
    #   cat(red(i), "\n")
    # }
  }
  unique(s@meta.data[["mutation_status"]])
  return(s)
}


# has the plus notation but overwrites any with braf as just braf. 
label_seurat_with_mutations2<-function(sce_subset, s, original_desired_variants = NA){ #original_desired_variants = c("BRAF", "NRAS", "KRAS", "PTPN11")
  ngt <- assay(sce_subset, "NGT")
  ngt_sparse <- as(ngt, "dgCMatrix")
  ngt_sparse@x[ngt_sparse@x == 2] <- 1# replace all 2s with 1s bc like just BRAF present or not.
  ngt <- as.matrix(ngt_sparse)
  ngt <- as.data.frame(ngt) 
  mutations<- rownames(ngt)
  label_key<- data.frame(og_desired_variants = original_desired_variants)
  label_key$columns_to_choose <- list(NA)
  if ("RAS" %in% label_key$og_desired_variants) {
    #ras_mutations <- c("PTPN11", "NRAS", "KRAS")
    ras_mutations <- rownames(ngt)[grepl("NRAS|PTPN11|KRAS", rownames(ngt))]
    ras_matches <- unlist(lapply(ras_mutations, function(gene) {
      mutations[grepl(gene, mutations)]
    }))
    ras_vars2 <- unique(ras_matches); print(ras_vars2)
    label_key$columns_to_choose[label_key$og_desired_variants == "RAS"] <- list(ras_vars2)
  }
  for (var in label_key$og_desired_variants){
    if(var != "RAS"){
      for (i in mutations){ 
        matches <- mutations[grepl(var, mutations)]
      }
      col3<- unique((matches))
      cat("var:", var, "\n", "mut", col3, "\n")
      label_key$columns_to_choose[label_key$og_desired_variants == var]<- list(col3)
    }
  }
  label_key
  cat(blue("Identified cells with these variant mutations:", paste(mutations, collapse = ", ")), "\n")
  cat(yellow("But we will identify", paste(label_key$og_desired_variants,collapse = ", "), "\n"))
  cells_with_mutations <- ngt[, colSums(ngt == 1) > 0] # select columns that have a 1 for any of their values
  cells_with_muts<-as.data.frame(t(cells_with_mutations))
  s@meta.data$mutation_status <- 'WT' # default is WT. 
  for(label_index in seq_along(label_key$og_desired_variants)){
    columns<- label_key$columns_to_choose[[label_index]]; columns
    cat("labeling",paste(columns, collapse=","), "as", label_key$og_desired_variants[[label_index]] , "\n")
    for(mutation_name in columns){
      cells_with_mut <- rownames(cells_with_muts[cells_with_muts[[mutation_name]] == 1, ])
      s@meta.data$mutation_status[rownames(s@meta.data) %in% cells_with_mut] <- label_key$og_desired_variants[[label_index]]
    }
  }
  # REmove the WT+
  s@meta.data$mutation_status <- ifelse(s@meta.data$mutation_status == "WT",
                                        "WT",gsub("^WT\\+", "", s@meta.data$mutation_status))
  
  as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
  cat(magenta("Overwritting mutation profiles with BRAF in them as just BRAF \n"))
  # make any braf win, relabel as braf if it contains one.
  for (idx in seq_along(s@meta.data$mutation_status)) {
    i <- s@meta.data$mutation_status[idx]
    if (grepl("BRAF", i)) {
      #cat(blue(i), "\n")
      braf_mut <- regmatches(i, regexpr("BRAF\\.[^+]+", i))
      # cat(green("----", braf_mut, "\n"))
      if (length(braf_mut) == 0) { # if there's nothing after BRAF so liek just BRAF instead of "BRAF.G469E"
        braf_mut <- "BRAF"
      }
      s@meta.data$mutation_status[idx] <- braf_mut
    }
  }
  unique(s@meta.data[["mutation_status"]])
  return(s)
}

#sce_subset <- readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/two_cohorts/two_cohorts_sce.rds")
#s<- readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/two_cohorts/two_cohorts_SEURAT.rds")
# all_mutations_plus_no_overwrite<-function(sce_subset, s, original_desired_variants = NA, general=FALSE){ #original_desired_variants = c("BRAF", "NRAS", "KRAS", "PTPN11")
#   ngt <- assay(sce_subset, "NGT")
#   ngt_sparse <- as(ngt, "dgCMatrix")
#   ngt_sparse@x[ngt_sparse@x == 2] <- 1# replace all 2s with 1s bc like just BRAF present or not.
#   ngt <- as.matrix(ngt_sparse)
#   ngt <- as.data.frame(ngt) 
#   
#   mutations<- rownames(ngt)
#   label_key<- data.frame(og_desired_variants = original_desired_variants)
#   label_key$columns_to_choose <- list(NA)
#   if ("RAS" %in% label_key$og_desired_variants) {
#     #ras_mutations <- c("PTPN11", "NRAS", "KRAS")
#     ras_mutations <- rownames(ngt)[grepl("NRAS|PTPN11|KRAS", rownames(ngt))]
#     ras_matches <- unlist(lapply(ras_mutations, function(gene) {
#       mutations[grepl(gene, mutations)]
#     }))
#     ras_vars2 <- unique(ras_matches); print(ras_vars2)
#     label_key$columns_to_choose[label_key$og_desired_variants == "RAS"] <- list(ras_vars2)
#   }
#   for (var in label_key$og_desired_variants){
#     if(var != "RAS"){
#       for (i in mutations){ 
#         matches <- mutations[grepl(var, mutations)]
#       }
#       col3<- unique((matches))
#       cat("var:", var, "\n", "mut", col3, "\n")
#       label_key$columns_to_choose[label_key$og_desired_variants == var]<- list(col3)
#     }
#   }
#   label_key
#   
#   cat(blue("Identified cells with these variant mutations:", paste(mutations, collapse = ", ")), "\n")
#   cat(yellow("But we will identify", paste(label_key$og_desired_variants,collapse = ", "), "\n"))
#   
#   cells_with_mutations <- ngt[, colSums(ngt == 1) > 0] # select columns that have a 1 for any of their values
#   cells_with_muts<-as.data.frame(t(cells_with_mutations))
#   s@meta.data$mutation_status <- 'WT' # default is WT. 
#   
#   for(label_index in seq_along(label_key$og_desired_variants)){
#     columns<- label_key$columns_to_choose[[label_index]]; columns
#     cat("labeling",paste(columns, collapse=","), "as", label_key$og_desired_variants[[label_index]] , "\n")
#     for (mutation_name in columns) {
#       # Identify cells with the mutation
#       cells_with_mut <- rownames(cells_with_muts[cells_with_muts[[mutation_name]] == 1, ])
#       target_rows <- rownames(s@meta.data) %in% cells_with_mut
#       # Get current mutation status for these cells
#       current_status <- s@meta.data[target_rows, "mutation_status"]
#       # New label to add
#       new_label <- label_key$og_desired_variants[[label_index]]
#       if(general == TRUE & (str_detect(new_label, "BRAF") == FALSE)){
#         new_label<-gsub("\\..*", "", new_label)
#       }
#       # Only add if new_label is NOT already a substring of current_status
#       updated_status <- ifelse(
#         grepl(new_label, current_status),
#         current_status,                         # if already present, keep current
#         paste0(current_status, "+", new_label)   # if not, append
#       )
#       # Update the metadata
#       s@meta.data[target_rows, "mutation_status"] <- updated_status
#     }
#   }
#   # Remove the WT+ from the mutation_status plus list bc every single one would have it.
#   s@meta.data$mutation_status <- ifelse(s@meta.data$mutation_status == "WT",
#                                         "WT",gsub("^WT\\+", "", s@meta.data$mutation_status))
#   
#   as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
#   unique(s@meta.data[["mutation_status"]])
#   return(s)
# }


all_mutations_plus_no_overwrite<-function(sce_subset, s, original_desired_variants = NA, general=FALSE,
                                          include_individual_braf_muts = FALSE){ #original_desired_variants = c("BRAF", "NRAS", "KRAS", "PTPN11")
  ngt <- assay(sce_subset, "NGT")
  ngt_sparse <- as(ngt, "dgCMatrix")
  ngt_sparse@x[ngt_sparse@x == 2] <- 1# replace all 2s with 1s bc like just BRAF present or not.
  ngt <- as.matrix(ngt_sparse)
  ngt <- as.data.frame(ngt) 
  
  mutations<- rownames(ngt)
  label_key<- data.frame(og_desired_variants = original_desired_variants)
  label_key$columns_to_choose <- list(NA)
  if ("RAS" %in% label_key$og_desired_variants) {
    #ras_mutations <- c("PTPN11", "NRAS", "KRAS")
    ras_mutations <- rownames(ngt)[grepl("NRAS|PTPN11|KRAS", rownames(ngt))]
    ras_matches <- unlist(lapply(ras_mutations, function(gene) {
      mutations[grepl(gene, mutations)]
    }))
    ras_vars2 <- unique(ras_matches); print(ras_vars2)
    label_key$columns_to_choose[label_key$og_desired_variants == "RAS"] <- list(ras_vars2)
  }
  for (var in label_key$og_desired_variants){
    if(var != "RAS"){
      for (i in mutations){ 
        matches <- mutations[grepl(var, mutations)]
      }
      col3<- unique((matches))
      cat("var:", var, "\n", "mut", col3, "\n")
      label_key$columns_to_choose[label_key$og_desired_variants == var]<- list(col3)
    }
  }
  label_key
  if(include_individual_braf_muts == TRUE){
    # add back individual BRAF mutations. 
    braf_muts<- unlist(label_key[label_key$og_desired_variants=="BRAF",]$columns_to_choose)
    label_key <- label_key[label_key$og_desired_variants != "BRAF", ]
    nrow(label_key)
    for (braf_i in seq_along(braf_muts)) {
      label_key[nrow(label_key) + 1, ] <- list(og_desired_variants = braf_muts[[braf_i]], 
                                               columns_to_choose = braf_muts[[braf_i]])
    }
  }
  
  cat(blue("Identified cells with these variant mutations:", paste(mutations, collapse = ", ")), "\n")
  cat(yellow("But we will identify", paste(label_key$og_desired_variants,collapse = ", "), "\n"))
  
  cells_with_mutations <- ngt[, colSums(ngt == 1) > 0] # select columns that have a 1 for any of their values
  cells_with_muts<-as.data.frame(t(cells_with_mutations))
  s@meta.data$mutation_status <- 'WT' # default is WT. 
  
  for(label_index in seq_along(label_key$og_desired_variants)){
    columns<- label_key$columns_to_choose[[label_index]]; columns
    cat("labeling",paste(columns, collapse=","), "as", label_key$og_desired_variants[[label_index]] , "\n")
    for (mutation_name in columns) {
      # Identify cells with the mutation
      cells_with_mut <- rownames(cells_with_muts[cells_with_muts[[mutation_name]] == 1, ])
      target_rows <- rownames(s@meta.data) %in% cells_with_mut
      # Get current mutation status for these cells
      current_status <- s@meta.data[target_rows, "mutation_status"]
      # New label to add
      new_label <- label_key$og_desired_variants[[label_index]]
      if(general == TRUE & (str_detect(new_label, "BRAF") == FALSE)){
        new_label<-gsub("\\..*", "", new_label)
      }
      # Only add if new_label is NOT already a substring of current_status
      updated_status <- ifelse(
        grepl(new_label, current_status),
        current_status,                         # if already present, keep current
        paste0(current_status, "+", new_label)   # if not, append
      )
      # Update the metadata
      s@meta.data[target_rows, "mutation_status"] <- updated_status
    }
  }
  # Remove the WT+ from the mutation_status plus list bc every single one would have it.
  s@meta.data$mutation_status <- ifelse(s@meta.data$mutation_status == "WT",
                                        "WT",gsub("^WT\\+", "", s@meta.data$mutation_status))
  
  #  replace exact "BRAF.L79Q" with "WT"
  s@meta.data$mutation_status <- ifelse( s@meta.data$mutation_status == "BRAF.L79Q", "WT",gsub("(^|\\+)BRAF\\.L79Q(\\+|$)", "\\1", s@meta.data$mutation_status))
  s@meta.data$mutation_status <- gsub("^\\+|\\+$", "", s@meta.data$mutation_status) # : Remove any leftover '+' at start or end
  s@meta.data$mutation_status[s@meta.data$mutation_status == ""] <- "WT"
  
  as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
  unique(s@meta.data[["mutation_status"]])
  cat(green("Final labels:", paste0(unique(s@meta.data[["mutation_status"]]), collapse = ", "), "\n"))
  return(s)
}



construct_map<- function(labeled_s, title = " ", facet_map = FALSE){
  color_palette<- set_color_pallette()
  if(facet_map == FALSE){
    umap_plot <- DimPlot(labeled_s, reduction = "umap", 
                         group.by = "mutation_status", 
                         alpha=1,
                         pt.size =0.15) +
      scale_color_manual(values = color_palette) +
      labs(title = title) + theme(legend.text = element_text(size = 10))
    umap_plot
  } else{
    labeled_s$Group <- labeled_s@meta.data$Group
    umap_plot1 <- DimPlot(labeled_s[, labeled_s$Group == "BRAF Cohort"], reduction = "umap", group.by = "mutation_status",
                          split.by = "Group", raster = TRUE, alpha = 0.8) +
      scale_color_manual(values = color_palette)  +
      theme(legend.text = element_text(size = 10), plot.title = element_blank())
    umap_plot2 <- DimPlot(labeled_s[, labeled_s$Group == "RAS Cohort"], reduction = "umap", group.by = "mutation_status",
                          split.by = "Group", raster = TRUE, alpha = 0.8) +
      scale_color_manual(values = color_palette) + labs(title = title) +
      theme(legend.text = element_text(size = 10))
    umap_plot <- plot_grid(umap_plot1, umap_plot2, ncol =2)
  }
  
  print(umap_plot)
  return(umap_plot)
}


set_color_pallette <- function(){ 
  kp_color_palette <- c()
  kp_color_palette["WT"] <-"#e3e3e3" 
  kp_color_palette["gray_theme"] <-"#e3e3e3" 
  kp_color_palette["NRAS"] <- "#BC271C"#E87D72"
  kp_color_palette["KRAS"] <- "#fa9d11" # fa9111
  
  kp_color_palette["PTPN11"] <- "black"  
  kp_color_palette["RAS"]<- "#AC2A26"
  kp_color_palette["BRAF"] <- "#426FB0"#4a84c2" ##709dcc
  kp_color_palette["theme1"] <-"#4a84c2" 
  kp_color_palette["theme2"]<- "#AC2A26"
  kp_color_palette["BRAF.V600E"] <- "#416A75"
  kp_color_palette["BRAF.G469"] <-"#9AC065"  
  kp_color_palette["BRAF.G469A"] <-"#9AC065" 
  kp_color_palette["BRAF.D594"] <- "#BE40EE"
  kp_color_palette["BRAF.D594G"] <- "#BE40EE"
  kp_color_palette["BRAF Cohort"] <-"#3170B5"
  kp_color_palette["RAS Cohort"] <-"#AC2A26"
  kp_color_palette[" "] <-"white"
  
  kp_color_palette["non-BRAF non-Ras Mutant WT"]<- "darkgray"

  kp_color_palette["TET2"] <- "peachpuff3"  
  kp_color_palette["TP53"] <- "forestgreen"  
  kp_color_palette["RUNX1"] <- "#F0E442" 
  kp_color_palette["NPM1"] <- "salmon"  
  kp_color_palette["FLT3"] <- "purple4"
  
  return(kp_color_palette)
}


get_most_intense_color<- function(seurat_obj, sample_cd = cds_of_interest[[1]], color_scheme = c("lightgrey", "blue")){ 
  exp_plot <- suppressWarnings(FeaturePlot(seurat_obj, features = sample_cd, raster = TRUE, cols= color_scheme))
  og_colors<- prep_data(exp_plot)
  rgb_vals <- col2rgb(og_colors)  
  diff_vals <- apply(rgb_vals, 2, function(x) max(x) - min(x))  # Compute the difference between max and min for each color: to measure how gray it is
  ordered_colors <- og_colors[order(diff_vals)] # order colors from smallest difference (most gray) to largest (least gray)
  last_color <- ordered_colors[length(ordered_colors)]
  ordered_colors <- og_colors[order(diff_vals)]
  n <- length(ordered_colors)
  last_percent_colors <- ordered_colors[ceiling(n * 0.99):n]
  rgb_last <- col2rgb(last_percent_colors)
  mean_rgb <- rowMeans(rgb_last)
  mean_color <- rgb(mean_rgb[1]/255, mean_rgb[2]/255, mean_rgb[3]/255)
  print_color <- make_style(mean_color, colors = 256)
  cat(print_color("Most intense colors:", mean_color, reset("\n")))
  #plot(1, 1, col = mean_color, pch = 16, cex = 5, xlim = c(0, 2), ylim = c(0, 2), main =mean_color)
  return(mean_color)
}

prep_data<- function(exp_plot){
  plot_build <- ggplot_build(exp_plot[[1]])
  og_point_data<- plot_build[["data"]][[1]]
  og_colors<- og_point_data$colour 
  return(og_colors)
}

adjust_alpha <- function(umap_plot, which_group_to_adjust_alpha_for, this_group = 1, the_rest=0.18, pt_size = 1.4){
  color_palette<- set_color_pallette()
  legend_to_use<- get_legend(umap_plot)
  plot_data <- umap_plot$data %>%
    mutate(alpha_val = ifelse(mutation_status %in% which_group_to_adjust_alpha_for, this_group, the_rest)) 
  
  custom_umap_plot <- ggplot(plot_data, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = mutation_status, alpha = alpha_val), size = pt_size) +
    scale_color_manual(values = color_palette) +
    scale_alpha_identity() + 
    labs(title = umap_plot[[1]][["labels"]][["title"]]) + theme(legend.text = element_text(size = 10))+
    theme_classic() 
  
  #custom_umap_plot <- custom_umap_plot  + NoLegend() + legend_to_use
  custom_umap_plot
  return(custom_umap_plot)
}


relab_cells_for_variant <- function(umap, cells_positive_for_variant, variant_name, new_colors=c("red", "green", "yellow"), dont_overwrite= c("BRAF")){
  color_palette<- set_color_pallette()
  names(new_colors) <- variant_name
  newcolor_palette <- list()
  for (i in variant_name) {
    cat(i, new_colors[i], "\n") 
    newcolor_palette[[i]] <- new_colors[i]
  }
  color_palette <- c(color_palette, newcolor_palette)
  print(color_palette)
  
  legend_to_use<- get_legend(umap_plot)
  umap_plot$data$mutation_status<- as.character(umap_plot$data$mutation_status)
  # change mutation_status to variant_name for cells_positive_for_variant
  plot_data <- umap_plot$data %>%
    dplyr::mutate(
      mutation_status1 = ifelse(
        rownames(.) %in% cells_positive_for_variant &
          !sapply(mutation_status, function(x) any(grepl(x, dont_overwrite))),
        variant_name,
        mutation_status
      )
    )
  
  
  custom_umap_plot <- ggplot(plot_data, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = mutation_status1, alpha = 0.5), size = 0.6) +
    scale_color_manual(values = color_palette) +
    scale_alpha_identity() + 
    labs(title = umap_plot[[1]][["labels"]][["title"]]) + theme(legend.text = element_text(size = 10))+
    theme_classic() 
  
  #custom_umap_plot <- custom_umap_plot  + NoLegend() + legend_to_use
  custom_umap_plot
  return(custom_umap_plot)
}


