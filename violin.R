library(dplyr)
library(ggpubr)


# the different mutations actually have different immunophenotypes  like G469 looks likes its the CD34 cellswhereas V600E and D594 are CD11b and CD14/16
#🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻🎻
violin<- function(s, cluster_of_differentiation, compare_cohort_violins=FALSE, most_intense_colors_g1 = "blue", most_intense_colors_g2 = "red" , mutations_to_compare = NA){
  color_palette<- set_color_pallette()
  
  #mutations_to_compare<- c( "WT", "NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594")
  if(compare_cohort_violins == TRUE){
    #braf_muts <- c("BRAF.G469", "BRAF.D594", "BRAF.V600E")
    x_axis_labels <- c("WT", "Non-RAS Mutant WT", "RAS", "BRAF.V600E", "BRAF.G469", "BRAF.D594")
    muts_order <- c("WT", "Non-RAS Mutant WT", "RAS", "BRAF.V600E", "BRAF.G469A", "BRAF.D594G") # order of violins
    s$mutation_status <- factor(s$mutation_status, levels = muts_order)
    
    original_violin_plot <- VlnPlot(
      s,#[, !(s$mutation_status %in% braf_muts)], 
      features = cluster_of_differentiation,
      pt.size = 0.0005, alpha = 0.02, 
      group.by = "mutation_status",
      cols = c(color_palette[["BRAF Cohort"]], color_palette[["RAS Cohort"]]),
      split.by = 'Group' # split.plot = TRUE
    )+ scale_x_discrete(labels = x_axis_labels) 
    original_violin_plot
    
    # violin_braf_comparison<- make_significance_brackets(s, mutations_to_compare, original_violin_plot, compare_cohorts= TRUE,
    #                                                     cluster_of_differentiation)
    
  } else{
    muts_order <- c("WT", "NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594", "RAS", "BRAF") # order of violins
    s$mutation_status <- factor(s$mutation_status, levels = muts_order)
    cat(blue("-------", cluster_of_differentiation, "----------\n"))
    original_violin_plot <- suppressWarnings(
      VlnPlot(s, features = cluster_of_differentiation,
              pt.size = 0.01, alpha = 0.17, 
              group.by = "mutation_status", cols = color_palette)) +
      theme(axis.title.x = element_blank());original_violin_plot
    
    #violin_braf_comparison<- make_significance_brackets(s, mutations_to_compare, original_violin_plot, compare_cohorts= FALSE)
  }
  # the stars
  #make_baby_legend(s, cluster_of_differentiation, mutations_to_compare, original_violin_plot)
  
  
  return(original_violin_plot)
}


compute_p_vals<- function(seurat_obj, cluster_of_differentiation, compare_what){
  comparing_groups <- FALSE # init
    if("BRAF Cohort" %in% compare_what){ # comparing groups.
      comparing_groups <- TRUE
      data <- FetchData(object = seurat_obj, vars = c(cluster_of_differentiation, "Group", "mutation_status"), layer = "data")
      fml <- reformulate("Group", response = cluster_of_differentiation)
    } else{
      
      data<- FetchData(object = seurat_obj, vars = c(cluster_of_differentiation, "mutation_status" ), layer = "data")
      fml <- reformulate("mutation_status", response = cluster_of_differentiation)
    }

    pvals <- compare_means(
      formula = fml,
      data = data,
      method = "wilcox.test",
      comparisons = compare_what, # mutations or groups
      paired = FALSE
    )
    pvals<-as.data.frame(pvals)
    pvals$Cohort <- NA
    pvals$mutation_status <- NA
    if(comparing_groups == TRUE){
      pvals$mutation_status <- unique(data$mutation_status)[1]
    } else {
      pvals$Cohort <- unique(seurat_obj@meta.data[["Group"]])[1]
    }
    return(pvals)
}
# how likely it is to observe a difference in the data (e.g., in CD34 expression) 
#as large or larger than what you observed, if there were actually no difference between 
#the groups.
get_p_vals_table<- function(s, cluster_of_differentiation, mutations_to_compare, compare_groups = FALSE){
  if (compare_groups==TRUE){
    Idents(s) <-s@meta.data[["Group"]]
    cat("Computing p vals for comparing groups \n")
    #pvals<- compute_p_vals_group(s, cluster_of_differentiation, mutations_to_compare)
    s_g1 <- subset(s, subset = Group == "BRAF Cohort")
    pvalss_g1<- compute_p_vals(s_g1, cluster_of_differentiation, compare_what = mutations_to_compare)
    
    s_wt <- subset(s, subset = mutation_status == "WT")
    pvalss_wt<- compute_p_vals(s_wt, cluster_of_differentiation, compare_what = c("BRAF Cohort", "RAS Cohort"))
    
    non_ras_mut <- subset(s, subset = mutation_status == "Non-RAS Mutant WT")
    pvals_non_ras_mut<- compute_p_vals(non_ras_mut, cluster_of_differentiation, compare_what = c("BRAF Cohort", "RAS Cohort"))
    
    ras_mut <- subset(s, subset = mutation_status == "RAS")
    pvalsras_mut<- compute_p_vals(ras_mut, cluster_of_differentiation, compare_what = c("BRAF Cohort", "RAS Cohort"))
    
    common_cols <- intersect(colnames(pvalss_g1), colnames(pvalsras_mut))
    pvalss_g1 <- pvalss_g1[, common_cols]
    pvalss_wt <- pvalss_wt[, common_cols]
    pvals_non_ras_mut <- pvals_non_ras_mut[, common_cols]
    
    pvals<-rbind(pvalss_g1, pvalss_wt, pvals_non_ras_mut, pvalsras_mut)
    excel_file_name<- paste0("two_cohorts/", cluster_of_differentiation, "_pval_results_two_cohorts.xlsx")
    write_xlsx(as.data.frame(pvals), excel_file_name) 
    
    pval_text <- paste(capture.output(print(pvals)), collapse = "\n")
    pval_text <- paste(cluster_of_differentiation, "\n", pval_text)
    pval_plot <- ggplot() +
      annotate("text", x = 0, y = 1, label = pval_text, hjust = 0, vjust = 1, size = 2, family = "mono") +
      xlim(0, 0.5) + ylim(0, 1) +
      theme_void()

  } else {
    pvals<- compute_p_vals(s, cluster_of_differentiation, compare_what = mutations_to_compare)
    excel_file_name<- paste0("just_braf/", cluster_of_differentiation, "_pval_results_just_braf.xlsx")
    write_xlsx(as.data.frame(pvals), excel_file_name) 
    pval_text <- paste(capture.output(print(pvals)), collapse = "\n")
    pval_text <- paste(cluster_of_differentiation, "\n", pval_text)
    pval_plot <- ggplot() +
      annotate("text", x = 0, y = 1, label = pval_text, hjust = 0, vjust = 1, size = 3, family = "mono") +
      xlim(0, 1) + ylim(0, 1) +
      theme_void()
  }
  
  #ggsave(file.path(path, "pval_table_plot.png"), plot = pval_plot, width = 8, height = 6, dpi = 300)
  return(pval_plot)
}






make_baby_legend<- function(s, cluster_of_differentiation, mutations_to_compare, original_violin_plot){
  pvals<- compute_p_vals(s, cluster_of_differentiation, mutations_to_compare)
  pvals_expanded<- pvals
  for (i in seq_len(nrow(pvals))){
    row<- pvals[i,]
    new_row<-row
    new_row$group1<-row$group2
    new_row$group2<-row$group1
    exists <- any( # to check if the reversed pair is already in pvals
      pvals$group1 == new_row$group1 & pvals$group2 == new_row$group2
    )
    if (!exists) {
      pvals_expanded <- rbind(pvals_expanded, new_row)
    }
  }
 
  color_palette <- set_color_pallette() 
  color_palette["WT"] <- "#4D4D4D" 
  
  pvals_expanded$color_of_group_2 <- color_palette[as.character(pvals_expanded$group2)]
  
  pvals_expanded$group2 <- factor(pvals_expanded$group2, levels = muts_order)
  pvals_expanded <- pvals_expanded %>%
    group_by(group1) %>%
    arrange(group2, .by_group = TRUE)
  
  violin_hat_p_vals <- data.frame(violin_label = unique(pvals$group1),violin_hat_stars = NA)
  violin_hat_plots<- list()
  
  for (group in unique(mutations_to_compare)) {
    pvals_group <- pvals_expanded[pvals_expanded$group1 == group, ];pvals_group
    pvals_group$y <- nrow(pvals_group):1  # y positions
    
    p <- ggplot(pvals_group, aes(x = 1, y = y, label = p.signif, color = color_of_group_2)) +
      geom_text(aes(size = ifelse(p.signif == "ns", 2, 2.5))) +  # conditional size
      scale_size_identity() +                                   # interpret as literal sizes
      scale_color_identity() +                                  # interpret hex codes directly
      xlim(0, 2) +
      ylim(0, nrow(pvals_group)) +
      theme_void() +
      theme(legend.position = "none");p
    violin_hat_plots[[group]] <- p
  }
  return(violin_hat_plots)
}


put_on_star_hats<- function(violin_hat_plots, original_violin_plot){
  expr_col <- colnames(original_violin_plot[[1]][["data"]])[1]
  violin_position_data <- original_violin_plot[[1]][["data"]] %>%
    group_by(ident) %>%
    summarize( x_pos = mean(as.numeric(ident)),
               max_y = max(.data[[expr_col]], na.rm = TRUE))
  accessorized_violin<- original_violin_plot
  for(i in seq_len(nrow(violin_position_data))){
    star_grob_plot <- ggplotGrob(violin_hat_plots[[i]])
    xmin_for_stars <- violin_position_data$x_pos[[i]] -1
    xmax_for_stars <- violin_position_data$x_pos[[i]] +1
    ymin_for_stars <- violin_position_data$max_y[[i]]
    ymax_for_stars <- violin_position_data$max_y[[i]] + 7
    accessorized_violin<-  accessorized_violin + annotation_custom(grob = star_grob_plot, 
                                             xmin = xmin_for_stars, 
                                             xmax = xmax_for_stars, 
                                             ymin = ymin_for_stars, 
                                             ymax = ymax_for_stars)
    
  }
  accessorized_violin
}




