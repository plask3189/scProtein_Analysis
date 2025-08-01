
source("R/protein_kp/violin.R")
make_significance_brackets<- function(s, mutations_to_compare, original_violin_plot, compare_cohorts= FALSE, cluster_of_differentiation){
  s_extended<- s
  
  compare_cohorts <- TRUE
  if(compare_cohorts == TRUE){
    #muts_order <- c(muts_order, c("BRAF", "NRAS"))
    muts_order <- c("WT", "Non-RAS Mutant WT", "RAS", "BRAF.V600E", "BRAF.G469A", "BRAF.D594G")
  } else {  
    muts_order <- c("WT", "NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594") 
  }
  s_extended$mutation_status <- factor(s$mutation_status, levels = muts_order)
  original_violin_plot2 <- VlnPlot(
    s_extended,#[, !(s$mutation_status %in% braf_muts)], 
    features = cluster_of_differentiation,
    pt.size = 0.0005, alpha = 0.02, 
    group.by = "mutation_status",
    cols = c(color_palette[["BRAF Cohort"]], color_palette[["RAS Cohort"]]),
    split.by = 'Group' # split.plot = TRUE
  )
  original_violin_plot2
  
  
  
  ## get the max y of the violins we're supposed to put the brackets over. 
  mutations_to_compare_data <- original_violin_plot[[1]]$data[
    original_violin_plot[[1]]$data$ident %in% mutations_to_compare,]
  
  og_max_y <- max(mutations_to_compare_data[[cluster_of_differentiation]], na.rm = TRUE);og_max_y
  each_bracket_separated_by = og_max_y/10 ; each_bracket_separated_by
  # TEEEEEMP:
  #new_max_y <- 60
  #if (compare_cohorts == FALSE){
    new_max_y <- og_max_y + (length(mutations_to_compare)*each_bracket_separated_by) + each_bracket_separated_by ;new_max_y
    new_max_y<-og_max_y+(length(mutations_to_compare)*each_bracket_separated_by)
    og_1_val<- s@assays[["Protein"]]@data[cluster_of_differentiation,][1]; og_1_val
    s_extended@assays[["Protein"]]@data[cluster_of_differentiation,][1]<-new_max_y #og_max_y +20 # to extend axis. replaced with original values later
    s_extended$mutation_status <- factor(s_extended$mutation_status, levels = muts_order)
    cat(blue("-------", cluster_of_differentiation, "---------"))
    if (compare_cohorts == TRUE){
      new_violin_plot <- VlnPlot(
        s_extended,
        features = cluster_of_differentiation,
        pt.size = 0.0005, alpha = 0.02, 
        group.by = "mutation_status",
        cols = c(color_palette[["BRAF Cohort"]], color_palette[["RAS Cohort"]]),
        split.by = 'Group', y.max = 35)
      new_violin_plot
    } else {
      new_violin_plot <- VlnPlot(s_extended, features = cluster_of_differentiation,
                                 pt.size = 0.0005, alpha = 0.04, group.by = "mutation_status", cols = color_palette, y.max = new_max_y, same.y.lims = FALSE)+
        theme(axis.title.x = element_blank())+theme(axis.title.x = element_blank())
    }
  
  # } else {
  #   new_max_y <- og_max_y+ each_bracket_separated_by +0.5;new_max_y
  #   og_1_val<- s@assays[["Protein"]]@data[cluster_of_differentiation,][1]; og_1_val
  #   s_extended@assays[["Protein"]]@data[cluster_of_differentiation,][1]<-new_max_y #og_max_y +20 # to extend axis. replaced with original values later
  #   s_extended$mutation_status <- factor(s_extended$mutation_status, levels = muts_order)
  #   
  #   new_violin_plot <- VlnPlot(
  #     s_extended,#[, !(s$mutation_status %in% braf_muts)], 
  #     features = cluster_of_differentiation,
  #     pt.size = 0.0005, alpha = 0.01,
  #     group.by = "mutation_status",
  #     split.by = 'Group',
  #     cols = c(most_intense_colors_g1, most_intense_colors_g2),y.max = new_max_y, same.y.lims = FALSE)+ 
  #     theme(axis.title.x = element_blank())+theme(axis.title.x = element_blank(),
  #                                                 inherit.aes = FALSE )
  #   new_violin_plot
  # }
  # 
  # return dummy point to original value:
  new_violin_plot[[1]][["data"]][[cluster_of_differentiation]][1]<-og_1_val
  
  if (compare_cohorts == TRUE){
    mutations_to_compare_combos <- lapply(mutations_to_compare, function(x) c(x, x))
  } else {
    mutations_to_compare_combos <- combn(mutations_to_compare , 2, simplify = FALSE); mutations_to_compare_combos
  }
  
  num_comparisons <- (length(mutations_to_compare_combos));num_comparisons
  
  
  start_level<- og_max_y # start level of lowest bracket
  mutations_to_compare_combos
  violin_with_brackets<- actually_make_brackets(mutations_to_compare, muts_order, start_level,
                                                mutations_to_compare_combos, cluster_of_differentiation,
                                                do_compare_cohorts=compare_cohorts, new_violin_plot)
  
  return(violin_with_brackets)
}


actually_make_brackets2<- function(mutations_to_compare, muts_order, start_level, mutations_to_compare_combos,
                                  cluster_of_differentiation,new_violin_plot){
  start_level <- 28
  new_violin_plot<- original_violin_plot2
  combo_lengths <- sapply(mutations_to_compare_combos, function(x) {
    abs(which(muts_order == x[1]) - which(muts_order == x[2]))
  }); combo_lengths
  
  if (sum(combo_lengths)==0){
    combo_lengths <- rep(1, length(combo_lengths))
  }
  ordered_indices <- order(combo_lengths)
  mutations_to_compare_combos <- mutations_to_compare_combos[ordered_indices]
  
  
  used_levels <- list()
  bracket_levels <- numeric(length(mutations_to_compare_combos))
  space_between_levels <- 1.2
  for (i in seq_along(mutations_to_compare_combos)) {
    combo <- mutations_to_compare_combos[[i]]
    left <- which(muts_order == combo[1])
    right <- which(muts_order == combo[2])
    range_i <- left:right
    level_found <- FALSE
    for (lvl in seq_along(used_levels)) {
      overlap <- any(range_i %in% used_levels[[lvl]])
      if (!overlap) {
        used_levels[[lvl]] <- c(used_levels[[lvl]], range_i)
        bracket_levels[i] <- start_level + (lvl - 1) * 1
        level_found <- TRUE
        break
      }
    }
    
    if (!level_found) {
      used_levels[[length(used_levels) + 1]] <- range_i
      bracket_levels[i] <- start_level + (length(used_levels) - 1) * space_between_levels
    }
  }
  
  bracket_df <- data.frame(
    Combos = sapply(mutations_to_compare_combos, function(x) paste(x, collapse = "_vs_")),
    bracket_level = bracket_levels+0.5,
    left_x_name = sapply(mutations_to_compare_combos, function(x) x[1]),
    left_x_pos = sapply(mutations_to_compare_combos, function(x) which(muts_order == x[1])),
    right_x_name = sapply(mutations_to_compare_combos, function(x) x[2]),
    right_x_pos = sapply(mutations_to_compare_combos, function(x) which(muts_order == x[2]))
  )
 
  size<-3

    violin_comparison <- new_violin_plot+ 
      stat_compare_means(
        comparisons = (mutations_to_compare_combos),
        paired = FALSE,
        bracket.size = 0.1,
        method = "wilcox.test",
        label = "p.signif",
        tip.length = 0.02,
        size = size,
        vjust =1,# space_between_levels, #space_between_levels+((size/2)-1), # 1.5,
        label.y = bracket_df$bracket_level) 
    
  print(violin_comparison)
  return(violin_comparison)
}


actually_make_brackets<- function(mutations_to_compare, muts_order, start_level, mutations_to_compare_combos,
                                  cluster_of_differentiation,
                                  do_compare_cohorts, new_violin_plot){
  
  combo_lengths <- sapply(mutations_to_compare_combos, function(x) {
    abs(which(muts_order == x[1]) - which(muts_order == x[2]))
  }); combo_lengths
  
  if (sum(combo_lengths)==0){
    combo_lengths <- rep(1, length(combo_lengths))
  }
  ordered_indices <- order(combo_lengths)
  mutations_to_compare_combos <- mutations_to_compare_combos[ordered_indices]
  
  
  used_levels <- list()
  bracket_levels <- numeric(length(mutations_to_compare_combos))
  space_between_levels <- 1.2
  for (i in seq_along(mutations_to_compare_combos)) {
    combo <- mutations_to_compare_combos[[i]]
    left <- which(muts_order == combo[1])
    right <- which(muts_order == combo[2])
    range_i <- left:right
    level_found <- FALSE
    for (lvl in seq_along(used_levels)) {
      overlap <- any(range_i %in% used_levels[[lvl]])
      if (!overlap) {
        used_levels[[lvl]] <- c(used_levels[[lvl]], range_i)
        bracket_levels[i] <- start_level + (lvl - 1) * 1
        level_found <- TRUE
        break
      }
    }
    
    if (!level_found) {
      used_levels[[length(used_levels) + 1]] <- range_i
      bracket_levels[i] <- start_level + (length(used_levels) - 1) * space_between_levels
    }
  }
  
  bracket_df <- data.frame(
    Combos = sapply(mutations_to_compare_combos, function(x) paste(x, collapse = "_vs_")),
    bracket_level = bracket_levels+0.5,
    left_x_name = sapply(mutations_to_compare_combos, function(x) x[1]),
    left_x_pos = sapply(mutations_to_compare_combos, function(x) which(muts_order == x[1])),
    right_x_name = sapply(mutations_to_compare_combos, function(x) x[2]),
    right_x_pos = sapply(mutations_to_compare_combos, function(x) which(muts_order == x[2]))
  )
  
  # If left_x_pos equals right_x_pos, add 1 to right_x_pos
  condition <- bracket_df$left_x_pos == bracket_df$right_x_pos
  bracket_df$left_x_pos[condition] <- bracket_df$left_x_pos[condition]-0.3
  bracket_df$right_x_pos[condition] <- bracket_df$right_x_pos[condition]+0.2
  bracket_df
  size<-3
  if(do_compare_cohorts == TRUE){
    p_vals_table<-compute_p_vals(s, cluster_of_differentiation, mutations_to_compare)
    p_val_data<- left_join(bracket_df, p_vals_table, by = c("left_x_name"= "mutation_status"))
    p_val_data$y.position <- p_val_data$bracket_level
    p_val_data$group1 <- p_val_data$left_x_pos
    p_val_data$group2 <- p_val_data$right_x_pos
    violin_comparison<- new_violin_plot + 
      stat_pvalue_manual(
        p_val_data,
        label = "p.signif", # show significance symbols, like *, **, etc
        tip.length = 0.01,
        bracket.size = 0.1,
        size = size,
        vjust = 0.7,
        inherit.aes = FALSE 
      ); violin_comparison
  } else {
    violin_comparison <- new_violin_plot+ 
      stat_compare_means(
        comparisons = (mutations_to_compare_combos),
        paired = FALSE,
        bracket.size = 0.1,
        method = "wilcox.test",
        label = "p.signif",
        tip.length = 0.02,
        size = size,
        vjust =1,# space_between_levels, #space_between_levels+((size/2)-1), # 1.5,
        label.y = bracket_df$bracket_level) 
    
    
  }
  
  print(violin_comparison)
  return(violin_comparison)
}



