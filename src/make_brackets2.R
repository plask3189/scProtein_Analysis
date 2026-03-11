
source("R/protein_kp/violin.R")
violins_list<- list()
for(plot_name in names(list_of_subfigs)){
  if (grepl("violin", plot_name, ignore.case = TRUE)) {
    violins_list[[plot_name]] <-  list_of_subfigs[[plot_name]] + theme(axis.title.x = element_blank())
  }
}
sample_violin<- violins_list[[3]]
print(names(violins_list))


# just for one cohort
make_combo_brackets<- function(sample_violin){
  rebuild_violin() 
  cluster_of_differentiation<- colnames(sample_violin[[1]][["data"]])[[1]]
  muts_order <- c("WT", "NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594") 
  
  ## get the max y of the violins we're supposed to put the brackets over. 
  mutations_to_compare_data <- sample_violin[[1]]$data[sample_violin[[1]]$data$ident %in% muts_order,]
  levels(mutations_to_compare_data$ident)
  max_data_val <- max(mutations_to_compare_data[[cluster_of_differentiation]], na.rm = TRUE);og_max_y
  y_axis_length<- sample_violin[[1]][["coordinates"]][["limits"]][["y"]][[2]]

  mutations_to_compare_combos <- combn(muts_order , 2, simplify = FALSE); mutations_to_compare_combos
  mutations_to_compare_combos_braf <- Filter(function(x) any(grepl("BRAF", x)), mutations_to_compare_combos)


  bracket_df<- make_distance_bracket_df(mutations_to_compare_combos_braf, max_data_val, y_axis_length)
  summary(bracket_df$bracket_level)
  summary(mutations_to_compare_data[[cluster_of_differentiation]])
  y_axis_length
  size<-3
  
  bracket_plot<- rebuild_violin(sample_violin)
  
  sample_violin + bracket_plot
    
  print(violin_comparison)
  
  return()
}

rebuild_violin<- function(sample_violin){
  muts_order <- c("WT", "NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594") 
  
  # stat.test <- compare_means(
  #   len ~ dose, data = sample_violin[[1]][["data"]],
  #   method = "t.test"
  # )
  # stat.test
  # 
  # 
  
  mutations_to_compare_data$ident <- factor(mutations_to_compare_data$ident, levels = muts_order)
  
  p <- ggplot(mutations_to_compare_data, aes(x = ident, y = !!sym(cluster_of_differentiation), fill = ident)) +
    geom_violin() +
    geom_jitter(width = 0.2, size = 0.01, alpha = 0.05) +
    stat_compare_means(
      comparisons = bracket_df$comparison,
      label.y = bracket_df$bracket_level,
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.02,
      size = 3
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, y_axis_length), clip = "off") +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 1, 3, 1), "lines"))

  built_p <- ggplot_build(p)
  which_signif <- sapply(p$layers, function(layer) inherits(layer$geom, "GeomSignif"))

  signif_data <- built_p$data[[which(which_signif)]]
  signif_data <- subset(signif_data, annotation != "ns")

  # Keep only the horizontal segments (x != xend & y == yend)
  annotation_segments <- subset(signif_data, x != xend & y == yend & annotation != "ns")
  
  bracket_plot <- ggplot(signif_data, aes(x = x, xend = xend, y = y, yend = yend)) +
    geom_segment(linewidth = 0.2) +  # all lines
    geom_text(
      data = annotation_segments,
      aes(x = (x + xend) / 2, y = y, label = annotation),
      vjust = 0, size = 3
    ) +
    scale_x_continuous(breaks = unique(c(signif_data$x, signif_data$xend))) +
    theme_minimal()+ theme(
      panel.grid = element_blank(),        # removes both major and minor grid lines
      axis.line = element_blank(),         # removes axes lines
      axis.ticks = element_blank(),        # removes tick marks
      axis.text = element_blank(),         # removes axis text
      axis.title = element_blank()         # removes axis titles
    )
  bracket_plot
  

  return(bracket_plot)
  
}




make_distance_bracket_df<- function(mutations_to_compare_combos_braf, max_data_val, y_axis_length){
  
  combo_lengths <- sapply(mutations_to_compare_combos_braf, function(x) {
    abs(which(muts_order == x[1]) - which(muts_order == x[2]))
  }); combo_lengths
  ordered_indices <- order(combo_lengths,  decreasing = TRUE) # indices of longest brackets
  mutations_to_compare_combos_braf <- mutations_to_compare_combos_braf[ordered_indices]
  used_levels <- list()
  bracket_levels <- numeric(length(mutations_to_compare_combos_braf))
  
  min_y_for_brackets<- max_data_val+0.2 ; min_y_for_brackets
  max_y_for_brackets <- y_axis_length-0.2 ;max_y_for_brackets
  num_brackets<- length(mutations_to_compare_combos_braf)
  space_for_each_bracket <- (max_y_for_brackets-min_y_for_brackets)/num_brackets
  
  space_between_levels <- space_for_each_bracket
  for (i in seq_along(mutations_to_compare_combos_braf)) {
    combo <- mutations_to_compare_combos[[i]]
    left <- which(muts_order == combo[1])
    right <- which(muts_order == combo[2])
    range_i <- left:right
    level_found <- FALSE
    for (lvl in seq_along(used_levels)) {
      overlap <- any(range_i %in% used_levels[[lvl]])
      if (!overlap) {
        used_levels[[lvl]] <- c(used_levels[[lvl]], range_i)
        bracket_levels[i] <- min_y_for_brackets + (lvl - 1) * 1
        level_found <- TRUE
        break
      }
    }
    if (!level_found) {
      used_levels[[length(used_levels) + 1]] <- range_i
      bracket_levels[i] <- max_data_val + (length(used_levels) - 1) * space_between_levels
    }
  }
  
  bracket_df <- data.frame(
    comparison = I(mutations_to_compare_combos_braf),  # list-column of pairs
    Combos = sapply(mutations_to_compare_combos_braf, function(x) paste(x, collapse = "_vs_")),
    bracket_level = bracket_levels,
    left_x_name = sapply(mutations_to_compare_combos_braf, function(x) x[1]),
    left_x_pos = sapply(mutations_to_compare_combos_braf, function(x) which(muts_order == x[1])),
    right_x_name = sapply(mutations_to_compare_combos_braf, function(x) x[2]),
    right_x_pos = sapply(mutations_to_compare_combos_braf, function(x) which(muts_order == x[2]))
  )
  
  return(bracket_df)
}

