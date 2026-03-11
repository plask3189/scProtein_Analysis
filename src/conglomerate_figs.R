library(ggpubr)
library(cowplot)
library(ggplot2)
# organize the subfigs. 
# three_mutaneers,       just_braf_cohort_braf_vs_ras, 
# just_braf_CD34_exp_umap,    just_braf_CD11b_exp_umap,
# just_braf_CD34violins, just_braf_CD11bviolins,        both_cohorts_ras_vs_braf_mutations_labeled_umap,
# BRAF_CD34_exp_umap,    RAS_CD34_exp_umap,             both_cohorts_CD34violins, 
# BRAF_CD11b_exp_umap,   RAS_CD11b_exp_umap,            both_cohorts_CD11bviolins
#                           dotplot
do_main_fig <- TRUE 
if (do_main_fig == FALSE){
  all_plots_dir<- "supplimentary"
  path<- "supplimentary"
} else {
  all_plots_dir<- "all_plots"
  path<- "all_plots"
}

rds_files <- list.files(all_plots_dir, pattern = "\\.rds$", full.names = TRUE); rds_files# List all .rds files in the directory
rds_names <- tools::file_path_sans_ext(basename(rds_files))
og_list_of_subfigs <- setNames(lapply(rds_files, readRDS), rds_names) 

if (do_main_fig == FALSE){
  # remove legends only of some plots, 
  list_of_subfigs <- strip_legends(og_list_of_subfigs, plots_to_strip_legend = c( "just_braf_CD123violins", 
                                                                                    "just_braf_CD14violins", 
                                                                                   "just_braf_FLT3violins",
                                                                                    "both_cohorts_CD123violins",
                                                                                   "both_cohorts_CD14violins")); names(list_of_subfigs) <- rds_names
} else {
  # remove legends only of some plots, 
  list_of_subfigs <- strip_legends(og_list_of_subfigs, plots_to_strip_legend = c("three_mutaneers","just_braf_cohort_braf_vs_ras","just_braf_CD34violins",
                                                                                 "just_braf_CD11bviolins", "both_cohorts_CD11bviolins","both_cohorts_CD34violins",
                                                                                 "BRAF_CD34_exp_umap", "RAS_CD34_exp_umap",
                                                                                 "BRAF_CD11b_exp_umap", "RAS_CD11b_exp_umap",
                                                                                 "both_cohorts_ras_vs_braf_mutations_labeled_umap")); names(list_of_subfigs) <- rds_names
}
legend_to_use <- make_main_legend(legend_text = c("WT", "BRAF.V600E", "BRAF.G469","BRAF.D594", "NRAS", "KRAS", "PTPN11", "BRAF", "RAS"))
legend_to_use <- make_main_legend(legend_text = c("WT", "BRAF", "RAS"))
list_of_subfigs <- lapply(list_of_subfigs, function(p) {p + theme(plot.title = element_text(size = 14))})

if (do_main_fig == FALSE){
  supp_fig <- make_supp_fig(list_of_subfigs)
  supp_fig_final_plot <- ggdraw() +
    draw_plot(supp_fig, x = 0, y = 0, width = 0.9, height = 1) +  
    draw_grob(legend_to_use, x = 0.9, y = 0.86, width = 0.1, height = 0.17) 
  ggsave(file.path(path, "supp_fig_final_plot.png"), plot = supp_fig_final_plot, width = 25, height =30, dpi = 150)
} else {
  main_fig<- make_main_fig(list_of_subfigs)
  final_plot <- ggdraw() +
    draw_plot(main_fig, x = 0, y = 0, width = 0.9, height = 1) +  
    draw_grob(legend_to_use, x = 0.9, y = 0.86, width = 0.1, height = 0.17)
  ggsave(file.path(path, "final_plot.png"), plot = final_plot, width = 21, height =36, dpi = 150)
}








#----------------------------------------------------
make_main_fig<- function(list_of_subfigs) {
  plot_x_val<- 7.5
  plot_y_val<-6
  list_of_subfigs <- lapply(list_of_subfigs, function(p) {
    p + theme(plot.title = element_text(size = 25),
      axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20),
      axis.title.x =  element_text(size = 18), axis.title.y =  element_text(size = 18),
      legend.text = element_text(size = 16), legend.title = element_text(size = 16),
      legend.key.size = unit(1.2, "lines")
    )})
  for(plot_name in names(list_of_subfigs)){
    if (grepl("violin", plot_name, ignore.case = TRUE)) {
      list_of_subfigs[[plot_name]] <-  list_of_subfigs[[plot_name]] + theme(axis.title.x = element_blank())
    }
    if (grepl("dot", plot_name, ignore.case = TRUE)) {
      list_of_subfigs[[plot_name]] <-  list_of_subfigs[[plot_name]] + theme(axis.title.x = element_blank(),
                                                                            axis.title.y = element_blank(),
                                                                            axis.text.x =  element_text(size = 17))
    }
  }

  two_figs_across_section_list <- list_of_subfigs[c("three_mutaneers","just_braf_cohort_braf_vs_ras","just_braf_CD34_exp_umap","just_braf_CD11b_exp_umap")]
  #two_figs_across_section <- plot_grid(plotlist = two_figs_across_section_list, ncol = 2, label_size = 26,align = "hv", axis = "tblr", labels = c("D", "E", "F", "G")) 
  two_figs_across_section_list_col1 <- list_of_subfigs[c("three_mutaneers","just_braf_CD34_exp_umap")] 
  two_figs_across_section_list_col2 <-list_of_subfigs[c("just_braf_cohort_braf_vs_ras", "just_braf_CD11b_exp_umap")]
  
  # Step 1: Create each column (vertical stack)
  col1 <- plot_grid(
    plotlist = two_figs_across_section_list_col1,
    ncol = 1,
    labels = c("D", "E"),
    label_size = 26,
    align = "v"
  )
  
  col2 <- plot_grid(
    plotlist = two_figs_across_section_list_col2,
    ncol = 1,
    labels = c("F", "G"),
    label_size = 26,
    align = "v"
  )
  
  # Step 2: Combine the two columns side by side
  two_figs_across_section <- plot_grid(
    col1,
    NULL,   # spacer column
    col2,
    ncol = 3,
    rel_widths = c(1, 0, 1),  # 0.1 controls space width
    align = "h",
    axis = "tblr"
  )
  #two_figs_across_section <- plot_grid(plotlist = two_figs_across_section_list, ncol = 2, label_size = 26,align = "hv", axis = "tblr", labels = c("D", "E", "F", "G")) 
  
  ggsave(file.path(path, "main_fig_part_1.pdf"), plot = two_figs_across_section, width = (plot_x_val *2), height = (plot_y_val *2))
  odd_line <- plot_grid(list_of_subfigs[["just_braf_CD34violins"]],
    list_of_subfigs[["just_braf_CD11bviolins"]],
    list_of_subfigs[["both_cohorts_ras_vs_braf_mutations_labeled_umap"]],
    ncol = 3, label_size = 26,align = "hv", axis = "tblr",
    rel_widths = c(0.55, 0.55, 1), labels = c("H", "I", "J")) + 
    theme(plot.margin = unit(c(5, 50, 5, 5), "pt"))
  #+ theme(plot.margin = margin( t = 5, r = 50, b = 5, l = 5, unit = "pt"))

  # odd_line <- plot_grid(list_of_subfigs[["just_braf_CD34violins"]],
  #                       list_of_subfigs[["just_braf_CD11bviolins"]],
  #                       list_of_subfigs[["both_cohorts_ras_vs_braf_mutations_labeled_umap"]], 
  #                       ncol = 3, label_size = 26,align = "hv", axis = "tblr", labels =c("H", "I", "J"))
  # 
  #ggsave(file.path(path, "main_fig_part_2.pdf"), plot = odd_line, width = (plot_x_val *2)+3, height = (plot_y_val *2))
  
  three_figs_across_section_list <- list_of_subfigs[c(
    "BRAF_CD34_exp_umap", "RAS_CD34_exp_umap", "both_cohorts_CD34violins", 
    "BRAF_CD11b_exp_umap", "RAS_CD11b_exp_umap", "both_cohorts_CD11bviolins")]
  # three_figs_across_section <- plot_grid(plotlist = three_figs_across_section, ncol = 3, label_size = 26,
  #                                        rel_widths = c(1, 1, 0.7), align = "hv", axis = "tblr", labels = c("K", "L", "M", "N", "O", "P"))
  three_figs_across_section <- plot_grid(plotlist = three_figs_across_section_list, ncol = 3, label_size = 26,
                                         rel_widths = c(1, 1, 0.2), align = "hv", axis = "tblr", labels = c("K", "L", "M", "N", "O", "P")) + theme(plot.margin = unit(c(0,0,0,0), "pt"))
  
  #ggsave(file.path(path, "main_fig_part_3.pdf"), plot = three_figs_across_section, width = (plot_x_val *3), height = (plot_y_val *2)+5)
  
  dotplot2<- plot_grid(list_of_subfigs[["dotplot"]], ncol = 1, label_size = 26,labels = c("Q")) + theme(axis.title.x = element_blank(),
                                                                                                      axis.title.y = element_blank())
  #ggsave(file.path(path, "dotplot.png"), plot = dotplot2, width =20, height = 6)
  main_fig <- plot_grid(
    two_figs_across_section, 
    odd_line, 
    three_figs_across_section, 
    dotplot2,
    ncol = 1, 
    label_size = 26,
    rel_heights = c(0.9, 0.48, 1.05, 0.25), # 0.25),
    rel_widths = c(0.9, 0.6, 1, 0.81)
  );cat(green("Made main fig. \n"))
  #ggsave(file.path(path, "main_fig.png"), plot = main_fig, width = (plot_x_val *2), height = (plot_y_val *6))
  
  
  final_plot <- ggdraw() +
    draw_plot(main_fig, x = 0, y = 0, width = 0.9, height = 1) +  
    draw_grob(legend_to_use, x = 0.9, y = 0.86, width = 0.1, height = 0.17)
  ggsave(file.path(path, "final_plot.png"), plot = final_plot, width = 21, height =36, dpi = 150)
  return(main_fig)
}




#----------------------------------------------------
make_supp_fig<- function(list_of_subfigs) {
  plot_x_val<- 7.5
  plot_y_val<-6
  
  # increase text size:
  list_of_subfigs <- lapply(list_of_subfigs, function(p) {
    p + theme(
      plot.title = element_text(size = 25),
      axis.text.x = element_text(size = 20),  # X-axis numbers
      axis.text.y = element_text(size = 20),
      axis.title.x =  element_text(size = 18),
      axis.title.y =  element_text(size = 18),
      legend.text = element_text(size = 16),      # Increase legend text
      legend.title = element_text(size = 16),
      legend.key.size = unit(1.2, "lines")
    )
  })
  for(plot_name in names(list_of_subfigs)){
    if (grepl("violin", plot_name, ignore.case = TRUE)) {
      list_of_subfigs[[plot_name]] <-  list_of_subfigs[[plot_name]] + theme(axis.title.x = element_blank())
    }
    if (grepl("dot", plot_name, ignore.case = TRUE)) {
      list_of_subfigs[[plot_name]] <-  list_of_subfigs[[plot_name]] + theme(axis.title.x = element_blank(),
                                                                            axis.title.y = element_blank(),
                                                                            axis.text.x =  element_text(size = 17))
    }
  }
  
  
  # ---------- three across ------------------------
  three_figs_across_section <- list_of_subfigs[c(
    "just_braf_CD14_exp_umap", "just_braf_CD123_exp_umap", "just_braf_FLT3_exp_umap", 
    "just_braf_CD14violins",    "just_braf_CD123violins",   "just_braf_FLT3violins",
    "BRAF_CD14_exp_umap", "RAS_CD14_exp_umap","both_cohorts_CD14violins",
    "BRAF_CD123_exp_umap", "RAS_CD123_exp_umap", "both_cohorts_CD123violins")]
  three_figs_across_section <- plot_grid(plotlist = three_figs_across_section, ncol = 3, label_size = 26,
                                         align = "hv", axis = "tblr", labels = "AUTO")
  ggsave(file.path(path, "three_figs_across_section.png"), plot = three_figs_across_section, width = (plot_x_val *2.8), height = (plot_y_val *4))
  
  supp_fig <- plot_grid(
    three_figs_across_section, 
    ncol = 1, 
    label_size = 26
    # rel_heights = c(1, 0.5, 1, 0.3),
    # rel_widths = c(1, 0.8, 1, 0.8)
  )
  ggsave(file.path(path, "main_fig.png"), plot = supp_fig, width = (plot_x_val *2), height = (plot_y_val *6))
  
  cat(green("Made supp fig. \n"))
  return(supp_fig)
}


strip_legends<- function(list_of_subfigs, plots_to_strip_legend ){
  list_of_subfigs <- lapply(names(list_of_subfigs), function(name) {
    p <- list_of_subfigs[[name]]
    if (name %in% plots_to_strip_legend) {
      p <- p + theme(legend.position = "none")
    }
    return(p)
  })
  return(list_of_subfigs)
}




# need legends list for the names of each of the colored boxes.
make_main_legend <- function(legend_text) {
  kp_color_palette<- set_color_pallette()
  kp_color_palette_selected <- kp_color_palette[legend_text]
  
  dummy_df <- data.frame(name_for_color = names(kp_color_palette_selected), value = 1)
  dummy_df$name_for_color <- factor(dummy_df$name_for_color, levels = legend_text)  # preserve order
  
  dummy_plot <- ggplot(dummy_df, aes(x = name_for_color, y = value, fill = name_for_color)) +
    geom_bar(stat = "identity") +  # <- this matters!
    scale_fill_manual(values = kp_color_palette_selected, name = NULL) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 21),
      legend.key.size = unit(2.05, "lines")
    ) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
  legend_for_fig <- cowplot::get_legend(dummy_plot)
  plot(legend_for_fig)
  return(legend_for_fig)
}























conglomerate_figs_one_cohort<- function(plots, legends){
  legends_list<- legends
  lapply(legends_list, function(legend) plot(legend))
  kp_color_palette<- set_color_pallette()
  mutation_plots <- lapply(mutations_plots_list, function(p) p + theme(legend.position = "none"))
  plots<- c(mutation_plots, exp_plots) 
  lapply(exp_plots, function(exp) plot(exp))
  plots_together <- plot_grid(plotlist = plots, ncol = 2, label_size = 26, labels = "AUTO")
  
  great_plot_with_right_padding <- ggdraw() + draw_plot(plots_together, x = 0, y = 0, width = 0.85, height = 1) 
  
  new_great_legend <- make_legends_from_pallete_and_legend_list(legends_list, kp_color_palette)
  
  great_plot_with_great_legend <- ggdraw() + draw_plot(great_plot_with_right_padding, 0, 0, 1, 1) + #  top, right, bottom, left
    draw_plot(new_great_legend, x = 0.83, y = 0.81, width = 0.2, height = 0.2)
  
  ggsave(file.path(path, "great_plot3.png"), plot = great_plot_with_great_legend, width = 17, height = 22)
  
}














conglomerate_figs_one_cohortv2<- function(plots, legends){
  legends_list<- legends
  lapply(legends_list, function(legend) plot(legend))
  kp_color_palette<- set_color_pallette()
  mutation_plots <- lapply(mutations_plots_list, function(p) p + theme(legend.position = "none"))
  plots<- c(mutation_plots, exp_plots) 
  
  violins<- plots[c(4,6)]
  violins <- lapply(violins, function(p) p + theme(plot.margin = unit(c(1, 1, 7, 1), "lines")))
  violins_together <- plot_grid(plotlist = violins, labels = c(" ", "F"),label_size = 26, ncol = 2); violins_together
  
  
  not_violins<- plots[c(1,2,3,5)]
  

  plots_again <- c(not_violins, list(violins_together))
  
  plots_together <- plot_grid(plotlist = plots_again, ncol = 2, label_size = 26, labels = "AUTO")
  
  great_plot_with_right_padding <- ggdraw() + draw_plot(plots_together, x = 0, y = 0, width = 0.85, height = 1) 
  
  new_great_legend <- make_legends_from_pallete_and_legend_list(legends_list, kp_color_palette)
  
  great_plot_with_great_legend <- ggdraw() + draw_plot(great_plot_with_right_padding, 0, 0, 1, 1) + #  top, right, bottom, left
    draw_plot(new_great_legend, x = 0.83, y = 0.81, width = 0.2, height = 0.2)
  
  ggsave(file.path(path, "great_plot4.pdf"), plot = great_plot_with_great_legend, width = 17, height = 22)
  
}


supp_conglomerate_figs_one_cohortv2<- function(plots, legends){
  names<- c("CD14", "FLT3", "CD123", "CD117")
  plots<- c(exp_plots) 
  #violin_with_title<- plot_grid(violin_plot, title, ncol = 2, rel_heights = c(1, 0.1), rel_widths= c(1, 0.08)); violin_with_title
  
 # mutations_plots_list <- lapply(plots, function(p) p + theme(plot.title = element_blank()))
  
  plots_together <- plot_grid(plotlist = plots, ncol = 2, label_size = 26, labels = "AUTO")
  
  great_plot_with_right_padding <- ggdraw() + draw_plot(plots_together, x = 0, y = 0, width = 0.85, height = 1) 
  
  ggsave(file.path(path, "supp.png"), plot = great_plot_with_right_padding, width = 17, height = 22)
  
}











# --------------------
legends_list<- legends
kp_color_palette<- set_color_pallette()

mutation_plots <- lapply(mutations_plots_list, function(p) p + theme(legend.position = "none"))

fav_cds<- c("FLT3", "CD11b", "CD14", "CD34")
which_plots_to_use<- match(fav_cds, cds_of_interest)
selected_expression_plots<- exp_plots[which_plots_to_use]


plots<- c(mutation_plots, selected_expression_plots) 
plots_together <- plot_grid(plotlist = plots, ncol = 3, label_size = 26, labels = "AUTO")

great_plot_with_right_padding <- ggdraw() + draw_plot(plots_together, x = 0, y = 0, width = 0.85, height = 1) 

new_great_legend <- make_legends_from_pallete_and_legend_list(legends_list, kp_color_palette)

great_plot_with_great_legend <- ggdraw() + draw_plot(great_plot_with_right_padding, 0, 0, 1, 1) + #  top, right, bottom, left
  draw_plot(new_great_legend, x = 0.83, y = 0.81, width = 0.2, height = 0.2)

ggsave(file.path(path, "great_plot2.png"), plot = great_plot_with_great_legend, width = 16, height = 28)









#------------mutations_plots_list---------
legends_list<- legends
kp_color_palette<- set_color_pallette()

mutations_plots_list <- lapply(mutations_plots_list, function(p) p + theme(legend.position = "none"))
mutations_plots <- plot_grid(plotlist = mutations_plots_list, ncol = 3, label_size = 26, labels = "AUTO")

new_great_legendmutations_plots <- small_make_legends_from_pallete_and_legend_list(legends_list, kp_color_palette)
great_plot_with_right_padding_mutations_plots <- ggdraw() +
  draw_plot(mutations_plots, x = 0, y = 0,
            width = 0.85, height = 1) # shrinks width to leave space on the right (15%)

mutations_plots_great_plot_with_great_legend <- ggdraw() +
  draw_plot(great_plot_with_right_padding_mutations_plots, 0, 0, 1, 1) + #  top, right, bottom, left
  draw_plot(new_great_legendmutations_plots, x = 0.89, y = 0.75, width = 0.01, height = 0.01)

#reat_plot_and_legend<- plot_grid(exp_plots_together, new_great_legend, ncol = 2, rel_widths = c(1, 0.1), labels = NULL)
ggsave(file.path(path, "great_plot3.png"), plot = mutations_plots_great_plot_with_great_legend, width = 21, height = 6)





make_cowplot_for_just_braf<- function(exp_plots, mutations_plots_list, legends_list, kp_color_palette){
  

  #selected_plots <- lapply(exp_plots, function(p) {p + theme(plot.title = element_text(size = 18))})
  #great_plot_list_nolegend <- lapply(selected_plots, function(p) p + theme(legend.position = "none"))
  
  great_plot <- plot_grid(plotlist = great_plot_list_nolegend, ncol = 2, label_size = 26, labels = "AUTO")
  
  great_plot_with_right_padding <- ggdraw() +
    draw_plot(great_plot, x = 0, y = 0,
              width = 0.85, height = 1) # shrinks width to leave space on the right (15%)
  
  new_great_legend <- make_legends_from_pallete_and_legend_list(legends_list, kp_color_palette)
  great_plot_with_great_legend <- ggdraw() +
    draw_plot(great_plot_with_right_padding, 0, 0, 1, 1) + #  top, right, bottom, left
    draw_plot(new_great_legend, x = 0.82, y = 0.81, width = 0.2, height = 0.2)
  
  ggsave(file.path(path, "great_plot.png"), plot = great_plot_with_great_legend, width = 22, height = 30)
  #ggsave(file.path(path, "supplementary.png"), plot = great_plot_with_great_legend, width = 22, height = 10)
}

small_make_legends_from_pallete_and_legend_list<- function(legends_list, kp_color_palette){
  legends_to_display<- get_legends_from_list(legends_list)
  kp_color_palette_df<- as.data.frame(kp_color_palette)
  palette_order <- rownames(kp_color_palette_df)
  legends_to_display <- palette_order[palette_order %in% legends_to_display];legends_to_display
  
  kp_color_palette_df_filtered <- kp_color_palette_df[rownames(kp_color_palette_df) %in% legends_to_display, , drop = FALSE]
  named_palette <- setNames(kp_color_palette_df_filtered$kp_color_palette,
                            rownames(kp_color_palette_df_filtered))
  
  dummy_df <- data.frame(name_for_color = rownames(kp_color_palette_df_filtered), value = 1)
  dummy_df$name_for_color <- factor(dummy_df$name_for_color, levels = legends_to_display) #preserve the order of the legend too
  dummy_plot <- ggplot(dummy_df, aes(x = value, fill = name_for_color)) +
    geom_bar() +
    scale_fill_manual(values = named_palette, name = NULL) +
    theme(legend.position = "right")+
    theme(
      legend.text = element_text(size = 16),
      legend.key.size = unit(1, "lines")  #  color boxes
    );dummy_plot
  
  new_great_legend <- get_legend(dummy_plot)
  return(new_great_legend)
}





# need legends list for the names of each of the colored boxes.
make_legends_from_pallete_and_legend_list<- function(legends_list, kp_color_palette){
  legends_to_display<- get_legends_from_list(legends_list)
  kp_color_palette_df<- as.data.frame(kp_color_palette)
  palette_order <- rownames(kp_color_palette_df)
  legends_to_display <- palette_order[palette_order %in% legends_to_display];legends_to_display
  
  kp_color_palette_df_filtered <- kp_color_palette_df[rownames(kp_color_palette_df) %in% legends_to_display, , drop = FALSE]
  named_palette <- setNames(kp_color_palette_df_filtered$kp_color_palette,
                            rownames(kp_color_palette_df_filtered))
  
  dummy_df <- data.frame(name_for_color = rownames(kp_color_palette_df_filtered), value = 1)
  dummy_df$name_for_color <- factor(dummy_df$name_for_color, levels = legends_to_display) #preserve the order of the legend too
  dummy_plot <- ggplot(dummy_df, aes(x = value, fill = name_for_color)) +
    geom_bar() +
    scale_fill_manual(values = named_palette, name = NULL) +
    theme(legend.position = "right")+
    theme(
      legend.text = element_text(size = 21),
      legend.key.size = unit(2.05, "lines")  #  color boxes
    );dummy_plot
  
  new_great_legend <- get_legend(dummy_plot)
  return(new_great_legend)
}


get_legends_from_list<- function(legends_list){
  get_labels_from_legend <- function(legend) {
    grobs <- legend[["grobs"]][[1]][["grobs"]]
    labels <- sapply(grobs, function(g) {
      children <- g[["children"]]
      text_grobs <- children[grep("GRID.text", names(children))]
      if (length(text_grobs) > 0) {
        return(text_grobs[[1]][["label"]])
      } else {return(NA)}
    })
    labels[!is.na(labels)]
  }
  all_labels <- lapply(legends_list, get_labels_from_legend)
  unique_labels <- unique(unlist(all_labels))
  cat(blue("Found", length(all_labels), " labels:", paste(all_labels, collapse = ", ") ,"\n"))
  cat(green("Found", length(unique_labels), "unique labels:", paste(unique_labels, collapse = ", ") ,"\n"))
  return(unique_labels)
  # combined_legend <- do.call(cowplot::plot_grid, c(legends_list, ncol = 1)); combined_plot
}





