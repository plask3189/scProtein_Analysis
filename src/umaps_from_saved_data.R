library(crayon)
library(RColorBrewer)
library(scales) 
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(gridExtra)
library(cowplot)
library(writexl)
library(ggpubr)
library(SingleCellExperiment)
source("R/protein_kp/violin.R")
source("R/protein_kp/visualization_helper_functions.R")
source("R/protein_kp/label_seurat_utilities.R")
source("R/protein_kp/significance_brackets.R")
just_BRAF<- FALSE
all_plots_dir = "all_plots"
#all_plots_dir = "supplimentary"
if(just_BRAF == TRUE){
  #🔵🔵🔵🔵🔵🔵🔵🔵🔵
  sce_subset  <- readRDS("just_braf_cohort_sce_subset_n7_v2.rds")
  sample_names<- unique(sub("_.*", "", colnames(sce_subset)))
  cat("Samples:", sample_names, "\n")
  cat("With variants:, ", rownames(sce_subset))
  comparing_cohorts <- FALSE
  protein_normalization_method <- "DSB_norm"
  path = "just_braf"
  seurat_save_name <- "just_braf_cohort_sce_subset_n7_v2_SEURAT.rds"
  s<- readRDS(seurat_save_name)
  if (grepl("supp", all_plots_dir)){
    cds_of_interest<- c("FLT3", "CD14", "CD123") # the supplementary markers
  } else {
    cds_of_interest<- c("CD11b", "CD34")
  }
  

  most_intense_colors<- get_most_intense_color(s, sample_cd = cds_of_interest[[1]])
} else{
  #🔵🔴🔵🔴🔵🔴🔵🔴🔵🔴
  #sce_subset <- readRDS("two_cohorts/two_cohorts_sce.rds")
  sce_subset<- readRDS(file.path("two_cohorts", "two_cohorts_sce_subset_v4.rds"))
  sample_names <- unique(sub("_[^_]*$", "", colnames(sce_subset))); cat("Samples:", sample_names, "\n"); cat("# of samples:", length(sample_names), "\n")
  comparing_cohorts <- TRUE
  protein_normalization_method <- "DSB_norm"
  path <- "two_cohorts"
  s<- readRDS("two_cohorts/two_cohorts_SEURAT.rds")
  s@meta.data[["Group"]][s@meta.data[["Group"]] == 1] <- "BRAF Cohort"
  s@meta.data[["Group"]][s@meta.data[["Group"]] == 2] <- "RAS Cohort"
  if (grepl("supp", all_plots_dir)){
    cds_of_interest<- c( "CD14", "CD123") # the supplementary markers
  } else {
    cds_of_interest<- c("CD11b", "CD34")
  }
  color_pallett <- set_color_pallette()
  most_intense_colors_g1<- color_pallett[["BRAF Cohort"]] # get_most_intense_color(subset(s, subset = Group == "BRAF Cohort"), sample_cd = cds_of_interest[[1]], color_scheme = c("lightgrey", "blue"))
  most_intense_colors_g2<-color_pallett[["RAS Cohort"]]  #get_most_intense_color(subset(s, subset = Group == "RAS Cohort"), sample_cd = cds_of_interest[[1]], color_scheme =  c("lightgrey", "red"))
}


#cds_of_interest <- rownames(sce_subset@int_colData@listData[["altExps"]]@listData[["Protein"]]@se@assays@data@listData[["DSB_norm"]]);cds_of_interest 
plot_protein_from_saved_data<- function(sce_subset,
                                        s,
                                        comparing_cohorts, 
                                        protein_normalization_method = "DSB_norm", 
                                        cds_of_interest = c("CD11b"),
                                        path = "protein_plots"){
  if (!dir.exists(path)) { dir.create(path, recursive = TRUE)}; full_path <- normalizePath(path, winslash = "/", mustWork = FALSE); cat(blue("saving plots to:",full_path))
  set.seed(1333)
  
  # -------------------- 🎨🖌️ Make UMAPs  🎨🖌-----------------------------
  umap_resolution<- 400 

  sample_names <- unique(sub("_[^_]*$", "", colnames(sce_subset))); cat("# of samples:", length(sample_names), "\n")
  text2<- paste(paste("----Samples---- \n", paste0(sample_names, collapse = "\n"), "\n"),"\n" ,"# of Samples:", length(sample_names), "\n")
  cohort_info<- ggplot() + annotate("text", x = 0.5, y = 0.5, label = text2, size = 3) + xlim(0, 1) + ylim(0, 1) +theme_void(); cohort_info
  ggsave(file.path(path, paste("cohort_info.png")), plot = cohort_info, width = 8, height = 6, dpi = umap_resolution)
  
  if (just_BRAF == TRUE){
    violin_variant_set<- c("NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469", "BRAF.D594")
    s<- label_seurat_with_mutations(sce_subset, s, original_desired_variants =violin_variant_set)
    as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
    mutations_plots_list<- list()
    legends<- list()
    three_mutaneers_lab <- label_seurat_with_mutations(sce_subset, s, original_desired_variants = c("BRAF.G469", "BRAF.D594", "BRAF.V600E"))
    three_mutaneers<- construct_map(three_mutaneers_lab, title = " "); three_mutaneers
    legends[[length(legends)+1]]<- get_legend(three_mutaneers)
    three_mutaneers<- adjust_alpha(three_mutaneers, c("BRAF.G469", "BRAF.D594", "BRAF.V600E"),
                                   this_group = 0.7, the_rest = 0.7, pt_size = 1)+ NoLegend();three_mutaneers
    mutations_plots_list[[length(mutations_plots_list)+1]]<- (three_mutaneers + NoLegend()+ theme(axis.title = element_text(size =11)))
    ggsave(file.path(path, paste("three_mutaneers.png")), plot = three_mutaneers, width = 8, height = 6, dpi = umap_resolution)
    saveRDS(three_mutaneers, file.path(all_plots_dir,paste0("three_mutaneers.rds")))
    
    braf_vs_ras_lab <- label_seurat_with_mutations(sce_subset, s, original_desired_variants = c("BRAF", "NRAS", "KRAS", "PTPN11"))
    braf_vs_ras<- construct_map(braf_vs_ras_lab, title = " ")+theme(plot.title = element_text(size = 15)) 
    legends[[length(legends)+1]]<- get_legend(braf_vs_ras)
    braf_vs_ras1<- adjust_alpha(braf_vs_ras, c("KRAS", "NRAS", "PTPN11"),
                                this_group = 1, the_rest = 0.15)+ NoLegend()+ theme(axis.title = element_text(size =11))
    print(braf_vs_ras1)
    ggplot_build(braf_vs_ras1)$data
    mutations_plots_list[[length(mutations_plots_list)+1]]<- braf_vs_ras1
    ggsave(file.path(path, paste("just_braf_cohort_braf_vs_ras.png")), plot = braf_vs_ras1, width = 8, height = 6, dpi = umap_resolution)
    saveRDS(braf_vs_ras1, file.path(all_plots_dir,paste0("just_braf_cohort_braf_vs_ras.rds")))
  } else {
   # muts_order <- c("WT", "Non-RAS Mutant WT", "RAS", "BRAF.V600E", "BRAF.G469A", "BRAF.D594G")
    mutations_plots_list<- list()
    legends<- list()
    color_palette<- set_color_pallette()
    braf_vs_ras_lab <- label_seurat_with_mutations(sce_subset, s, original_desired_variants = c("BRAF", "RAS"))
    both_cohorts_ras_vs_braf <- DimPlot(braf_vs_ras_lab, reduction = "umap", group.by = "mutation_status", pt.size = 0.4, 
                                        alpha=0.6, raster = FALSE) +
      scale_color_manual(values = color_palette) + labs(title = " ")+
      theme(legend.text = element_text(size = 10),plot.title = element_text(size = 26)+ theme(axis.title = element_text(size =11))); both_cohorts_ras_vs_braf
    ggsave(file.path(path, paste("both_cohorts_ras_vs_braf_together.png")), plot = both_cohorts_ras_vs_braf, width = 6, height = 6, dpi = umap_resolution)
     # ------------------------------------- seurat object label edits------------------------------------- 
    s<- relabel_seurat_for_two_cohort_run(sce, s)
    saveRDS(both_cohorts_ras_vs_braf, file.path(all_plots_dir,paste0("both_cohorts_ras_vs_braf_mutations_labeled_umap.rds")))
    
  }
  
  #🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊
  for(cluster_of_differentiation in cds_of_interest){ 
    cat(blue("Generating plots for", cluster_of_differentiation, "\n"))
    plot_name <- paste0(cluster_of_differentiation, "_umap"); file_plot_name <- paste0(plot_name, ".png")
    if (protein_normalization_method == "DSB_norm"){ # we make negative DSB values 0 bc they are background. https://www.nature.com/articles/s41467-022-29356-8
       s@meta.data[[cluster_of_differentiation]] <- s@assays$Protein@data[cluster_of_differentiation, ]
      s@meta.data[[cluster_of_differentiation]][s@meta.data[[cluster_of_differentiation]] < 0] <- 0
     # range(sce_subset@int_colData@listData[["altExps"]]@listData[["Protein"]]@se@assays@data@listData[["DSB_norm"]])
    }
    title <- ggdraw() + draw_label(cluster_of_differentiation, fontface = 'bold', size = 25, angle = 270, x=0, hjust = 0.5, vjust = 0.5)
    if(comparing_cohorts == TRUE){
      p1 <- FeaturePlot(s[, s$Group == "BRAF Cohort"], features = cluster_of_differentiation, pt.size = 0.4, label.size = 5.5, raster = FALSE, cols = c("lightgrey", "#548ede")) + theme(axis.title = element_text(size =11))#+ ggtitle(" ") ; p1 #        ggtitle(paste("BRAF Cohort", cluster_of_differentiation)) 
      p2 <- FeaturePlot(s[, s$Group == "RAS Cohort"], features = cluster_of_differentiation, pt.size = 0.4, label.size = 5.5, raster = FALSE, cols = c("lightgrey", "#9E352D"))+ theme(axis.title = element_text(size =11)) #+ ggtitle(" ") ;p2
      saveRDS(p1, file.path(all_plots_dir,paste0("BRAF_", cluster_of_differentiation, "_exp_umap.rds")))
      saveRDS(p2, file.path(all_plots_dir,paste0("RAS_", cluster_of_differentiation, "_exp_umap.rds")))
      
      # make the p values between groups
      p_val_table_plot<- get_p_vals_table(s, cluster_of_differentiation,
                                          mutations_to_compare= c("RAS", "BRAF.V600E", "BRAF.G469A", "BRAF.D594G"),
                                          compare_groups = TRUE); p_val_table_plot

      ggsave(file.path(path, paste0(cluster_of_differentiation,"_p_values.png")), plot = p_val_table_plot, width = 9.5, height = 5, dpi = umap_resolution)

      main_violin_plot<- violin(s, cluster_of_differentiation, compare_cohort_violins = comparing_cohorts) + 
        theme(axis.title.x = element_blank()); main_violin_plot
                            #mutations_to_compare = c("WT", violin_variant_set2)) 
      # add p value brackets:
      # main_violin_plot<- make_significance_brackets(s, mutations_to_compare = c("BRAF.V600E", "BRAF.G469", "BRAF.D594"), main_violin_plot, compare_cohorts= TRUE,
      #                                                     cluster_of_differentiation); print(main_violin_plot)

      max_height_of_violin<- max(main_violin_plot[[1]][["data"]][[cluster_of_differentiation]]);max_height_of_violin
      violin_plot_with_space<- main_violin_plot + NoLegend() + 
        coord_cartesian(ylim = c(0, (max_height_of_violin+4)),clip = "off")+ theme(axis.title = element_text(size =11));violin_plot_with_space
      
      saveRDS(violin_plot_with_space, file.path(all_plots_dir,paste0("both_cohorts_", cluster_of_differentiation, "violins.rds")))
      
      violin_plot_no_legend<- main_violin_plot + NoLegend() + 
        theme(axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 17),
              plot.margin = unit(c(8, 1, 2.77, 1), "lines")) + ggtitle(" ")+
        coord_cartesian(ylim = c(0, (max_height_of_violin+1)),clip = "off");violin_plot_no_legend
      
      main_violin_plot<- main_violin_plot+ theme(legend.position = "top", plot.title =element_blank())
      p<- p1 +p2+ main_violin_plot+ plot_layout(ncol = 3) 
      p2_with_title<- plot_grid(p, title, ncol = 2, rel_heights = c(1, 0.1), rel_widths= c(1, 0.08)); p2_with_title
      ggsave(file.path(path, paste0(cluster_of_differentiation,"_umap.png")), plot = p2_with_title, width = 18, height = 6, dpi = umap_resolution)
      
      } else { # comparing_cohorts == FALSE ... JUST BRAF:
      p <- FeaturePlot(s, features = cluster_of_differentiation,label.size = 5.5,  pt.size = 0.4, cols = c("lightgrey", "#548ede"))+
        ggtitle(cluster_of_differentiation)  + theme(axis.title = element_text(size =11))
      saveRDS(p, file.path(all_plots_dir,paste0("just_braf_", cluster_of_differentiation, "_exp_umap.rds")))
      
      p_val_table_plot<- get_p_vals_table(s, cluster_of_differentiation, mutations_to_compare= c("BRAF.G469", "BRAF.D594", "BRAF.V600E"), compare_groups = FALSE)
      ggsave(file.path(path, paste0(cluster_of_differentiation,"_p_values.png")), plot = p_val_table_plot, width = 12, height = 8, dpi = umap_resolution)
      
      violin_plot<- violin(s, cluster_of_differentiation, compare_cohort_violins = comparing_cohorts, 
                            mutations_to_compare = to_display_p_vals_on) + theme(axis.title = element_text(size =11))
      max_height_of_violin<- max(violin_plot[[1]][["data"]][[cluster_of_differentiation]]);max_height_of_violin
      violin_plot<- violin_plot +coord_cartesian(ylim = c(0, (max_height_of_violin+7)),clip = "off") # 5/6
      saveRDS(violin_plot, file.path(all_plots_dir,paste0("just_braf_", cluster_of_differentiation, "violins.rds")))
      
      max_height_of_violin <- max(violin_plot[[1]][["data"]][[cluster_of_differentiation]], na.rm = TRUE) + 2
      og_max_y<- max(violin_plot[[1]][["data"]][[cluster_of_differentiation]])
      violin_plot<- violin_plot + NoLegend() + 
        ggtitle(cluster_of_differentiation) +
        theme(#axis.text.x = element_blank(),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 17)#, plot.margin = unit(c(3.1, 1, 2.77, 1), "lines")
              ) #+ggtitle(" ") #+
        #coord_cartesian(ylim = c(0, (og_max_y+1)),clip = "off");violin_plot # to trim y axis
      violin_with_title<- plot_grid(violin_plot, title, ncol = 2, rel_heights = c(1, 0.1), rel_widths= c(1, 0.08)); violin_with_title
      
      p<- p + ggtitle(cluster_of_differentiation)
      violin_to_save<- violin_plot + ggtitle(" ") + theme(axis.title.y = element_text(size = 8),axis.text.x  = element_blank(), axis.title.x = element_blank())+
        coord_cartesian(ylim = c(0, (max_height_of_violin+5)),clip = "off")
      for_individual_exp_umap_saving <- plot_grid(p , violin_to_save, ncol = 2); for_individual_exp_umap_saving
      ggsave(file.path(path, paste0(cluster_of_differentiation,"_umap.png")), plot = for_individual_exp_umap_saving, width = 10, height = 4, dpi = umap_resolution)
    }
  }

}





relabel_seurat_for_two_cohort_run<- function(sce, s){
  sce <- readRDS("two_cohorts/two_cohorts_sce_subset_v4.rds")
  s<- readRDS("two_cohorts/two_cohorts_SEURAT.rds")
  all_variants<- rownames(sce_subset); all_variants
  s<- all_mutations_plus_no_overwrite(sce_subset, s, original_desired_variants = all_variants, general=TRUE)
  mutation_status_counts <- as.data.frame(s@meta.data) %>% dplyr::count(mutation_status, sort = TRUE); mutation_status_counts; sum(mutation_status_counts$n)
  
  s<- exclude_other_variants_from_seurat_label(s, variants_to_keep = c("NRAS", "KRAS", "PTPN11", "BRAF.V600E", "BRAF.G469A", "BRAF.D594G"))
  s<- make_gene_dominant(s) # braf
  mutation_status_counts <- as.data.frame(s@meta.data) %>% dplyr::count(mutation_status, sort = TRUE); mutation_status_counts; sum(mutation_status_counts$n)
  #s<- relabel_cells_with_multiple_variants(s)
  s<- rename_label_if_contains(s, to_replace = c("NRAS", "PTPN11", "KRAS"), new_label = "RAS")
  
  s@meta.data[["Group"]][s@meta.data[["Group"]] == 1] <- "BRAF Cohort"
  s@meta.data[["Group"]][s@meta.data[["Group"]] == 2] <- "RAS Cohort"
  
  mutation_status_counts_by_group <- as.data.frame(s@meta.data) %>%
    dplyr::count(Group, mutation_status, sort = TRUE)
  cell_count<- sum(mutation_status_counts_by_group$n)
  print(mutation_status_counts_by_group)
  cat(green("Our",cell_count,"cells are labeled as any of:",paste0(unique(s@meta.data[["mutation_status"]]), collapse = ", ")), "\n")
  return(s)
}


