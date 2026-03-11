library(scales)
source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/R/protein_kp/visualization_helper_functions.R", echo=TRUE)
library(crayon)
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(stringr)
library(cowplot)
# initialize data ----------------------
basic_seurat<- initialize_seurat_object_for_heatmap(seurat = readRDS("two_cohorts/two_cohorts_SEURAT.rds"),
                                                      sce =readRDS("two_cohorts/two_cohorts_sce_subset_v4.rds"))
mutation_status_counts <- as.data.frame(basic_seurat@meta.data) %>% dplyr::count(mutation_status, sort = TRUE); mutation_status_counts; sum(mutation_status_counts$n)
mutation_status_counts_by_group <- as.data.frame(basic_seurat@meta.data) %>%
  dplyr::count(Group, mutation_status, sort = TRUE)
as.data.frame(basic_seurat@meta.data) %>%dplyr::count(Group, mutation_status, sort = TRUE)

# differential expression ----------------
table(Idents(basic_seurat))# view the current idents
markers <- FindMarkers(basic_seurat, ident.1 = "BRAF Cohort", ident.2 = "RAS Cohort") %>% arrange(p_val_adj)
top_cds <- head(rownames(markers), 10)
top_cds <- top_cds[!grepl("Ig", top_cds)]
top_cds<- c(top_cds, "CD14", "CD123", "CD11b", "CD34")

# top_cds<- c("CD3","CD19","CD14","CD33","CD11b",# lin- (e.g., CD3 for T cells, CD19 for B cells, CD14 for monocytes).
#             "CD34","CD117") # Lin +

#top_cds<-  c("CD14", "CD123", "CD11b", "CD34")

# make heatmap ----------------
basic_heatmap_plot<- basic_heatmap(basic_seurat, top_cds)

# ---------------- The heatmap separated by original ident group ------------------
#color_palette_kp <- set_color_pallette()
identity_levels <- unique(basic_heatmap_plot[["data"]]$Identity); identity_levels
#identity_colors <- hue_pal(h = c(0, 180))(length(identity_levels))
manual_colors <- c("#4a84c2" , "#AC2A26")
identity_colors <- setNames(manual_colors[seq_along(identity_levels)], identity_levels)
names(identity_colors) <- identity_levels
identity_color_vector <- identity_colors[basic_heatmap_plot[["data"]]$Identity]

pbuild <- ggplot_build(basic_heatmap_plot) #  get plot layout
x.min <- min(pbuild$layout$panel_params[[1]]$x.range)
x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
y.range <- diff(pbuild$layout$panel_params[[1]]$y.range)
y.top <- max(pbuild$layout$panel_params[[1]]$y.range)

y_of_bar_bottom<- y.top #(y.top + 0.06 * y.range)-1.5; y_of_bar_bottom
bar_girth <- 1
y_of_bar_top = y_of_bar_bottom+bar_girth; y_of_bar_top

labeled_mutations <- unique(unlist(strsplit(unique(basic_seurat@meta.data[["mutation_status"]]), split = "\\+"))) ; labeled_mutations
space_for_the_new_mutation_bars<- ((1+bar_girth)*length(labeled_mutations)); space_for_the_new_mutation_bars # 1 is the separations space. 

basic_seurat@meta.data[["sample_id"]]  <- substr(basic_seurat@meta.data[["Cell"]], 1, nchar(basic_seurat@meta.data[["Cell"]]) - 19)
sample_id_color_vector<- make_sample_id_bar(basic_seurat, basic_heatmap_plot)

basic_heatmap_plot_with_base_split <- basic_heatmap_plot +
  coord_cartesian(ylim = c(0, y_of_bar_top+11), clip = "off") +
  annotation_raster(raster = t(identity_color_vector),
                    xmin = x.min, xmax = x.max, ymin =  y_of_bar_bottom, ymax =  y_of_bar_top)#+
  # annotation_raster(raster = t(sample_id_color_vector),
  #                   xmin = x.min, xmax = x.max, ymin =  y_of_bar_top, ymax =  y_of_bar_top+1)

ggsave(file.path("two_cohorts", "plot1.png"), plot = basic_heatmap_plot_with_base_split, width = 8, height = 6, dpi = 100)

# -------------------------------------------------------------

mutations_to_show<- c("TET2","TP53","BRAF.V600E", "BRAF.G469A", "BRAF.D594G", "RUNX1", "NPM1", "NRAS","FLT3","PTPN11", "KRAS")
mutation_status_bar_data <- list()
plot<- basic_heatmap_plot
for (variant in mutations_to_show) {
  mutation_color_vector <- make_mutation_subset_bar_data(basic_seurat, plot, variant)
  mutation_status_bar_data[[variant]] <- mutation_color_vector
}


dynamic_y_max <- y.top#+1  # initialize
heatmap_plot_with_variants <- basic_heatmap_plot_with_base_split  # init
mutation_bars_list <- list()

for (variant in rev(mutations_to_show)) {
  cat(blue("Variant:", variant, "---------------\n"))
  tip_of_bar_y <- dynamic_y_max + 1
  mutation_bar <- annotation_raster(
    raster = t(mutation_status_bar_data[[variant]]),
    xmin = x.min, xmax = x.max,
    ymin = tip_of_bar_y, ymax = tip_of_bar_y + 1
  )
  
  mutation_bars_list[[variant]] <- mutation_bar
  label_layer <- geom_text(
    data = data.frame(x = x.min - 0.8, y = tip_of_bar_y + 0.36, label = variant),
    aes(x = x, y = y, label = label),
    hjust = 1.2, size = 3
  )
  heatmap_plot_with_variants <- heatmap_plot_with_variants + mutation_bar + label_layer
  dynamic_y_max <- tip_of_bar_y
}

ggsave(file.path("two_cohorts", "plot2.png"), plot = heatmap_plot_with_variants, width = 12, height = 6, dpi = 300)



legend_text_order <- c("BRAF Cohort", "RAS Cohort")
legend_to_use <- make_heatmap_legend(legend_text = legend_text_order) 
heatmap_plot_with_variants<- heatmap_plot_with_variants + guides(colour = "none") + # remove the old group guide legend
  theme(plot.margin = margin(t = 0, r = 40, b = 0, l = 0, unit = "pt")) 

heatmap_plot_with_variants
final_plot <- ggdraw() +
  draw_plot(heatmap_plot_with_variants, x = 0, y = 0, width = 1, height = 1)  +
  draw_plot(legend_to_use, x = 0.92, y = 0.25, width = 0.05, height = 0.12)
#final_plot
ggsave(file.path("two_cohorts", "plot3.png"), plot = final_plot, width = 15, height = 6, dpi = 300)




# need legends list for the names of each of the colored boxes.
make_heatmap_legend <- function(legend_text) {
  kp_color_palette<- set_color_pallette()
  kp_color_palette_selected <- kp_color_palette[legend_text]
  
  dummy_df <- data.frame(
    name_for_color = factor(legend_text, levels = legend_text), 
    value = 1)
  dummy_df$name_for_color <- factor(dummy_df$name_for_color, 
                                    levels = legend_text)  # preserve order
  
  dummy_plot <- ggplot(dummy_df, aes(x = name_for_color, y = value, fill = name_for_color)) +
    geom_bar(stat = "identity") +  
    scale_fill_manual(values = kp_color_palette_selected, name = NULL) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.size = unit(1.2, "lines"),
      plot.margin = margin(0, 0, 0, 0)
    )
  legend_for_fig <- cowplot::get_legend(dummy_plot)
  plot(legend_for_fig)
  return(legend_for_fig)
}


make_sample_id_bar<- function(basic_seurat, plot){
  sample_id_df <- data.frame(Cell = rownames(basic_seurat@meta.data),
                            sample_id = basic_seurat@meta.data[["sample_id"]],
                            stringsAsFactors = FALSE)
  print(table(sample_id_df$sample_id))
  plot_sample_id_data <- left_join(plot[["data"]], sample_id_df, by = "Cell") %>% 
    dplyr::select(c("Cell", "Identity", "sample_id")) %>%  distinct()
  
  sample_id_levels <- unique(plot_sample_id_data$sample_id)
  n_samples <- length(sample_id_levels)
  my_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_samples)
  sample_id_colors <- setNames(my_colors, sample_id_levels)
  names(sample_id_colors) <- sample_id_levels
  sample_id_color_vector <- sample_id_colors[plot_sample_id_data$sample_id]
  return(sample_id_color_vector)
}





make_mutation_subset_bar_data<- function(basic_seurat, plot, variant){
  
  color_palette_kp <- set_color_pallette()
  color_palette_kp[["WT"]] <- "#f5f5f5"
  mutation_df <- data.frame(Cell = rownames(basic_seurat@meta.data),
                            mutation_status = basic_seurat@meta.data[["mutation_status"]],
                            stringsAsFactors = FALSE)
  mutation_df_subset<- mutation_df
  cat(blue("Variant:", variant, "---------------\n"))
  mutation_df_subset$mutation_status <- ifelse(
    grepl(variant, mutation_df$mutation_status), 
    variant, 
    "WT")
  print(table(mutation_df_subset$mutation_status))
  plot_with_mut_data <- left_join(plot[["data"]], mutation_df_subset, by = "Cell") %>% 
    dplyr::select(c("Cell", "Identity", "mutation_status")) %>%  distinct()

  mutation_levels <- unique(plot_with_mut_data$mutation_status)
  mutation_colors <- color_palette_kp[mutation_levels]
  names(mutation_colors) <- mutation_levels
  mutation_color_vector <- mutation_colors[plot_with_mut_data$mutation_status]
  return(mutation_color_vector)
}




basic_heatmap<- function (basic_seurat, top_cds){
  object<- basic_seurat
  features<- top_cds
  group.by = "Group"
  cells <- NULL
  group.colors <- NULL
  disp.min <- -2.5
  disp.max <- NULL
  slot <- "scale.data"
  assay <- "Protein"
  label <- TRUE
  hjust <- 0
  vjust <- 0
  raster <- TRUE
  
  second_grouping<- "mutation_status"
  
  # ----------------------- get features 
  assay <- assay %||% DefaultAssay(object = object); assay
  DefaultAssay(object = object) <- assay
  cells <- cells %||% colnames(x = object[[assay]])
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
  } ; cat(features)
  # -----------------------
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, slot = slot)[features, cells, drop = FALSE])))
  object <- suppressMessages(expr = StashIdent(object = object, save.name = "ident"))
  Idents(object) <-object@meta.data[["ident"]]
  group.by <- "ident"
  groups.use <- object[[group.by]][cells, , drop = FALSE]; as.data.frame(object@meta.data) %>%dplyr::count(ident, sort = TRUE)

  plots <- vector(mode = "list", length = ncol(x = groups.use)) # we just want one plot so fine.
  i<-1
  data.group <- data
  group.use <- groups.use[, i, drop = TRUE] ; cat("groups use:", unique(group.use), "\n")
  if (!is.factor(x = group.use)) {
    group.use <- factor(x = group.use)
  }
  names(x = group.use) <- cells ; cat("cell names:", head(names(x = group.use)), "\n")
  lgroup <- length(levels(group.use)) # number of groups
  
  # ********** THE BASIC HEATMAP PLOT: ***********************
  basic_heatmap_plot <- SingleRasterMap(data = data.group, raster = raster, feature.order = features,
                          cell.order = names(x = sort(x = group.use)), group.by = group.use) 
  
  return(basic_heatmap_plot)
}



initialize_seurat_object_for_heatmap<- function(seurat, sce){ 
  cat(blue("Using given seurat object and sce subset \n"))
  if(! "BRAF Cohort" %in% unique(seurat@meta.data[["Group"]]) ){
    seurat@meta.data[["Group"]][seurat@meta.data[["Group"]] == 1] <- "BRAF Cohort"
    seurat@meta.data[["Group"]][seurat@meta.data[["Group"]] == 2] <- "RAS Cohort"
  }


  for(cluster_of_differentiation in rownames(seurat@assays$Protein@data)){ # for all the markers: 
    seurat@meta.data[[cluster_of_differentiation]] <- seurat@assays$Protein@data[cluster_of_differentiation, ] # makes the marker a column in the seurat object
    seurat@meta.data[[cluster_of_differentiation]][seurat@meta.data[[cluster_of_differentiation]] < 0] <- 0 # https://www.nature.com/articles/s41467-022-29356-8
  }

  variants_in_sce_clean <- unique(sub("\\..*", "", rownames(sce)))
  seurat<- all_mutations_plus_no_overwrite(sce, seurat, original_desired_variants =variants_in_sce_clean,
                                           include_individual_braf_muts = TRUE)
  #seurat<- all_mutations_plus_no_overwrite(sce, seurat, original_desired_variants =variants_in_sce_clean,
                                           #include_individual_braf_muts = TRUE)
  Idents(seurat) <-seurat@meta.data[["Group"]]
  
  cat(blue("Labeled the cells with variants:", paste0(variants_in_sce_clean, collapse = ", "), "\n"))
  
  return(seurat)
}


add_mutation_legend_to_heatmap<- function(mutation_levels){
  mutation_legend_plot <- ggplot(data.frame(x = mutation_levels, y = 1), aes(x = x, y = y, color = x)) +
    geom_point(size = 5) +
    scale_color_manual(values = mutation_colors, name = "Mutation Status") +
    theme_void() +
    theme(legend.position = "right")
  return(mutation_legend_plot)
}





