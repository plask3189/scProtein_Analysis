library(ggalluvial)
library(dplyr)
library(circlize)
library(scales)
library(RColorBrewer)
library(cowplot)
library(png)
library(grid)

source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/R/protein_kp/visualization_helper_functions.R", echo=TRUE)
# basic_seurat<- initialize_seurat_object_for_wheel(seurat = readRDS("two_cohorts/two_cohorts_SEURAT.rds"),
#                                                     sce =readRDS("two_cohorts/two_cohorts_sce_subset_v4.rds"))
basic_seurat_just_braf_cohort<- initialize_seurat_object_for_wheel(seurat = readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/just_braf_cohort_sce_subset_n7_v2_SEURAT.rds"),
                                                  sce =readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/just_braf_cohort_sce_subset_n7_v2.rds"))


mutation_status_counts <- as.data.frame(basic_seurat_just_braf_cohort@meta.data) %>% dplyr::count(mutation_status, sort = TRUE); mutation_status_counts; sum(mutation_status_counts$n)

meta <- basic_seurat_just_braf_cohort@meta.data
filtered_meta <- meta %>% filter(!grepl("BRAF.L79Q", mutation_status))
mutation_lists <- strsplit(filtered_meta$mutation_status, "\\+")
mutation_vector <- unlist(mutation_lists)
mutation_summary <- sort(table(mutation_vector), decreasing = TRUE)
mutation_summary
co_mutation_patterns <- table(sapply(mutation_lists, function(x) paste(sort(x), collapse = "+")))
head(sort(co_mutation_patterns, decreasing = TRUE), 10)

cell_names <- rownames(filtered_meta)

# Get the list of mutation pairs with the corresponding cell name
mutation_pair_cells <- do.call(rbind, lapply(seq_along(mutation_lists), function(i) {
  muts <- mutation_lists[[i]]
  if (length(muts) < 2) return(NULL)
  pairs <- combn(muts, 2, simplify = FALSE)
  # Create data.frame for each pair with the corresponding cell name
  do.call(rbind, lapply(pairs, function(p) {
    data.frame(from = p[1], to = p[2], cell = cell_names[i])
  }))
}))


# Make sure `from < to` for consistent pair naming
mutation_pair_cells <- mutation_pair_cells %>%
  rowwise() %>%
  mutate(pair_from = min(c(from, to)), pair_to = max(c(from, to))) %>%
  ungroup() %>%
  select(from = pair_from, to = pair_to, cell)

cooccurrence_df_with_cells <- mutation_pair_cells %>%
  group_by(from, to) %>%
  summarise(
    weight = n(),
    cells = list(unique(cell)),
    .groups = "drop"
  )




#chordDiagram(cooccurrence_df, transparency = 0.3,directional =1)
#chord<- make_chord(cooccurrence_df_with_cells[,c("from", "to", "weight")])


markers<- basic_seurat_just_braf_cohort@assays[["Protein"]]@counts@Dimnames[[1]]
# Keep only markers that do NOT contain "Ig"
markers <- markers[!grepl("Ig", markers)]
#markers<- c("FLT3", "CD14", "CD123", "CD11b", "CD34")
for(marker in markers){
  column_name<- paste0(marker, "_avg_exp")
  cooccurrence_df_with_cells[[column_name]]<- NA
  basic_seurat_just_braf_cohort@meta.data[[marker]] <- basic_seurat_just_braf_cohort@assays$Protein@data[marker, ]
  basic_seurat_just_braf_cohort@meta.data[[marker]][basic_seurat_just_braf_cohort@meta.data[[marker]] < 0] <- 0

  # get cell ids for each co-occurance
  for(i in seq_len(nrow(cooccurrence_df_with_cells))){ # get avg exp of these cells for this marker
    cells_with_these_muts <- cooccurrence_df_with_cells$cells[[i]]
    avg<- mean(basic_seurat_just_braf_cohort@meta.data[cells_with_these_muts, marker])
    cooccurrence_df_with_cells[[column_name]][i]<- avg
  }
}

cooccurrence_df_scaled <- cooccurrence_df_with_cells %>%
  mutate(across(ends_with("_avg_exp"), ~ rescale(., to = c(0, 1))))

plot_list<- list()
cols<-colnames(cooccurrence_df_scaled[,5:(length(markers)+4)])
for (col in cols) {
  cooccurrence_df_scaled_marker <- cooccurrence_df_scaled[, c("from", "to", "weight", col)]
  make_chord2(cooccurrence_df_scaled_marker, col)
}

png_files <- list.files(path = "chords", pattern = "\\.png$", full.names = TRUE)
grob_list <- lapply(png_files, function(png_path) {
  img <- readPNG(png_path)
  grid::rasterGrob(img, interpolate = TRUE)
})

chords<- plot_grid(plotlist = grob_list, ncol = 4)
ggsave("chords/chords.pdf", chords, width = 3, height = 6.2)







make_chord2 <- function(cooccurrence_df_scaled_marker, col) {
  clean_col <- sub("_avg_exp$", "", col)
  filename <- paste0("chords/chord_plot_", clean_col, ".png")
  
  png(filename, width = 1500, height = 1450, res = 150)
  circos.clear()  

  pair_ids <- paste(cooccurrence_df_scaled_marker$from, cooccurrence_df_scaled_marker$to, sep = "|")
  alphas <- cooccurrence_df_scaled_marker[[col]]
  
  get_pair_color <- function(from, to) {
    key1 <- paste(from, to, sep = "|")
    key2 <- paste(to, from, sep = "|")  # in case direction is reversed
    name1 <- ifelse(from %in% names(color_pallett), color_pallett[[from]], "grey")
    name2 <- ifelse(to %in% names(color_pallett), color_pallett[[to]], "grey")
    # Mix the two colors if needed, or choose one
    return(colorRampPalette(c(name1, name2))(3)[2])  # middle blend
  }
  
  pair_colors <- mapply(get_pair_color,
                        cooccurrence_df_scaled_marker$from,
                        cooccurrence_df_scaled_marker$to)
  colored_links <- alpha(pair_colors, alphas)
  names(colored_links) <- pair_ids
  
  circos.par(canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2))  # Shrink the plot space
  
  chordDiagram(
    x = cooccurrence_df_scaled_marker[, c("from", "to", "weight")],
    grid.col = "grey",
    col = colored_links,
    transparency = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.05)
  )
  
  # Add upright sector labels
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        x = mean(xlim),
        y = ylim[1] + mm_y(2),
        labels = sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.8
      )
    },
    bg.border = NA
  )
  clean_col <- sub("_avg_exp$", "", col)
  title(main = clean_col, cex.main = 3, line = -3.4)
  dev.off()
}


make_chord<- function(cooccurrence_df){
  chordDiagram(
    cooccurrence_df,
    transparency = 0.3,
    annotationTrack = "grid",     # Adds sector labels
    preAllocateTracks = list(track.height = 0.05)
  )
  
  # Make labels upright
  chord<- circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        x = mean(xlim),
        y = ylim[1] + mm_y(2),
        labels = sector.name,
        facing = "clockwise",      # always clockwise
        niceFacing = TRUE,         # makes it upright
        adj = c(0, 0.5),cex = 0.8)},
    bg.border = NA
  )
  return(chord)
}



# df<- basic_seurat@meta.data
# # df is your DE result: columns = gene, log2FC, p_val_adj
# ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = (p_val_adj < 0.05))) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "Volcano Plot: BRAF vs RAS", x = "log2 Fold Change", y = "-log10(adj. p-value)")



initialize_seurat_object_for_wheel<- function(seurat, sce){ 
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
  Idents(seurat) <-seurat@meta.data[["Group"]]
  
  cat(blue("Labeled the cells with variants:", paste0(variants_in_sce_clean, collapse = ", "), "\n"))
  
  return(seurat)
}

