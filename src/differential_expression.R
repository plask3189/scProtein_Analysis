# FOR FUN -----------------------
sce_subset <- readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/two_cohorts/two_cohorts_sce.rds")
sample_names <- unique(sub("_[^_]*$", "", colnames(sce_subset))); cat("Samples:", sample_names, "\n"); cat("# of samples:", length(sample_names), "\n")
comparing_cohorts <- TRUE
protein_normalization_method <- "DSB_norm"
path <- "two_cohorts"
s<- readRDS("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/two_cohorts/two_cohorts_SEURAT.rds")
s@meta.data[["Group"]][s@meta.data[["Group"]] == 1] <- "BRAF Cohort"
s@meta.data[["Group"]][s@meta.data[["Group"]] == 2] <- "RAS Cohort"

variants_in_sce_clean <- sub("\\..*", "", rownames(sce_subset))
s<- label_seurat_with_mutations(sce_subset, s, original_desired_variants =variants_in_sce_clean)
as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
if (protein_normalization_method == "DSB_norm"){ 
  s@meta.data[[cluster_of_differentiation]] <- s@assays$Protein@data[cluster_of_differentiation, ]
  s@meta.data[[cluster_of_differentiation]][s@meta.data[[cluster_of_differentiation]] < 0] <- 0 # https://www.nature.com/articles/s41467-022-29356-8
}
# remove wt cells
mutations_of_interest <- c("BRAF", "NRAS", "KRAS", "PTPN11")
pattern <- paste(mutations_of_interest, collapse = "|")
seurat_obj_subsetted <- subset(s, subset = mutation_status != "WT" & grepl(pattern, mutation_status))
as.data.frame(seurat_obj_subsetted@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)

# Find DE features between groups
Idents(seurat_obj_subsetted) <-seurat_obj_subsetted@meta.data[["Group"]]
markers <- FindMarkers(seurat_obj_subsetted, ident.1 = "BRAF Cohort", ident.2 = "RAS Cohort")
markers <- markers %>% arrange(p_val_adj) ; head(markers, 10) 
top_cds <- head(rownames(markers), 30)  # select top 20 genes
two_cohorts_heatmap<- DoHeatmap(seurat_obj_subsetted, features = top_cds, group.by = c("Group", "mutation_status")) + NoLegend()



ggsave(file.path(path, paste("two_cohorts_heatmap.png")), plot = two_cohorts_heatmap, width = 8, height = 6, dpi = umap_resolution)


