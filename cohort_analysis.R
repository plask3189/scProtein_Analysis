library(crayon)
library(scDNA)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
source("R/protein_kp/braf_protein_helper_functions.R")
source("R/variant_ID.R")
source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/R/clonograph.R")
library(ggplot2)


sample_files_group1 <-c("BRAF/M7456braf.dna+protein.h5",
                        "BRAF/M1912braf.dna+protein.h5", 
                        "BRAF/M0292braf.dna+protein.h5",
                        "BRAF/2459_Braf.dna+protein.h5", #"BRAF/M0626braf.dna+protein.h5", "BRAF/6232_braf.dna+protein.h5",
                        "BRAF/4629_braf.dna+protein.h5",
                        "BRAF/A5330braf.dna+protein.h5",
                        "BRAF/A0634braf.dna+protein.h5")
sample_files_group2 = c("./data/A0290.dna+protein.h5",
                        "./data/A1107.dna+protein.h5",
                        "./data/A0259.dna+protein.h5",
                        "./data/Sample17020.dna+protein.h5",
                        "./data/LAM_1962.dna+protein.h5",
                        "./data/LAM_3307.dna+protein.h5")
sample_files<-c(sample_files_group1,sample_files_group2)
variants_to_select<- c("NRAS", "KRAS", "PTPN11")
variants_to_select<- c("TP53", "BRAF", "NRAS", "PTPN11", "KRAS", "NPM1", "TET2", "FLT3","IDH1","IDH2", "DNMT3A", "RUNX1", "ASXL1", "FLT3", "EZH2")

# 🚧🚧🚧🚧🚧🚧🚧🚧 big run 🚧🚧🚧🚧🚧🚧🚧🚧🚧
variant_output_list<- variant_output(sample_files_group1, variants_to_select)
saveRDS(variant_output_list, "b1_variants_output_list.rds")
# 🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧🚧

variant_output_list<- readRDS("b1_variants_output_list.rds")
print(names(variant_output_list))
# variant_output_list to show the types of BRAF Mutations
variant_output_list2 <- variant_output_list

variant_output_list_top<- identify_top_variants_to_examine(variant_output_list)


variant_output_list3 <- lapply(variant_output_list, function(df) { # top for each sample
  df[1:min(5, nrow(df)), ]
})

variant_data <- mapply(function(df1, df2) {
  merged_df <- rbind(df1, df2)
  merged_df <- unique(merged_df)  # remove duplicate rows
  return(merged_df)
}, variant_output_list_top, variant_output_list3, SIMPLIFY = FALSE)

all_stats<- make_stats(variant_data)
all_stats_summarized <- all_stats %>%
  group_by(sample_id, final_annot) %>%
  summarize(VAF = max(VAF, na.rm = TRUE), .groups = "drop")


vaf_heat_mat <- all_stats_summarized %>%
  pivot_wider(
    names_from = sample_id,
    values_from = VAF,
    values_fill = list(VAF = 0)
  ) %>%
  column_to_rownames("final_annot") %>%
  as.matrix()


ordered_samples <- all_stats %>%distinct(sample_id, group_id) %>%arrange(group_id, sample_id) %>%pull(sample_id)
all_stats$sample_id <- factor(all_stats$sample_id, levels = ordered_samples)

vaf_ggplot <- ggplot() + geom_tile(data = all_stats[all_stats$group_id == "1", ],
            aes(x = sample_id, y = final_annot, fill = VAF), color = "white") +
  scale_fill_gradient(low = "white", high = "blue", name = "VAFs for BRAF Cohort") +
  new_scale_fill() +  
  geom_tile(data = all_stats[all_stats$group_id == "2", ],
            aes(x = sample_id, y = final_annot, fill = VAF), color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "VAFs for RAS Cohort") +
  labs(x = "Sample ID", y = "Mutation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle("VAFs")
vaf_ggplot
# -----------------------------

all_stats_summarized <- all_stats %>%
  group_by(sample_id, final_annot) %>%
  summarize(genotyping_rate = mean(genotyping_rate, na.rm = TRUE), .groups = "drop")

genotyping_rate_heat_mat <- all_stats_summarized %>%
  pivot_wider(
    names_from = sample_id,
    values_from = genotyping_rate,
    values_fill = list(genotyping_rate = 0)  # fix: make values_fill a named list
  ) %>%
  column_to_rownames("final_annot") %>%
  as.matrix()

genotyping_rate_ggplot <- ggplot() + geom_tile(data = all_stats[all_stats$group_id == "1", ],
                                   aes(x = sample_id, y = final_annot, fill = genotyping_rate), color = "white") +
  scale_fill_gradient(low = "white", high = "blue", name = "genotyping_rate for BRAF Cohort") +
  new_scale_fill() +  
  geom_tile(data = all_stats[all_stats$group_id == "2", ],
            aes(x = sample_id, y = final_annot, fill = genotyping_rate), color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "genotyping_rate for RAS Cohort") +
  labs(x = "Sample ID", y = "Mutation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle("Genotyping Rates")
genotyping_rate_ggplot

combined_plot <- vaf_ggplot / genotyping_rate_ggplot
ggsave("data.png", combined_plot, width = 8, height =16)

valid_variants<- identify_valid_variants(genotyping_rate_heat_mat, vaf_heat_mat, vaf_thresh=5, gt_thresh = 93)
valid_variants


identify_valid_variants<- function(genotyping_rate_heat_mat,vaf_heat_mat, vaf_thresh=5, gt_thresh = 93){
  head(genotyping_rate_heat_mat)
  head(vaf_heat_mat)
  rows_with_high_vaf <- rownames(vaf_heat_mat)[apply(vaf_heat_mat, 1, function(x) any(x > vaf_thresh))]
  rows_with_high_gt <- rownames(genotyping_rate_heat_mat)[apply(genotyping_rate_heat_mat, 1, function(x) any(x > gt_thresh))]
  rows_with_both <- intersect(rows_with_high_vaf, rows_with_high_gt)
  vaf_heat_mat[rows_with_both, ]
  genotyping_rate_heat_mat[rows_with_both, ]
  return(rows_with_both)
}




#were there BRAF mutations that were not: V600, G469, or D594 in the samples
variant_output<- function(sample_files, variants_to_select){
  variant_output_list <- vector("list", length(sample_files))
  names(variant_output_list) <- sample_files
  for(sample_file in sample_files){
    cat(blue("----------", sample_file, "--------------------\n"))
    variant_output<-variant_ID(file=sample_file,panel="MSK_RL",GT_cutoff=90, VAF_cutoff=0.01)
    variant_output <- variant_output %>%
      distinct()%>%
      dplyr::filter(!CONSEQUENCE == "synonymous" ) %>%
      dplyr::filter(SYMBOL %in% c("NRAS", "KRAS", "PTPN11")) %>%
      arrange(desc(VAF))%>% 
      dplyr::slice(1:30)
    variant_output_list[[sample_file]]<- variant_output
  }
  return(variant_output_list)
}

make_stats<- function(variant_output_list2){
  all_stats <- do.call(rbind, lapply(seq_along(variant_output_list2), function(i) {
    df <- variant_output_list2[[i]] %>%
      dplyr::select(final_annot, VAF, genotyping_rate) %>%
      mutate(
        file =  names(variant_output_list2)[i],
        sample_id =  names(variant_output_list2)[i],
        num_variants = nrow(variant_output_list2[[i]])
      )
    return(df)
  }))
  all_stats$sample_id <- sub("\\..*", "", basename(as.character(all_stats$sample_id)))
  all_stats$group_id <- ifelse(
    all_stats$file %in% sample_files_group1, "1",
    ifelse(all_stats$file %in% sample_files_group2, "2", NA)
  )
  all_stats
  library(ggplot2)
  library(ggnewscale)
  library(reshape2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(patchwork)

  return(all_stats)
}

identify_top_variants_to_examine<- function(variant_output_list){

  all_variants_df <- do.call(rbind, lapply(names(variant_output_list), function(name) {
    df <- variant_output_list[[name]]
    df$sample_id <- name
    return(df)
  }))
  
  vaf_sums <- aggregate(VAF ~ AA_change, data = all_variants_df, sum)
  
  top_vafs <- vaf_sums[order(vaf_sums$VAF, decreasing = TRUE), ]
  print(head(top_vafs))
  variants_to_select <- head(top_vafs$AA_change, 5)
  
  variant_output_list_top <- lapply(variant_output_list2, function(df) {
    df[df$AA_change %in% variants_to_select, ]
  })
  return(variant_output_list_top)
}

