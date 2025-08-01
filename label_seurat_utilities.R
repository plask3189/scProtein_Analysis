


# give me some variants. the cells with entirely not these are put into "non-BRAF non-Ras Mutant WT". IF a cell does have
# one of the variants to exclude but another variant that is to keep, just delete the variant to exclude
exclude_other_variants_from_seurat_label<- function(s, variants_to_keep){
  # HANDLE NOT ACTUALLY WT CELLS. remove the variant not in the violins set. 
  pattern <- paste(c(variants_to_keep,"WT"), collapse = "|")
  no_violin_variant <- !grepl(pattern, s@meta.data[["mutation_status"]]); 
  variants_to_remove<- unique(s@meta.data[no_violin_variant, "mutation_status"])
  cat(blue("Relabeling variants:", paste0(variants_to_remove, collapse=",  "), " \n")); cat(magenta("as 'Non-RAS Mutant WT'\n"))
  s@meta.data[no_violin_variant, "mutation_status"] <- "Non-RAS Mutant WT"# "non-BRAF non-Ras Mutant WT"
  as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE)
  # remove any non-violin plot variants. 
  printed_messages <- new.env()  # Create an empty environment to store printed messages
  s@meta.data$mutation_status <- sapply(s@meta.data$mutation_status, function(x) {
    key <- x  # Create a key that summarizes this specific case
    if (x %in% c("WT", "Non-RAS Mutant WT")) { 
      if (!exists(key, envir = printed_messages)) {
        cat("OG mut:", x, " -> WT or Non-RAS Mutant WT. Keep.\n")
        assign(key, TRUE, envir = printed_messages)
      }
      return(x)
    }
    
    variants <- strsplit(x, "\\+")[[1]]
    kept_variants <- variants[variants %in% variants_to_keep]
    if (length(kept_variants) == 0) { 
      if (!exists(key, envir = printed_messages)) {
        cat("OG mut:", x, " -> No variants kept. Assigning Non-RAS Mutant WT.\n")
        assign(key, TRUE, envir = printed_messages)
      }
      return("Non-RAS Mutant WT")
    } else {
      final_status <- paste(kept_variants, collapse = "+")
      if (!exists(key, envir = printed_messages)) {
        cat("OG mut:", x)
        if (final_status != x) {
          #cat(blue(" -> Variants kept:", paste(kept_variants, collapse = ", "), "\n"))
          cat(blue(" -> New mutation_status assigned:", final_status, "\n"))
        } else {
          cat(blue(" -> No change to mutation_status.\n"))
        }
        assign(key, TRUE, envir = printed_messages)
      }
      return(final_status)
    }
  })
  return(s)
}


# Make braf dominant
make_gene_dominant<- function(s){
  s@meta.data$mutation_status <- sapply(s@meta.data$mutation_status, function(x) {
    if (grepl("BRAF\\.", x)) {# If it contains "BRAF."
      match <- regmatches(x, regexpr("BRAF\\.[^+]*", x))# Extract and keep the full BRAF.something
      return(match)
    } else {# Otherwise keep the original
      return(x)
    }
  })
  print(as.data.frame(s@meta.data) %>%dplyr::count(mutation_status, sort = TRUE))
  return(s)
}

rename_label_if_contains <- function(s, to_replace = c("NRAS", "PTPN11", "KRAS"), new_label = "RAS") {
  pattern <- paste(to_replace, collapse = "|") # Create regex pattern: "NRAS|PTPN11|KRAS"
  
  s@meta.data$mutation_status <- sapply(s@meta.data$mutation_status, function(x) {
    if (grepl(pattern, x)) {  # Check if any of the to_replace patterns are found
      return(new_label)       # If yes, replace with new_label
    } else {
      return(x)               # Otherwise keep the original
    }
  })
  s@meta.data$mutation_status <- as.character(s@meta.data$mutation_status) # ensure vector stays character
  print(as.data.frame(s@meta.data) %>% dplyr::count(mutation_status, sort = TRUE))# preview counts
  return(s)
}

relabel_cells_with_multiple_variants<- function(s){
  # Relabel cells with multiple variants to just the first variant
  s@meta.data$mutation_status <- sapply(s@meta.data$mutation_status, function(x) {
    strsplit(x, "\\+")[[1]][1]  # Split by "+" and take the first part
  })
  
  mutation_status_counts <- as.data.frame(s@meta.data) %>% 
    dplyr::count(mutation_status, sort = TRUE)
  
  print(mutation_status_counts)
  cat("Total cells:", sum(mutation_status_counts$n), "\n")
  print(as.data.frame(s@meta.data) %>% dplyr::count(mutation_status, sort = TRUE)) # preview counts
  return(s)
}
