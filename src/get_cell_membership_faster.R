
GetMembershipFaster <- function(droplet_metadata) { # kp added
  
  library(data.table)
  library(mvtnorm)

  setDT(droplet_metadata)  # convert to data.table
  
  # separates Cell and Empty droplet types
  cell_data <- droplet_metadata[Droplet_type == "Cell", .(dna_size, protein_size)]
  empty_data <- droplet_metadata[Droplet_type == "Empty", .(dna_size, protein_size)]
  
  Cell_cov_mat <- cov(cell_data)
  Empty_cov_mat <- cov(empty_data)
  
  cell_mu_vec <- colMeans(cell_data)
  empty_mu_vec <- colMeans(empty_data)
  
  comp1.prod <- dmvnorm(droplet_metadata[, .(dna_size, protein_size)], mean = cell_mu_vec, sigma = Cell_cov_mat)
  comp2.prod <- dmvnorm(droplet_metadata[, .(dna_size, protein_size)], mean = empty_mu_vec, sigma = Empty_cov_mat)
  
  sum_of_comps <- comp1.prod + comp2.prod
  droplet_metadata[, `:=`(
    comp1_prod = comp1.prod,
    comp2_prod = comp2.prod,
    sum_of_comps = sum_of_comps,
    comp1_post = comp1.prod / sum_of_comps,
    comp2_post = comp2.prod / sum_of_comps,
    positive_cell_calls = ifelse(comp1.prod > comp2.prod, "Cell", "Empty")
  )]
  
  return(droplet_metadata)
}
