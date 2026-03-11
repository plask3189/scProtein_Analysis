

# new_colors<- adjust_hue(og_colors)
# point_data$colour <- new_colors
# p_new <- ggplot(point_data, aes(x = x, y = y, color = colour)) +
#   geom_point(size = 0.0005, alpha = 0.7) +
#   labs(x = "umap_1", y = "umap_2", title = cluster_of_differentiation) +
#   theme_classic() +
#   scale_color_identity() +
#   theme(plot.title = element_text(hjust = 0.5)); p_new
# 
prep_data<- function(p){
  p <- FeaturePlot(s, features = cluster_of_differentiation) 
  plot_build <- ggplot_build(p[[1]])
  og_point_data<- plot_build[["data"]][[1]]
  og_colors<- og_point_data$colour 
  cat("Num of og colors:", length(og_colors))
  return(og_colors)
}

make_spectrum_legend <- function(p){
  og_colors<- prep_data(p)
  rgb_vals <- col2rgb(og_colors)  
  diff_vals <- apply(rgb_vals, 2, function(x) max(x) - min(x))  # Compute the difference between max and min for each color: to measure how gray it is
  # order colors from smallest difference (most gray) to largest (least gray)
  ordered_colors <- og_colors[order(diff_vals)]
  last_color <- ordered_colors[length(ordered_colors)]
  ordered_colors <- og_colors[order(diff_vals)]
  color_df <- data.frame(
    x = 1,
    y =seq_along(ordered_colors),
    col = ordered_colors
  )
  
  ggplot(color_df, aes(x = x, y = y, fill = col)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    theme(axis.title.x = element_blank())
  
  
}

get_most_intense_color1<- function(p){
  og_colors<- prep_data(p)
  rgb_vals <- col2rgb(og_colors)  
  diff_vals <- apply(rgb_vals, 2, function(x) max(x) - min(x))  # Compute the difference between max and min for each color: to measure how gray it is
  # order colors from smallest difference (most gray) to largest (least gray)
  ordered_colors <- og_colors[order(diff_vals)]
  last_color <- ordered_colors[length(ordered_colors)]
  ordered_colors <- og_colors[order(diff_vals)]
  
  n <- length(ordered_colors)
  last_percent_colors <- ordered_colors[ceiling(n * 0.999):n]
  rgb_last <- col2rgb(last_percent_colors)
  mean_rgb <- rowMeans(rgb_last)
  mean_color <- rgb(mean_rgb[1]/255, mean_rgb[2]/255, mean_rgb[3]/255)
  print(mean_color)
  plot(1, 1, col = mean_color, pch = 16, cex = 5, xlim = c(0, 2), ylim = c(0, 2), main =mean_color)
  return(mean_color)
}


adjust_hue<-function(og_colors_){
  hsv_matrix <- rgb2hsv(col2rgb(og_colors_))
  # Adjust the hue by adding a shift value (here, 0.1) and wrap around using modulo 1
  hsv_matrix["h", ] <- (hsv_matrix["h", ] + 0.1) %% 1
  # Convert the adjusted HSV values back to hex colors
  adjusted_colors <- hsv(hsv_matrix["h", ], hsv_matrix["s", ], hsv_matrix["v", ])
  cat("Num of new colors:", length(adjusted_colors))
  
  par(mfrow = c(2, 1))
  og <- barplot(rep(1, length(og_colors)),
                col = og_colors,
                border = NA,
                space = 0,
                main = "OG Color Spectrum",
                axes = FALSE)
  new <- barplot(rep(1, length(adjusted_colors)),
                 col = adjusted_colors,
                 border = NA,
                 space = 0,
                 main = "New Color Spectrum",
                 axes = FALSE)
  
  return(adjusted_colors)
}


modify_point_data<- function(point_data){
  library(colorspace)
  idx <- og_point_data$indent != 1
  # Get the colors that need to be transformed
  cols_to_transform <- og_point_data$colour[idx]
  rgb_mat <- col2rgb(cols_to_transform)
  
  # Set a tolerance for detecting gray (colors with similar R, G, B values)
  tol <- 60
  max_diff <- pmax(abs(rgb_mat[1,] - rgb_mat[2,]),
                   abs(rgb_mat[1,] - rgb_mat[3,]),
                   abs(rgb_mat[2,] - rgb_mat[3,]))
  is_gray <- max_diff < tol
  
  # Convert these colors to an RGB object (note: pass values positionally)
  rgb_obj <- RGB(rgb_mat[1,] / 255, rgb_mat[2,] / 255, rgb_mat[3,] / 255)
  
  # Convert the RGB object to HCL using polarLUV
  hcl_obj <- as(rgb_obj, "polarLUV")
  
  # Extract the H, C, and L values
  h_vals <- hcl_obj@coords[, "H"]
  c_vals <- hcl_obj@coords[, "C"]
  l_vals <- hcl_obj@coords[, "L"]
  
  # For non-gray colors, force the hue to 0 (red)
  h_vals[!is_gray] <- 0
  
  # Reconstruct the new HCL colors and convert back to hex codes
  new_hcl <- polarLUV(H = h_vals, C = c_vals, L = l_vals)
  transformed_colors <- hex(new_hcl)
  
  # Create a new full color vector, replacing only the transformed rows
  new_colors <- og_point_data$colour
  new_colors[idx] <- transformed_colors
  
  # Optionally, update your data frame with the new colors:
  og_point_data$colour <- new_colors
  
  print(new_colors)
  return(new_colors)
}
