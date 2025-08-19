## Packages ----
#--------------------------------------------------------------------------#
library(tidyverse)
library(here)
library(reshape2)
library(parallel)
library(pbmcapply)
library(hms)
library(ecotaxarapi)
library(unix)
library(castr)
library(furrr)
library(vegan)
library(MLmetrics)


# Null hypothesis
#library(MASS)
library(twosamples)
library(quantreg)
library(broom)

# Reading data
library(arrow)

# Spatial point pattern
library(spatstat)
library(bayestestR)

# Plots
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(cmocean)
library(scales)
library(ggtext)
library(paletteer)
library(khroma)
library(chroma)
library(ggrepel)
library(geomtextpath)
library(ggpattern)



theme_set(theme_minimal())

## Subsampling ----
#--------------------------------------------------------------------------#
sub_sample <- FALSE # whether to subsample or not
n_img <- 10000 # number of images to consider for subsampling
n_dist <- 10000 # number of distances to retain for comparisons


## Distance threshold ----
#--------------------------------------------------------------------------#
# Minimum number of distances to include a group/pair in analyses
n_dist_min <- 10000 

# Distance threshold
# NB: computed at step 3
dist_thr <- 11 # in cm


## Quantiles probabilities ----
#--------------------------------------------------------------------------#
# Prepare probabilities for 10000-quantiles
probs <- seq(0, 1, length.out = 10000)


## Image volume ----
#--------------------------------------------------------------------------#
# image volume in pixels
vol <- c()
vol$x <- 10240
vol$y <- 2048
vol$z <- 9572

## Directories ----
#--------------------------------------------------------------------------#
data_dir <- here("data")

## Parallel ----
#--------------------------------------------------------------------------#
n_cores <- 8

## Seed ----
#--------------------------------------------------------------------------#
seed <- 1
set.seed(seed)

## Coast data ----
#--------------------------------------------------------------------------#
coast <- read_csv(file.path(data_dir, "raw/coast.csv"), show_col_types = FALSE)


## Functions ----
#--------------------------------------------------------------------------#
# Floor and ceil with given precision
floor_dec <- function(x, level = 1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level = 1) round(x + 5*10^(-level-1), level)


# Function to increment file numbers in batch
increment_files <- function(to_increment){
  # List files with this number
  to_rename <- list.files(pattern = paste0(to_increment, ".*"))
  
  # Compute new number and generate new file names
  new_nb <- as.character(as.numeric(to_increment) + 1) %>% str_pad(width = 2, pad = "0")
  renamed <- str_replace_all(to_rename, pattern = to_increment, replacement = new_nb)
  
  # Rename files
  file.rename(from = to_rename, to = renamed)
}


#' Compute distances between individuals of two taxonomic groups within an ensemble of images
#' 
#' In each image, the distance between all unique pairs of points is computed.
#'
#' @param points dataframe with columns `x` and `y` for point position, `img_name` for image reference and `taxon` for taxonomic groups
#' @param n_cores number of cores to use for parallel processing
#' @param z_dim whether to consider 3rd dimension (`z`)
#'
#' @return a dataframe with all distances for each image
#' 
compute_all_dist_inter <- function(points, n_cores = 12, z_dim = FALSE) {
  dist_all <- pbmclapply(unique(points$img_name), function(name) {
    
    # Get points within image, split dataframe by taxon
    if (z_dim){ # 3D case
      points_img <-  points %>% 
        filter(img_name == name) %>% 
        dplyr::select(taxon, x, y, z) %>% 
        group_by(taxon) %>% 
        group_split(.keep = F)
    } else { # 2D case
      points_img <- points %>% 
        filter(img_name == name) %>% 
        dplyr::select(taxon, x, y) %>% 
        group_by(taxon) %>% 
        group_split(.keep = F)
    }
    
    # Use proxy::dist to compute distances between two matrices
      dist_all <- proxy::dist(points_img[[1]], points_img[[2]])
      dist_all <- `dim<-`(c(dist_all), dim(dist_all))
      dist_all <- dist_all %>% 
        melt(varnames = c("p1", "p2")) %>% 
        as_tibble() %>% 
        rename(dist = value) %>% 
        mutate(img_name = name)
    
    return(dist_all)
   }, mc.cores = n_cores, ignore.interactive = TRUE)
  
  # Combine results in a dataframe
  dist_all <- do.call(bind_rows, dist_all)

  # Check that we have the number of expected distances within each image
  # Number of distances in each image
  #dist_ok <- dist_all %>% 
  #  count(img_name, name = "n_dist") %>% 
  #  # Join with theoretical number of distances for each image
  #  left_join(
  #    points %>% 
  #      count(img_name, taxon) %>% 
  #      group_by(img_name) %>% 
  #      summarise(n_dist_th = prod(n)),
  #    by = join_by(img_name)
  #  ) %>% 
  #  mutate(ok = n_dist == n_dist_th)
    
  #if(!all(dist_ok$ok)){stop("Number of computed distances differs from what is expected based on the number of points and their type.")} # should return TRUE
  
  return(dist_all)
}



#' Compute distances between all points within an ensemble of images
#' 
#' In each image, the distance between all unique pairs of points is computed.
#'
#' @param points dataframe with columns `x` and `y` for point position and `img_name` for image reference
#' @param n_cores number of cores to use for parallel processing
#' @param z_dim whether to consider 3rd dimension (`z`)
#'
#' @return a dataframe with all distances for each image
#' 
compute_all_dist <- function(points, n_cores = 12, z_dim = FALSE) {
  dist_all <- pbmclapply(unique(points$img_name), function(name) {
    
    # Get points within image
    if (z_dim){ # 3D case
      points_img <-  points %>% 
        filter(img_name == name) %>% 
        dplyr::select(x, y, z) %>% 
        as.matrix()
    } else { # 2D case
      points_img <-  points %>% 
        filter(img_name == name) %>% 
        dplyr::select(x, y) %>% 
        as.matrix()
    }
    
    # Compute distances between points
    melt(as.matrix(dist(points_img)), varnames = c("p1", "p2")) %>% 
      as_tibble() %>% 
      filter(p1 != p2) %>% 
      rename(dist = value) %>% 
      mutate(
        img_name = name,
        # Keep only one of distances compute between A and B (A to B and B to A were computed)
        pair = ifelse(p1 < p2, paste(p1, p2), paste(p2, p1))
      ) %>% 
      distinct(pair, dist, .keep_all = TRUE) %>% 
      dplyr::select(-pair)
  }, mc.cores = n_cores, ignore.interactive = TRUE
  ) %>% 
    bind_rows()
  
  # Check that we have the number of expected distances within each image
  # For a set of n points, the number of unique distances is n(n-1)/2
  #dist_ok <- left_join(
  #  points %>% count(img_name, name = "n_obj"),
  #  dist_all %>% count(img_name, name = "n_dist"),
  #  by = join_by(img_name)
  #) %>% 
  #  mutate(ok = n_dist == (n_obj * (n_obj-1))/2)
  
  #if(!all(dist_ok$ok)){stop("Number of computed distances differs from what is expected based on the number of points.")} # should return TRUE
  
  return(dist_all)
}


#' Get Correlation Between Two Taxa
#'
#' This function computes the correlation between the abundances of two taxa, with options to choose
#' the correlation method (Pearson or Spearman) and log-transform the abundances.
#'
#' @param t1 A character string representing the first taxon.
#' @param t2 A character string representing the second taxon.
#' @param abundance_mat A matrix where rows represent different taxa and columns represent samples. The values in the matrix correspond to the abundances of each taxon in each sample.
#' @param taxa_indices A named list where names correspond to taxon names and values are the row indices of those taxa in `abundance_mat`.
#' @param method A character string specifying the correlation method to use. Options are `"pearson"` (default) or `"spearman"`.
#' @param log_transform A logical value indicating whether to log-transform the abundances before calculating correlation. Default is `FALSE`. The log transformation applied is `log(abundance + 1)` to avoid issues with zeros.
#'
#' @return A numeric value representing the correlation between the two taxa, or `NA` if the taxa have zero variance or are identical.
#'
#' @details
#' The function first extracts the abundances for the specified taxa (`t1` and `t2`) from `abundance_mat`. If `log_transform = TRUE`, the abundances are log-transformed using the natural logarithm (with a pseudocount of 1). The function then checks for zero variance in either taxon. If either taxon has zero variance, it returns `NA`. If the taxa are different, the function computes the correlation using the specified method (Pearson or Spearman). If `t1` and `t2` are the same taxon, the function returns `NA`.
#'
#' @examples
#' # Example usage:
#' abundance_matrix <- matrix(abs(rnorm(100)), nrow = 10)  # Example abundance matrix (10 taxa, 10 samples)
#' taxa_idx <- list("taxon1" = 1, "taxon2" = 2)  # Example taxa indices
#' get_corr("taxon1", "taxon2", abundance_matrix, taxa_idx, method = "pearson", log_transform = TRUE)
#'
get_corr <- function(t1, t2, abundance_mat, taxa_indices, method = "spearman", log_transform = FALSE) {
  idx1 <- taxa_indices[[t1]]
  idx2 <- taxa_indices[[t2]]
  ab_t1 <- abundance_mat[idx1, ]
  ab_t2 <- abundance_mat[idx2, ]
  
  # Log-transform the abundances if specified
  if (log_transform) {
    ab_t1 <- log1p(ab_t1)  # log1p to avoid log(0) issues
    ab_t2 <- log1p(ab_t2)
  }
  
  # Check for zero variance
  if (sd(ab_t1) == 0 || sd(ab_t2) == 0) {
    return(NA)  # Return NA if variance is zero
  }
  
  # Compute correlation if taxa are different
  if (t1 != t2) {
    cor(ab_t1, ab_t2, method = method)
  } else {
    NA
  }
}

#' Get Correlation and P-Value Between Two Taxa
#'
#' This function computes the correlation between the abundances of two taxa and returns both
#' the correlation value and the p-value. It provides options to choose the correlation method
#' (Pearson or Spearman) and log-transform the abundances.
#'
#' @param t1 A character string representing the first taxon.
#' @param t2 A character string representing the second taxon.
#' @param abundance_mat A matrix where rows represent different taxa and columns represent samples. The values in the matrix correspond to the abundances of each taxon in each sample.
#' @param taxa_indices A named list where names correspond to taxon names and values are the row indices of those taxa in `abundance_mat`.
#' @param method A character string specifying the correlation method to use. Options are `"pearson"` or `"spearman"` (default) .
#' @param log_transform A logical value indicating whether to log-transform the abundances before calculating correlation. Default is `FALSE`. The log transformation applied is `log(abundance + 1)` to avoid issues with zeros.
#'
#' @return A list containing:
#'   - `corr`: A numeric value representing the correlation between the two taxa, or `NA` if the taxa have zero variance or are identical.
#'   - `p_val`: The associated p-value of the correlation test, or `NA` if variance is zero.
#'
get_corr_pval <- function(t1, t2, abundance_mat, taxa_indices, method = "spearman", log_transform = FALSE) {
  idx1 <- taxa_indices[[t1]]
  idx2 <- taxa_indices[[t2]]
  ab_t1 <- abundance_mat[idx1, ]
  ab_t2 <- abundance_mat[idx2, ]
  
  # Log-transform the abundances if specified
  if (log_transform) {
    ab_t1 <- log1p(ab_t1)  # log1p to avoid log(0) issues
    ab_t2 <- log1p(ab_t2)
  }
  
  # Check for zero variance
  if (sd(ab_t1) == 0 || sd(ab_t2) == 0) {
    return(list(correlation = NA, p_value = NA))  # Return NA if variance is zero
  }
  
  # Compute correlation if taxa are different
  if (t1 != t2) {
    cor_test <- cor.test(ab_t1, ab_t2, method = method, exact = FALSE)
    return(list(corr = cor_test$estimate, p_val = cor_test$p.value))
  } else {
    return(list(corr = NA, p_val = NA))
  }
}


# Function to read vignette
create_image_plot <- function(img_path) {
  tryCatch({
    # Read image
    img <- imager::load.image(img_path)
    # Drop 29 pixels at the bottom to remove scale bar
    img <- imager::imsub(img, y < imager::height(img) - 29)
    # Create a plot object instead of displaying directly
    p <- ggplot() +
      annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
      theme_void() +
      theme(plot.margin = unit(c(0,0,0,0), "mm"))
    return(p)
  }, error = function(e) {
    # Return empty plot if image fails
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Failed", size = 2) +
      theme_void() +
      theme(plot.margin = unit(c(0,0,0,0), "mm"))
  })
}



# Function to create image mosaic and return both plot and sampled data
create_image_mosaic <- function(object_data, mosaic_size = 10) {
  n_images <- mosaic_size^2
  
  # Sample images
  if (nrow(object_data) < n_images) {
    warning(paste("Only", nrow(object_data), "images available"))
    sampled_objects <- object_data
    n_missing <- n_images - nrow(object_data)
  } else {
    sampled_objects <- object_data |> slice_sample(n = n_images)
    n_missing <- 0
  }
  
  # Create list of plots
  plot_list <- list()
  
  # Add actual image plots
  for (i in 1:nrow(sampled_objects)) {
    plot_list[[i]] <- create_image_plot(sampled_objects$path_to_img[i])
  }
  
  # Add empty plots if needed
  if (n_missing > 0) {
    for (i in (nrow(sampled_objects) + 1):n_images) {
      plot_list[[i]] <- ggplot() + theme_void()
    }
  }
  
  # Arrange in grid
  mosaic_plot <- do.call(arrangeGrob, c(plot_list, list(ncol = mosaic_size)))
  
  # Convert to ggplot object for patchwork compatibility
  mosaic_ggplot <- ggplot() +
    annotation_custom(mosaic_plot) +
    theme_void() +
    coord_fixed()
  
  # Return both plot and data
  return(mosaic_ggplot)
}
