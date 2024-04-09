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


# Null hypothesis
#library(MASS)
library(twosamples)

# Reading data
library(arrow)

# Spatial point pattern
library(spatstat)

# Plots
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(cmocean)
library(scales)
library(ggtext)


theme_set(theme_minimal())

## Subsampling ----
#--------------------------------------------------------------------------#
sub_sample <- TRUE # whether to subsample or not
n_img <- 10000 # number of images to consider for subsampling
n_dist <- 10000 # number of distances to retain for comparisons


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
n_cores <- 2

## Seed ----
#--------------------------------------------------------------------------#
seed <- 1
set.seed(seed)

## Coast data ----
#--------------------------------------------------------------------------#
coast <- read_csv(file.path(data_dir, "raw/coast.csv"), show_col_types = FALSE)


## Functions ----
#--------------------------------------------------------------------------#
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
