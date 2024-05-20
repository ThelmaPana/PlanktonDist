#/usr/local/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Generate 3 big null datasets and compute distances.
# Date: 17/05/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")

## Load data ----
#--------------------------------------------------------------------------#
set.seed(seed)

message("Reading data")

# Correct image volume for x-axis
load("data/01.corr_factor.Rdata")
vol$x <- vol$x * med_corr

sub_sample <- FALSE
## Subsampling
if (sub_sample){
  load("data/00.images_sub.Rdata")
  load("data/01.x_corrected_plankton_sub.Rdata")
  images <- images_sub
  plankton <- plankton_sub
} else {
  ## All data
  images <- read_parquet("data/00.images_clean.parquet")
  plankton <- read_parquet("data/01.x_corrected_plankton_clean.parquet")
}

# Count objects per image
counts <- plankton %>% count(img_name)


## Generate null datasets ----
#--------------------------------------------------------------------------#
message("Generating null data")
n_img <- nrow(images) # Number of images to generate
n_sets <- 3 # Number of datasets to generate

# Representative number of objets per images
n_pts <- counts %>% slice_sample(n = n_img, replace = TRUE) %>% pull(n)

# Generate sets of random images
rand_points <- lapply(1:n_sets, function(i_set) {
  message(paste("Generating set", i_set))
  # Pick random points within image volumes
  mclapply(1:n_img, function(i){
    # Number of points to sample within image
    n <- n_pts[i]
    # Draw points
    d_points <- tibble(
      x = runif(n = n, min = 1, max = vol$x),
      y = runif(n = n, min = 1, max = vol$y),
      z = runif(n = n, min = 1, max = vol$z)
    ) %>% # Add information for img name
      mutate(img_name = paste0("img_", str_pad(i, nchar(n_img), pad = "0")))
  }, mc.cores = n_cores) %>% 
    bind_rows()
})

# Names of new images
images <- rand_points[[1]] %>% count(img_name)


## Break into chunks ----
#--------------------------------------------------------------------------#
message("Generating chunks")

# Break into n pieces of images to process 
n_chunk <- 10 # number of chunks
chunk_size <- ceiling(nrow(images) / n_chunk) # size of chunks
# compute chunks
images <- images %>% 
  mutate(
    chunk = as.character(floor((row_number() - 1) / chunk_size) + 1),
    .before = img_name
  )

# propagate chunks to points
rand_points <- lapply(rand_points, function(el) {
  el %>% left_join(images %>% select(-n), by = join_by(img_name))
}) 




## Loop over datasets and compute distances ----
#--------------------------------------------------------------------------#
lapply(1:n_sets, function(i) {
  message(paste("Processing dataset", i, "out of", n_sets))
  # Get dataset of interest
  el <- rand_points[[i]]
  
  # Split by chunk
  el <- el %>% 
    group_by(chunk) %>% 
    group_split()
  
  # Loop over chunks and compute distances
  dist_all_el <- lapply(1:n_chunk, function(j) {
    message(paste("Processing chunk", j, "out of", n_chunk))
    
    # Get chunk data
    df <- el[[j]]
    # Compute distances
    dist_all <- compute_all_dist(df, n_cores = n_cores)
    # Sleep
    Sys.sleep(30)
    
    
    # Return results
    return(dist_all)
  }) %>% 
    bind_rows() %>% 
    mutate(set = as.character(i))
  
  # Save it 
  filename <- paste0("data/08a.distances_set_", i,".Rdata")
  write_parquet(dist_all_el, sink = filename)
  
  gc(verbose = FALSE)
})


## Read saved distances and compute quantiles ----
#--------------------------------------------------------------------------#
# List of saved files
files <- list.files(path = "data", pattern = "08a.d", full.names = TRUE)
# Number of distances, initially set to NULL
n_dist <- NULL

# Quantiles to compute
probs <- seq(0, 1, length.out = 10000)

# Read files and get quantiles
df_quant <- lapply(files, function(file) {
  message(paste("Processing", file))
  
  # Read file
  df_dist <- read_parquet(file)
  
  # Number of distances, compute it only once
  if(is.null(n_dist)) {n_dist <- nrow(df_dist)}
  
  # Compute quantiles
  df_dist <- quantile(df_dist$dist, probs = probs, names = FALSE)
  gc(verbose = FALSE) # Memory cleaning

  # Return quantiles
  df_dist
  
})


## Compare quantiles ----
#--------------------------------------------------------------------------#
# Pair of sets to compare
set_pairs <- crossing(set_a = 1:n_sets, set_b = 1:n_sets) %>% filter(set_a < set_b)

# Loop over pairs
res <- lapply(1:nrow(set_pairs), function(i) {
  # Get sets of interest
  i_set_a <- set_pairs %>% slice(i) %>% pull(set_a)
  i_set_b <- set_pairs %>% slice(i) %>% pull(set_b)
  dist_set_a <- df_quant[[i_set_a]]
  dist_set_b <- df_quant[[i_set_b]]
  
  # Perform kuiper test between sets
  kt <- kuiper_test(dist_set_a, dist_set_b)
  
  # Return results
  tibble(
    n_dist = n_dist,
    test_stat = kt[1],
    set_a = as.factor(i_set_a),
    set_b = as.factor(i_set_b)
  )
}) %>% bind_rows()


## Save test results ----
#--------------------------------------------------------------------------#
save(res, file = "data/08a.big_null_data.Rdata")

