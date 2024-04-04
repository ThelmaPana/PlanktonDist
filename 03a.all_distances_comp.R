#/usr/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute all distances between objects
# Date: 11/03/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Set-up and read data ----
#--------------------------------------------------------------------------#
source("utils.R")
#set.seed(seed)

message("Reading data")

# Correct image volume for x-axis
load("data/02.corr_factor.Rdata")
vol$x <- vol$x * med_corr


sub_sample <- TRUE

## Subsampling
if (sub_sample){
  load("data/00.images_sub.Rdata")
  load("data/02.x_corrected_plankton_sub.Rdata")
  images <- images_sub
  plankton <- plankton_sub
} else {
  ## All data
  images <- read_parquet("data/00.images_clean.parquet")
  plankton <- read_parquet("data/02.x_corrected_plankton_clean.parquet")
}


## Break into chunks ----
#--------------------------------------------------------------------------#
message("Generating chunks")

# Break into n pieces of images to process 
n_chunk <- 6 # number of chunks
chunk_size <- ceiling(nrow(images) / 7) # size of chunks
# compute chunks
images <- images %>% 
  mutate(
    chunk = as.character(floor((row_number() - 1) / chunk_size) + 1),
    .before = transect
  )
# propagate chunks to plankton
plankton <- plankton %>% left_join(images %>% select(img_name, chunk), by = join_by(img_name))

# split by chunk
plankton <- plankton %>% 
  group_by(chunk) %>% 
  group_split()

gc(verbose = FALSE)


## Distances between all organisms ----
#--------------------------------------------------------------------------#
# Loop over chunks and compute distances between all points within each image

message("Computing plankton distances")

dist_all <- lapply(1:n_chunk, function(i) {
  message(paste("Processing chunk", i, "out of", n_chunk))
  
  # Get chunk data
  df <- plankton[[i]]
  # Compute distances
  dist_all <- compute_all_dist(df, n_cores = n_cores)
  # Sleep
  Sys.sleep(30)
  
  
  # Return results
  return(dist_all)
}) %>% 
  bind_rows()

message("Done with computing plankton distances")


## Generate null data and compute distances ----
#--------------------------------------------------------------------------#
# Loop over chunks, generate null data and computed distances

message("Computing null distances")

dist_all_rand <- lapply(1:n_chunk, function(i) {
  message(paste("Processing chunk", i, "out of", n_chunk))
  
  # Get chunk data
  df <- plankton[[i]]
  
  # Count points per image
  n_pts <- df %>% count(img_name)
  
  # Generate representative null data
  # Pick random points within image volumes
  rand_points <- pbmclapply(1:nrow(n_pts), function(j){
    # Number of points to sample within image
    n <- n_pts$n[j]
    # Draw points
    d_points <- tibble(
      x = runif(n = n, min = 1, max = vol$x),
      y = runif(n = n, min = 1, max = vol$y),
      z = runif(n = n, min = 1, max = vol$z)
    ) %>% # Add information for img name
      mutate(img_name = paste0("img_", str_pad(j, nchar(nrow(n_pts)), pad = "0")))
  }, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
    bind_rows()
  
  # Compute distances
  dist_all_rand <- compute_all_dist(rand_points, n_cores = n_cores)
  
  # Sleep
  Sys.sleep(30)
  
  # Return results
  return(dist_all_rand)
}) %>% 
  bind_rows()

message("Done with computing null distances")


## Compare to null data ----
#--------------------------------------------------------------------------#
message("Comparing to null data")

# First, we extract 10000-quantiles
probs <- seq(0, 1, length.out = 10000)
dist_all_rand <- quantile(dist_all_rand$dist, probs = probs, names = FALSE)
dist_all <- quantile(dist_all$dist, probs = probs, names = FALSE)

# Perform kuiper test
out <- kuiper_test(dist_all, dist_all_rand)
out


## Save results ----
#--------------------------------------------------------------------------#
message("Saving")
# Store results
df_all <- tibble(
  taxon = "all",
  n_obj = nrow(plankton),
  n_img = nrow(images),
  test_stat = out[1],
  p_value = out[2],
  dist = list(dist_all),
  dist_rand = list(dist_all_rand),
)

# Save
save(df_all, file = "data/03.all_distances.Rdata")

