#/usr/local/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute intrataxa distances for all taxa
# Date: 12/02/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Set-up and read data ----
#--------------------------------------------------------------------------#
source("utils.R")
set.seed(seed)

message("Reading data")

# Correct image volume for x-axis
load("data/02.corr_factor.Rdata")
vol$x <- vol$x * med_corr

sub_sample <- F
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


## Select images ----
#--------------------------------------------------------------------------#
message("Selecting images")

## Keep only images containing several organisms of each type
plankton <- plankton %>% 
  add_count(img_name, taxon) %>% 
  filter(n > 1) %>% 
  select(-n)


## Work only with taxa that are present in at least 10 images
taxa <- plankton %>% 
  count(img_name, taxon) %>% 
  count(taxon) %>% 
  filter(n > 10) %>% 
  pull(taxon)
plankton <- plankton %>% filter(taxon %in% taxa)

# List retained images
images <- images %>% filter(img_name %in% plankton$img_name)
img_names <- sort(unique(images$img_name))

gc(verbose = FALSE)


## Break into chunks ----
#--------------------------------------------------------------------------#
message("Generating chunks")

# Need to do chunks of 100000 images within taxa
chunk_size <- 100000

## For plankton data
img_chunks <- plankton %>% 
  # Generate combinations of images and taxa
  count(taxon, img_name) %>% 
  select(-n) %>% 
  # Compute chunks within taxa
  group_by(taxon) %>% 
  mutate(w_chunk = floor((row_number() - 1) / chunk_size) + 1) %>% 
  ungroup()

# Get chunk numbers and taxa
chunks_info <- img_chunks %>% select(w_chunk, taxon) %>% unique()

# Each image is assigned to a chunk
# Propagate this to plankton data
plankton <- plankton %>% left_join(img_chunks, by = join_by(img_name, taxon))


## Compute distances ----
#--------------------------------------------------------------------------#
# Loop over plankton groups, and then chunk
message("Computing plankton distances")


df_intra <- lapply(taxa, function(my_taxon){

  # Inform user
  message(paste("Processing",  my_taxon))
  
  ## Get data for this taxon
  # Plankton data
  t_plankton <- plankton %>% filter(taxon == my_taxon) 
  
  # List of chunks
  t_img_chunks <- img_chunks %>% filter(taxon == my_taxon)
  chunks <- unique(t_img_chunks$w_chunk)
  
  ## Generate null data 
  message("Generating null data")
  
  # Images names
  t_img_names <- t_plankton %>% pull(img_name) %>% unique()
  # Counts per image
  t_counts <- t_plankton %>% count(img_name)
  # The number of images to generate is now limited by the number of images with the given taxon
  t_n_img <- nrow(t_counts)
  # Pick random points within image volumes
  t_null_plankton <- mclapply(1:t_n_img, function(i) {
    # Number of points to sample within image
    n <- t_counts %>% slice(i) %>% pull(n)
    
    # Draw points
    d_points <- tibble(
      x = runif(n = n, min = 1, max = vol$x),
      y = runif(n = n, min = 1, max = vol$y),
      z = runif(n = n, min = 1, max = vol$z)
    ) %>% # Add information for img name: represented image and null image name
      mutate(
        img_name = t_img_names[i],
        null_img_name = paste0("img_", str_pad(i, nchar(t_n_img), pad = "0"))
      )
  }, mc.cores = n_cores) %>% 
    bind_rows() %>% 
    mutate(taxon = my_taxon)
  
  # Propagate chunks to null data
  t_null_plankton <- t_null_plankton %>% left_join(t_img_chunks, by = join_by(img_name, taxon))
  
  ## Loop over chunks to compute distances
  dist_all <- lapply(chunks, function(j) {
    
    message(paste("Processing chunk",  j, "out of", length(chunks)))
    
    ## Get data for this chunk
    # Plankton data
    c_plankton <- t_plankton %>% filter(w_chunk == j)
    # Null data
    c_null_plankton <- t_null_plankton %>% filter(w_chunk == j)
    
    ## Compute distances
    # Plankton
    message("Plankton distances")
    dist_plank <- compute_all_dist(c_plankton, n_cores = n_cores)
    # Null data
    message("Null distances")
    dist_null <- compute_all_dist(c_null_plankton, n_cores = n_cores) %>% rename(rand_dist = dist)
    
    # Join together
    dist_all <- dist_plank %>% left_join(dist_null, by = join_by(p1, p2, img_name))
    
    # Return the result
    return(dist_all)
  }) %>% 
    bind_rows()
  
  # Sleep
  Sys.sleep(30)
  
  ## Compute 10000-quantiles
  # Save X-quantiles, depending on the number of observations
  if (nrow(dist_all) > 10000){ # if more than 10,000 distances, save 10000-quantiles
    probs <- seq(0, 1, length.out = 10000)
    dist_plank <- quantile(dist_all$dist, probs = probs, names = FALSE)
    dist_null <- quantile(dist_all$rand_dist, probs = probs, names = FALSE)
  } else { # otherwise, keep all values
    dist_plank <- dist_all$dist
    dist_null <- dist_all$rand_dist
  }
  
  ## Perform kuiper test and return results
  # Kuiper-test
  out <-  kuiper_test(dist_plank, dist_null)
  
  # Return outputs
  res <- tibble(
    taxon = my_taxon,
    n_obj = nrow(t_plankton),
    n_img = t_n_img,
    test_stat = out[1],
    p_value = out[2],
    dist = list(dist_plank),
    dist_rand = list(dist_null),
  )
  return(res)

}) %>% 
  bind_rows()



## Reformat results ----
#--------------------------------------------------------------------------#
message("Reformating")

# One big df with distances
df_intra_dist <- df_intra %>% 
  select(taxon, dist, dist_rand) %>% 
  unnest(c(dist, dist_rand))

# One small df with summary
df_intra <- df_intra %>% select(-c(dist, dist_rand))


## Save results ----
#--------------------------------------------------------------------------#
message("Saving")

save(df_intra, df_intra_dist, file = "data/05.intra_distances.Rdata")


