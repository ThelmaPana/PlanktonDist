#/usr/local/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Generate null VS null for Acantharea-like data
# Date: 14/05/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Set-up and read data ----
#--------------------------------------------------------------------------#
source("utils.R")
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


## Select images ----
#--------------------------------------------------------------------------#
message("Selecting images")

## Keep only images with several Acantharea
plankton <- plankton %>% 
  filter(taxon == "Acantharea") %>% 
  add_count(img_name, taxon) %>% 
  filter(n > 1) %>% 
  select(-n)

# List retained images
images <- images %>% filter(img_name %in% plankton$img_name)
img_names <- sort(unique(images$img_name))

# Count objects per image
counts <- plankton %>% count(img_name)

gc(verbose = FALSE)


## Generate null data ----
#--------------------------------------------------------------------------#
message("Generating null data")
n_img <- nrow(images) # Number of images to generate
n_sets <- 5 # Number of datasets to generate

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

toto <- rand_points[[1]]


## Compute distances and compare ----
#--------------------------------------------------------------------------#
message("Computing distances")
dist_rand <- lapply(rand_points, compute_all_dist)


message("Performing comparison of subsets")
# Size of subsets
dist_count <- dist_rand[[1]] %>% count(img_name) %>% mutate(n_cum = cumsum(n))

## List images to retain for each target
#img_lists <- lapply(n_tar_dist, function(n_tar) {
#  dist_count %>% 
#    # keep image to reach the number of target distances
#    mutate(keep = lag(n_cum < n_tar, default = TRUE)) %>%
#    filter(keep) %>% 
#    pull(img_name)
#})


# Pair of sets to compare
set_pairs <- crossing(set_a = 1:n_sets, set_b = 1:n_sets) %>% filter(set_a < set_b)

# Loop on pairs of sets and perform kuiper test on each subset
f_val_acant_dist <- pbmclapply(1:nrow(set_pairs), function(i) {
  message(paste0("Processing pair ", i, " out of ", nrow(set_pairs)))
  
  # Get sets of interest
  i_set_a <- set_pairs %>% slice(i) %>% pull(set_a)
  i_set_b <- set_pairs %>% slice(i) %>% pull(set_b)
  set_a <- dist_rand[[i_set_a]]
  set_b <- dist_rand[[i_set_b]]
  
  ## Get list of images of interest
  #tar_img_list <- img_lists[[j]]
  #
  ## Get objects in these images
  #tar_set_a <- set_a %>% filter(img_name %in% tar_img_list)
  #tar_set_b <- set_b %>% filter(img_name %in% tar_img_list)
    
  # Subsample to 10-000 quantiles if needed
  if (nrow(set_a) > 10000) {
    probs <- seq(0, 1, length.out = 10000)
    dist_set_a <- quantile(set_a$dist, probs = probs, names = FALSE)
    dist_set_b <- quantile(set_b$dist, probs = probs, names = FALSE)
  } else {
    dist_set_a <- set_a$dist
    dist_set_b <- set_b$dist
  }
    
  # Perform kuiper test between sets
  kt <- kuiper_test(dist_set_a, dist_set_b)
  
  # Return results
  tibble(
    n_dist = nrow(set_a),
    test_stat = kt[1],
    set_a = as.factor(i_set_a),
    set_b = as.factor(i_set_b)
  )
}, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
  bind_rows()


# Compute log10 of n_dist and test_stat
f_val_acant_dist <- f_val_acant_dist %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_test_stat = log10(test_stat)
  )
  

## Save ----
#--------------------------------------------------------------------------#
save(f_val_acant_dist, file = "data/07.f_val_acant_dist.Rdata")
