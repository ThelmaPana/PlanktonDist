#/usr/local/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Generate multiple null datasets and compute distances.
# Date: 09/04/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")

## Load data ----
#--------------------------------------------------------------------------#
## Load data
images <- read_parquet("data/00.images_clean.parquet")
plankton <- read_parquet("data/01.x_corrected_plankton_clean.parquet")
load("data/01.corr_factor.Rdata")

# list img names
img_names <- sort(unique(images$img_name))

## Apply correction factor in x
vol$x <- vol$x * med_corr

## Count number of objects per image
counts <- plankton %>% count(img_name)

# We target the following number of distances
n_tar_dist <- c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)


## Generate null datasets ----
#--------------------------------------------------------------------------#
message("Generating null data")
n_img <- 55000 # Number of images to generate
n_sets <- 1 # Number of datasets to generate

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


## Compute distances and compare ----
#--------------------------------------------------------------------------#
message("Computing distances")
dist_rand <- lapply(rand_points, compute_all_dist)

# Keep only distances smaller than threshold
dist_rand <- lapply(dist_rand, function(el) {
  el %>% filter(dist < dist_thr_px)
})

el <- dist_rand[[1]]
nrow(el) > n_tar_dist

message("Performing comparison of subsets")
# Size of subsets
dist_count <- dist_rand[[1]] %>% count(img_name) %>% mutate(n_cum = cumsum(n))

# List images to retain for each target
img_lists <- lapply(n_tar_dist, function(n_tar) {
  dist_count %>% 
    # keep image to reach the number of target distances
    mutate(keep = lag(n_cum < n_tar, default = TRUE)) %>%
    filter(keep) %>% 
    pull(img_name)
})


# Pair of sets to compare
set_pairs <- crossing(set_a = 1:n_sets, set_b = 1:n_sets) %>% filter(set_a < set_b)

# Loop on pairs of sets and perform kuiper test on each subset
f_val_dist <- pbmclapply(1:nrow(set_pairs), function(i) {
  message(paste0("Processing pair ", i, " out of ", nrow(set_pairs)))
  
  # Get sets of interest
  i_set_a <- set_pairs %>% slice(i) %>% pull(set_a)
  i_set_b <- set_pairs %>% slice(i) %>% pull(set_b)
  set_a <- dist_rand[[i_set_a]]
  set_b <- dist_rand[[i_set_b]]
  
  # Loop on targets to reach
  lapply(1:length(n_tar_dist), function(j) {
    
    # Get list of images of interest
    tar_img_list <- img_lists[[j]]
    
    # Get objects in these images
    tar_set_a <- set_a %>% filter(img_name %in% tar_img_list)
    tar_set_b <- set_b %>% filter(img_name %in% tar_img_list)
    
    # Subsample to 10-000 quantiles if needed
    if (nrow(tar_set_a) > 10000) {
      probs <- seq(0, 1, length.out = 10000)
      dist_set_a <- quantile(tar_set_a$dist, probs = probs, names = FALSE)
      dist_set_b <- quantile(tar_set_b$dist, probs = probs, names = FALSE)
    } else {
      dist_set_a <- tar_set_a$dist
      dist_set_b <- tar_set_b$dist
    }
    
    # Perform kuiper test between sets
    kt <- kuiper_test(dist_set_a, dist_set_b)
    
    # Return results
    tibble(
      n_tar_dist = n_tar_dist[j],
      n_dist = nrow(tar_set_a),
      test_stat = kt[1],
      set_a = as.factor(i_set_a),
      set_b = as.factor(i_set_b)
    )
  }) %>% 
    bind_rows()
}, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
  bind_rows()

# Check that we have the correct number of distances
ggplot(f_val_dist) + 
  geom_point(aes(x = n_tar_dist, y = n_dist)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  scale_x_log10() +
  scale_y_log10()

# The number of distances varies across sets because of filtering on small distances
# Replace it by the targeted number of distances
f_val_dist <- f_val_dist %>% 
  select(-n_dist) %>% 
  rename(n_dist = n_tar_dist)


## Plot results ----
#--------------------------------------------------------------------------#
ggplot(f_val_dist) + 
  geom_point(aes(x = n_dist, y = test_stat)) +
  scale_x_log10() + scale_y_log10() 

ggplot(f_val_dist) + 
  geom_boxplot(aes(x = n_dist, y = test_stat, group = n_dist)) +
  scale_x_log10() + scale_y_log10() 


## Save results ----
#--------------------------------------------------------------------------#
# Kuiper statistic VS number of distances
save(f_val_dist, file = "data/02a.f_val_dist_small.Rdata")

# Save on set of distances
dist_rand <- dist_rand[[1]]
save(dist_rand, file = "data/02a.dist_rand_small.Rdata")
