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
load("data/01a.corr_factor.Rdata")
vol$x <- vol$x * med_corr

sub_sample <- FALSE
## Subsampling
if (sub_sample){
  load("data/00.images_sub.Rdata")
  load("data/01a.x_corrected_plankton_sub.Rdata")
  images <- images_sub
  plankton <- plankton_sub
} else {
  ## All data
  images <- read_parquet("data/00.images_clean.parquet")
  plankton <- read_parquet("data/01a.x_corrected_plankton_clean.parquet")
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


## Compute distances and compare ----
#--------------------------------------------------------------------------#
message("Computing distances")

dist_rand <- lapply(1:n_sets, function(i) {
  compute_all_dist(rand_points[[i]]) %>% 
    mutate(dist = dist * 51 / 10000) %>% # convert from px to cm
    mutate(set = as.character(i)) # flag set
}) %>% 
  bind_rows()

# Apply distance threshold
dist_rand <- dist_rand |> filter(dist < dist_thr)

message("Performing comparison of subsets")
# Size of subsets
dist_count <- dist_rand |> filter(set == "1") %>% count(img_name) %>% mutate(n_cum = cumsum(n))


# Pair of sets to compare
set_pairs <- crossing(set_a = 1:n_sets, set_b = 1:n_sets) %>% filter(set_a < set_b)

# Loop on pairs of sets and perform kuiper test on each subset
null_ks_n_dist_acant <- pbmclapply(1:nrow(set_pairs), function(i) {
  message(paste0("Processing pair ", i, " out of ", nrow(set_pairs)))
  
  # Get sets of interest
  i_set_a <- set_pairs %>% slice(i) %>% pull(set_a)
  i_set_b <- set_pairs %>% slice(i) %>% pull(set_b)
  set_a <- dist_rand |> filter(set == as.character(i_set_a))
  set_b <- dist_rand |> filter(set == as.character(i_set_b))
  
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
    kuiper_stat = kt[1],
    set_a = as.factor(i_set_a),
    set_b = as.factor(i_set_b)
  )
}, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
  bind_rows()

# Compute log10 of n_dist and kuiper_stat
null_ks_n_dist_acant <- null_ks_n_dist_acant %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_kuiper_stat = log10(kuiper_stat)
  )


## Save ----
#--------------------------------------------------------------------------#
save(null_ks_n_dist_acant, file = "data/10a.null_ks_n_dist_acant.Rdata")
