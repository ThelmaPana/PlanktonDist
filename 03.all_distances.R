#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute all distances between objects
# Date: 11/03/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute intrataxa distances for all taxa
# Date: 12/02/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Set-up and read data ----
#--------------------------------------------------------------------------#
source("utils.R")


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

## Null hypothesis data
load("data/02.x_corrected_null_data.Rdata")



## Distances between all organisms ----
#--------------------------------------------------------------------------#
# Loop over images and compute distances between all points within each image
dist_all <- compute_all_dist(plankton, n_cores = n_cores)



## Compare to null data ----
#--------------------------------------------------------------------------#
# First, we extract 10000-quantiles
probs <- seq(0, 1, length.out = 10000)
dist_all_rand <- quantile(dist_all_rand$dist, probs = probs, names = FALSE)
dist_all <- quantile(dist_all$dist, probs = probs, names = FALSE)

# Perform kuiper test
out <- kuiper_test(dist_all, dist_all_rand)
out


## Save results ----
#--------------------------------------------------------------------------#
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

