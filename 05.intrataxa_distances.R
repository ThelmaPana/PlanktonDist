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

# list img names
img_names <- sort(unique(images$img_name))


# list taxa
taxa <- plankton %>% pull(taxon) %>% unique() %>% sort()


## Loop over taxa ----
#--------------------------------------------------------------------------#
start_time = Sys.time()
df_intra <- mclapply(taxa, function(my_taxon){
#for (my_taxon in taxa){
  print(my_taxon)
  
  ## Taxon set-up, number of objects per image
  # Keep organisms of given taxa, only in images with > 1 organism
  t_plankton <- plankton %>%
    filter(taxon == my_taxon) %>%
    add_count(img_name) %>%
    filter(n > 1) %>%
    select(-n)
  
  # If at least 10 images with 2 or more organisms
  if (nrow(t_plankton) > 10){
    
    # Images names
    t_img_names <- t_plankton %>% pull(img_name) %>% unique()
    
    # Counts per image
    t_counts <- t_plankton %>% count(img_name)
    
    
    # The number of images to generate is now limited by the number of images with the given taxon
    t_n_img <- min(n_img, nrow(t_counts))
    n_pts <- t_counts %>% slice_sample(n = n_img, replace = TRUE) %>% pull(n)
    
    
    ## Generate random data for given taxon
    # Pick random points within image volumes
    rand_points <- mclapply(1:t_n_img, function(i) {
      # Number of points to sample within image
      n <- n_pts[i]
      # Draw points
      d_points <- tibble(
        x = runif(n = n, min = 1, max = vol$x),
        y = runif(n = n, min = 1, max = vol$y),
        z = runif(n = n, min = 1, max = vol$z)
      ) %>% # Add information for img name
        mutate(img_name = paste0("img_", str_pad(i, nchar(t_n_img), pad = "0")))
    }, mc.cores = n_cores)
    
    # Store this in a df
    rand_points <- do.call(bind_rows, rand_points)
    
    
    # Loop over images and compute distances between all points within each image
    dist_all_rand <- compute_all_dist(rand_points, n_cores = n_cores)
    # Save X-quantiles, depending on the number of observations
    if (nrow(dist_all_rand) > 10000){ # if more than 10,000 distances, save 10000-quantiles
      probs <- seq(0, 1, length.out = 10000)
      dist_all_rand <- quantile(dist_all_rand$dist, probs = probs, names = FALSE)
    } else { # otherwise, keep all values
      dist_all_rand <- dist_all_rand$dist
    }
    
    ## Compute distances between all organisms
    # Loop over images and compute distances between all points within each image
    dist_all <- compute_all_dist(t_plankton, n_cores = n_cores)
    # Save X-quantiles, depending on the number of observations
    if (nrow(dist_all) > 10000){ # if more than 10,000 distances, save 10000-quantiles
      probs <- seq(0, 1, length.out = 10000)
      dist_all <- quantile(dist_all$dist, probs = probs, names = FALSE)
    } else { # otherwise, keep all values
      dist_all <- dist_all$dist
    }

    # Kuiper-test
    s1 <-  dist_all
    s2 <-  dist_all_rand
    out <-  kuiper_test(s1, s2)
    
    # Return outputs
    tibble(
      taxon = my_taxon,
      n_obj = nrow(t_plankton),
      n_img = t_n_img,
      test_stat = out[1],
      p_value = out[2],
      dist = list(dist_all),
      dist_rand = list(dist_all_rand),
    )

  }  
}, mc.cores = n_cores) %>% 
  bind_rows()

end_time = Sys.time()
end_time - start_time

## Reformat results ----
#--------------------------------------------------------------------------#
# One big df with distances
df_intra_dist <- df_intra %>% 
  pivot_longer(c(dist, dist_rand)) %>% 
  mutate(value = map(value, `length<-`, max(lengths(value)))) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  unnest(c(dist, dist_rand)) %>% 
  select(taxon, dist, dist_rand)

# One small df with summary
df_intra <- df_intra %>% select(-c(dist, dist_rand))


## Save results ----
#--------------------------------------------------------------------------#
save(df_intra, df_intra_dist, file = "data/05.intra_distances.Rdata")


