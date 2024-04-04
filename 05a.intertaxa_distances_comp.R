#/usr/local/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute intertaxa distances for all pairs of taxa
# Date: 13/02/2024
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


## Subsampling
sub_sample <- FALSE
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

# list img names
img_names <- sort(unique(images$img_name))


# list taxa
taxa <- plankton %>% pull(taxon) %>% unique() %>% sort()

# list pairs of taxa
pairs <- crossing(t1 = taxa, t2 = taxa) %>% filter(t1 != t2)

## Remove taxa already processed
# list save files
if (sub_sample){ # with subsampling
  saved <- list.files("data/", pattern = "05.sub")  
} else { # withouth subsampling
  saved <- list.files("data/", pattern = "05.inter")
}


# whether a taxon was processed or not
processed <- sapply(taxa, function(i) {any(str_detect(saved, i))}, USE.NAMES = FALSE)

# remove already processed from the list
taxa <- taxa[!processed]


## Compute distances ----
#--------------------------------------------------------------------------#
walk(taxa, function(my_taxon) {
  message(paste0("Processing ", my_taxon))

  
  # List pairs with given taxon
  my_pairs <- pairs %>% filter(t1 == my_taxon)
  
  # Loop over pairs
  res_pair <- lapply(1:nrow(my_pairs), function(i) {
    
    # Get pair
    my_pair <- my_pairs %>% slice(i)
    t1 <- my_pair %>% pull(t1)
    t2 <- my_pair %>% pull(t2)
    my_pair_str <- paste(t1, t2, sep = " - ")
    
    message(paste0("Processing ", my_pair_str, ", pair ", i, " out of ", nrow(my_pairs), " for ", my_taxon))
    
    ## Taxon set-up, number of objects per image
    # Keep images where both organisms are present
    img_both <- plankton %>%
      filter(taxon %in% c(t1, t2)) %>% 
      count(taxon, img_name) %>% 
      mutate(taxon = ifelse(taxon == t1, "t1", "t2")) %>% 
      pivot_wider(names_from = taxon, values_from = n, values_fill = 0) %>% 
      filter(t1 > 0 & t2 > 0)
    
    # Keep organisms within these images
    t_plankton <- plankton %>% 
      filter(img_name %in% img_both$img_name) %>% 
      filter(taxon %in% c(t1, t2))
    
    # If at least 10 images with 2 or more organisms
    if (nrow(img_both) > 10){
      
      # Images names
      t_img_names <- img_both$img_name
      
      # The number of images to generate is now limited by the number of images with both given taxon
      t_n_img <- nrow(img_both)
      
      message("Generating null data")
      
      ## Generate random data for given taxon
      # Pick random points within image volumes
      #set.seed(seed)
      rand_points <- pbmclapply(1:t_n_img, function(j) {
        # Number of points to sample within image
        n1 <- img_both %>% slice(j) %>% pull(t1) # first taxon
        n2 <- img_both %>% slice(j) %>% pull(t2) # second taxon
        
        n_tot <- n1 + n2 # total number of points
        # Draw points
        d_points <- tibble(
          x = runif(n = n_tot, min = 1, max = vol$x),
          y = runif(n = n_tot, min = 1, max = vol$y),
          z = runif(n = n_tot, min = 1, max = vol$z)
        ) %>% 
          # Simulate taxa with representative proportions
          mutate(taxon = c(rep(t1, times = n1), rep(t2, times = n2))) %>%
          # Shuffle taxon for randomness
          transform(taxon = sample(taxon)) %>% 
          as_tibble() %>% 
          # Add information for img name
          mutate(img_name = paste0("img_", str_pad(j, nchar(t_n_img), pad = "0")))
      }, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
        bind_rows()
      
      Sys.sleep(10)
      
      # Loop over images and compute distances between all points within each image
      message("Computing null distances")
      dist_all_rand <- compute_all_dist_inter(rand_points, n_cores = n_cores)
      n_dist_rand <- nrow(dist_all_rand)
      # Save X-quantiles, depending on the number of observations
      if (n_dist_rand > 10000){ # if more than 10,000 distances, save 10000-quantiles
        probs <- seq(0, 1, length.out = 10000)
        dist_all_rand <- quantile(dist_all_rand$dist, probs = probs, names = FALSE)
      } else { # otherwise, keep all values
        dist_all_rand <- dist_all_rand$dist
      }
      
      Sys.sleep(10)
      
      ## Compute distances between all organisms
      message("Computing plankton distances")
      # Loop over images and compute distances between all points within each image
      dist_all <- compute_all_dist_inter(t_plankton, n_cores = n_cores)
      n_dist <- nrow(dist_all)
      # Save X-quantiles, depending on the number of observations
      if (n_dist > 10000){ # if more than 10,000 distances, save 10000-quantiles
        probs <- seq(0, 1, length.out = 10000)
        dist_all <- quantile(dist_all$dist, probs = probs, names = FALSE)
      } else { # otherwise, keep all values
        dist_all <- dist_all$dist
      }
      
      Sys.sleep(10)
      
      ## Comparison with null data
      message("Performing kuiper test")
      # Kuiper-test
      s1 <-  dist_all
      s2 <-  dist_all_rand
      out <-  kuiper_test(s1, s2)
      
      # Return outputs
      tibble(
        pair = my_pair_str,
        n_obj = nrow(t_plankton),
        n_img = t_n_img,
        n_dist = n_dist,
        test_stat = out[1],
        p_value = out[2],
        dist = list(dist_all),
        dist_rand = list(dist_all_rand),
      )
    }  
    
  }) %>% 
    bind_rows()
  
  ## Reformat and save output for given taxon
  message("Reformating & saving")
  # One big df with distances
  res_pair_dist <- res_pair %>% 
    pivot_longer(c(dist, dist_rand)) %>% 
    mutate(value = map(value, `length<-`, max(lengths(value)))) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 
    unnest(c(dist, dist_rand)) %>% 
    select(pair, dist, dist_rand)
  
  # One small df with summary
  res_pair <- res_pair %>% select(-c(dist, dist_rand))
  
  # Save
  if (sub_sample){ # add a subsampling tag to file name
    file_name <- paste0("data/05.sub_inter_distances_", my_taxon, ".Rdata")  
  } else {
    file_name <- paste0("data/05.sub_inter_distances_", my_taxon, ".Rdata")  
  }
  save(res_pair, res_pair_dist, file = file_name)
  
  # Clean
  rm(res_pair, res_pair_dist)
  gc(verbose = FALSE)
  
  Sys.sleep(10)
})

