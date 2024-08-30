#/usr/local/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute intertaxa distances for all pairs of taxa
# Date: 29/08/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Set-up and read data ----
#--------------------------------------------------------------------------#
source("utils.R")
set.seed(seed)
sleep_time <- 1

message("Reading data")

# Correct image volume for x-axis
load("data/01a.corr_factor.Rdata")
vol$x <- vol$x * med_corr


## Subsampling
sub_sample <- FALSE
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


## Prepare pairs ----
#--------------------------------------------------------------------------#
# List img names
img_names <- sort(unique(images$img_name))

# List taxa
taxa <- plankton %>% pull(taxon) %>% unique() %>% sort()

# Generate all pairs of taxa
pairs <- crossing(t1 = taxa, t2 = taxa) %>% 
  filter(t1 != t2) %>% 
  mutate(name = paste(t1, t2, sep = "-"))

## Remove taxa already processed
# list save files
if (sub_sample){ # with subsampling
  saved <- list.files("data/distances", pattern = "^(02c.sub).*(.parquet)$")
} else { # without subsampling
  saved <- list.files("data/distances", pattern = "^(02c.inter).*(.parquet)$")
}

# Detect whether a pair was processed or not
processed <- sapply(pairs$name, function(i) {any(str_detect(saved, i))}, USE.NAMES = FALSE)

# Remove already processed from the list
if (any(processed)) {pairs <- pairs[!processed, ]}

gc(verbose = FALSE)


## Compute distances ----
#--------------------------------------------------------------------------#
# Loop on pairs
# - generate null data
# - compute plankton distances
# - compute null distances

walk(1:nrow(pairs), function(i) {
  
  # Get pair name
  my_pair <- pairs %>% slice(i)
  message(paste0("Processing ", my_pair$name, " pair ", i, " out of ", nrow(pairs)))
    
  # Get involved groups
  t1 <- my_pair %>% pull(t1)
  t2 <- my_pair %>% pull(t2)
    
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
    
    Sys.sleep(sleep_time)
    
    # Loop over images and compute distances between all points within each image
    message("Computing null distances")
    dist_null <- compute_all_dist_inter(rand_points, n_cores = n_cores) %>% 
      mutate(dist = dist * 51 / 10000) 
    
    Sys.sleep(sleep_time)
    
    ## Compute distances between all organisms
    message("Computing plankton distances")
    # Loop over images and compute distances between all points within each image
    dist_plank <- compute_all_dist_inter(t_plankton, n_cores = n_cores) %>% 
      mutate(dist = dist * 51 / 10000)
    
    Sys.sleep(sleep_time)
    
    # Join together
    dist_all <- bind_cols(
      dist_plank %>% select(img_name, dist),
      dist_null %>% select(rand_dist = dist)
    )
    # Store information regarding pair
    dist_all <- crossing(my_pair, dist_all)
    
    # Save
    if (sub_sample){ # add a subsampling tag to file name
      file_name <- paste0("data/distances/02c.sub_inter_distances_", my_pair$name, ".parquet")  
    } else {
      file_name <- paste0("data/distances/02c.inter_distances_", my_pair$name, ".parquet")  
    }
    write_parquet(dist_all, sink = file_name)
  }  
})

