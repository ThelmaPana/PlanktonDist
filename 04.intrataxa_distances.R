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
  load("data/00.subsample.Rdata")
  images <- images_sub
  plankton <- plankton_sub
} else {
  ## All data
  images <- read_parquet("data/00.images_clean.parquet")
  plankton <- read_parquet("data/00.plankton_clean.parquet")
}

# list img names
img_names <- sort(unique(images$img_name))


# list taxa
taxa <- plankton %>% pull(taxon) %>% unique() %>% sort()


## Loop over taxa ----
#--------------------------------------------------------------------------#

res <- mclapply(taxa, function(my_taxon){
  
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
    
    
    ## Compute distances between all organisms
    # Loop over images and compute distances between all points within each image
    dist_all <- compute_all_dist(t_plankton, n_cores = n_cores)
    
    
    ## Comparison with null data
    # Store distances and null distances together
    if (sub_sample){
      set.seed(seed)
      df_dist <- bind_rows(
        dist_all %>% select(dist) %>% mutate(data = my_taxon) %>% slice_sample(n = n_dist),
        dist_all_rand %>% select(dist) %>% mutate(data = "null") %>% slice_sample(n = n_dist)
      )
    } else {
      df_dist <- bind_rows(
        dist_all %>% select(dist) %>% mutate(data = my_taxon),
        dist_all_rand %>% select(dist) %>% mutate(data = "null"))
    }
    
    # Kuiper-test
    s1 <-  df_dist %>% filter(data == my_taxon) %>% pull(dist)
    s2 <-  df_dist %>% filter(data == "null") %>% pull(dist)
    out <-  kuiper_test(s1, s2)
    
    # Return outputs
    tibble(
      taxon = my_taxon,
      n_obj = nrow(t_plankton),
      n_img = t_n_img,
      test_stat = out[1],
      p_value = out[2],
      dist = list(df_dist %>% filter(data != "null") %>% pull(dist)),
      dist_rand = list(df_dist %>% filter(data == "null") %>% pull(dist)),
    )

  }  
}, mc.cores = n_cores)


# Combine into a tibble
df_intra <- do.call(bind_rows, res)


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

save(df_intra, df_intra_dist, file = "data/03.intra_distances.Rdata")





### ----quadrats--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 1 image = 5 square frames
## we can do 10 quadrats in x and 2 in y: square quadraits, 20 per image
#nx <- 10 # quadrats in x dimension
#ny <- 2  # quadrats in y dimension
#
## p-value threshold for randomness
#thres <- 0.01
#
## Loop over images and perform the quadrat test in each one
#qt_all <- mclapply(img_names, function(name) {
#  # Get points within image
#  points <- t_plankton %>% 
#    filter(img_name == name) %>% 
#    select(x, y)
#  
#  # Convert to ppp
#  points_ppp <- ppp(points$x, points$y, window = owin(xrange = c(1, vol$x), yrange = c(1, vol$y)))
#  
#  # Perform quadrat test and extract p-value
#  qt <- quadrat.test(points_ppp, nx = nx, ny = ny, method = "MonteCarlo", conditional = TRUE)
#  
#  # Store results in a tibble and return it
#  tibble(img_name = name, n_obj = nrow(points), p_value = qt$p.value)
#}, mc.cores = n_cores)
#
## Transform list to one tibble
#df_qt <- do.call(bind_rows, qt_all) %>% 
#  # distribution is random if p_value > thres
#  mutate(random = ifelse(p_value < thres, FALSE, TRUE))
#summary(df_qt)
#
#
### ----plot_p_val------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ggplot(df_qt) + geom_histogram(aes(x = p_value), bins = 100)
#ggplot(df_qt) + geom_density_2d(aes(x = n_obj, y = p_value))

