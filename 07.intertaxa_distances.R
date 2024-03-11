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

# list pairs of taxa
pairs <- crossing(t1 = taxa, t2 = taxa) %>% filter(t1 != t2)
#pairs <- pairs %>% slice_head(n = 100)
#pairs <- pairs %>% slice_sample(n = 100)


## Loop over pairs ----
#--------------------------------------------------------------------------#

#j = 3
#j = 85
start_time = Sys.time()
# Use indices instead of taxa
df_inter <- mclapply(1:nrow(pairs), function(j){
  
  # Get pair
  my_pair <- pairs %>% slice(j)
  t1 <- my_pair %>% pull(t1)
  t2 <- my_pair %>% pull(t2)
  my_pair_str <- paste(t1, t2, sep = " - ")
  
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
  if (nrow(t_plankton) > 10){
    
    # Images names
    t_img_names <- img_both$img_name
    
    # The number of images to generate is now limited by the number of images with both given taxon
    t_n_img <- min(n_img, nrow(img_both))
    
    ## Generate random data for given taxon
    # Pick random points within image volumes
    rand_points <- mclapply(1:t_n_img, function(i) {
      set.seed(seed)
      # Number of points to sample within image
      n1 <- img_both %>% slice(i) %>% pull(t1) # first taxon
      n2 <- img_both %>% slice(i) %>% pull(t2) # second taxon
      
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
    
    ## Comparison with null data
    # Kuiper-test
    s1 <-  dist_all
    s2 <-  dist_all_rand
    out <-  kuiper_test(s1, s2)
    
    # Return outputs
    tibble(
      pair = my_pair_str,
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

df_inter %>% 
  ggplot() + 
  geom_jitter(aes(x = 0, y = p_value, colour = str_detect(pair, "Acant")))


## Reformat results ----
#--------------------------------------------------------------------------#
# One big df with distances
df_inter_dist <- df_inter %>% 
  pivot_longer(c(dist, dist_rand)) %>% 
  mutate(value = map(value, `length<-`, max(lengths(value)))) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  unnest(c(dist, dist_rand)) %>% 
  select(pair, dist, dist_rand)

# One small df with summary
df_inter <- df_inter %>% select(-c(dist, dist_rand))


## Save results ----
#--------------------------------------------------------------------------#
save(df_inter, df_inter_dist, file = "data/06.inter_distances.Rdata")



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

