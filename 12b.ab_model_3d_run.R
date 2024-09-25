#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Build a 3D agent-based model to simulate attraction between organisms
# Date: 21/06/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")
source("utils_ab_model.R")


## Input data ----
#--------------------------------------------------------------------------#
## Load image info
images <- read_parquet("data/00.images_clean.parquet")
plankton <- read_parquet("data/01a.x_corrected_plankton_clean.parquet")
# Count objects per image
counts <- plankton %>% count(img_name)
# Generate a representative sample of number of objects per image, with more than 2 points per image
n_pts <- counts %>% filter(n > 2) %>%  slice_sample(n = n_img) %>% pull(n)

# Input data
n_img <- 1000 # number of images

## Generate images
#n_pts <- rep(n_pts, times = n_img)
d_points <- mclapply(1:n_img, function(i){
  # Number of points to sample within image
  n <- n_pts[i]
  # Draw points
  tibble(
    x = round(runif(n = n, min = 1, max = vol$x)),
    y = round(runif(n = n, min = 1, max = vol$y)),
    z = round(runif(n = n, min = 1, max = vol$z))
  ) %>% # Add information for img name
    #arrange(x, y) %>% 
    arrange(desc(y), x) %>% # sort point for rasterize
    mutate(
      img_name = paste0("img_", str_pad(i, nchar(n_img), pad = "0")),
      id = paste0(str_pad(row_number(), 3, pad = "0"))
      #id = row_number() %>% as.character()
    )
}, mc.cores = n_cores) %>% 
  bind_rows()

# Separate tibbles
d_points <- d_points %>% 
  group_by(img_name) %>% 
  group_split()



## Settings ----
#--------------------------------------------------------------------------#
# Define parameter grid
param_grid <- crossing(
  h = c(200000, 250000, 300000),   # Density bandwidth
  d_length = c(50, 100, 150),      # Displacement length (in px)
  prop_mv = 0.96                   # Proportion of points to move, representative of ISIIS dataset
)


## Run ----
#--------------------------------------------------------------------------#
# Process each parameter combination in parallel
res <- lapply(1:nrow(param_grid), function(i) {
  
  message(paste0("Processing parameter set ", i, " out of ", nrow(param_grid)))
  
  # Get parameters from grid
  params <- param_grid %>% slice(i)
  h <- params$h
  #gridsize <- unlist(params$gridsize)
  d_length <- params$d_length
  prop_mv <- params$prop_mv
  
  # Process all images with current parameter combination
  processed_images <- pbmclapply(d_points, function(dfi) {
    process_dataframe(
      dfi, 
      h = h,
      vol = vol,
      d_length = d_length, 
      prop_mv = prop_mv
    ) %>% 
      mutate(dist = dist * 51 / 10000)
  }, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
    bind_rows()
  
  # Apply distance threshold
  
  # Extract 10-000 quantiles
  #processed_images <- processed_images %>% 
  #  select(when, dist) %>% 
  #  group_by(when) %>% 
  #  reframe(dist = quantile(dist, probs = seq(0, 1, length.out = 10000), names = FALSE))
  
  # Return results with parameter info
  processed_images %>% 
    mutate(
      d_length = d_length,
      prop_mv = prop_mv,
      h = h
    )
}) %>% 
  bind_rows()


# Set before/after as a factor
sim_res <- res %>% 
  mutate(when = factor(when, levels = c("before", "after")))


## Save results ----
#--------------------------------------------------------------------------#
save(sim_res, file = "data/fine_simulation_results.Rdata")



#
### Recall tests ----
##--------------------------------------------------------------------------#
## Input data
#n_img <- 10000 # number of images
#
### Generate images
##n_pts <- rep(n_pts, times = n_img)
#d_points <- mclapply(1:n_img, function(i){
#  # Number of points to sample within image
#  n <- n_pts[i]
#  # Draw points
#  tibble(
#    x = round(runif(n = n, min = 1, max = vol$x)),
#    y = round(runif(n = n, min = 1, max = vol$y)),
#    z = round(runif(n = n, min = 1, max = vol$z))
#  ) %>% # Add information for img name
#    #arrange(x, y) %>% 
#    arrange(desc(y), x) %>% # sort point for rasterize
#    mutate(
#      img_name = paste0("img_", str_pad(i, nchar(n_img), pad = "0")),
#      id = paste0(str_pad(row_number(), 3, pad = "0"))
#      #id = row_number() %>% as.character()
#    )
#}, mc.cores = n_cores) %>% 
#  bind_rows()
#
## Separate tibbles
#d_points <- d_points %>% 
#  group_by(img_name) %>% 
#  group_split()
#
#dfi <- d_points[[1]]
#dfi
#
#
#
## Compute density
#h <- 400000
#d_length <- 300
#prop_mv <- 1
#
#
#
#new_points <- pbmclapply(d_points, function(dfi) {
#  # Compute density
#  dens_3d <- calculate_density_3d(dfi, h = h, vol = vol)  
#  
#  # Compute density gradient
#  grad <- calculate_density_gradient_3d(dens_3d)
#  
#  # And extract it at points of interest
#  grad_points <- extract_gradient(dfi, gradient = grad, kde = dens_3d)
#  
#  # Move points
#  new_d_points <- move_points(
#    dfi,
#    gradient_values = grad_points,
#    d_length = d_length,
#    prop_mv = prop_mv
#  )
#  
#  return(new_d_points)
#}, mc.cores = n_cores) %>% 
#  bind_rows()
#
#ori_points <- d_points %>% bind_rows()
#
#
## Compute distances between new points
## Convert from px to cm
## Keep only distances below threshold
## Recall rates to try
#rec_rates <- seq(0.1, 1, by = 0.1)
## Loop over recall rates
#rate <- rec_rates[1]
#foo <- lapply(rec_rates, function(rate) {
#  message(paste0("Processing recall of ", rate))
#  
#  ## Plankton points
#  # In each image, we only retain this proportion of points
#  sub_plankton_points <- new_points %>% 
#    group_by(img_name) %>% 
#    slice_sample(prop = rate) %>% 
#    ungroup()
#  
#  # Compute distances between subsampled points
#  plankton_dist_sub <- compute_all_dist(sub_plankton_points, n_cores = n_cores) %>% 
#    mutate(dist = dist * 51 / 10000) %>% 
#    filter(dist < dist_thr) %>% 
#    mutate(
#      recall = rate,
#      type = "plankton"
#    )
#  
#  ## Random points
#  # Subsample
#  sub_rand_points <- ori_points %>% 
#    group_by(img_name) %>% 
#    slice_sample(prop = rate) %>% 
#    ungroup()
#  
#  # Compute distances between random points
#  rand_dist_sub <- compute_all_dist(sub_rand_points, n_cores = n_cores) %>% 
#    mutate(dist = dist * 51 / 10000) %>% 
#    filter(dist < dist_thr) %>% 
#    mutate(
#      recall = rate,
#      type = "rand"
#    )
#  
#  return(bind_rows(plankton_dist_sub, rand_dist_sub))
#}) %>% 
#  bind_rows()
#
#
#ggplot(foo) +
#  geom_density(aes(x = dist, colour = recall, group = recall)) +
#  scale_colour_viridis_c() +
#  labs(x = "Distance (cm)", y = "Density", colour = "Recall") +
#  facet_wrap(~type) +
#  theme_classic()
#
## Kuiper stat between plankton and rand distances
#bar <- foo %>% 
#  select(recall, dist, type) %>% 
#  pivot_wider(names_from = type, values_from = dist, values_fn = list)
#
## Loop over recall values
#i <- 1
#baz <- lapply(1:nrow(bar), function(i) {
#  # Extract row
#  r <- bar %>% slice(i)
#  
#  # Get distances
#  plankton_dist <- unlist(r$plankton)
#  rand_dist <- unlist(r$rand)
#  
#  # Get number of distances
#  n_dist_plankton <- length(plankton_dist)
#  n_dist_rand <- length(rand_dist)
#  n_dist <- (n_dist_plankton + n_dist_rand) / 2
#  
#  # If needed, get quantiles
#  if (n_dist > 10000) {
#    plank_qt <- quantile(plankton_dist, probs = probs, names = FALSE)
#    rand_qt <- quantile(rand_dist, probs = probs, names = FALSE)
#  } else {
#    plank_qt <- plankton_dist
#    rand_qt <- rand_dist
#  }
#  
#  # Compute Kuiper stat
#  ks <- kuiper_stat(plank_qt, rand_qt)
#  
#  # Return result
#  res <- tibble(
#    recall = r$recall,
#    n_dist = n_dist,
#    log_n_dist = log10(n_dist),
#    kuiper_stat = ks,
#    log_kuiper_stat = log10(kuiper_stat)
#  )
#  
#  return(res)
#}) %>% 
#  bind_rows()
#
#
#ggplot(baz) +
#  geom_point(aes(x = log_n_dist, y = log_kuiper_stat, colour = recall)) +
#  scale_colour_viridis_c() +
#  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = seq(3, 6, by = 1)) +
#  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
#  labs(x = "N distances", y = "Kuiper statistic", colour = "Recall") 
