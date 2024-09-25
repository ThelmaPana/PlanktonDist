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



## Investigate density bandwidth ----
#--------------------------------------------------------------------------#

df <- d_points[[1]]
ggplot(df) +
  geom_point(aes(x = x, y = y)) +
  coord_fixed()


# Compute distances before
dist_before <- compute_pairwise_distances(df) %>% mutate(when = "before")

# Compute 3D density
dens_3d <- calculate_density_3d(df, h = h, vol = vol)

# Compute density gradient
grad_3d <- calculate_density_gradient_3d(dens_3d)

# Extract gradient at points
d_grad <- extract_gradient(
  df,
  gradient = grad_3d,
  kde = dens_3d
)

# Move points
new_d_points <- move_points(
  df,
  gradient_values = d_grad,
  d_length = d_length,
  prop_mv = prop_mv
)

# Compute distances after moving
dist_after <- compute_pairwise_distances(new_d_points) %>% mutate(when = "after")
