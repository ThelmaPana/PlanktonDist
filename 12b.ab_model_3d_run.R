#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Build a 3D agent-based model to simulate attraction between organisms
# Date: 21/06/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")
source("utils_ab_model.R")


## Settings ----
#--------------------------------------------------------------------------#
# Whether to try multiple parameter values
gridsearch <- TRUE

# Read parameters grid
load("data/12a.parameters_grid.Rdata")


if (!gridsearch) {
  # If no gridsearch, use selected parameters
  sel_param <- list()
  sel_param$sr <- 10 # sensory radius in cm
  sel_param$d_length_cm <- 0.4 # displacement length in cm
  
  # Get corresponding row in grid
  sel_row <- param_grid %>% filter(sr == sel_param$sr & d_length_cm == sel_param$d_length_cm)
  
  # Extract parameters
  h <- sel_row$h
  d_length_px <- sel_row$d_length_px
  prop_mv <- sel_row$prop_mv
}


## Input data ----
#--------------------------------------------------------------------------#
# Input data
n_img <- 1000 # number of images

## Load image info
images <- read_parquet("data/00.images_clean.parquet")
plankton <- read_parquet("data/01a.x_corrected_plankton_clean.parquet")
# Count objects per image
counts <- plankton %>% count(img_name)
# Generate a representative sample of number of objects per image, with more than 2 points per image
n_pts <- counts %>% filter(n > 2) %>%  slice_sample(n = n_img) %>% pull(n)


message("Generating images")

## Generate images
#n_pts <- rep(n_pts, times = n_img)
d_points <- pbmclapply(1:n_img, function(i){
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
}, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
  bind_rows()

# Separate tibbles
d_points <- d_points %>% 
  group_by(img_name) %>% 
  group_split()





## Gridsearch ----
#--------------------------------------------------------------------------#
if (gridsearch) {
  
  
  # Process each parameter combination in parallel
  res <- lapply(1:nrow(param_grid), function(i) {
    
    message(paste0("Processing parameter set ", i, " out of ", nrow(param_grid)))
    
    # Get parameters from grid
    params <- param_grid %>% slice(i)
    h <- params$h
    d_length_px <- params$d_length_px
    prop_mv <- params$prop_mv
    
    # Process all images with current parameter combination
    processed_images <- pbmclapply(d_points, function(dfi) {
      process_dataframe(
        dfi, 
        h = h,
        vol = vol,
        d_length = d_length_px, 
        prop_mv = prop_mv
      ) %>% 
        mutate(dist = dist * 51 / 10000)
    }, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
      bind_rows()
    
    # Return results with parameter info
    processed_images %>% 
      mutate(
        d_length_px = d_length_px,
        prop_mv = prop_mv,
        h = h
      )
  }) %>% 
    bind_rows()
  
  
  # Set before/after as a factor
  sim_res_grid <- res %>% 
    mutate(when = factor(when, levels = c("before", "after")))
  
  # Join with meaningful parameters of the grid
  sim_res_grid <- sim_res_grid %>% left_join(param_grid, by = join_by(d_length_px, prop_mv, h))
  
  ## Save results 
  write_parquet(sim_res_grid, sink = "data/12b.simulation_results_gridsearch.parquet")
}


## Single run ----
#--------------------------------------------------------------------------#

if (!gridsearch) {
  # Process all images with current parameter combination
  sim_res <- pbmclapply(d_points, function(dfi) {
    process_dataframe(
      dfi, 
      h = h,
      vol = vol,
      d_length = d_length_px, 
      prop_mv = prop_mv
    ) %>% 
      mutate(dist = dist * 51 / 10000)
  }, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
    bind_rows() %>% 
    mutate(when = factor(when, levels = c("before", "after")))
  
  ## Save results
  write_parquet(sim_res, sink = "data/12b.simulation_results.parquet")
  
}