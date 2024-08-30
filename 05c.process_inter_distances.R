#/usr/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Filter intertaxonomic distances and perform Kuiper test.
# Date: 30/08/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")


## Read distances ----
#--------------------------------------------------------------------------#
saved <- list.files("data/distances", pattern = "^(02c.inter).*(.parquet)$", full.names = TRUE)


## Loop over taxa and process ----
#--------------------------------------------------------------------------#
# - filter distances
# - compute quantiles
# - compute Kuiper statistic

df_inter <- lapply(saved, function(el) {
  
  # Get pair name
  pair <- str_split_fixed(el, "_", n = 3)[,3] %>% str_remove_all(".parquet")
  message(paste0("Processing ", pair))
  
  # Read data
  df <- read_parquet(el)
  
  # Extract plankton and random distances
  plank_dist <- df %>% select(dist) %>% filter(dist < dist_thr)
  rand_dist <- df %>% select(dist = rand_dist) %>% filter(dist < dist_thr)
  
  # Compute 10000-quantiles
  plank_qt <- quantile(plank_dist$dist, probs = probs, names = FALSE)
  rand_qt <- quantile(rand_dist$dist, probs = probs, names = FALSE)
  
  # Compute Kuiper statistic
  ks <- kuiper_stat(plank_qt, rand_qt)
  
  # Format results
  res <- tibble(
    pair = pair,
    dist_thres = dist_thr, # distance threshold
    n_dist_plank = nrow(plank_dist), # number of plankton distances
    n_dist_rand = nrow(rand_dist), # number of random distances
    n_dist = (nrow(plank_dist) + nrow(rand_dist)) / 2, # average between these
    kuiper_stat = ks, # Kuiper statistic
    plank_qt = plank_qt, # plankton quantiles
    rand_qt = rand_qt # random quantiles
  ) %>% 
    nest(qt = c(plank_qt, rand_qt))
  
  # Return results
  return(res)
  
}) %>% 
  bind_rows()


## Format and save ----
#--------------------------------------------------------------------------#
# We save two tables
# - the results of the Kuiper test, 1 row per pair
# - the 10000 quantiles for both plankton and null distances

df_inter_dist <- df_inter %>% 
  select(pair, qt) %>% 
  unnest(qt)

df_inter <- df_inter %>% select(-qt)

save(df_inter, df_inter_dist, file = "data/05c.inter_distances_ks.Rdata")
