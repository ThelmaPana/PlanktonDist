#/usr/bin/R
#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Filter intra distances and perform Kuiper test.
# Date: 29/08/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")


## Read distances ----
#--------------------------------------------------------------------------#
plank_dist <- read_parquet("data/distances/02a.all_distances_plankton.parquet")
rand_dist <- read_parquet("data/distances/02a.all_distances_random.parquet")


## Filter distances using threshold ----
#--------------------------------------------------------------------------#
plank_dist <- plank_dist %>% filter(dist < dist_thr)
rand_dist <- rand_dist %>% filter(dist < dist_thr)


## Kuiper statistic ----
#--------------------------------------------------------------------------#
# First we need to extract 10000-quantiles.
plank_qt <- quantile(plank_dist$dist, probs = probs, names = FALSE)
rand_qt <- quantile(rand_dist$dist, probs = probs, names = FALSE)

# Compute Kuiper statistic
ks <- kuiper_stat(plank_qt, rand_qt)


## Format and save ----
#--------------------------------------------------------------------------#
# We save two tables
# - the results of the Kuiper test, 1 row
# - the 10000 quantiles for both plankton and null distances

df_all <- tibble(
  taxon = "all",
  dist_thres = dist_thr, # distance threshold
  n_dist_plank = nrow(plank_dist), # number of plankton distances
  n_dist_rand = nrow(rand_dist), # number of random distances
  n_dist = (nrow(plank_dist) + nrow(rand_dist)) / 2, # average between these
  kuiper_stat = ks
)

df_all_dist <- tibble(
  taxon = "all",
  plank = plank_qt,
  rand = rand_qt
)

save(df_all, df_all_dist, file = "data/05a.all_distances_ks.Rdata")
