#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Filter null datasets and perform pairwise comparisons.
# Date: 29/08/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")

## Set-up and load data ----
#--------------------------------------------------------------------------#
# Read distances from null datasets
null_df <- read_parquet("data/04a.null_datasets_distances.parquet")

# We target the following number of distances
n_tar_dist <- c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)


## Filter distances ----
#--------------------------------------------------------------------------#
# Group and split by set
null_df <- null_df %>% 
  group_by(set) %>% 
  group_split()
# Extract number of sets
n_sets <- length(null_df)

# This returns a list, apply a filter on each element of the list
null_df <- lapply(null_df, function(el) {
  el %>% filter(dist < dist_thr)
})


## Prepare subsets for various number of distances ----
#--------------------------------------------------------------------------#
# Count distances per image
dist_count <- null_df[[1]] %>% count(img_name) %>% mutate(n_cum = cumsum(n))

# List images to retain for each target
img_lists <- lapply(n_tar_dist, function(n_tar) {
  dist_count %>% 
    # keep image to reach the number of target distances
    mutate(keep = lag(n_cum < n_tar, default = TRUE)) %>%
    filter(keep) %>% 
    pull(img_name)
})


## Perform pairwise comparison ----
#--------------------------------------------------------------------------#
# Generate all combinations of pairs between sets
set_pairs <- crossing(set_a = 1:n_sets, set_b = 1:n_sets) %>% filter(set_a < set_b)

# Loop on pairs of sets and perform kuiper test on each subset
null_ks_n_dist <- pbmclapply(1:nrow(set_pairs), function(i) {
  message(paste0("Processing pair ", i, " out of ", nrow(set_pairs)))
  
  # Get sets of interest
  i_set_a <- set_pairs %>% slice(i) %>% pull(set_a)
  i_set_b <- set_pairs %>% slice(i) %>% pull(set_b)
  set_a <- null_df[[i_set_a]]
  set_b <- null_df[[i_set_b]]
  
  # Loop on targets to reach
  lapply(1:length(n_tar_dist), function(j) {
    
    # Get list of images of interest
    tar_img_list <- img_lists[[j]]
    
    # Get objects in these images
    tar_set_a <- set_a %>% filter(img_name %in% tar_img_list)
    tar_set_b <- set_b %>% filter(img_name %in% tar_img_list)
    
    # Subsample to 10-000 quantiles if needed
    if (nrow(tar_set_a) > 10000) {
      probs <- seq(0, 1, length.out = 10000)
      dist_set_a <- quantile(tar_set_a$dist, probs = probs, names = FALSE)
      dist_set_b <- quantile(tar_set_b$dist, probs = probs, names = FALSE)
    } else {
      dist_set_a <- tar_set_a$dist
      dist_set_b <- tar_set_b$dist
    }
    
    # Perform kuiper test between sets
    ks <- kuiper_stat(dist_set_a, dist_set_b)
    
    # Return results
    tibble(
      n_tar_dist = n_tar_dist[j], # target number of distances
      n_dist_a = nrow(tar_set_a), # number of distances in set a
      n_dist_b = nrow(tar_set_b), # number of distances in set b
      n_dist = (n_dist_a + n_dist_b) / 2, # average between these
      kuiper_stat = ks, # kuiper statistic
      set_a = as.factor(i_set_a),
      set_b = as.factor(i_set_b)
    )
  }) %>% 
    bind_rows()
}, mc.cores = n_cores, ignore.interactive = TRUE) %>% 
  bind_rows()

null_ks_n_dist <- null_ks_n_dist %>% select(n_tar_dist, n_dist_a, n_dist_b, n_dist, kuiper_stat, set_a, set_b)

# Check that we have the correct number of distances
ggplot(null_ks_n_dist) + 
  geom_point(aes(x = n_tar_dist, y = n_dist)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  scale_x_log10() +
  scale_y_log10()
# OK!

# The number of distances varies across sets because of filtering on small distances
# Replace it by the targeted number of distances
null_ks_n_dist <- null_ks_n_dist %>% 
  select(-c(n_dist, n_dist_a, n_dist_b)) %>% 
  rename(n_dist = n_tar_dist)


## Perform quantile regression to estimate null Kuiper statistic for larger number of distances ----
#--------------------------------------------------------------------------#
ggplot(null_ks_n_dist) + 
  geom_point(aes(x = n_dist, y = kuiper_stat)) +
  scale_x_log10() +
  scale_y_log10()

# 190 values for each targeted number of distances

#null_ks_n_dist
# Apply log-transformation
null_ks_n_dist <- null_ks_n_dist %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_kuiper_stat = log10(kuiper_stat)
  )

# Perform quantile regression on 0.05 and 0.95 quantiles
rq_fit <- rq(log_kuiper_stat ~ log_n_dist, data = null_ks_n_dist, tau = c(0.05, 0.95))
# Perform linear regression -> mean
lm_fit <- lm(log_kuiper_stat ~ log_n_dist, data = null_ks_n_dist)

#  Extract coefficients
rq_coef <- tidy(rq_fit) %>% 
  select(term, estimate, tau) %>% 
  mutate(tau = as.character(tau)) %>% 
  mutate(term = case_when(
    term == "(Intercept)" ~ "intercept",
    term == "log_n_dist" ~ "slope",
  )) %>% 
  pivot_wider(names_from = "term", values_from = estimate)
lm_coef <- tidy(lm_fit) %>% 
  select(term, estimate) %>% 
  mutate(term = case_when(
    term == "(Intercept)" ~ "intercept",
    term == "log_n_dist" ~ "slope",
  )) %>% 
  pivot_wider(names_from = "term", values_from = estimate) %>% 
  mutate(tau = "mean", .before = intercept)

# Store all coef together
null_ks_rq_coef <- rq_coef %>% 
  bind_rows(lm_coef) %>% 
  mutate(tau = as.factor(tau))


## Check that our regression is ok
# Generate data to plot a ribbon between the regression lines
lim_dist <- c(5e1, 2e9) # extend of the ribbon in x-axis
null_ks_rib <- tibble(n_dist = lim_dist) %>% 
  mutate(
    # apply log-transformation
    log_n_dist = log10(n_dist),
    # compute estimated kuiper-stat from slope and intercept
    ymin = null_ks_rq_coef %>% filter(tau == 0.05) %>% pull(slope) * log_n_dist + null_ks_rq_coef %>% filter(tau == 0.05) %>% pull(intercept),
    ymax = null_ks_rq_coef %>% filter(tau == 0.95) %>% pull(slope) * log_n_dist + null_ks_rq_coef %>% filter(tau == 0.95) %>% pull(intercept)
  ) %>% 
  # reformat and reorder to plot a polygon
  pivot_longer(ymin:ymax, values_to = "y") %>% 
  mutate(order = c(1, 2, 4, 3)) %>% 
  arrange(order)

# Plot
ggplot(null_ks_n_dist) +
  geom_boxplot(aes(x = log_n_dist, y = log_kuiper_stat, group = log_n_dist), colour = "gray") +
  geom_polygon(data = rib_data, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_abline(data = null_ks_rq_coef %>% filter(tau != "mean"), aes(slope = slope, intercept = intercept, group = tau), linetype = 3) +
  #geom_abline(data = rq_coef %>% filter(tau == "mean"), aes(slope = slope, intercept = intercept, group = tau), linetype = 2) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), expand = c(0, 0)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
  labs(x = "N distances", y = "Kuiper statistic")
# OK!


## Save results ----
#--------------------------------------------------------------------------#
save(null_ks_n_dist, null_ks_rq_coef, null_ks_rib, file = "data/04b.null_ks.Rdata")
