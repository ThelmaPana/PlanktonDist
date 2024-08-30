## Investigate intra distances for Doliolida

source("utils.R")

## Load data ----
#--------------------------------------------------------------------------#
load("data/02a.f_val_dist_small.Rdata")
load("data/02b.rq_coef_small.Rdata")

## Intra distances
# list processed files
processed <- list.files("data", pattern = "04\\.", full.names = TRUE)

# load data
res <- sapply(processed, function(x) mget(load(x)), simplify = TRUE)

# get summary data (1st line)
df_intra <- res[1,] %>% bind_rows()
# get distances (2nd line)
df_intra_dist <- res[2,] %>% bind_rows()

## Drop colonial Collodaria because of segmentation artifact
df_intra <- df_intra %>% filter(!str_detect(taxon, "Collodaria_colonial"))
df_intra_dist <- df_intra_dist %>% filter(!str_detect(taxon, "Collodaria_colonial"))

## Assign colour and shape to each taxon
df_intra <- df_intra %>% 
  mutate(
    colour = as.character(paletteer_d(`"khroma::discreterainbow"`, n = 27)),
    shape = rep(21:25, 6)[1:27]
  )

## Process null ----
#--------------------------------------------------------------------------#
# Apply log transformation for plotting
f_val_dist <- f_val_dist %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_test_stat = log10(test_stat)
  )

## Generate data to plot a ribbon between the regression lines
# limits for n_dist
lim_dist <- c(50, 2e9)
# Generate ribbon
rib_data <- tibble(n_dist = lim_dist) %>% 
  mutate(
    # apply log-transformation
    log_n_dist = log10(n_dist),
    # compute estimated kuiper-stat from slope and intercept
    ymin = rq_coef %>% filter(tau == 0.05) %>% pull(slope) * log_n_dist + rq_coef %>% filter(tau == 0.05) %>% pull(intercept),
    ymax = rq_coef %>% filter(tau == 0.95) %>% pull(slope) * log_n_dist + rq_coef %>% filter(tau == 0.95) %>% pull(intercept)
  ) %>% 
  # reformat and reorder to plot a polygon
  pivot_longer(ymin:ymax, values_to = "y") %>% 
  mutate(order = c(1, 2, 4, 3)) %>% 
  arrange(order)

## Process plankton ----
#--------------------------------------------------------------------------#
## Keep only small distances
#df_intra_dist_small <- df_intra_dist %>% 
#  pivot_longer(dist:dist_rand) %>% 
#  mutate(value = value * 51 / 10000) %>% # from pixel to cm
#  filter(value < dist_thr_cm)

df_intra_dist_small <- df_intra_dist %>% 
  pivot_longer(dist:dist_rand) %>% 
  mutate(value = value * 51 / 10000) %>% 
  left_join(plankton_esd, by = join_by(taxon)) #%>% 
  #filter(value < dist_thr)


# Perform Kuiper test using only small distances
df_intra_small <- mclapply(df_intra$taxon, function (ta) {
  # Get distances
  d1 <- df_intra_dist_small %>% filter(taxon == ta) %>% filter(name == "dist") %>% pull(value) # plankton
  d2 <- df_intra_dist_small %>% filter(taxon == ta) %>% filter(name == "dist_rand") %>% pull(value) # null
  
  # Number of distances
  n_dist_small <- (length(d1) + length(d2)) / 2
  
  # Kuiper test
  kt <- kuiper_test(d1, d2)
  
  # Return results
  return(
    tibble(
      taxon = ta,
      prop_dist_small = n_dist_small / 10000, # proportion of small distances based on the quantiles
      test_stat = kt[1],
      p_value = kt[2]
    ))
  
}, mc.cores = n_cores) %>% 
  bind_rows()

# compute number of small distances based on proportion of small distances and the total number of distances
df_intra_small <- df_intra %>% 
  select(taxon, n_dist, colour, shape) %>% 
  left_join(df_intra_small, by = join_by(taxon)) %>% 
  mutate(
    n_dist_small = round(n_dist * prop_dist_small), # total number of small distances
    log_n_dist_small = log10(n_dist_small), # log transform number of small distances
    log_test_stat = log10(test_stat) # log transform test stat
  ) %>% 
  filter(n_dist_small > 100) # keep only taxons with at least 100 small distances

# Detect taxa which are above the ribbon and have the minimum number of distances required
df_intra_small <- df_intra_small %>% 
  mutate(
    above = log_test_stat > log_n_dist_small * (rq_coef %>% filter(tau == 0.95) %>% pull(slope)) + (rq_coef %>% filter(tau == 0.95) %>% pull(intercept)), # above the polygon
    keep = n_dist_small > n_dist_min, # enough distances
    above = above & keep # combination of both
  ) %>% 
  select(-keep)
df_intra_small_above <- df_intra_small %>% filter(above) # needed to keep the same colour and shape across plots


## Plot distributions of distances ----
#--------------------------------------------------------------------------#
df_intra_dist_small <- df_intra_dist_small %>% 
  mutate(name = ifelse(name == "dist", "plankton", "null"))

p1a <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  geom_vline(xintercept = c(5, 10, 20), alpha = 0.5, linewidth = 0.5) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", title = "All distances") +
  theme_classic()

p2a <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 20) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", title = "Below 20") +
  theme_classic()

p3a <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 10) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", title = "Below 10") +
  theme_classic()

p4a <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 5) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", title = "Below 5") +
  theme_classic()

p1a + p2a + p3a + p4a


# Try clipping x axis for all distances
df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  geom_vline(xintercept = c(5, 10, 20), alpha = 0.5, linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 20)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", title = "All distances") +
  theme_classic()
## Same result as for filtering distances


## Plot ECDF of distances ----
#-------------------------------------------------------------------------#
p1b <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name)) +
  geom_vline(xintercept = c(5, 10, 20), alpha = 0.5, linewidth = 0.5) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "All distances") +
  theme_classic()

p2b <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 20) %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 20") +
  theme_classic()

p3b <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 10) %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 10") +
  theme_classic()

p4b <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 5) %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 5") +
  theme_classic()


p1b + p2b + p3b + p4


(p2a + geom_vline(xintercept = 7.35, alpha = 0.5, linewidth = 0.5)) + (p3a + geom_vline(xintercept = 6.15, alpha = 0.5, linewidth = 0.5)) + p2b + p3b

(p3a + geom_vline(xintercept = 6.15, alpha = 0.5, linewidth = 0.5)) + (p4a) + p3b + p4b

# Try clipping x axis for all distances
df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name)) +
  #geom_vline(xintercept = c(5, 10, 20), alpha = 0.5, linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 20)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "All distances") +
  theme_classic()
p2b
## Same result as for filtering distances


## Plot differences in ECDF ----
#--------------------------------------------------------------------------#
# Reshape df to have one row per taxon
df_wide <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  pivot_wider(names_from = "name", values_from = "value", values_fn = list)

# Compute ecdf
# Needs a nested call: ecdf(x) returns a function, x being the values of interest
# Call this function with argument x2 to get values along x2 (here, plankton distances so that we have a common x axis) 
#x_axis_seq <- 
summary(df_intra_dist_small$value)
x_axis_seq <- seq(0, 70, length.out = 1000)
plank_ecdf <- ecdf(unlist(df_wide$plankton))(x_axis_seq)
null_ecdf <- ecdf(unlist(df_wide$null))(x_axis_seq)
#TODO check this

# Store ECDFs and their difference
res <- tibble(
  x_axis_seq = x_axis_seq, # we need this as the x axis
  plank_ecdf = plank_ecdf,
  null_ecdf = null_ecdf,
  diff = plank_ecdf - null_ecdf
) %>% 
  mutate(rank = row_number()) # not necessarly needed
summary(res)

# Check that generated ECDFs are the same as the plotted ones
p1c <- res %>% 
  select(x_axis_seq, plank_ecdf, null_ecdf) %>% 
  pivot_longer(plank_ecdf:null_ecdf) %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "All distances") +
  theme_classic()

p2c <- res %>% 
  select(x_axis_seq, plank_ecdf, null_ecdf) %>% 
  pivot_longer(plank_ecdf:null_ecdf) %>% 
  filter(x_axis_seq < 20) %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 20") +
  theme_classic()

p3c <- res %>% 
  select(x_axis_seq, plank_ecdf, null_ecdf) %>% 
  pivot_longer(plank_ecdf:null_ecdf) %>% 
  filter(x_axis_seq < 10) %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 10") +
  theme_classic()

p4c <- res %>% 
  select(x_axis_seq, plank_ecdf, null_ecdf) %>% 
  pivot_longer(plank_ecdf:null_ecdf) %>% 
  filter(x_axis_seq < 5) %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 5") +
  theme_classic()

p1c + p2c + p3c + p4c
## These are different, ECDF should be generated from truncated distribution of distances?
## It’s probably bad to cut ECDF as they should go from 0 to 1, or at least we need to rescale them

p2b / p2c



## OK for consistency
## But there is a difference if we compute ECDF on truncated distributions.
# Let’s try it

df_wide_20 <- df_intra_dist_small %>% 
  filter(taxon == "Doliolida") %>% 
  filter(value < 20) %>% 
  pivot_wider(names_from = "name", values_from = "value", values_fn = list)
df_wide_20

# Compute ecdf
# Needs a nested call: ecdf(x) returns a function, x being the values of interest
# Call this function with argument x2 to get values along x2 (here, plankton distances so that we have a common x axis) 
x_axis_seq_20 <- seq(0, 20, length.out = 1000)
plank_ecdf_20 <- ecdf(unlist(df_wide_20$plankton))(x_axis_seq_20)
null_ecdf_20 <- ecdf(unlist(df_wide_20$null))(x_axis_seq_20)

res_20 <- tibble(
  x_axis_seq = x_axis_seq_20, # we need this as the x axis
  plank_ecdf = plank_ecdf_20,
  null_ecdf = null_ecdf_20
) %>% 
  mutate(diff = plank_ecdf - null_ecdf)

# Plot ECDFs
res_20 %>% 
  select(x_axis_seq, plank_ecdf, null_ecdf) %>% 
  pivot_longer(plank_ecdf:null_ecdf) %>% 
  mutate(name = ifelse(name == "null_ecdf", "null", "plankton")) %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = value, colour = name)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type", title = "Below 20") +
  theme_classic()

# Plot ECDF diff
res_20 %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = diff), colour = "#00BFC4") +
  geom_hline(yintercept = 0, colour = "#F8766D") + 
  geom_vline(xintercept = c(5, 10, 20), alpha = 0.5, linewidth = 0.5) +
  labs(x = "Distance (cm)", y = "Plankton ECDF - Null ECDF", title = "Below 20") +
  theme_classic()

## Generate ECDFs for all distance thresholds
dist_thres <- tibble(
    dist_rthr_cm = c(70, 20, 10, 5),
    name = c("All distances", "Below 20", "Below 10", "Below 5")
  )

i <- 1
res <- lapply(1:nrow(dist_thres), function(i) {
  # Get distance threshold and plot name
  r <- dist_thres %>% slice(i)
  
  # Generate the wide dataset and apply threshold
  df_wide_r <- df_intra_dist_small %>% 
    filter(taxon == "Doliolida") %>% 
    filter(value < r$dist_thr_cm) %>% 
    pivot_wider(names_from = "name", values_from = "value", values_fn = list)
  
  # Compute ECDF on a common x axis sequence
  x_axis_seq_r <- seq(0, r$dist_thr_cm, length.out = 1000)
  plank_ecdf_r <- ecdf(unlist(df_wide_r$plankton))(x_axis_seq_r)
  null_ecdf_r <- ecdf(unlist(df_wide_r$null))(x_axis_seq_r)
  
  # Store results
  res <- tibble(
    x_axis_seq = x_axis_seq_r,
    plank_ecdf = plank_ecdf_r,
    null_ecdf = null_ecdf_r
  ) %>% 
    mutate(
      diff = plank_ecdf - null_ecdf,
      dist_thr_cm = r$dist_thr_cm,
      name = r$name
    )
  return(res)
}) %>% 
  bind_rows()


res %>% 
  select(-diff) %>% 
  pivot_longer(plank_ecdf:null_ecdf, names_to = "type", values_to = "value") %>% 
  ggplot() +
  geom_path(aes(x = x_axis_seq, y = value, colour = type)) +
  labs(x = "Distance (cm)", y = "ECDF") +
  facet_wrap(~name, scales = "free") +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white"))


ggplot(res) +
  geom_path(aes(x = x_axis_seq, y = diff), colour = "#00BFC4") +
  geom_hline(yintercept = 0, colour = "#F8766D") +
  labs(x = "Distance (cm)", y = "Plankton ECDF - Null ECDF") +
  facet_wrap(~name, scales = "free") +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white"))



#p1d <- res %>% 
#  ggplot() +
#  geom_path(aes(x = plank_dist, y = diff), colour = "#00BFC4") +
#  geom_hline(yintercept = 0, colour = "#F8766D") + 
#  geom_vline(xintercept = c(5, 10, 20), alpha = 0.5, linewidth = 0.5) +
#  labs(x = "Distance (cm)", y = "Plankton ECDF - Null ECDF", title = "All distances") +
#  theme_classic()
#
#p2d <- res %>% 
#  filter(plank_dist < 20) %>% 
#  ggplot() +
#  geom_path(aes(x = plank_dist, y = diff), colour = "#00BFC4") +
#  geom_hline(yintercept = 0, colour = "#F8766D") + 
#  labs(x = "Distance (cm)", y = "Plankton ECDF - Null ECDF", title = "Below 20") +
#  theme_classic()
#
#p3d <- res %>% 
#  filter(plank_dist < 10) %>% 
#  ggplot() +
#  geom_path(aes(x = plank_dist, y = diff), colour = "#00BFC4") +
#  geom_hline(yintercept = 0, colour = "#F8766D") + 
#  labs(x = "Distance (cm)", y = "Plankton ECDF - Null ECDF", title = "Below 10") +
#  theme_classic()
#
#p4d <- res %>% 
#  filter(plank_dist < 5) %>% 
#  ggplot() +
#  geom_path(aes(x = plank_dist, y = diff), colour = "#00BFC4") +
#  geom_hline(yintercept = 0, colour = "#F8766D") + 
#  labs(x = "Distance (cm)", y = "Plankton ECDF - Null ECDF", title = "Below 5") +
#  theme_classic()
#
#p1d + p2d + p3d + p4d