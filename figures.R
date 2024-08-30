#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Generate all figures for presentation & paper
# Date: 01/08/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#


## Loading and set-up ----
#--------------------------------------------------------------------------#
source("utils.R")
source("utils_ab_model.R")
# All distances
load("data/03c.all_distances_processed.Rdata")
# Intra distances
load("data/04c.df_intra_scores_small.Rdata")
# Inter distances
load("data/05c.df_inter_scores_small.Rdata")
# Simulation results
load("data/simulation_results.Rdata")



## Illustration of distribution of distances ----
#--------------------------------------------------------------------------#
# Plankton distances all
df_all_dist_small %>% 
  filter(type == "data") %>% 
  ggplot() +
  geom_density(aes(x = dist), colour = "#00B2FF") +
  labs(x = "Distance (cm)", y = "Density") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/plankton_dist.png", width = 60, height = 35, units = "mm")

# Plankton distances small
df_all_dist_small %>% 
  filter(type == "data") %>% 
  ggplot() +
  geom_density(aes(x = dist)) +
  labs(x = "Distance (cm)", y = "Density")

# Random distances
df_all_dist_small %>% 
  filter(type == "rand") %>% 
  ggplot() +
  geom_density(aes(x = dist), colour = "grey") +
  labs(x = "Distance (cm)", y = "Density") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/null_dist.png", width = 60, height = 35, units = "mm")



## Kuiper stat nor null data ----
#--------------------------------------------------------------------------#
ggplot() +
  geom_boxplot(data = f_val_dist, aes(x = log_n_dist, y = log_test_stat, group = log_n_dist), colour = "gray", linewidth = 0.2) +
  geom_polygon(data = rib_data, aes(x = log_n_dist, y = y), alpha = 0.1) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = seq(2, 8, by = 2)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
  labs(x = "Number of distances", y = "Measurement of randomness / non randomness") +
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/kt_null.png", width = 110, height = 70, units = "mm")



## Kuiper stat for all distances ----
#--------------------------------------------------------------------------#
ggplot() +
  geom_boxplot(data = f_val_dist, aes(x = log_n_dist, y = log_test_stat, group = log_n_dist), colour = "gray", linewidth = 0.2) +
  geom_polygon(data = rib_data, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_all, aes(x = log_n_dist_small, y = log10(test_stat_small)), colour = "#00B2FF", size = 0.5) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = seq(2, 8, by = 2)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
  labs(x = "Number of distances", y = "Measurement of randomness / non randomness") +
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/kt_all_dist.png", width = 110, height = 70, units = "mm")



## Distribution & ECDF for all distances ----
#--------------------------------------------------------------------------#
df_all_dist_small %>% 
  ggplot() +
  geom_density(aes(x = dist, colour = type)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  theme_classic() +
  theme(text = element_text(size = 8))
#ggsave(filename = "figures/all_dist.png", width = 120, height = 60, units = "mm")

# Difference in distribution
tibble(
  x = density(df_all_dist_small %>% filter(type == "data") %>% pull(dist), cut = 0)$x,
  data = density(df_all_dist_small %>% filter(type == "data") %>% pull(dist), cut = 0)$y,
  rand = density(df_all_dist_small %>% filter(type == "rand") %>% pull(dist), cut = 0)$y
) %>% 
  mutate(diff = data - rand) %>% 
  ggplot() +
  geom_path(aes(x = x, y = diff), colour = "#00B2FF") +
  geom_hline(yintercept = 0, colour = "grey") +
  labs(x = "Distance (cm)", y = "Plankton density - null density") +
  theme_classic()

df_all_dist_small %>% 
  ggplot() +
  stat_ecdf(aes(x = dist, colour = type)) +
  labs(x = "Distance (cm)", y = "ECDF", colour = "Type") +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001)) +
  theme_classic()


## Kuiper stat for intra distances ----
#--------------------------------------------------------------------------#
# Number of colour and shapes we need
n_col_intra <- sum(df_intra_scores$above)
cols_intra <- tibble(
  colour = color("discrete rainbow")(n_col_intra),
  shape = rep(21:25, 6)[1:n_col_intra]
)

ggplot() +
  geom_boxplot(data = f_val_dist, aes(x = log_n_dist, y = log_test_stat, group = log_n_dist), colour = "gray", outlier.shape = NA, linewidth = 0.2) +
  geom_polygon(data = rib_data, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_intra_scores %>% filter(above), aes(x = log_n_dist_small, y = log_test_stat, fill = taxon, colour = taxon, shape = taxon), size = 0.5) +
  geom_point(data = df_intra_scores %>% filter(!above), aes(x = log_n_dist_small, y = log_test_stat), colour = "gray", size = 0.5) +
  scale_fill_manual(values = cols_intra$colour) +
  scale_colour_manual(values = cols_intra$colour) +
  scale_shape_manual(values = cols_intra$shape) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = seq(2, 8, by = 2)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
  labs(x = "Number of distances", y = "Measurement of randomness / non randomness", colour = "Taxon", fill = "Taxon", shape = "Taxon") +
  theme_classic() +
  guides(
    fill = guide_legend(byrow = TRUE, ncol = 2), 
    colour = guide_legend(byrow = TRUE, ncol = 2), 
    shape = guide_legend(byrow = TRUE, ncol = 2)
  ) +
  theme(text = element_text(size = 8), legend.spacing.y = unit(0.001, "lines"), legend.text = element_text(size = 6))
ggsave(filename = "figures/kt_intra_dist.png", width = 120, height = 80, units = "mm")




## Distribution & ECDF for intra distances ----
#--------------------------------------------------------------------------#
# Distribution for Acantharea, Cop_small, Cop_calanoida and Cop_other
df_intra_dist_small %>% 
  filter(taxon %in% c("Acantharea", "Cop_small", "Cop_Calanoida", "Cop_other")) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  facet_wrap(~taxon) +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white")) +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/intra_dist.png", width = 130, height = 70, units = "mm")

# Distribution for all above
df_intra_dist_small %>% 
  filter(taxon %in% (df_intra_scores %>% filter(above) %>% pull(taxon))) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name), linewidth = 0.3) +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  facet_wrap(~taxon, scales = "free") +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white"))

# ECDF for all above
df_intra_dist_small %>% 
  filter(taxon %in% (df_intra_scores %>% filter(above) %>% pull(taxon))) %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name), linewidth = 0.3) +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  facet_wrap(~taxon, scales = "free") + 
  theme_classic() +
  theme(strip.background = element_rect(colour = "white"))


## Direction of variation for intra distances ----
#--------------------------------------------------------------------------#


## Kuiper stat for inter distances ----
#--------------------------------------------------------------------------#

# Number of colour and shapes we need
n_col_inter <- sum(df_inter_scores$above)
cols_inter <- tibble(
  colour = color("smooth rainbow")(n_col_inter, range = c(0.25, 1)),
  shape = rep(21:25, 20)[1:n_col_inter]
) %>% 
  slice_sample(prop = 1)

# With all taxa
ggplot() + 
  geom_boxplot(data = f_val_dist, aes(x = log_n_dist, y = log_test_stat, group = log_n_dist), colour = "gray", linewidth = 0.2) +
  geom_polygon(data = rib_data, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_inter_scores %>% filter(!above), aes(x = log_n_dist_small, y = log_test_stat), colour = "grey", size = 0.5) +
  geom_point(data = df_inter_scores %>% filter(above), aes(x = log_n_dist_small, y = log_test_stat, colour = pair, fill = pair, shape = pair), size = 0.5, show.legend = F) +
  scale_fill_manual(values = cols_inter$colour) +
  scale_colour_manual(values = cols_inter$colour) +
  scale_shape_manual(values = cols_inter$shape) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = seq(2, 8, by = 2)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
  labs(x = "Number of distances", y = "Measurement of randomness / non randomness") +
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/kt_inter_dist.png", width = 110, height = 70, units = "mm")

  
## Distribution & ECDF for inter distances ----
#--------------------------------------------------------------------------#
# Plot density for 4 pairs
df_inter_dist_small %>% 
  filter(pair %in% c("Acantharea - Appendicularia", "Cop_Oithona - Ctenophora", "Cop_Harpacticoida - Cop_other", "Cop_Calanoida - Doliolida")) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name)) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  facet_wrap(~pair) +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white")) +
  theme(text = element_text(size = 8))
ggsave(filename = "figures/inter_dist.png", width = 130, height = 70, units = "mm")

# Plot density for 20 random pairs
pairs_to_plot <- df_inter_scores %>% filter(above) %>% slice_sample(n = 20) %>% pull(pair)
df_inter_dist_small %>% 
  filter(pair %in% pairs_to_plot) %>% 
  ggplot() +
  geom_density(aes(x = value, colour = name), linewidth = 0.3) +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  facet_wrap(~pair, scales = "free") +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white"), strip.text = element_text(size = 6))

# Plot ECDF for the same pairs
df_inter_dist_small %>% 
  filter(pair %in% pairs_to_plot) %>% 
  ggplot() +
  stat_ecdf(aes(x = value, colour = name), linewidth = 0.3) +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type") +
  facet_wrap(~pair, scales = "free") +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white"), strip.text = element_text(size = 6))



## Direction of variation for inter distances ----
#--------------------------------------------------------------------------#

## Interaction matrix ----
#--------------------------------------------------------------------------#
# Intra
df_intra_mat <- df_intra_scores %>% 
  select(t1 = taxon, t2 = taxon, test_stat, z_score, dir, above)

# Inter
df_inter_mat <- df_inter_scores %>% 
  select(pair, test_stat, z_score, dir, above) %>% 
  separate_wider_delim(pair, delim = " - ", names = c("t1", "t2"))
# and generate reverse of t1 / t2
df_inter_mat_bis <- df_inter_mat %>% rename(t2 = t1, t1 = t2)

# Store all together
df_pairs <- df_intra_mat %>% 
  bind_rows(df_inter_mat) %>% 
  bind_rows(df_inter_mat_bis)

# Create all combination of taxa and join scores
taxa <- sort(unique(df_pairs$t1))
df_mat <- crossing(t1 = taxa, t2 = taxa) %>% 
  left_join(df_pairs %>% filter(above), by = join_by(t1, t2))

# Plot Kuiper stat and Z-score
ggplot(df_mat) +
  geom_raster(aes(x = t1, y = t2, fill = test_stat, alpha = z_score)) +
  scale_fill_viridis_c(na.value = NA) +
  labs(fill = "Interaction\nstrength", alpha = "Interaction\nsignificance") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + 
  coord_fixed() +
  theme(
    text = element_text(size = 6),
    legend.key.width = unit(1, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )
ggsave(filename = "figures/matrix.png", width = 100, height = 100, units = "mm")



## Model set-up ----
#--------------------------------------------------------------------------#
set.seed(13)
# Generate an image with organisms in 3D and plot it in 2D
df <- tibble(
  x = round(runif(n = 15, min = 1, max = vol$x)),
  y = round(runif(n = 15, min = 1, max = vol$y)),
  z = round(runif(n = 15, min = 1, max = vol$z))
) %>% # Add information for img name
  mutate(
    img_name = "img_001",
    id = paste0(str_pad(row_number(), 3, pad = "0"))
  )

library(plotly)

# Show points in 3D
plot_ly(df, x = ~x, y = ~y, z = ~z) %>% 
  add_markers(  
    marker = list(
      color = 'black',
      size = 3)
  ) %>%
  layout(scene = list(
    aspectmode = 'data', 
    camera = list(
      eye = list(x = 0, y = 1, z = -3),
      center = list(x = 0, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1)
    ))
    )


# Show points in 2D
ggplot(df) +
  geom_point(aes(x, y), show.legend = F) +
  #scale_colour_gradient(high = "grey90", low = "black") +
  scale_x_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$x), expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$y), expand = c(0, 0), breaks = NULL, name = NULL) +
  theme_classic() +
  coord_fixed()

# Show points with arrows
arrows <- crossing(
  # generate all points combination
  id = df$id,
  id_end = df$id
) %>% 
  # remove self distance
  filter(id != id_end) %>% 
  # join with starting points
  left_join(df %>% select(id, x, y), by = join_by(id)) %>% 
  # join with end points
  left_join(df %>% select(id, x, y) %>% rename_with(~paste0(.x, "_end")), by = join_by(id_end)) %>% 
  # keep only small distances
  mutate(
    dist_px = sqrt((x_end - x)^2 + (y_end - y)^2),
    dist_cm = dist_px * 51 / 10000
  ) %>% 
  filter(dist_cm < 10)


ggplot(df) +
  geom_point(aes(x, y), show.legend = F) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = x_end, yend = y_end), colour = "grey", alpha = 0.5) +
  scale_x_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$x), expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$y), expand = c(0, 0), breaks = NULL, name = NULL) +
  theme_classic() +
  coord_fixed()




## Density
# Compute density
h <- 400000
dens_3d <- calculate_density_3d(df, h = h, vol = vol)

# Extract results
dens_to_plot <- dens_3d$eval.points %>% 
  as_tibble() %>% 
  mutate(dens = dens_3d$estimate) %>% 
  group_by(x, y) %>% 
  summarise(dens = mean(dens), .groups = "drop") %>% 
  mutate(h = h)

# Plot density
ggplot(dens_to_plot) +
  geom_raster(aes(x, y, fill = dens), show.legend = F) +
  geom_point(data = df, aes(x, y)) +
  #scale_fill_distiller(palette = "Oranges", direction = 1) +
  scale_fill_gradient(low = "#ffffff", high = "#FFEF68") +
  scale_x_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$x), expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$y), expand = c(0, 0), breaks = NULL, name = NULL) +
  theme_classic() +
  coord_fixed()


## Plot displacement
# First we need to compute gradient
grad <- calculate_density_gradient_3d(dens_3d)
# And extract it at points of interest
grad_points <- extract_gradient(df, gradient = grad, kde = dens_3d)
grad_points

d_length <- 300
grad_points <- tibble(
  dx = grad_points$dx,
  dy = grad_points$dy,
  dz = grad_points$dz
) %>% 
  # compute magnitude
  mutate(mag = sqrt(dx^2 + dy^2 + dz^2)) %>% 
  # normalise
  mutate(
    dx_norm = (dx / mag) * d_length,
    dy_norm = (dy / mag) * d_length,
    dz_norm = (dz / mag) * d_length
  )
# Store gradients with points
df <- df %>% bind_cols(grad_points)

# Plot displacement arrows
ggplot(df) +
  geom_raster(data = dens_to_plot, aes(x, y, fill = dens), show.legend = F) +
  geom_point(aes(x, y)) +
  geom_segment(aes(x, y, xend = x + dx_norm, yend = y + dy_norm), colour = "#00B2FF", arrow = arrow(length = unit(0.03, "npc"))) +
  scale_fill_gradient(low = "#ffffff", high = "#FFEF68") +
  scale_x_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$x), expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$y), expand = c(0, 0), breaks = NULL, name = NULL) +
  theme_classic() +
  coord_fixed()
  
# Plot new points
ggplot(df) +
  #geom_raster(data = dens_to_plot, aes(x, y, fill = dens), show.legend = F) +
  geom_point(aes(x + dx_norm, y + dy_norm), colour = "#00B2FF") +
  #geom_segment(aes(x, y, xend = x + dx_norm, yend = y + dy_norm), colour = "red", arrow = arrow(length = unit(0.03, "npc"))) +
  #scale_fill_gradient(low = "#ffffff", high = "#FFEF68") +
  scale_x_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$x), expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$y), expand = c(0, 0), breaks = NULL, name = NULL) +
  theme_classic() +
  coord_fixed()

# Plot new_points with arrows
arrows_after <- crossing(
  # generate all points combination
  id = df$id,
  id_end = df$id
) %>% 
  # remove self distance
  filter(id != id_end) %>% 
  # join with starting points
  left_join(df %>% mutate(x = x + dx_norm, y = y + dy_norm) %>% select(id, x, y), by = join_by(id)) %>% 
  # join with end points
  left_join(df %>% mutate(x = x + dx_norm, y = y + dy_norm) %>% select(id, x, y) %>% rename_with(~paste0(.x, "_end")), by = join_by(id_end)) %>% 
  # keep only small distances
  mutate(
    dist_px = sqrt((x_end - x)^2 + (y_end - y)^2),
    dist_cm = dist_px * 51 / 10000
  ) %>% 
  filter(dist_cm < 10)

ggplot(df) +
  geom_point(aes(x + dx_norm, y + dy_norm), colour = "#00B2FF") +
  geom_segment(data = arrows_after, aes(x = x, y = y, xend = x_end, yend = y_end), colour = "#67D0FD", alpha = 0.5) +
  scale_x_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$x), expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_continuous(sec.axis = dup_axis(breaks = NULL, name = NULL), limits = c(0, vol$y), expand = c(0, 0), breaks = NULL, name = NULL) +
  theme_classic() +
  coord_fixed()


## Model VS observations ----
#--------------------------------------------------------------------------#
# The parameter set we want
param_set <- list(
  d_length = 100,
  h = 250000,
  prop_mv = 0.8
)


# Get corresponding model simulation
one_sim_res <- sim_res %>% 
  # keep only distances below threshold
  filter(dist < dist_thr_cm) %>% 
  # Keep only results for this parameter set
  filter(
    d_length == param_set$d_length & 
      h == param_set$h & 
      prop_mv == param_set$prop_mv
  ) %>% 
  select(type = when, dist) %>% 
  mutate(type = ifelse(type == "before", "rand", "data")) %>% 
  mutate(from = "Model")

# Assemble and plot, model only
bind_rows(df_all_dist_small %>% mutate(from = "Obs"), one_sim_res) %>% 
  filter(from == "Model") %>% 
  ggplot() +
  geom_density(aes(x = dist, colour = type, linetype = from), linewidth = 0.3) +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", linetype = "From") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "top")
ggsave(filename = "figures/model_dist.png", width = 80, height = 60, units = "mm")


# Assemble and plot
bind_rows(df_all_dist_small %>% mutate(from = "Obs"), one_sim_res) %>% 
  ggplot() +
  geom_density(aes(x = dist, colour = type, linetype = from), linewidth = 0.3) +
  scale_colour_manual(values = c("#00B2FF", "grey"), labels = c("Plankton", "Null")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Distance (cm)", y = "Density", colour = "Type", linetype = "From") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.001), breaks = c(0, 0.05, 0.10)) +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "right", legend.key.size = unit(0.5, "line"))
ggsave(filename = "figures/model_obs_dist.png", width = 80, height = 60, units = "mm")
