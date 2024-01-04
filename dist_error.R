#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Inspect effect of projecting 3D volume onto 2D image
# Date: 04/01/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")

# points per img
pt <- 20
# number of images
n_img <- 1


## Generate random points within an image ----
#--------------------------------------------------------------------------#
# Pick random points within an image
r_points <- tibble(
  x = runif(n = pt*n_img, min = 1, max = 10240),
  y = runif(n = pt*n_img, min = 1, max = 2048),
  z = runif(n = pt*n_img, min = 1, max = 10240)
  ) %>% 
  mutate(img = ceiling(row_number() / pt))


# Plot one image
r_points %>% 
  filter(img == 1) %>% 
  ggplot() +
  geom_point(aes(x = x, y = y, colour = z)) +
  coord_fixed()


## NN distance ----
#--------------------------------------------------------------------------#
x2D <- ppp(r_points$x, r_points$y, window = owin(xrange = c(1, 10240), yrange = c(1, 2048)))
plot(x2D)
nn2D <- nndist(x2D, k = 5)
summary(nn2D)

x3D <- pp3(r_points$x, r_points$y, r_points$z, box3(xrange = c(1, 10240), yrange = c(1, 2048), zrange = c(1, 10240)))
plot(x3D)
nn3D <- nndist(x3D, k = 5)
summary(nn3D)

tibble(nn2D, nn3D) %>% 
  ggplot() +
  geom_point(aes(x = nn3D, y = nn2D)) +
  geom_abline(slope = 1, color = "red") + 
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_fixed()

cor(nn2D, nn3D, method = "spearman")
# Correlation of NN distance is very bad
# We need to go beyond NN distance

# Check how the number of considered NN affect correlation between true and perceived distance
k_nn_max <- min(100, pt)
# Initiate tibble to sort result
res <- tibble(k_nn = c(1:k_nn_max), rho = NA)
for (k_nn in c(1:k_nn_max)){
  nn2D <- nndist(x2D, k = k_nn)
  nn3D <- nndist(x3D, k = k_nn)
  res$rho[k_nn] = cor(nn2D, nn3D, method = "spearman")  
}
res %>% 
  ggplot() +
  geom_path(aes(x = k_nn, y = rho))
# For a large number of points, there seems to be an optimal number of NN to consider
# This is better than using all distances
# When there are fewer points, the correlation is much better

# How does this optimum varies according to the number of points?
res %>% arrange(desc(rho))

# Draw distributions from 1 to 50 points
# For each distribution, find tho optimal number of NN to consider





# Pick one reference point per image
ref <- tibble(
  x_ref = runif(n = n_img, min = 4096, max = 6144),
  y_ref = runif(n = n_img, min = 0, max = 2048),
  z_ref = runif(n = n_img, max = 10240),
) %>% 
  mutate(img = row_number())


# Join random and reference points
r_points <- r_points %>% left_join(ref)


# Compute true (3D) distance and perceived (2D) distance
r_points <- r_points %>% 
  mutate(
    true_dist = sqrt((x - x_ref)^2 + (y - y_ref)^2 + (z - z_ref)^2),
    per_dist = sqrt((x - x_ref)^2 + (y - y_ref)^2),
    error = abs((true_dist-per_dist)/per_dist)
  ) 


# Plot distances
r_points %>% 
  pivot_longer(true_dist:per_dist, names_to = "type", values_to = "dist") %>% 
  ggplot() +
  geom_density(aes(x = dist, color = type)) +
  facet_wrap(~img) + 
  theme_minimal() 

r_points %>% 
  pivot_longer(true_dist:per_dist, names_to = "type", values_to = "dist") %>% 
  ggplot() +
  geom_point(aes(x = type, y = dist, color = type)) +
  geom_segment(aes(x = "per_dist", y = per_dist, xend = "true_dist", yend = true_dist), alpha = 0.1, data = r_points) +
  facet_wrap(~img) + 
  theme_minimal()

# Plot error distribution
r_points %>% 
  ggplot() +
  geom_density(aes(x = error, color = factor(img))) +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal()

# Plot error distribution
r_points %>% 
  ggplot() +
  geom_density(aes(x = error, color = factor(img))) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 5)) +
  theme_minimal()

# Plot error distribution
r_points %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(img), y = error, color = factor(img)))


# Average error per image
r_img <- r_points %>% 
  group_by(img) %>% 
  summarise(
    true_dist = mean(true_dist),
    per_dist  = mean(per_dist),
    avg_error = mean(error), # average all errors
    med_error = median(error) # median of all errors
  ) %>% 
  ungroup() %>% 
  mutate(error_avg = abs((true_dist-per_dist)/per_dist)) # error from average distance
r_img



r_points %>% 
  ggplot(aes(x = per_dist, y = true_dist)) +
  geom_point(alpha = 0.2) + 
  geom_smooth()


formula <- y ~ poly(x, 2, raw = TRUE)

r_points %>%   
  ggplot(aes(x = per_dist, y = true_dist)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm", formula=formula) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula=formula) +
  labs(y = "True distance (px)", x = "Perceived distance (px)") +
  theme_minimal()

r_points %>%   
  ggplot(aes(x = true_dist, y = per_dist)) +
  geom_point(alpha = 0.2) +
  #stat_smooth(method = "lm") +
  labs(x = "True distance (px)", y = "Perceived distance (px)") +
  coord_fixed() +
  theme_minimal()

r_points %>%   
  ggplot(aes(x = true_dist, y = per_dist)) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) + 
  scale_fill_viridis_c() +
  #stat_smooth(method = "lm") +
  labs(x = "True distance (px)", y = "Perceived distance (px)") +
  coord_fixed() +
  theme_minimal()

r_points %>%   
  ggplot(aes(x = true_dist, y = per_dist)) +
  geom_density_2d_filled(alpha = 0.8) +
  geom_point(alpha = 0.05, size = 0.5) +
  #stat_smooth(method = "lm") +
  labs(x = "True distance (px)", y = "Perceived distance (px)") +
  coord_fixed() +
  theme_minimal()


r_points <- r_points %>% 
  mutate(
    true_dist_mm = true_dist * 51 / 1000,
    per_dist_mm = per_dist * 51 / 1000
  )


c(min(r_points$true_dist_mm), max(r_points$true_dist_mm))

# Test normality
r_points_sample <- r_points %>% slice_sample(n = 5000)
shapiro.test(r_points_sample$per_dist)
shapiro.test(r_points_sample$true_dist)
# normality not verified

rho = cor(r_points$per_dist, r_points$true_dist, method = "spearman")
?format
format(rho, digits = 2)

po <- r_points %>%   
  ggplot(aes(x = true_dist_mm, y = per_dist_mm)) +
  geom_point(alpha = 0.05, size = 0.5) +
  geom_density_2d(color = rgb(0.06,0.41,0.55)) +
  #geom_abline(slope = 1, intercept = 0) +
  #stat_smooth(method = "lm") +
  labs(x = "True distance (mm)", y = "Perceived distance (mm)") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  annotate("text", x = 60, y = 275, label = paste0("ρ = ", format(rho, digits = 2))) +
  #coord_fixed() +
  theme_minimal() +
  theme(text = element_text(size = 10))

d1 <- r_points %>% 
  ggplot() +
  geom_histogram(aes(x = true_dist_mm), binwidth = 10, fill = "white", color = "black") +
  theme_void() +
  #coord_fixed() +
  scale_x_continuous(expand = c(0,0))

d2 <- r_points %>% 
  ggplot() +
  geom_histogram(aes(x = per_dist_mm), binwidth = 10, fill = "white", color = "black") +
  scale_x_continuous(expand = c(0,0)) + 
  theme_void() +
  #coord_fixed() +
  coord_flip()

p <- d1 + plot_spacer() + po + d2 + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(1, 4),
    byrow = TRUE
  )
p
#ggsave(plot = p, filename = "figures/true_vs_per_dist.png", width = 160, height = 100, unit = "mm")

