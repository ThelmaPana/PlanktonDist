#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Generate figures for the paper
# Date: 31/01/2025
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")

# PNAS figures
# Format eps or pdf.
# Images must be provided at final size, preferably 1 column
# width (8.7cm). Figures wider than 1 column should be sized
# to 11.4cm or 17.8cm wide. Numbers, letters, and symbols
# should be no smaller than 6 points (2mm) and no larger than
# 12 points (6mm) after reduction and must be consistent.


## Figure 1: Distance-based associations ----
#--------------------------------------------------------------------------#
# Load distances
load("data/06a.df_all_qt.Rdata")
load("data/06b.df_intra_qt.Rdata")
load("data/06c.df_inter_qt.Rdata")


# Replace Hydrozoa by Cnidaria
df_intra_qt <- df_intra_qt %>% 
  mutate(taxon = ifelse(taxon == "Hydrozoa", "Cnidaria", taxon)) %>% 
  arrange(taxon)

# All distances + intra
df_1a <- df_all_qt %>% 
  mutate(taxon = "all", .before = null) %>% 
  bind_rows(df_intra_qt) %>% 
  mutate(
    taxon = str_replace_all(taxon, "_", " "),
    taxon = fct_inorder(taxon)
  )

unique(df_1a$taxon)

# Define colours for each taxonomic group
cols <- df_1a %>% 
  select(taxon) %>% 
  distinct() %>% 
  mutate(
    colour = c(
      "black", # all groups
      "#004586", 
      "#FF420E", 
      "#FFD320", 
      "#579D1C", 
      "#7E0021", 
      "#83CAFF", 
      "#314004", 
      "#AECF00", 
      "#4B1F6F", 
      "#FF950E", 
      "#C5000B"
  ))


# Plot
p1a <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, colour = "grey80") +
  geom_path(data = df_1a, aes(x = x_axis, y = diff, group = taxon, colour = taxon), alpha = 0.5, linewidth = 0.4) +
  geom_path(data = df_1a %>% filter(taxon == "all"), aes(x = x_axis, y = diff)) +
  labs(x = "Distance quantiles", y = "Null - plankton quant. (cm)", colour = "Plankton group") +
  scale_colour_manual(values = cols$colour) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  guides(colour = guide_legend(
    ncol = 5, 
    override.aes = list(alpha = 1, linewidth = 0.8),
    position = "top",
    title.position = "top"
  )) +
  theme(
    legend.key.spacing.x = unit(3, "pt"),
    legend.key.spacing.y = unit(0, "pt"),
    text = element_text(size = 8),
    axis.ticks.length = unit(0.15, "cm")
  )

## Inter distances
p1b <- ggplot(df_inter_qt) +
  geom_hline(yintercept = 0, colour = "grey", linewidth = 1) +
  geom_path(aes(x = x_axis, y = diff, group = pair), alpha = 0.2, linewidth = 0.3) +
  labs(x = "Distance quantiles", y = "Null - plankton quant. (cm)") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    axis.ticks.length = unit(0.15, "cm")
  )

p1 <- (p1a + p1b) + plot_layout(nrow = 2, axis_titles = "collect_x") + plot_annotation(tag_levels = "a")
p1
ggsave(p1, filename = "figures/figure_1.pdf", width = 11.4, height = 13, units = "cm", bg = "white")
ggsave(p1, filename = "figures/figure_1.png", width = 11.4, height = 13, units = "cm", bg = "white")




# Do smaller organisms have smaller distances?
## Read plankton and images
#plankton <- read_parquet("data/00.plankton_clean.parquet") # no need to use the X correction here
#
## List taxonomic groups
#taxa <- plankton %>% select(taxon) %>% distinct() %>% pull(taxon) %>% sort()
## Drop unwanted groups
#taxa <- setdiff(taxa, c("Collodaria_colonial", "Rhizaria"))
#plankton <- plankton %>% filter(taxon %in% taxa)
#
## Convert ESD from px to mm
#plankton <- plankton %>% mutate(esd = esd * 51 / 1000)
#plankton_esd <- plankton %>% 
#  group_by(taxon) %>% 
#  summarise(median_esd = median(esd)) %>% 
#  ungroup() %>% 
#  mutate(taxon = ifelse(taxon == "Hydrozoa", "Cnidaria", taxon)) # Rename hydrozoa to cnidaria



#toto <- df_1a %>% 
#  left_join(plankton_esd, by = join_by(taxon))
#
#ggplot() +
#  geom_hline(yintercept = 0, linewidth = 1, colour = "grey80") +
#  geom_path(data = toto, aes(x = x_axis, y = diff, group = taxon, colour = median_esd), alpha = 0.5, linewidth = 0.4) +
#  geom_path(data = toto %>% filter(taxon == "all"), aes(x = x_axis, y = diff)) +
#  labs(x = "Distance quantiles", y = "Null - plankton quant. (cm)", colour = "ESD") +
#  #scale_colour_manual(values = cols$colour) +
#  scale_colour_viridis_c() +
#  scale_x_continuous(expand = c(0, 0)) +
#  theme_classic() +
#  #guides(colour = guide_legend(
#  #  ncol = 5, 
#  #  override.aes = list(alpha = 1, linewidth = 0.8),
#  #  position = "top",
#  #  title.position = "top"
#  #)) +
#  theme(
#    legend.key.spacing.x = unit(3, "pt"),
#    legend.key.spacing.y = unit(0, "pt"),
#    text = element_text(size = 8),
#    axis.ticks.length = unit(0.15, "cm")
#  )
#

## Figure 2: Agent-based model ----
#--------------------------------------------------------------------------#
load("data/12d.df_ab_qt.Rdata")
df_2 <- df_ab_qt %>% 
  group_by(type) %>% 
  reframe(dist = quantile(dist, probs = seq(0, 1, by = 0.01))) %>% # compute percentiles
  pivot_wider(names_from = "type", values_from = "dist", values_fn = list) %>% 
  unnest(c(mod, obs, null)) %>% 
  mutate(
    diff_null_mod = null - mod, # compute difference between null and model quantiles
    diff_null_obs = null - obs, # compute difference between null and observed quantiles
    diff_obs_mod = obs - mod, # compute difference between observed and model quantiles
    x_axis = rep(seq(0, 1, by = 0.01), length.out = n()) # add percentiles for X axis sequence, 
  )  %>% 
  pivot_longer(diff_null_mod:diff_obs_mod) %>% 
  mutate(name = case_when(
    name == "diff_null_mod" ~ "null - model",
    name == "diff_null_obs" ~ "null - obs",
    name == "diff_obs_mod" ~ "obs - model"
  ))

p2 <- ggplot(df_2) +
  geom_hline(yintercept = 0, linewidth = 1, colour = "grey80") +
  geom_path(aes(x = x_axis, y = value, colour = name, linetype = name)) +
  scale_color_manual(values = c("null - model" = "#9ecae1", "null - obs" = "#08519c", "obs - model" = "grey")) +
  scale_linetype_manual(values = c("null - model" = "solid", "null - obs" = "solid", "obs - model" = "dashed")) +
  #scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  labs(x = "Distance quantiles", y = "Quantile diff. (cm)", colour = "", linetype = "") +
  theme_classic() +
  guides(colour = guide_legend(position = "inside")) +
  theme(
    legend.position.inside = c(0.5, 0.6),
    legend.background = element_blank(),
    text = element_text(size = 8),
    axis.ticks.length = unit(0.15, "cm")
  )
p2
ggsave(p2, filename = "figures/figure_2.pdf", width = 11.4, height = 5, units = "cm", bg = "white")
ggsave(p2, filename = "figures/figure_2.png", width = 11.4, height = 5, units = "cm", bg = "white")


## Figure 3: Metrics comparison ----
#--------------------------------------------------------------------------#
# Read metrics values
load("data/07.distance_matrix.Rdata")
load("data/15b.co_occurrence_matrix.Rdata")
load("data/16.size_matrix.Rdata")

# Select relevant columns
df_dist <- df_mat_ks %>% select(t1, t2, dist_int = int_dist)
df_cooc <- df_cooc %>% select(t1, t2, cooc_int)

# Assemble matrices
df_mat <- df_cooc %>% 
  left_join(df_size, by = join_by(t1, t2)) %>% 
  left_join(df_dist, by = join_by(t1, t2)) %>% 
  mutate(pair = paste(t1, t2, sep = " - "), .after = t2)

# Reain pairs for which at least one value is missing but flag them
df_retain <- df_mat %>% 
  rowwise() %>%
  # if mean(cooc_int, dist_int, size_int) is NA, then one is missing
  mutate(missing = is.na(sum(c_across(cooc_int:dist_int)))) %>% 
  ungroup()

# Reshape to longer df
df_retain_long <- df_retain %>% pivot_longer(cooc_int:dist_int, names_to = "method")

# Little trick for polar plot
# Duplicate the first variable to close the circle
df_retain_1 <- df_retain_long %>% 
  arrange(method) %>% 
  mutate(
    var = case_when(
      method == "cooc_int" ~ "1",
      method == "size_int" ~ "2",
      method == "dist_int" ~ "3",
    )
  )
df_retain_2 <- df_retain_1 %>% 
  filter(var == "1") %>% 
  mutate(var = "6") %>% 
  bind_rows(df_retain_1)

# Compute Spearman correlations
# Compute correlations for each pair of metrics
corr_size_cooc <- cor.test(
  df_retain$size_int, 
  df_retain$cooc_int, 
  method = "spearman", 
  use = "complete.obs"
)
corr_size_dist <- cor.test(
  df_retain$size_int, 
  df_retain$dist_int, 
  method = "spearman", 
  use = "complete.obs"
)
corr_dist_cooc <- cor.test(
  df_retain$dist_int, 
  df_retain$cooc_int, 
  method = "spearman", 
  use = "complete.obs"
)

# Circular plot with correlations
p3 <- ggplot(df_retain_2, aes(x = var, y = value)) +
  # Ticks annotation
  annotate("text", x = c(1.03, 2.03, 2.97), y = 1, label = "1", size = 3) +
  annotate("text", x = 1, y = 1, label = "-", colour = "grey90", size = 5) +
  geom_textpath(aes(x = 2, y = 1, label = "-"), linetype = 0, colour = "grey90", size = 5) +
  geom_textpath(aes(x = 3, y = 1, label = "-"), linetype = 0, colour = "grey90", size = 5) +
  #geom_point(aes(x = var, y = 0), colour = "red") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  # Lines
  geom_line(aes(group = pair, colour = !missing, alpha = missing), linewidth = 0.1) +
  # Points on axes 
  geom_point(size = 0.5) +
  scale_x_discrete(expand = c(0, 0), breaks = c("1", "2", "3"), labels = c("1" = "Co-oc.", "2" = "Size", "3" = "Dist.")) +
  labs(colour = "All metrics", alpha = "All metrics") +
  scale_colour_manual(values = c("TRUE" = "darkblue", "FALSE" = "orange")) +
  scale_alpha_manual(values = c(0.5, 0.5)) +
  # Correlation values
  annotate("text", x = 1.5, y = 1,   label = paste0("ρ = ", round(corr_size_cooc$estimate, digits = 2)), size = 3) +
  annotate("text", x = 2.5, y = 0.8, label = paste0("ρ = ", round(corr_size_dist$estimate, digits = 2)), size = 3) +
  annotate("text", x = 3.5, y = 1,   label = paste0("ρ = ", round(corr_dist_cooc$estimate, digits = 2)), size = 3) +
  coord_radial() +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1))) 
p3
#ggsave(p3, filename = "figures/figure_3.pdf", width = 11.4, height = 9, units = "cm", bg = "white")
# ρ does not render in the pdf, convert png to pdf outside of R
ggsave(p3, filename = "figures/figure_3.png", width = 11.4, height = 7, units = "cm", bg = "white")


## Figure S2: Distance thresholds ----
#--------------------------------------------------------------------------#
## Overall threshold
load("data/03a.all_thr.Rdata")
ps2a <- ggplot(all_thr) +
  geom_vline(xintercept = all_thr %>% arrange(desc(kuiper_stat)) %>% slice(1) %>% pull(dist_thres), colour = "grey") +
  geom_point(aes(x = dist_thres, y = kuiper_stat), size = 0.8) +
  labs(x = "Distance threshold (cm)", y = "Kuiper statistic") +
  theme_classic()
ps2a


## Intra-group thresholds
load("data/03b.intra_opt.Rdata")
intra_opt <- intra_opt %>% 
  filter(n_dist > n_dist_min) %>% 
  mutate(
    # Correct Hydrozoa to Cnidaria
    taxon = ifelse(taxon == "Hydrozoa", "Cnidaria", taxon),
    # Replace underscore by space
    taxon = str_replace_all(taxon, "_", " ")
  ) %>% 
  arrange(taxon) %>% 
  # Join with colours
  left_join(cols, by = join_by(taxon))

ps2b <- ggplot(intra_opt) +
  geom_hline(yintercept = 11, colour = "grey") +
  geom_point(aes(x = log10(n_dist), y = dist_thres, colour = taxon), size = 0.8) +
  scale_colour_manual(values = intra_opt$colour) +
  labs(x = "Number of distances", y = "Identified threshold (cm)", colour = "Plankton group") +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(4, 5, 6, 7, 8), limits = c(3.9, 8.1)) +
  scale_y_continuous(breaks = c(9, 10, 11)) +
  theme_classic() +
  guides(colour = guide_legend(
    ncol = 5, 
    override.aes = list(alpha = 1, linewidth = 0.8),
    position = "top",
    title.position = "top"
  )) +
  theme(
    legend.key.spacing.x = unit(3, "pt"),
    legend.key.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = 8)
  )
ps2b
# Crustacea other?
# Because optimal threshold for Crustacea_other is 5 cm, thus with a number of retained distances < 10,000

## Inter-group thresholds
load("data/03c.inter_opt.Rdata")
inter_opt <- inter_opt %>% 
  filter(n_dist > n_dist_min)

ps2c <- ggplot(inter_opt) +
  geom_hline(yintercept = 11, colour = "grey") +
  geom_point(aes(x = log10(n_dist), y = dist_thres), size = 0.8) +
  labs(x = "Number of distances", y = "Identified threshold (cm)") +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(4, 5, 6, 7, 8), limits = c(3.9, 8.1)) +
  scale_y_continuous(breaks = seq(5, 17, by = 2)) +
  theme_classic()
ps2c

# Assemble and save
ps2 <- (ps2a + (ps2b + ps2c + plot_layout(nrow = 1, axis_titles = "collect", guides = "collect"))) + 
  plot_layout(nrow = 2) + plot_annotation(tag_levels = "a") & 
  theme(legend.position = "bottom")
ps2
ggsave(ps2, filename = "figures/figure_s2.png", width = 17.8, height = 15, units = "cm", bg = "white")



## Figure S3: NNKS explanations ----
#--------------------------------------------------------------------------#
# A: NNKS only, with null Acantharea and big null
load("data/04b.null_ks.Rdata")
load("data/10a.null_ks_n_dist_acant.Rdata")
load("data/11a.null_ks_n_dist_big.Rdata")

# Prepare data for regression lines
reg_lines <- tibble(log_n_dist = c(3.8, 9.3)) %>% 
  crossing(null_ks_rq_coef %>% filter(tau != "mean")) %>% 
  mutate(log_ks = slope * log_n_dist + intercept) %>% 
  arrange(tau)

# Redo the ribbon for the same range as for regression lines
reg_poly <- reg_lines %>% 
  mutate(order = c(1,4,2,3)) %>% 
  arrange(order) %>% 
  select(log_n_dist, y = log_ks) %>% 
  mutate(name = c("ymin", "ymax", "ymin", "ymax"))

# For Acantharea null, round number of distances for boxplot
null_ks_n_dist_acant <- null_ks_n_dist_acant %>% 
  mutate(log_n_dist = round(log_n_dist, digits = 1))

# For big null, log-transform number of distances and Kuiper statistic
null_ks_n_dist_big <- null_ks_n_dist_big %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_kuiper_stat = log10(kuiper_stat)
  )

# B: NNKS and PNKS for all distances
load("data/05a.all_distances_ks.Rdata")

# For plankton distances, log-transform number of distances and Kuiper statistic
df_all <- df_all %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_kuiper_stat = log10(kuiper_stat)
  )

# C: NNKS and PNKS for intra distances
load("data/05b.intra_distances_ks.Rdata")

# Drop Rhizaria and colonial Collodaria, replace Hydroza by Cnidaria
df_intra <- df_intra %>% 
  filter(!taxon %in% c("Rhizaria", "Collodaria_colonial")) %>% 
  mutate(taxon = ifelse(taxon == "Hydrozoa", "Cnidaria", taxon))

# For plankton distances, log-transform number of distances and Kuiper statistic
df_intra <- df_intra %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_kuiper_stat = log10(kuiper_stat)
  ) %>% 
  filter(n_dist > n_dist_min)

# Flag significant points
# Get slope and intercept of the 0.95 quantile reg
intercept <- null_ks_rq_coef %>% filter(tau == 0.95) %>% pull(intercept)
slope <- null_ks_rq_coef %>% filter(tau == 0.95) %>% pull(slope)
df_intra <- df_intra %>% 
  mutate(sig = log_kuiper_stat > slope * log_n_dist + intercept) # flag points above the 0.95 regression line

# D: NNKS and PNKS for inter distances
load("data/05c.inter_distances_ks.Rdata")

# Remove Rhizaria and colonial Collodaria
df_inter <- df_inter %>% 
  filter(!str_detect(pair, "Rhizaria")) %>% 
  filter(!str_detect(pair, "Collodaria_colonial"))

# For plankton distances, log-transform number of distances and Kuiper statistic
df_inter <- df_inter %>% 
  mutate(
    log_n_dist = log10(n_dist),
    log_kuiper_stat = log10(kuiper_stat)
  ) %>% 
  filter(n_dist > n_dist_min)

# Flag significant points
df_inter <- df_inter %>% 
  mutate(sig = log_kuiper_stat > slope * log_n_dist + intercept) # flag points above the 0.95 regression line


## Plot all
ps3a <- ggplot() + 
  geom_path(data = reg_lines, aes(x = log_n_dist, y = log_ks, group = tau), linewidth = 0.3, linetype = 2, colour = "grey20") +
  geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, group = log_n_dist, y = log_kuiper_stat), colour = "grey", width = 0.4, linewidth = 0.3, outlier.size = 0.3) +
  geom_boxplot(data = null_ks_n_dist_acant, aes(x = log_n_dist, y = log_kuiper_stat), width = 0.4, linewidth = 0.3, colour = "black") +
  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = null_ks_n_dist_big, aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "dodgerblue") +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(2, 4, 6, 8), limits = c(3.8, 9.3)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(-4, -3, -2, -1, 0), limits = c(-4.5, -0.9)) +
  labs(x = "Number of distances", y = "Kuiper statistic") +
  theme_classic()

ps3b <- ggplot() + 
  #geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, group = log_n_dist, y = log_kuiper_stat), colour = "grey", width = 0.4, linewidth = 0.3, outlier.size = 0.3) +
  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_all, aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "black") +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(2, 4, 6, 8), limits = c(3.8, 9.3)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(-4, -3, -2, -1, 0), limits = c(-4.5, -0.9)) +
  labs(x = "Number of distances", y = "Kuiper statistic") +
  theme_classic()

ps3c <- ggplot() + 
  #geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, group = log_n_dist, y = log_kuiper_stat), colour = "grey", width = 0.4, linewidth = 0.3, outlier.size = 0.3) +
  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_intra %>% filter(!sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "grey") +
  #geom_point(data = df_intra %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat, colour = taxon), size = 0.5) +
  geom_point(data = df_intra %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "black") +
  #scale_colour_manual(values = cols$colour[-1]) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(2, 4, 6, 8), limits = c(3.8, 9.3)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(-4, -3, -2, -1, 0), limits = c(-4.5, -0.9)) +
  labs(x = "Number of distances", y = "Kuiper statistic", colour = "Taxon") +
  theme_classic() 
ps3c

ps3d <- ggplot() + 
  #geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, group = log_n_dist, y = log_kuiper_stat), colour = "grey", width = 0.4, linewidth = 0.3, outlier.size = 0.3) +
  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_inter %>% filter(!sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.1, colour = "grey") +
  geom_point(data = df_inter %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.1, colour = "black") +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(2, 4, 6, 8), limits = c(3.8, 9.3)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(-4, -3, -2, -1, 0), limits = c(-4.5, -0.9)) +
  labs(x = "Number of distances", y = "Kuiper statistic") +
  theme_classic()


## Assemble
ps3 <- ps3a + ps3b + ps3c + ps3d + plot_layout(nrow = 2, ncol = 2, axis_titles = "collect") + 
  plot_annotation(tag_levels = "a") & theme(plot.tag.location = "panel", plot.tag.position = c(-0.05, 1.05), legend.position = "bottom") 
ps3
ggsave(ps3, filename = "figures/figure_s3.png", width = 17.8, height = 12, units = "cm", bg = "white")

## Second version with coloured taxa

design <- "
  12
  34
  55
"
ps3c <- ggplot() + 
  #geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, group = log_n_dist, y = log_kuiper_stat), colour = "grey", width = 0.4, linewidth = 0.3, outlier.size = 0.3) +
  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = df_intra %>% filter(!sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "grey") +
  geom_point(data = df_intra %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat, colour = taxon), size = 0.5) +
  #geom_point(data = df_intra %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "black") +
  scale_colour_manual(values = cols$colour[-1]) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(2, 4, 6, 8), limits = c(3.8, 9.3)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(-4, -3, -2, -1, 0), limits = c(-4.5, -0.9)) +
  labs(x = "Number of distances", y = "Kuiper statistic", colour = "Plankton group") +
  theme_classic() 

ps3 <- ps3a + ps3b + ps3c + ps3d + guide_area() + 
  plot_layout(design = design, guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") & 
  theme(legend.direction = "horizontal")
ps3
ggsave(ps3, filename = "figures/figure_s3.png", width = 17.8, height = 15, units = "cm", bg = "white")


# Is there a link with organisms size?

## Read plankton and images
#plankton <- read_parquet("data/00.plankton_clean.parquet") # no need to use the X correction here
#
## List taxonomic groups
#taxa <- plankton %>% select(taxon) %>% distinct() %>% pull(taxon) %>% sort()
## Drop unwanted groups
#taxa <- setdiff(taxa, c("Collodaria_colonial", "Rhizaria"))
#plankton <- plankton %>% filter(taxon %in% taxa)
#
## Convert ESD from px to mm
#plankton <- plankton %>% mutate(esd = esd * 51 / 1000)
#plankton_esd <- plankton %>% 
#  group_by(taxon) %>% 
#  summarise(median_esd = median(esd)) %>% 
#  ungroup() %>% 
#  mutate(taxon = ifelse(taxon == "Hydrozoa", "Cnidaria", taxon)) # Rename hydrozoa to cnidaria
#df_intra <- df_intra %>% 
#  left_join(plankton_esd, by = join_by(taxon))
#
#ggplot() + 
#  #geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, group = log_n_dist, y = log_kuiper_stat), colour = "grey", width = 0.4, linewidth = 0.3, outlier.size = 0.3) +
#  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
#  geom_point(data = df_intra %>% filter(!sig), aes(x = log_n_dist, y = log_kuiper_stat), size = 0.5, colour = "grey") +
#  #geom_point(data = df_intra %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat, colour = taxon), size = 0.5) +
#  geom_point(data = df_intra %>% filter(sig), aes(x = log_n_dist, y = log_kuiper_stat, colour = median_esd), size = 0.5) +
#  #scale_colour_manual(values = cols$colour[-1]) +
#  scale_colour_viridis_c() +
#  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(2, 4, 6, 8), limits = c(3.8, 9.3)) +
#  scale_y_continuous(labels = label_math(expr = 10^.x, format = force), breaks = c(-4, -3, -2, -1, 0), limits = c(-4.5, -0.9)) +
#  labs(x = "Number of distances", y = "Kuiper statistic", colour = "ESD") +
#  theme_classic() 
#
## No link between non-randomness and ESD.


## Figure S4: Agent-based model parameter choice ----
#--------------------------------------------------------------------------#
load("data/12c.best_mods.Rdata")

ps4 <- ggplot(best_mods) +
  geom_count(aes(x = sr, y = d_length_cm, color = after_stat(n), size = after_stat(n))) +
  scale_colour_viridis_c(breaks = c(1, 10, 20, 30)) +
  scale_y_discrete(drop = F) +
  scale_size_continuous(breaks = c(1, 10, 20, 30)) +
  labs(x = "Sensory radius (cm)", y = "Displacement length (cm)", colour = "Number of\nselected\nmodels", size = "Number of\nselected\nmodels") +
  guides(color = 'legend') +
  theme_classic()
ps4
ggsave(ps4, filename = "figures/figure_s4.png", width = 17.8, height = 8, units = "cm", bg = "white")


## Figure S5: Association matrices ----
#--------------------------------------------------------------------------#

# Read matrices
load("data/07.distance_matrix.Rdata")
df_dist <- df_mat_ks %>% select(t1, t2, assoc = int_dist)
load("data/15b.co_occurrence_matrix.Rdata")
df_cooc <- df_cooc %>% select(t1, t2, assoc = cooc_int)
load("data/16.size_matrix.Rdata")
df_size <- df_size %>% select(t1, t2, assoc = size_int)

# Replace Hydrozoa by Cnidaria
df_dist <- df_dist %>% 
  mutate(
    t1 = ifelse(t1 == "Hydrozoa", "Cnidaria", t1),
    t2 = ifelse(t2 == "Hydrozoa", "Cnidaria", t2)
  )
df_cooc <- df_cooc %>% 
  mutate(
    t1 = ifelse(t1 == "Hydrozoa", "Cnidaria", t1),
    t2 = ifelse(t2 == "Hydrozoa", "Cnidaria", t2)
  )
df_size <- df_size %>% 
  mutate(
    t1 = ifelse(t1 == "Hydrozoa", "Cnidaria", t1),
    t2 = ifelse(t2 == "Hydrozoa", "Cnidaria", t2)
  )

# Assemble all for common colour scale
df_all <- df_dist %>% 
  bind_rows(df_cooc) %>% 
  bind_rows(df_size)

# Get list of all taxa
taxa <- df_all %>% select(t1) %>% distinct() %>% pull(t1) %>% sort()

# Generate all combinations of taxonomic groups to plot complete matrices + taxa as factor
df_dist <- crossing(t1 = taxa, t2 = taxa) %>% 
  filter(t1 <= t2) %>% 
  left_join(df_dist, by = join_by(t1, t2)) %>% 
  mutate(
    t1 = factor(t1, levels = taxa),
    t2 = factor(t2, levels = taxa)
  )

df_cooc <- crossing(t1 = taxa, t2 = taxa) %>% 
  filter(t1 <= t2) %>% 
  left_join(df_cooc, by = join_by(t1, t2)) %>% 
  mutate(
    t1 = factor(t1, levels = taxa),
    t2 = factor(t2, levels = taxa)
  )

df_size <- crossing(t1 = taxa, t2 = taxa) %>% 
  filter(t1 <= t2) %>% 
  left_join(df_size, by = join_by(t1, t2)) %>% 
  mutate(
    t1 = factor(t1, levels = taxa),
    t2 = factor(t2, levels = taxa)
  )

# Subplots
ps5a <- ggplot(df_dist) +
  geom_raster(aes(x = t1, y = t2, fill = assoc)) +
  labs(x = "", y = "", fill = "Association") +
  scale_fill_gradient2(na.value = NA, low = "#ca0020", high = "#0571b0", limits = c(min(df_all$assoc, na.rm = TRUE), max(df_all$assoc, na.rm = TRUE))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 8)
  )

ps5b <- ggplot(df_cooc) +
  geom_raster(aes(x = t1, y = t2, fill = assoc)) +
  labs(x = "", y = "", fill = "Association") +
  scale_fill_gradient2(na.value = NA, low = "#ca0020", high = "#0571b0", limits = c(min(df_all$assoc, na.rm = TRUE), max(df_all$assoc, na.rm = TRUE))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 8)
  )

ps5c <- ggplot(df_size) +
  geom_raster(aes(x = t1, y = t2, fill = assoc)) +
  labs(x = "", y = "", fill = "Association") +
  scale_fill_gradient2(na.value = NA, low = "#ca0020", high = "#0571b0", limits = c(min(df_all$assoc, na.rm = TRUE), max(df_all$assoc, na.rm = TRUE))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 8)
  )

# Assemble all
ps5 <- ps5a + guide_area() + ps5b + ps5c +
  plot_layout(ncol = 2, guides = "collect", axes = "collect", axis_titles = "collect") + 
  plot_annotation(tag_levels = "a")
ps5
ggsave(ps5, filename = "figures/figure_s5.png", width = 17.8, height = 18, units = "cm", bg = "white")


## Figure S7: Effect of recall ----
#--------------------------------------------------------------------------#
# Recall values for taxonomic groups: "data/raw/classif_report_aft_thres.csv"
load("data/04b.null_ks.Rdata")
load("data/13.rec_dist.Rdata")

# Prepare data for regression lines
reg_lines <- tibble(log_n_dist = c(3.8, 9.3)) %>% 
  crossing(null_ks_rq_coef %>% filter(tau != "mean")) %>% 
  mutate(log_ks = slope * log_n_dist + intercept) %>% 
  arrange(tau)

# Redo the ribbon for the same range as for regression lines
reg_poly <- reg_lines %>% 
  mutate(order = c(1,4,2,3)) %>% 
  arrange(order) %>% 
  select(log_n_dist, y = log_ks) %>% 
  mutate(name = c("ymin", "ymax", "ymin", "ymax"))

# Plot
ps7 <- ggplot() +
  #geom_boxplot(data = null_ks_n_dist, aes(x = log_n_dist, y = log_kuiper_stat, group = log_n_dist), colour = "gray", outlier.shape = NA) +
  geom_polygon(data = reg_poly, aes(x = log_n_dist, y = y), alpha = 0.1) +
  geom_point(data = rec_dist, aes(x = log_n_dist, y = log_kuiper_stat, colour = rec_val), size = 0.8) +
  labs(x = "Number of distances", y = "Kuiper statistic", colour = "Points\nrecall") +
  scale_colour_viridis_c(limits = c(0, 1)) +
  scale_x_continuous(labels = label_math(expr = 10^.x, format = force), breaks = seq(2, 8, by = 2)) +
  scale_y_continuous(labels = label_math(expr = 10^.x, format = force)) +
  theme_classic()
ps7
ggsave(ps7, filename = "figures/figure_s7.png", width = 17.8, height = 8, units = "cm", bg = "white")


