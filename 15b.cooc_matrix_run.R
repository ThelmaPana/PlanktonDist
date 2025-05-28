#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Compute co-occurrence matrix.
# Date: 06/02/2025
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")

# Pearson correlation on log-transformed data


## Read data ----
#--------------------------------------------------------------------------#
# Read plankton and images
plankton <- read_parquet("data/00.plankton_clean.parquet") # no need to use the X correction here
images <- read_parquet("data/00.images_clean.parquet")

# List taxonomic groups
taxa <- plankton %>% select(taxon) %>% distinct() %>% pull(taxon) %>% sort()
# Drop unwanted groups
taxa <- setdiff(taxa, c("Collodaria_colonial", "Rhizaria"))
plankton <- plankton %>% filter(taxon %in% taxa)




## Compute correlation ----
#--------------------------------------------------------------------------#
# Number of repetitions
n_rep <-  100

df_corr <- mclapply(1:n_rep, function(i) {
  
  # Number of images to select
  n_sub <- 10000
  
  # Select images
  images_ssub <- images %>% 
    slice_sample(n = n_sub) %>% 
    mutate(img = paste0("img_", str_pad(row_number(), width = nchar(n_sub), pad = "0")))
  
  # Get plankton in selected images
  plankton_ssub <- plankton %>% 
    filter(img_name %in% images_ssub$img_name) %>% 
    left_join(images_ssub %>% select(img_name, img), by = join_by(img_name))
  
  ## Abundance matrix
  # Create all combination of imgs by taxon
  mat <- crossing(taxon = taxa, img = images_ssub$img) %>% 
    left_join(count(plankton_ssub, img, taxon), by = join_by(taxon, img)) %>% 
    replace_na(list(n = 0)) %>% 
    pivot_wider(names_from = img, values_from = n, values_fill = 0) %>% 
    arrange(taxon) %>% 
    select(-taxon) %>% 
    as.matrix() ## Convert abundance data to matrix form for faster access
  
  
  # Prepare storage for correlations using only taxa that are present
  df_corr <- crossing(t1 = taxa, t2 = taxa) %>% 
    filter(t1 <= t2) 
  
  # Create a named vector for taxa to index the abundance matrix
  taxa_indices <- setNames(1:nrow(mat), taxa)
  
  # Compute correlations, with pearson on log-transformed abundances
  df_corr <- df_corr %>%
    mutate(
      res = map2(t1, t2, ~get_corr_pval(.x, .y, mat, taxa_indices, method = "pearson", log_transform = TRUE)),
      corr = map_dbl(res, "corr"),
      p_val = map_dbl(res, "p_val")
    ) %>% 
    select(-res) %>% 
    mutate(idx = i) # index label, will be used to recover block / repetition
  
  
  return(df_corr)
}, mc.cores = n_cores) %>% 
  bind_rows()



# Set non significant correlations to NA and compute mean corr across blocks.
df_cooc <- df_corr %>% 
  mutate(corr = ifelse(p_val < 0.01, corr, NA)) %>% 
  group_by(t1, t2) %>% 
  summarise(
    mean_corr = mean(corr, na.rm = TRUE),
    sd_corr = sd(corr, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(df_cooc) +
  geom_tile(aes(x = t1, y = t2, fill = mean_corr)) +
  scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf", limits = c(-1, 1), na.value = "grey90") +
  labs(x = "Taxon 1", y = "Taxon 2", fill = "Mean\ncorr.") +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(df_cooc) +
  geom_tile(aes(x = t1, y = t2, fill = sd_corr)) +
  scale_fill_viridis_c(na.value = "grey90") +
  labs(x = "Taxon 1", y = "Taxon 2", fill = "Sd\ncorr.") +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Save co-occurrence matrix ----
#--------------------------------------------------------------------------#
df_cooc <- df_cooc %>% select(t1, t2, cooc_int = mean_corr)

# Save results
save(df_cooc, file = "data/15b.co_occurrence_matrix.Rdata")
