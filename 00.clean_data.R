#--------------------------------------------------------------------------#
# Project: PlanktonDist
# Script purpose: Clean data
# Date: 26/01/2024
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

source("utils.R")


## Load raw data ----
#--------------------------------------------------------------------------#
plankton <- read_parquet("data/raw/plankton_predictions.parquet")
cops <- read_parquet("data/raw/copepods.parquet")
images <- read_parquet("data/raw/predicted_images.parquet")


## Refine copepods taxonomy ----
#--------------------------------------------------------------------------#
plankton <- plankton %>% 
  # keep only relevant columns
  select(
    transect, img_name, object_id = id, taxon,
    datetime, lat, lon, dist, depth,
    x = centroid_1, y = centroid_0,
  ) %>% 
  # join with cops
  left_join(cops %>% select(object_id, new_taxon = taxon), by = join_by(object_id)) %>% 
  # keep only new taxonomy if present
  mutate(taxon = ifelse(!is.na(new_taxon), new_taxon, taxon)) %>% 
  select(-new_taxon)

# Read clean taxa names
clean_taxa <- read_csv("data/raw/taxa_list.csv", show_col_types = FALSE) %>% select(taxon = orig_name, new_taxon = new_name)
taxa <- clean_taxa %>% drop_na(new_taxon) %>% pull(new_taxon) %>% sort()

# Add clean taxa names to plankton table
plankton <- plankton %>% 
  left_join(clean_taxa, by = join_by(taxon)) %>% 
  select(-taxon) %>% 
  rename(taxon = new_taxon) %>% 
  # drop objects that should be ignored
  drop_na(taxon)


## Keep only images with more than 1 object ----
#--------------------------------------------------------------------------#
# names of images with more than 1 objects
img_names <- count(plankton, img_name) %>% filter(n > 1) %>% pull(img_name)

# keep only objects in these images
plankton <- plankton %>% filter(img_name %in% img_names)
images <- images %>% filter(img_name %in% img_names)


## Save clean data ----
#--------------------------------------------------------------------------#
write_parquet(plankton, sink = "data/00.plankton_clean.parquet")
write_parquet(images, sink = "data/00.images_clean.parquet")


## Save a subsample ----
#--------------------------------------------------------------------------#
# Generate a subsample of 100 images for tests
n_img <- 100
# subsample images
images_sub <- images %>% slice_sample(n = n_img) 
# keep plankton in sampled images
plankton_sub <- plankton %>% filter(img_name %in% images_sub$img_name)

# Save
save(plankton_sub, images_sub, file = "data/00.subsample.Rdata")
