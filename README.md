# PlanktonDist

Analyze distances between plankton organisms in situ at the centimeter scale and investigate potential ecological interactions.

## Overview

PlanktonDist is an R-based pipeline to study the spatial distribution of plankton. It compares observed distances to those expected under a random distribution of such organisms, infers interactions using multiple metrics, and models observed distances with a simple agent-based attraction model.

### Hypothesis

Distances between planktonic organisms carry ecological information that can be used to infer potential interactions.

### Workflow

1.  **Compare observed distances to random expectations**\
    Generate null datasets with randomly positioned organisms, compute distances under random distribution, and compare them to observed distances (overall, intra-group, and inter-group).

2.  **Compare inferred interactions to other metrics**\
    Use non-randomness of observed distances as an interaction metric and compare it to size-based (i.e. empirical) and co-occurrence metrics.

3.  **Build an agent-based model**\
    Implement a simple attraction model and optimize its parameters to reproduce observed distances.

------------------------------------------------------------------------

## Repository Structure

### Data

-   `raw/`: raw data files
-   `distances/`: saved distance calculations as Parquet files

### Scripts & Quarto Documents

**Note:**
- Files are to be run sequentially.
- `.R` files are R scripts.
- Files without an extension in this list are Quarto documents (`.qmd`) that generate `.html` reports.

#### Data cleaning & correction

-   `00.clean_data.R`: clean plankton data
-   `01a.x_axis_correction`: correct pixel deformation on x axis
-   `01b.z_axis_test`: check that the lack of information on the z axis is not a problem

#### Distances computation

-   `02a.compute_all_distances.R`: compute all distances for all images
-   `02b.compute_intra_distances.R`: compute intra-group distances for all images
-   `02c.compute_inter_distances.R`: compute inter-group distances for all images

#### Thresholding

-   `03a.threshold_all_distances`: find the threshold (max distance) for all distances
-   `03b.threshold_intra_distances`: find threshold for intra distances
-   `03c.threshold_inter_distances`: find threshold for inter distances

#### Null datasets

-   `04a.generate_null_datasets.R`: generate multiple null datasets and compute random distances
-   `04b.process_null_datasets.R`: filter random distances and perform pairwise comparisons

#### Distance analyses

-   `05a.process_all_distances.R`: compare all distances to random distances
-   `05b.process_intra_distances.R`: compare intra distances to random distances
-   `05c.process_inter_distances.R`: compare inter distances to random distances
-   `06a.process_all_distances.R`: analyze results for all distances
-   `06b.process_intra_distances.R`: analyze results for intra distances
-   `06c.process_inter_distances.R`: analyze results for inter distances
-   `07.distance_matrix`: build an association matrix based on non-random distances

#### Specific null datasets

-   `10a.null_acantharea_prep.R`: generate null datasets for Acantharea only and compute distances
-   `10b.null_acantharea_analysis.R`: analyze results for null Acantharea
-   `11a.big_null_datasets.R`: generate 3 big (10\^8 distances) null datasets
-   `11b.big_null_datasets_analysis.R`: analyze results for big null datasets

#### Agent-based model

-   `12a.ab_model_3d_set-up.R`: set-up the attraction agent-based model
-   `12b.ab_model_3d_run.R`: run the agent-based model for various parameters
-   `12c.ab_model_3d_gridsearch_analysis.R`: analyze results of gridsearch for optimal parameters
-   `12d.ab_model_3d_analysis.R`: analyze results of final agent-based model

#### Additional analyses

-   `13.recall_effect.R`: investigate effect of recall on plankton distances
-   `14.day_vs_night.R`: check whether plankton distances differ between day and night

#### Other association matrices

-   `15a.cooc_matrix_explore.R`: explore co-occurrence methods
-   `15b.cooc_matrix_run.R`: compute co-occurrence metrics
-   `16.size_matrix.R`: compute size-based (empirical) metrics
-   `17.compare_matrices.R`: compare distance-based, co-occurrence, and size-based metrics

#### Additional checks on Acantharea

-   `18.check_acantharea.R`: additional checks on Acantharea distances

#### Figures

-   `figures.R`: generate figures for presentations
-   `figures_paper.R`: generate figures for the paper

------------------------------------------------------------------------

## Setup & Usage

1.  **Install R and Quarto** if not already installed.

2.  **Install dependencies with renv**:

    ``` r
    renv::restore()
    ```

    This installs all required R packages for the project.

3.  Place raw data in `data/raw/`.

4.  Run scripts in numerical order to reproduce analyses and generate figures.
