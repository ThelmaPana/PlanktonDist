# PlanktonDist

Investigate distances between plankton organisms in situ at the centimetre scale.

### Analyses

#### Steps

To detect patterns in distances between planktonic organisms, we need a null hypothesis: "planktonic organisms are randomly distributed". Null data is generated for this hypothesis, by generating random points within a set of images, so that the number of random points in images is representative of the number of planktonic organisms in images. Distances within ISIIS  data are then compared to null data to investigate the following points.

-   **overall** distances: if overall distances differ from what expected by the null data, this might reveal global interactions / non-randomness in plankton distribution.

-   **intrataxa** distances: if distances within a taxonomic group differ from what expected by the null data, this might reveal intrataxa interactions.

-   **intertaxa** distances: if distances between a pair of taxonomic groups differ from what expected by the null data, this might reveal intertaxa interactions.

#### Elements of discussion

Things to account for:

-   biological

    -   motility & trophic status

    -   division

    -   size of organisms

-   non-biological

    -   segmentation recall: 92% overall but vary between classes

    -   classification precision & recall + interclass confusions

    -   variations in environmental conditions

    -   variations in towing speed, i.e. pixel size in `x` axis can vary. We can use round solitary blacks to detect these changes.

    -   turbulence

    -   approximation of 3D distances (imaged volume) to 2D distances (image)

## Repo organisation

### Data

Where data lives, includes:

-   `raw`: raw data

-   `distances`: saved distances as parquet files

### Scripts

#### Main steps

-   00 − Clean data

    -   `00.clean_data`: Read raw data, clean it, save it.

-   01 − Perform tests and corrections

    -   `01a.x_axis_correction`: Generate data to correct distortion in `x` axis using Solitary Collodaria and perform correction.

    -   `01b.z_axis_test`: Check the effect of the lack of information on the `z` axis.

-   02 − Compute distances

    -   `02a.compute_all_distances`: Compute all plankton distances regardless of taxonomy.

    -   `02b.compute_intra_distances`: Compute intrataxa distances.

    -   `02c.compute_inter_distances`: Compute intertaxa distances.

-   03 − Identify distance threshold

    -   `03a.threshold_all_distances`: Find optimal distance threshold using all distances.

    -   `03b.threshold_intra_distances`: Find optimal distance threshold using intrataxa distances.

    -   `03c.threshold_inter_distances`: Find optimal distance threshold using intertaxa distances.

-   04 − Prepare null datasets

    -   `04a.generate_null_datasets`: Generate multiple null datasets with same properties as ISIIS dataset.

    -   `04b.process_null_datasets`: Filter null datasets distances using the previously identified threshold and perform pairwise comparisons for various number of distances.

-   05 − Filter and process plankton distances

    -   `05a.process_all_distances`: Filter all distances using the previously identified threshold and compute Kuiper statistic.

    -   `05b.process_intra_distances`: Filter intrataxa distances using the previously identified threshold and compute Kuiper statistic.

    -   `05c.process_inter_distances`: Filter intertaxa distances using the previously identified threshold and compute Kuiper statistic.

-   06 − Analyse results

    -   `06a.analyse_all_distances`: Analyse all distances.

    -   `06b.analyse_intra_distances`: Analyse intrataxa distances.

    -   `06c.analyse_inter_distances`: Analyse intertaxa distances.

-   07 − Interaction matrix

    -   `07.interaction_matrix`: Build the interaction matrix.

#### Other stuff & checks

-   [ ] Agent-based model

    -\> Develop a 3D agent based model to reproduce observations for all distances regardless of taxonomy.

    `12a.ab_model_3d_run.R`

    `12b.ab_model_3d_analysis`

-   [ ] Recall effect

    -\> Make sure that the recall has no effect on detected patterns.

    `13.recall_effect`

-   [ ] Big null datasets

    -\> Make sure the Kuiper statistic values can be extrapolated for very large number of distances.

    `11a.big_null_datasets.R`

    `11b.big_null_datasets.analysis`

-   [ ] Null Acantharea dataset

    -\> Make sure that the global null dataset is representative of a null dataset generated for one taxonomic group only.

    `10a.null_acantharea_prep.R`

    `10b.null_acantharea_analysis`

## Loose ideas

-   H0: distances between individuals are independent from taxonomy

-   <https://www.zoology.ubc.ca/~krebs/downloads/krebs_chapter_06_2017.pdf>

-   Multitype point pattern

-   <https://www.emilyburchfield.org/courses/gsa/point_pattern_lab>

-   <https://geographicdata.science/book/notebooks/08_point_pattern_analysis.html>
