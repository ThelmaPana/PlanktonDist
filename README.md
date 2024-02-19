# PlanktonDist

Investigate what the distances between planktonic organisms can tell us about their ecology.

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

Where data lives.

### Notes

Some ideas / files to present the dataset.

### Scripts

-   `00.clean_data.R`: Read raw data, clean it, save it.
-   `01.null_hypothesis.qmd`: Generate data for the null hypothesis.
-   `02.x_axis_correction.qmd`: Generate data to correct distortion in `x` axis using Solitary Collodaria.
-   `03.overall_randomness.qmd`: Test for overall randomness of plankton distribution.
-   `04.intrataxa_distances.R`: Compute intrataxa distances for each taxa.
-   `05.intrataxa_distances_analysis.R`: Analize intrataxa distances.
-   `06.intertaxa_distances.R`: Compute intertaxa distances for all pairs of taxa.
-   `07.intertaxa_distances_analysis.R`: Analize intertaxa distances.

## Loose ideas

-   H0: distances between individuals are independent from taxonomy

-   <https://www.zoology.ubc.ca/~krebs/downloads/krebs_chapter_06_2017.pdf>

-   Multitype point pattern

-   <https://www.emilyburchfield.org/courses/gsa/point_pattern_lab>

-   <https://geographicdata.science/book/notebooks/08_point_pattern_analysis.html>
