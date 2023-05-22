# Early-melting snowpatch plant communities are transitioning into novel ecosystems
John Morgan and Zac Walker

## About
This repository contains the code and data associated with an unpublished manuscript. The paper assesses vegetation change over a ~40 yr period in early-melting snowpatches in southern Australia. These analyses were done using R version 4.2.3.

## Abstract
Snowpatch plant community distribution and composition is strongly tied to the duration of long-lasting snow cover in alpine areas; they are vulnerable to global climatic changes that result in warmer temperatures and longer growing seasons. We used a multi-decadal revisitation study to quantify early-melting snowpatch floristic and functional diversity change in southern Australia, and document shrub size-class distributions over time to detect evidence for their encroachment patterns into snowpatches, a key prediction with climatic change. Early-melting snowpatches have declined in areal extent, changed in species composition, shrub and tussock grass cover has increased, but the effects of these changes on functional trait diversity were less clear. Species gains accounted for most of the floristic change observed over time. Shrub recruitment into early-melting snowpatches is ongoing by most shrub species. Biotic differentiation is occurring; many early-melting snowpatches are transitioning to a novel community with changed composition and taller vegetation structure, but there is little evidence for species loss having occurred. Given enough time, however, the long-term loss of snowpatch species is likely (i.e. biotic homogenisation) if taller shrubs outcompete short-statured snowpatch species. We believe that the changes observed provide evidence that this alpine ecosystem is on the verge of collapse.


## Description of the data and file structure.
There are six csv data files and four R scripts.

The 'snowpatch_diversity.R' code is used for the diversity analyses assessing differences between 1982 snowpatch floristics and 2022 low shrub cover plots and 2022 high shrub cover plots. This includes calculating and plotting hill numbers (q = 0,1,2) and measures of beta-diversity. This code uses snowpatch_df.csv.
The 'snowpatch_floristics.R' code is used for producing the NMDS ordination. This code uses snowp_floristics.csv and snowp_additional.csv.
The 'snowapatch_functrait.R' code is used for the functional trait analyses, which includes producing functional trait diversity metrics. The community weighted means of each trait are also calculated for 1982, 2022 low shrub cover and 2022 high shrub cover plots. This code uses snowp_floristics.csv, snowp_additional.csv, snowp_traits.csv.
The 'snowpatch_shrubs.R' code is used to create size-class distribution histograms of main shrub species and produce condit curves. This code uses snowpatch_shrubdata.csv and snowpatch_condit.csv.

