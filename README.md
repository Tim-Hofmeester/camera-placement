# Camera-placement
Repository for data and code for the publication: Tim R. Hofmeester, Neri H. Thorsen, Joris P. G. M. Cromsigt, Jonas Kindberg, Henrik Andrén, John D. C. Linnell & John Odden. 2021. Effects of camera-trap placement and number on detection of members of a mammalian assemblage. Ecosphere 12: e03662. https://doi.org/10.1002/ecs2.3662

<b>Abstract</b>

A central goal in camera trapping (CT) studies is to maximize detection probability and precision of occupancy estimates while minimizing the number of CTs to reduce equipment and labor costs. Few studies, however, have examined the effect of CT number on detection probability. Moreover, historically, most studies focused on a specific species and the design could be tailored towards maximizing detection of this target species. Increasingly, however, such studies use data for all captured, non-target, species (“by-catch data”) for animal community-level analyses. It remains unclear if, and how, the targeting of CTs towards one species affects the detection of non-target species. We paired CTs from a permanent camera-trapping grid (with 38 CTs) targeted at monitoring Eurasian lynx (<i>Lynx lynx</i>) in Innlandet County, Norway, with additional randomly-placed CTs at two spatial scales (38 CTs within the same habitat patch and 38 CTs within the same 50 km2 grid cell as the lynx-targeted CTs) for three months. We combined multi-scale occupancy models that enable the separation of large-scale occupancy, CT-scale site use and detection probability with single-scale occupancy models. This allowed us to study the effects of targeted placement and CT number on the detection probability of the target species (lynx) and seven non-target mammal species (4 carnivores, 3 herbivores and 1 rodent). We found that all species, except moose, had the highest detection probability at lynx-targeted CTs. Moose had equal detection probabilities at all three placement types. Adding extra CTs generally increased detection probabilities. Consequently, for all species, combining a lynx-targeted CT with one or more randomly placed CTs, increased the accuracy and precision of occupancy estimates for 50 km2-grid cells compared to single CT estimates. The placement of single CTs underestimated grid-cell occupancy compared to known minimum occupancy and were similar to site-use probability estimates of multi-scale models. It is, however, uncertain to which spatial extent these site-use probabilities refer. We therefore recommend the use of multiple targeted CTs to estimate occupancy in large grid cells and to interpret occupancy estimates from single CTs as site-use of an, as of yet undefined, area surrounding the CT.

<b>Data</b>

The MSocc- files contain the detection history and covariate data as used in the manuscript. Detection histories are given in an array of three dimensions: 1) grid cells (i), 2) camera sites within grid cells in the order Targeted camera, Habitat-patch random camera, Landscape random camera, (j) and 3) detection periods of 5 days (k). See manuscript for more information.

The R-script contains the description of the models and code to run the models using the jagsUI package and different MSocc- files.

<b>Repository</b>
A secured version of the data and code as used in the published article can be found on Zenodo:

<a href="https://zenodo.org/badge/latestdoi/216810454"><img src="https://zenodo.org/badge/216810454.svg" alt="DOI"></a>
