# JARA: Just Another Red-list Assessment


### Henning Winker, Nathan Pacoureau & Richard Sherley

JARA (*Just Another Red-list Assessment*) is a Bayesian state-space trend analysis tool that is designed to objectively incorporate uncertainty into the IUCN Red Listing evaluation process. To ensure to a high degree of transparency and reproducibility, JARA has now implemented as an R package that hosted online on the global open-source platform GitHub (https://github.com/henning-winker/JARA). 
Installing JARA requires the `librabry(devtools)`, which can be install by 'install.packages('devtools')' and a R version >= 3.5. Then simply install JARA from github with the command:

`install_github("henning-winker/JARA")`

`library(JARA)`

The name ‘Just Another Red List Assessment’ is a reference to [JAGS](https://sourceforge.net/projects/mcmc-jags/) (Just Another Gibbs Sampler, Plummer, 2003), which is the software called from R to run the Bayesian state-space model application. The name reference, together with user-friendly R interface and modulated coding structure of JARA follows the example of the new open source fisheries stock assessment R Package ‘Just Another Bayesian Biomass Assessment‘ ([`JABBA`](https://github.com/jabbamodel/JABBA); [Winker et al. 2018](https://www.sciencedirect.com/science/article/pii/S0165783618300845))
  
The fully commented R JARA functions enable to easy adjustments of the input data and parameters by conservation practitioners to apply JARA to their count or relative abundance data. We provide a detailed Technical Documentation of the JARA R package in the form of a preprint draft ([Winker, Pacoureau & Sherley, 2020](https://www.biorxiv.org/content/biorxiv/early/2019/08/19/672899.full.pdf)) available on [bioRxiv](https://www.biorxiv.org/content/10.1101/672899v3.article-info). Code to reproduce the main figures in the preprint is available as [`jara_MSplots_v1.2.R`](https://github.com/Henning-Winker/JARA/blob/master/Paper/jara_MSplots_v1.2.R). The following seven datasets are made available with the JARA library and upon loading the library can be accessed via the 'jrdat' command.  They are also compiled in  [`jaradata`](https://github.com/Henning-Winker/JARA/tree/master/data) and worked examples of each can be reproduced by running [`JARA_7examples.R`](https://github.com/Henning-Winker/JARA/blob/master/JARA_7examples.R)

- [Cape gannet](https://www.tandfonline.com/doi/abs/10.2989/00306525.2019.1684396) <i>Morus capensis</i> (Breeding counts, in thousands of pairs, from the species' six breeding colonies)
- [African Penguin](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.6554) <i>Spheniscus demersus</i> (Breeding pair counts from 11 major colonies in South Africa and Namibia)
- [Cape Mountain Zebra](https://www.researchgate.net/publication/338510670_Population_trends_and_management_strategy_tools_for_Cape_Mountain_Zebra) <i>Equus zebra zebra</i> (Census data with Carrying Capacity option. Counts of individuals from 9 protected areas in South Africa)
- [Yellowspotted skate](https://github.com/Henning-Winker/JARA/tree/master/YellowspottedSkate) <i>Leucoraja wallacei</i> (South African demersal survey abundance indices)
- [Smoothhound shark](https://www.researchgate.net/publication/338491221_Assessment_of_smoothhound_shark_Mustelus_mustelus_in_South_Africa) <i>Mustelus mustelus</i> (Combining South African demersal survey indices and recreational angling survey indices from the De Hoop MPA)
- [Striped Marlin](https://www.iotc.org/documents/WPB/16/16-MLS_JABBA) <i>Tetrapturus audax</i> (Fitting multiple fisheries Catch-Per-Unit Effort (CPUE) dataset from the Indian Ocean from the 2018 IOTC Striped Marlin Assessment)
- [Atlantic Blue Marlin](https://www.iccat.int/Documents/SCRS/DetRep/BUM_SA_ENG.pdf) <i>Makaira nigricans</i> (Comparing estimated Biomass trends and fitted CPUE indices from the 2018 ICCAT stock assessment) 

The JARA framework provides the option to simultaniously fit multiple relative abundance indices and estimate the underlying mean trend (Figure 1). 

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/Indices2x1.png" width = "800" >

<i> <b> Figure 1.</b> (a) Relative abundance indices with 95% CIs for smoothhound shark <i>Mustelus mustelus</i>, depecting for abundance indices from demersal trawl surveys along the South African South Coast and one from research angling surveys conducted in the De Hoop MPA and (b) Standardized Catch-Per-Unit-Effort indices for Indian Striped marlin.
</i>.
<br />
<br />

To evaluate model fit, JARA provides the user with three plots. The first shows the observed and predicted abundance values  for each time series together with the 95% posterior predictive credibility intervals (Figure 2). 
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/trj2x1.png" width = "800" >

<i> <b> Figure 2.</b>  Time-series of observed (color coded dots) relative indices and a joint predicted (balck solid line) trend for (a) smoothhound shark and (b) Indian Ocean striped marlin. Shaded grey area indicates 95% credibility intervals.
</i>.
<br />
<br />

The second shows individual fits, as well as the 95% credible intervals (CI) derived from the observation variance  (Figure 3).
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/Fits_MLS_s1.png" width = "800" >

<i> <b> Figure 3.</b> JARA fits  to six standardized Catch-Per-Unit-Effort indices for Indian Ocean striped marlin. The solid black line is the model predicted CPUE and the circles are observed CPUE values. Grey-shaded areas denote the 95% Credibility Intervals.
</i>.
<br />
<br />

The third plot illustrates the fits on log scale. 
<br />
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/logFits_MLS_s1.png" width = "800" >

<i> <b> Figure 4.</b> JARA fits (on log-scale) to six standardized Catch-Per-Unit-Effort indices for Indian Ocean striped marlin. The solid blue line is the model predicted CPUEand the circles are observed CPUE values. Error bars denote the assumed observation variance (in 95% CIs) for the observed CPUE values.
</i>
<br />
<br />

If absolute abundance estimates are available (e.g. census data from different breeding colonies) JARA can also produce a summed population trend for the total population, where each subpopulations' trend is modelled seperately within JARA given a common process variance (Figure 5).  
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/AP2x1.png" width = "1000" >

<i> <b> Figure 5.</b> (a) Trends in estimated numbers (coloured points) and JARA fits
(coloured lines) with 95% credible intervals (grey polygons) of African penguin breeding pairs at its eleven colonies and (b) total estimated numbers for the 'global' population obtained by summing the 11 posteriors from each colony. The colored dashed lines denote 1 x generation length (GL), 2 x GL and 3 X GL, going backwards from the last year.</i>
<br />
<br />  
  
If the time series is longer than than three generation lengths (GL) the %decline is automatically estimated from the last assessment year minus 3 x GL. If the time series is shorter than 3 x GT, JARA projects forward to by passing the number of desired future years without observations to the state-space model to achive a time horizon of 3 x GL + 2 (Figure 6).

<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/SHSH_trend.png" width = "1000" >
<i> <b> Figure 6.</b> (a) Alligned Trends in relative abundance indices (coloured points) and JARA fits
(black line) with 95% credible intervals (grey polygons) of smoothhound shark and (b) projected abundance trends over 3 x GL. </i>
<br />
<br />  
  
To facilitate Red List assessment decision-making, JARA automatically produces two key plots for each input dataset showing: (1) the median and posterior probabilities for the percentage annual population change calculated from all the observed data, and from each of the most recent 1 GL, 2 GL, and 3 GL (depending on the length of the observed time-series), shown relative to a stable population (%C = 0) (Figure 7a); and (2) how the posterior distribution for the percentage change in abundance (%C) over 3 GL aligns against the thresholds for the Red List categories (Least Concern LC, Vulnerable VU, Endangered EN or Critically Endangered CR) under criteria A2–A4 (by default, Figure 7b) or under A1 (if specified by the user using criteria = "A1" in the call to jrplot_iucn). <strong>NOTE</strong> that, by default, Near Threatened NT is not included in the output because in the current version of the [IUCN Red List Categories and Criteria](https://www.iucnredlist.org/resources/categories-and-criteria) the category NT is not specified by its own criteria, but instead by the proximity of a species to the criteria for the category VU. However, as outlined in the [Guidelines for Using the IUCN Red List Categories and Criteria](https://www.iucnredlist.org/resources/redlistguidelines), a species can qualify for listing in the NT category if it is close to qualifying for the VU category. For example, a species could reasonably be listed as NT if the percentage decline in abundance over 3 GL is close to the 30% threshold for VU or if, given data uncertainties, the range of plausible categories include both LC and VU (or EN) then the taxon can be classified as NT (unless the best estimate is VU or EN). For this reason, users can opt to include NT in the JARA output (by specifying NT_opt = TRUE in the call to jrplot_iucn). This will add NT (declines of >10 and <30%) to the plot and NT will be selected as the recommended status if this category contains the highest percentage of the posterior distribution for %C. <strong>NOTE</strong> that even if NT_opt = FALSE (the default) in the call to jrplot_iucn, JARA will return NT as the recommended status in the console output if the combined percentage of the %C posterior distribution that overlap with the threatened categories (CR to VU) is more than 50% (i.e. LC contains <50% of the %C posterior distribution), but LC contains the single highest percentage of the %C posterior distribution. This is intended to indicate that, given data uncertainties, it is at least as likley that the species meets a threshold to be considered in a threatened category as it is to be LC, and to act as a flag to assessors that they may wish to consider listing as NT (or VU) (see the [Guidelines for Using the IUCN Red List Categories and Criteria](https://www.iucnredlist.org/resources/redlistguidelines) for more information on the use of NT).
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/APstatus2x1.png" width = "1000" >

<i> <b> Figure 7.</b> (a) Posterior probability densities for the annual rate of change of the combined population African Penguin over three generation lengths (3 x GL), 2 x GL, the most recent 1 x GL and for all available years (All.yrs) and (b) posterior probability density of the percentage change (C%) over 3 x GL. Note that the probabilities have become more negative in more recent 1 x GL and 2 x GL and the median rate of change over 3 x GL is estiamted −54.4%.
 </i>
<br />
<br />   

To assist with presenting the JARA outputs in scientific reports, we also provide a plotting R script [`JARAplotting.R`](https://github.com/Henning-Winker/JARA/blob/master/JARAplotting.R) that entails examples of the here presented plots . We also include R code for reproducing the JARA output format used for production of the [IUCN Redlist Supplemetary Information](https://www.iucnredlist.org/species/pdf/2903170/attachment) for recent and on-going global shark assessments by the Shark Specialist Group of the IUCN (see e.g. Figure 8).
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/YellowspottedSkate_status.png" width = "800" >
<i> <b> Figure 8.</b> JARA results for <i> Leucoraja wallacei </i> based on data from demersal trawl surveys off South Africa showing (a) the JARA fit to the observed time-series, (b) the posterior probability for the percentage annual population change calculated from all the observed data (in black), from the last 1 generation length of data (in blue) and from the last 2 generation lengths of data (in green), with the mean (solid lines) shown relative to a stable population (% change = 0, black dashed line), (c) the observed (black line) and predicted (red line) population trajectory over three generations (75 years, dashed grey lines) and (d) the median decline over three generation lengths (dashed line) and corresponding probabilities for rates of population decline falling within the IUCN threat criteria.
 </i>
<br />
<br />  


