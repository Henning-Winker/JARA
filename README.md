---
title: "JARA: Just Another Redlist Assessment"
author: "Henning Winker, Nathan Pacoureau & Richard Sherley 
date: "Cape Town, 2020"
output: html_document
---

# JARA: Just Another Redlist Assessment

### Henning Winker & Richard Sherley

JARA (*Just Another Redlist Assessment*) is a Bayesian state-space trend analysis tool that is designed to objectively incorporate uncertainty into the IUCN Red Listing evaluation process. To ensure to a high degree of transparency and reproducibility, JARA has now implemented as an R package that hosted online on the global open-source platform GitHub (https://github.com/henning-winker/JARA). 
Installing JARA requires the `librabry(devtools)`, which can be install by 'install.packages('devtools')' and a R version >= 3.5. Then simply install JARA from github with the command:

`install_github("henning-winker/JARA")`

`library(JARA)`

The name ‘Just Another Red List Assessment’ is a reference to [JAGS](https://sourceforge.net/projects/mcmc-jags/) (Just Another Gibbs Sampler, Plummer, 2003), which is the software called from R to run the Bayesian state-space model application. The name reference, together with user-friendly R interface and modulated coding structure of JARA follows the example of the new open source fisheries stock assessment R Package ‘Just Another Bayesian Biomass Assessment‘ ([`JABBA`](https://github.com/jabbamodel/JABBA); [Winker et al. 2018](https://www.sciencedirect.com/science/article/pii/S0165783618300845))
  
The fully commented R JARA functions enable to easy adjustments of the input data and parameters by conservation practitioners to apply JARA to their count or relative abundance data. We provide a detailed Technical Documentation of JARA in the form of a Preprint Draft ([Winker and Sherley, 2019](https://www.biorxiv.org/content/biorxiv/early/2019/08/19/672899.full.pdf)) available on [bioRxiv](https://www.biorxiv.org/content/10.1101/672899v2.article-info). We provide an easy "test-drive" of the JARA options `JARApkg_Test.R`[blah]. Example of the following example applications can be reproduced by running [`JARA_Examples.R`](https://github.com/Henning-Winker/JARA/blob/master/JARA_Prime.v1.1.R)

- [Cape gannet](https://github.com/Henning-Winker/JARA/tree/master/Cape_gannet/output1) (Breeding pair counts from multiple colonies)  
- [African Penguin](https://github.com/Henning-Winker/JARA/tree/master/Afr_penguin/output1) (Breeding pair counts from multiple colonies)  
- [Mountain Zebra](https://github.com/Henning-Winker/JARA/tree/master/Mountain_Zebra/output1) (Census data, Carrying Capacity option)
- [Yellowspotted skate](https://github.com/Henning-Winker/JARA/tree/master/Yellowspot_skate/output1) <i>Leucoraja wallacei</i> (South African demersal survey abudance indices)
- [Smoothhound shark](https://github.com/Henning-Winker/JARA/tree/master/SmoothhoundShark/output1) <i>Mustelus mustelus</i> (Combining South African demersal survey indices and recreational angling survey indices from the De Hoop MPA)
- [Indian Ocean Striped Marlin](https://github.com/Henning-Winker/JARA/tree/master/StripedMarlin_IO_CPUE/output1) (Fitting multiple fisheries Catch-Per-Unit Effort (CPUE) from the [2018 IOTC Striped Marlin Assessment](https://www.iotc.org/documents/WPB/16/16-MLS_JABBA))  
- [2018 ICCAT Atlantic Blue Marlin stock assessment](https://www.iccat.int/Documents/SCRS/DetRep/BUM_SA_ENG.pdf) (Comparing estimated [Biomass trends](https://github.com/Henning-Winker/JARA/tree/master/Bluemarlin_ICCAT/output1) and fitted [CPUE indices](https://github.com/Henning-Winker/JARA/tree/master/Bluemarlin_Atl_CPUE/output1)) 

The JARA framework provides the option to simultaniously fit multiple relative abundance indices and estimate the underlying mean trend (Figure 1). 

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/output1/AbundanceData_SmoothhoundShark.png" width = "500" >

<i> <b> Figure 1.</b> Relative abudance abundance indices with 95% CIs for smoothhound shark <i>Mustelus mustelus</i>, depecting for abudance indices from demersal trawl surveys along the South African South Coast and one from research angling surveys conducted in the De Hoop MPA.
</i>.
<br />
<br />

To evaluate model fit, JARA provides the user with three plots. The first shows the observed and predicted abundance values  for each time series together with the 95% posterior predictive credibility intervals (Figure 2). 
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/output1/Fits_SmoothhoundShark.png" width = "500" >

<i> <b> Figure 2.</b> Five time-series of observed (color coded dots) relative indices and a joint predicted (balck solid line) trend for smoothhound shark. Shaded grey area indicates 95% credibility intervals.
</i>.
<br />
<br />

The second shows individual fits, as well as the 95% credible intervals (CI) derived from the observation variance  (Figure 3).
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/StripedMarlin_IO_CPUE/output1/logFits_StripedMarlin_IO_CPUE.png" width = "800" >

<i> <b> Figure 3.</b> JARA fits (on log-scale) to six standardized Catch-Per-Unit-Effort indices for Indian Ocean striped marlin. The solid blue line is the model predicted CPUEand the circles are observed CPUE values. Error bars denote the assumed observation variance (in 95% CIs) for the observed CPUE values.
</i>.
<br />
<br />

The third is a residual plot, which illustrates potential data conflict when fitting multiple time series (Winker et al. 2018) and includes: (i) colour-coded lognormal residuals of observed versus predicted abundance indices i, (ii) boxplots indicating the median and quantiles of all residuals available for any given year; the area of each box indicates the strength of the discrepancy between the abundance indices (larger box means higher degree of conflicting information), and (iii) a loess smoother through all residuals which can aid to identify  systematically auto-correlated residual patterns. 
<br />
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/StripedMarlin_IO_CPUE/output1/Residuals_StripedMarlin_IO_CPUE.png" width = "500" >

<i> <b> Figure 4.</b> Residual diagnostic plot for six abudance indices for Indian Ocean striped marlin. Boxplots indicate the median and quantiles of all residuals available for any given year, and solid black lines indicate a loess smoother through all residuals.
</i>.
<br />
<br />

If absolute abundance estimates are available (e.g. census data from different breeding colonies) JARA can also produce a summed population trend for the total population, where each subpopulations' trend is modelled seperately within JARA given a common process variance (Figure 5).  
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/Afr_penguin/PenguinTrends_JARA.png" width = "1000" >

<i> <b> Figure 5.</b> (a) Trends in estimated numbers (coloured points) and JARA fits
(coloured lines) with 95% credible intervals (grey polygons) of African penguin breeding pairs at its eleven colonies and (b) total estimated numbers for the 'global' population obtained by summing the 11 posteriors from each colony. The colored dashed lines denote 1 x generation length (GL), 2 x GL and 3 X GL, going backwards from the last year.</i>
<br />
<br />  
  
If the time series is longer than than three generation lengths (GL) the %decline is automatically estimated from the last assessment year minus 3 x GL. If the time series is shorter than 3 x GT, JARA projects forward to by passing the number of desired future years without observations to the state-space model to achive a time horizon of 3 x GL + 2 (Figure 6).

<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/SmoothhoundPrj_JARA.png" width = "1000" >
<i> <b> Figure 6.</b> (a) Alligned Trends in relative abundance indices (coloured points) and JARA fits
(black line) with 95% credible intervals (grey polygons) of smoothhound shark and (b) projected abundance trends over 3 x GL. </i>
<br />
<br />  
  
To facilitate Red List assessment decision-making, JARA automatically produces two key plots for each input dataset showing: (1) the median and posterior probabilities for the percentage annual population change calculated from all the observed data, and from each of the most recent 1 GL, 2 GL, and 3 GL (depending on the length of the observed time-series), shown relative to a stable population (%C = 0) (Figure 7a); and (2) how the posterior distribution for the percentage change in abundance (%C) over 3 GL aligns against the thresholds for the Red List categories (Least Concern LC, Near Threatened NT, Vulnerable VU, Endangered EN or Critically Endangered CR) under criteria A2–A4 (Figure 7b) or under A1.
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/Afr_penguin/PenguinStatus_JARA.png" width = "1000" >

<i> <b> Figure 7.</b> (a) Posterior probability densities for the annual rate of change of the combined population African Penguin over three generation lengths (3 x GL), 2 x GL, the most recent 1 x GL and for all available years (All.yrs) and (b) posterior probability density of the percentage change (C%) over 3 x GL. Note that the probabilities have become more negative in more recent 1 x GL and 2 x GL and the median rate of change over 3 x GL is estiamted −54.4%.
 </i>
<br />
<br />   

To assist with presenting the JARA outputs in scientific reports, we also provide a plotting R script [`JARA_multiplots.R`](https://github.com/Henning-Winker/JARA/blob/master/JARA_multiplots.R) that entails examples of how to combine individual JARA figures into labeled multiplots (e.g Figures 5-7). We also include R code for reproducing the JARA output format used for production of the [IUCN Redlist Supplemetary Information](https://www.iucnredlist.org/species/pdf/2903170/attachment) for recent and on-going global shark assessments by the Shark Specialist Group of the IUCN (see e.g. Figure 8).
<br />

<img src="https://github.com/Henning-Winker/JARA/blob/master/Yellowspot_skate/IUCN_SoupfinShark_JARA.png" width = "800" >
<i> <b> Figure 8.</b> JARA results for <i> Leucoraja wallacei </i> based on data from demersal trawl surveys off South Africa showing (a) the JARA fit to the observed time-series, (b) the posterior probability for the percentage annual population change calculated from all the observed data (in black), from the last 1 generation length of data (in blue) and from the last 2 generation lengths of data (in green), with the mean (solid lines) shown relative to a stable population (% change = 0, black dashed line), (c) the observed (black line) and predicted (red line) population trajectory over three generations (75 years, dashed grey lines) and (d) the median decline over three generation lengths (dashed line) and corresponding probabilities for rates of population decline falling within the IUCN threat criteria.
 </i>
<br />
<br />  


