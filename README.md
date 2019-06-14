---
title: "JARA: Just Another Redlist Assessment"
author: "Henning Winker & Richard Sherley"
date: "Cape Town/Exeter, 2019"
output: html_document
---


JARA (*Just Another Redlist Assessment*) is a Bayesian state-space trend analysis tool that is designed to to objectively incorporate uncertainty into the IUCN Red Listing evaluation process. 

To ensure to a high degree of transparency and reproducibility, JARA is hosted online on the global open-source platform GitHub (https://github.com/henning-winker/JARA). The fully commented R code can be easily modified by conservation practitioners to apply JARA to their count or relative abundance data. We include a detailed [Technical Documentation](https://github.com/Henning-Winker/JARA/blob/master/JARAmodel.pdf) of JARA and a   [Tutorial](https://github.com/Henning-Winker/JABBAbeta/blob/master/Tutorial_Vignette.md) that describes how to set up the [JARA 'Prime' file](https://github.com/Henning-Winker/JARA/blob/master/JARA_Prime.v1.1.R) fro running the following example JARA assessments:


- [Cape gannet](https://github.com/Henning-Winker/JARA/tree/master/Cape_gannet/output1) (Breeding pair counts from multiple colonies)  
- [African Penguin](https://github.com/Henning-Winker/JARA/tree/master/Afr_penguin/output1) (Breeding pair counts from multiple colonies)  
- [Mountain Zebra](https://github.com/Henning-Winker/JARA/tree/master/Mountain_Zebra/output1) (Census data, Carrying Capacity option)
- [Soupfin Shark](https://github.com/Henning-Winker/JARA/tree/master/SoupfinShark/output1) <i>Galeorhinus galeus</i> (South African demersal survey abudance indices)
- [Smoothhound shark](https://github.com/Henning-Winker/JARA/tree/master/SmoothhoundShark/output1) <i>Mustelus mustelus</i> (Combining South African demersal survey indices and recreational angling survey indices from the De Hoop MPA)
- [Indian Ocean Striped Marlin](https://github.com/Henning-Winker/JARA/tree/master/StripedMarlin_IO_CPUE/output1) (Fitting multiple fisheries Catch-Per-Unit Effort (CPUE) from the [2018 IOTC Striped Marlin Assessment](https://www.iotc.org/documents/WPB/16/16-MLS_JABBA))  
- [2018 ICCAT Atlantic Blue Marlin stock assessment](https://www.iccat.int/Documents/SCRS/DetRep/BUM_SA_ENG.pdf) (Comparing estimated [Biomass trends](https://github.com/Henning-Winker/JARA/tree/master/Bluemarlin_ICCAT/output1) and fitted [CPUE indices](https://github.com/Henning-Winker/JARA/tree/master/Bluemarlin_Atl_CPUE/output1) 


The name ‘Just Another Red List Assessment’ is a reference to JAGS (Just Another Gibbs Sampler, Plummer, 2003), which is the software called from R to run the Bayesian state-space model application. The name reference, together with user-friendly R interface and modulated coding structure of JARA follows the example of the new open source fisheries stock assessment software ‘Just Another Bayesian Biomass Assessment‘ ([JABBA](https://github.com/Henning-Winker/JABBAbeta); [Winker et al. 2018](https://www.sciencedirect.com/science/article/pii/S0165783618300845))


The JARA framework provides the option to simultaniously fit multiple relative abundance indices and estimate the underlying mean trend (Fig. 1). 

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/output1/AbundanceData_SmoothhoundShark.png" width = "500" >

<i> <b> Figure 1.</b> Relative abudance abundance indices with 95% CIs for smoothhound shark <i>Mustelus mustelus</i>, depecting for abudance indices from demersal trawl surveys along the South African South Coast and one from research angling surveys conducted in the De Hoop MPA.
</i>.
<br />

To evaluate model fit, JARA provides the user with three plots. The first shows the observed and predicted abundance values  for each time series together with the 95% posterior predictive credibility intervals (Figure 2). 

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/output1/Fits_SmoothhoundShark.png" width = "500" >

<i> <b> Figure 2.</b> Relative abudance abundance indices with 95% CIs for smoothhound shark <i>Mustelus mustelus</i>, depecting for abudance indices from demersal trawl surveys along the South African South Coast and one from research angling surveys conducted in the De Hoop MPA.
</i>.
<br />

The second shows individuals fits , as well as the 95% credible intervals (CI) derived from the observation variance  (Figure 3).

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/output1/logFits_SmoothhoundShark.png" width = "500" >

<i> <b> Figure 2.</b> JARA fits (on log-scale) to five relative abudance indices for smoothhound shark. The solid blue line is the model predicted CPUEand the circles are observed CPUE values. Error bars denote the assumed observation variance (in 95% CIs) for the observed CPUE values.
</i>.
<br />

The third is a residual plot, which illustrates potential data conflict when fitting multiple time series (Winker et al. 2018) and includes: (i) colour-coded lognormal residuals of observed versus predicted abundance indices i, (ii) boxplots indicating the median and quantiles of all residuals available for any given year; the area of each box indicates the strength of the discrepancy between the abundance indices (larger box means higher degree of conflicting information), and (iii) a loess smoother through all residuals which can aid to identify  systematically auto-correlated residual patterns. 

<img src="https://github.com/Henning-Winker/JARA/blob/master/SmoothhoundShark/output1/Fits_SmoothhoundShark.png" width = "500" >

<i> <b> Figure 4.</b> Residual diagnostic plot for five abudance indices for smoothhound shark. Boxplots indicate the median and quantiles of all residuals available for any given year, and solid black lines indicate a loess smoother through all residuals.
</i>.
<br />


If absolute abundance estimates are available (e.g. census data from different breeding colonies) JARA can also produce a summed population trend for the total population, where each subpopulations' trend is modelled seperately within JARA given a common process variance.  

<br />
<img src="https://github.com/Henning-Winker/JARA/blob/master/Afr_penguin/output1/Fits_Afr_penguin.png" width = "500" >

If the time series is longer than than three generation times (GT) the %decline is automatically estimated from the last assessment year minus 3 x GT. If the time series is shorter than 3 x GT, JARA projects forward.

<img src="https://github.com/Henning-Winker/JARA/blob/master/SMA_NAtl/output1/PopTrend_SMA_Natl.png" width = "500" >



As key ouput, JARA provides a posterior distribution for the estimated population change and produces easy to interpret graphics in which the probability of decline is displayed against the thresholds specified by the IUCN Red List criterion A2 (or optionally A1). 


<img src="https://github.com/Henning-Winker/JARA/blob/master/Afr_penguin/output1/IUCNplot_Afr_penguin.png" width = "500" >

<i> <b> JARA IUCN plot </b> showing the median decline reduction over three generation lengths (dashed line) and corresponding probabilities for rates of population decline falling within the IUCN threat criteria.Red List categories. 
</i>.
<br />





As an advanced feature, JARA provides an easy to implement means to conduct retrospective analysis by sequentially removing observations from the last year and refitting JARA. The idea is that learning from the past can aid to understand the influence of new data on our inference about a species threat status. 

<img src="https://github.com/Henning-Winker/JARA/blob/master/Cape_gannet/RetroTrends.Cape_gannet.png" width = "500" >
<img src="https://github.com/Henning-Winker/JARA/blob/master/Cape_gannet/RetroPosteriors.Cape_gannet.png" width = "500" >


All examples can be run using the development [`JARA_Prime.v1.3beta.R file`](https://github.com/Henning-Winker/JARA/blob/master/JARA_Prime.v1.3beta.R) from which the JARA model [`JARA.v1.3beta.R`](https://github.com/Henning-Winker/JARA/blob/master/JARA.v1.3beta.R) is executed.

*Still needs a lot of beautification!*..... 
