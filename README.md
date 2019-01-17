# JARA: Just Another Redlist Assessment

JARA (*Just Another Redlist Assessment*) builds on a Bayesian state-space framework with the aim to objectively incorporate uncertainty into the IUCN Red Listing evaluation process. As key ouput, JARA provides a posterior distribution for the estimated population change and produces easy to interpret graphics in which the probability of decline is displayed against the thresholds specified by the IUCN Red List criterion A2 (A1 option still to be implemented). 

To ensure to a high degree of transparency and reproducibility, JARA is hosted online on the global open-source platform GitHub (https://github.com/henning-winker/JARA), where we provide fully commented code that can be easily modified by conservation practitioners to apply JARA to their count or relative abundance data. The name ‘Just Another Red List Assessment’ is a reference to JAGS (Just Another Gibbs Sampler, Plummer, 2003), which is the software called from R to run the Bayesian state-space model application. The name reference, together with user-friendly R interface and modulated coding structure of JARA follows the example of the new open source fisheries stock assessment software ‘Just Another Bayesian Biomass Assessment‘ ([JABBA](https://github.com/Henning-Winker/JABBAbeta); [Winker et al. 2018](https://www.sciencedirect.com/science/article/pii/S0165783618300845))


The JARA framework provides a user options for either fitting multiple relative abundance indices, to estimate a mean trend, or multiple absolute abundance indices (e.g. census data from different breeding colonies) to produce a summed population trend for the total population.If the time series is longer than than three generation times (GT) the %decline is automatically estimated from the last assessment year minus 3 x GT. If the time series is shorter than 3 x GT, JARA projects forward.

Current examples include:  
- [Cape gannet](https://github.com/Henning-Winker/JARA/tree/master/Cape_gannet/output1) (census)  
- [African Penguin](https://github.com/Henning-Winker/JARA/tree/master/Afr_penguin/output1) (Census)  
- [Chinstrap penguin](https://github.com/Henning-Winker/JARA/tree/master/Chin_penguin/output1) (Census)  
- [Mountain Zebra](https://github.com/Henning-Winker/JARA/tree/master/Mountain_Zebra/output1) (Census, Carrying Capacity option) 
- [North Atlantic Shortfin mako shark](https://github.com/Henning-Winker/JARA/tree/master/SMA_NAtl/output1) (relative indices)  
- [South Atlantic blue shark](https://github.com/Henning-Winker/JARA/tree/master/BSH_SAtl/output1) (relative indices)  
- [Indian Ocean striped marlin](https://github.com/Henning-Winker/JARA/tree/master/StripedMarlin.IOTC/output1) (relative indices)  


All examples can be run using the development [`JARA_Prime.v1.3beta.R file`](https://github.com/Henning-Winker/JARA/blob/master/JARA_Prime.v1.3beta.R) from which the JARA model [`JARA.v1.3beta.R`](https://github.com/Henning-Winker/JARA/blob/master/JARA.v1.3beta.R) is executed.

*Still needs a lot of beautification!*..... 
