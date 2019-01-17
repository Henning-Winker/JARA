# JARA: Just Another Redlist Assessment

JARA (*Just Another Redlist Assessment*) builds on a Bayesian state-space framework with the aim to objectively incorporate uncertainty into the Red Listing evaluation process . In general, state-space models provide a powerful tool for time series analysis, as they allow accounting for both process error (environmental year to year variation) and observation (or reporting) error simultaneously. The posterior for the estimated population reduction provides a natural way to assign probabilities of the size reduction falling within each of the Red Listing categories. For this purpose, we developed an easy to interpret graph, in which the posterior of the reduction estimates is plotted against the IUCN Red List criteria. One of my main motivations was that there are currently no specific guidelines for dealing with uncertainty in the Red Listing process, which can result in inconsistent, subjective assessments, especially when dealing with multiple abundance indices. 


The  current framework can deal with both census and relative abundance indices. JARA fits all available time series data simultaneously. If the time series is longer than than three generation times (GT) the %decline is automatically estimated from the last assessment year minus 3 x GT. If the time series is shorter than 3 x GT, JARA projects forward.

Current examples include:  
- [Cape gannet](https://github.com/Henning-Winker/JARA/tree/master/Cape_gannet/output1) (census)  
- [African Penguin](https://github.com/Henning-Winker/JARA/tree/master/Afr_penguin/output1) (Census)  
- [Chinstrap penguin](https://github.com/Henning-Winker/JARA/tree/master/Chin_penguin/output1) (Census)  
- [Mountain_Zebra](https://github.com/Henning-Winker/JARA/tree/master/Mountain_zebra/output1) (Census, Carrying Capacity option) 
- [North Atlantic Shortfin mako shark](https://github.com/Henning-Winker/JARA/tree/master/SMA_NAtl/output1) (relative indices)  
- [South Atlantic blue shark](https://github.com/Henning-Winker/JARA/tree/master/BSH_SAtl/output1) (relative indices)  
- [Indian Ocean striped marlin](https://github.com/Henning-Winker/JARA/tree/master/StripedMarlin.IOTC/output1) (relative indices)  


All examples can be run using the development [`JARA_Prime.v1.3beta.R file`](https://github.com/Henning-Winker/JARA/blob/master/JARA_Prime.v1.3beta.R) from which the JARA model [`JARA.v1.3beta.R`](https://github.com/Henning-Winker/JARA/blob/master/JARA.v1.3beta.R) is executed.

*Still needs a lot of beautification!*..... 
