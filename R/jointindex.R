
#' jointindex()
#'
#' Generates a joint index from the population and joint residuals     
#' @param jara output list from fit_jara
#' @param credibility Confidence interval credibility interval with default of 0.95 
#' @return data.frame
#' @export
#' @author Henning Winker 

jointindex <- function(jara,credibility=0.95){
  
  yr = jara$yr
  trj = jara$pop.posterior[,1:length(yr)]
  mures = (apply(jara$residuals,2,mean,na.rm=T))
  index.posterior =exp(t(t(log(trj))+mures))
  
  
  
  Index = apply(index.posterior,2,median)
  trend = apply(trj,2,median)
  mu = mean(trend)
  CIs = HDInterval::hdi(index.posterior,credMass=credibility)
  log.se = apply(log(index.posterior),2,sd)
  out = data.frame(year=yr,trend=trend/mu,index=Index/mu,log.se,lci=CIs[1,]/mu,uci=CIs[2,]/mu)
  rownames(out) = 1:nrow(out)
  return(out)
}
