#' rqpois()
#'
#' Random generator for quasi-poisson variables    
#' @param n number of random observations
#' @param lambda vector of (non-negative) expected means
#' @param phi over-dispersion (variance) parameter 
#' @return random quasi-poisson observations
#' @export
rqpois = function(n, lambda, phi) {
  mu = lambda
  k = mu/phi/(1-1/phi)
  r = stats::rnbinom(n, mu = mu, size = k)
  return(r)
}

#' assign_iucn()
#'
#' Function to assign IUCN status based on decline and criteria     
#' @param change measured change over 3 generations 
#' @param criteria A1 or A2 for decline
#' @return IUCN status
#' @export
assign_iucn <- function(change,criteria=c("A2","A1")[1]){
  A = ifelse(criteria=="A1",1,2)
  status = ifelse(change< -80,"CR",ifelse(change>= c(-90,-80)[A] & change< c(-70,-50)[A],"EN",ifelse(change>= c(-70,-50)[A] & change< c(-50,-30)[A],"VU","LC"))) 
  return(status)
}  

#' iucn_estimators()
#'
#' Current estimators based on two points (p2) and log-linear regression (reg)      
#' @param I obserseved abundance index data.frame(yr,y)  
#' @param GL generation length
#' @param criteria A1 or A2 for decline
#' @return change over 3 x GL for methods p2 and reg
#' @export

iucn_estimators = function(I=NULL,GL=10,criteria=c("A2","A1")[1]){
  d. = data.frame(I)[,1:2]
  colnames(d.) = c("yr","y")
  p2obs = d.[nrow(d.),2]/d.[1,2]
  fit.reg = lm(log(y)~yr,d.)
  pr.reg = exp(predict(fit.reg))
  regobs = as.numeric(pr.reg[nrow(d.)]/pr.reg[1])
  p2change = (as.numeric(p2obs^(3*GL/length(d.$yr)))-1)*100
  regchange = (as.numeric(regobs^(3*GL/length(d.$yr)))-1)*100
  
  status = data.frame()
  out = list()
  out$reg.r = data.frame(est=coef(fit.reg)[2],se=summary(fit.reg)$coef[2,2])
  out$reg.pr = data.frame(yr=d.$yr,yobs=d.$y,yhat=pr.reg)
  out$change = data.frame(p2=p2change,reg=reg.change)
  out$status = data.frame(p2=assign_iucn(p2change,criteria),reg=assign_iucn(regchange,criteria))
  return(out)
}  
  




