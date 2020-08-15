#' jrdata_simulator()
#'
#' simulates observed and "true" abundance time series    
#' @param n0 start population size
#' @param yrs number of simulated years with observations  
#' @param r.mu mean rate of change (default random uniform number)
#' @param proc.pop lognormal process error for natural population variation  
#' @param AR1.pop 1st order autoregressive error on poplation size 
#' @param proc.imp lognormal process error on (anthropogenic) impact
#' @param AR1.imp 1st order autoregressive error on impact
#' @param CV.obs observation error (~ log.sd)
#' @param GL Generation length to measure decline over 3 x GL 
#' @param criteria A1 or A2 for decline
#' @param Plot Shows plot with with simulated data and processes  
#' @param add if TRUE par is not called to enable manual multiplots
#' @return simulated data and 'true' decline metrics
#' @export
#' @author Henning Winker (JRC-EC, Ispra)
#' @examples
#' simdat <- jrdat_simulator(yrs=40,GL=15)
#' inpsim <- build_jara(I = simdat$I,GL=simdat$GL)
#' fitsim <- fit_jara(inpsim)
#' jrpar(mfrow=c(1,2))
#' jrplot_poptrj(fitsim,add=TRUE)
#' lines(simdat$N,col="cyan3",lwd=2)
#' jrplot_iucn(fitsim,add=TRUE)
#' abline(v=simdat$perc.change,col="cyan3",lwd=2)

jrdat_simulator = function(n0 = 1000 ,yrs = 40,r.mu = runif(1,-0.05,0),proc.pop = 0.1,AR1.pop=0.3,proc.imp=0.05,AR1.imp=0.8,CV.obs=0.15,GL=10,criteria=c("A2","A1")[1],Plot=TRUE,add=FALSE){ 
  if(max(yrs)>=round(3*GL)+1){
    simyrs = yrs
  } else {
    simyrs = round(3*GL)+1
  }
  # Process error
  # (1) random variation in popdyn 
  pop_var = stats::rnorm(simyrs,0,proc.pop) 
  # Add AR1 (serial auto-correlation)
  pop_dev = pop_var[1]
  for(y in 2:simyrs) pop_dev[y] =  AR1.pop * pop_dev[y-1]+pop_var[y]*sqrt(1-AR1.pop^2)
  # (2) random variation in impact (e.g. offtake) 
  imp_var = stats::rnorm(simyrs,0,proc.imp) 
  # Add AR1 (serial auto-correlation)
  imp_dev = imp_var[1]
  for(y in 2:simyrs) imp_dev[y] =  AR1.imp * imp_dev[y-1]+imp_var[y]*sqrt(1-AR1.imp^2)
  
  # Combine all as latent effect in r vector given r.mu
  r = r.mu+pop_dev-imp_dev
  
  # Add random observation error
  obs_dev = stats::rnorm(simyrs,0,CV.obs)
  
  # True Population
  n_t = n0*exp(r[1]) 
  for(t in 2:simyrs){
    n_t[t] =  n_t[t-1]*exp(r[t-1]) 
    if(n_t[t]>n0*1.5) n_t[t] = n0*1.5}
  
  # Deterministic trend
  nbar = n0 
  for(t in 2:simyrs){
    nbar[t] =  nbar[t-1]*exp(r.mu)}
  
  # Imperfect observation
  y_t =exp(log(n_t)+obs_dev)
  
  prjyrs = simyrs-yrs
  obsyrs = yrs
  change=(n_t[simyrs]/n_t[simyrs-round(3*GL)+1]-1)*100
  y_t = y_t[1:yrs]
  if(Plot==TRUE){
    jrpar(mfrow=c(1,2))
    # Error 
    plot(n_t,ylim=c(0,max(n_t,y_t,n0*1.05)),type="l",col=2,lwd=2,ylab="Abundance",xlab="Year")
    points(y_t,pch=21,bg="white",type="p")
    lines(nbar,lty=2)
    legend("topright",c("True","Observed","Determistic"),pch=c(-1,1,-1),lwd=c(2,-1,1),col=c(2,1,1),lty=c(1,-1,2),bty="n",cex=0.8)
    legend("topleft","a)",bty="n",cex=1.1,x.intersp = -0.5,y.intersp = -0.2)
    plot(pop_dev,ylim=range(c(pop_dev,imp_dev,obs_dev)),type="l",col=3,lwd=2,ylab="Deviations",xlab="Year")
    lines(imp_dev,col=4,lwd=2)
    legend("topright",c("Stochasticity","Impact","Observation"),pch=c(-1,-1,1),lwd=c(2,2,1),col=c(3,4,1),bty="n",cex=0.8)
    points(obs_dev[1:obsyrs],type="b")
    abline(h=0,lty=2)
    legend("topleft","b)",bty="n",cex=1.1,x.intersp = -0.5,y.intersp = -0.2)
    
  }
  
  
  return(list(I=data.frame(yr=1:obsyrs,y=y_t),
              N=data.frame(yr=1:simyrs,N=n_t),r.mu=r.mu,
              change=change,status=assign_iucn(change,criteria),
              simyrs=simyrs,obsyrs=obsyrs,prjyrs=prjyrs,GL=GL))

  
}
