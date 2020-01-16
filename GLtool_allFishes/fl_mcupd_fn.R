#----------------------------------------------------------------------------------------
# Proto-type MVN Monte-Carlo tool to condition stock specific life-history traits,
# generation length (GL) and  intrinsic rate of population increase (r) for all fishes 
# by sampling of predictive distributions from FishLife2.0 (Thorson 2020).
# https://github.com/James-Thorson-NOAA/FishLife
# 
# Developed by Henning Winker, Cape Town 2020
# henning.winker@gmail.com
#  
# Note: 
# The initial motivation was to generate objective: 
# (1) Generation Length for JARA
# https://github.com/Henning-Winker/JARA
# (2) r priors for Schaefer model runs in JABBA, AMSY and CMSY+
# https://github.com/jabbamodel/JABBA;https://github.com/SISTA16/amsy; https://github.com/SISTA16/cmsy 
# (3) life history input parameters for JABBA-SELECT
# https://github.com/jabbamodel/JABBA-Select
#
#-------------------------------------------------------------------------

#installation
#install.packages("devtools")
#install.packages("mvtnorm")
#devtools::install_github("james-thorson/FishLife")

library(FishLife)
library(mvtnorm)

setwd("C:/Work/Research/Github/JARA/GLtool_allFishes")
#--------------------------------------------------------
# Monte-Carlo MVN FishLife2.0 update function for tuning 
#--------------------------------------------------------
fl_mcupd <- function(Genus="Rhabdosargus",Species="globiceps",parms=c("Loo","K","Lm","tm","tmax","M","logitbound_h","ln_r","ln_G"),Loo = NULL,K=NULL,tmax=NULL,tm=NULL,M=NULL,Lm=NULL,h=c(0.3,0.9),nrep = 5000,nmc = 2*10^5,Plot=TRUE,PlotPath=getwd()){
  from_logith <- function(logit_h){
    out=  0.2001 + 0.7998*1/(1+exp(-logit_h))
    out}
  # Search Taxa
  taxa = Search_species(Genus=Genus,Species = Species,add_ancestors=TRUE)$match_taxonomy
  
  # Predict LifeLife Traits  
  png(file = paste0(PlotPath,"/",Genus,".",Species,".traits.png"), width = 6, height = 7, 
      res = 200, units = "in")
  predfl =Plot_taxa(taxa,mfrow=c(3,2))
  dev.off()
  mu= predfl[[1]]$Mean_pred
  covs = predfl[[1]]$Cov_pred
  
  priors = parms
  priors[parms=="logitbound_h"] = "h"
  priors[parms=="ln_r"] = "r"
  priors[parms=="ln_G"] = "G"
  pars= priors
  priors = priors[1:6]
  
  mcparms = data.frame(rmvnorm(nmc,mean =  mu,sigma = covs,  method=c("eigen")))
  
  nprior = 1
  prd = NULL
  for(i in 1:6){
    getpr = get(priors[i])  
    
    # pdf of prior
    if(is.null(getpr)==FALSE){
      nprior = nprior+1 
      prd = cbind(prd,dnorm(mcparms[,parms[i]],log(getpr[1]),getpr[2]))
    }
  } # end of loop
  
  if(is.null(prd)==FALSE){
    prp = apply(prd,1,sum)/max(apply(prd,1,sum)) # Convert to probability  
    # subsample conditioned given prp
    rand = runif(nrow(mcparms))
    sel = which(rand < prp)
    subparms = mcparms[sel,]
  } else {
    subparms = mcparms
  }
  # subset h range
  subparms = subparms[subparms$h>h[1] & subparms$h<h[2],]
  # Get means and vars from subset
  mu.upd = apply(subparms,2,mean)
  covs.upd = covs 
  diag(covs.upd) = apply(subparms,2,var)
  updfl = predfl
  updfl[[1]]$Mean_pred = mu.upd # Here only updating on species (stock level) 
  updfl[[1]]$Cov_pred = covs.upd # Here only updating on species (stock level)
  
  # Plot
  png(file = paste0(PlotPath,"/",Genus,".",Species,".flupd.png"), width = 7, height = 6, 
      res = 200, units = "in")
  Par = list(mfrow=c(3,3),mai=c(0.5,0.5,0,.1),omi = c(0.2,0.2,0.1,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.7)
  par(Par)
  for(i in 1:9){
    if(pars[i]!="h"){
      xi = sort(rlnorm(1000,mu[parms[i]],sqrt(diag(covs)[parms[i]])))  
      pdf.upd =dlnorm(xi,mu.upd[parms[i]],sqrt(diag(covs.upd)[parms[i]]))
      pdf.fl = dlnorm(xi,mu[parms[i]],sqrt(diag(covs)[parms[i]])) 
    } else {
      xi = sort((rnorm(1000,mu[parms[i]],sqrt(diag(covs)[parms[i]]))))  
      pdf.upd =dnorm(xi,mu.upd[parms[i]],sqrt(diag(covs.upd)[parms[i]]))
      pdf.fl = dnorm(xi,mu[parms[i]],sqrt(diag(covs)[parms[i]]))  
      xi = from_logith(xi)
    }  
    plot(0,0,type="h",xlab=paste(pars[i]),ylab="Density",xlim=quantile(xi,c(0.01,0.99)),ylim=range(c(pdf.fl,pdf.upd),c(0.01)))
    polygon(c(xi,rev(xi)),c(pdf.fl,rep(0,1000)),col="grey")  
    polygon(c(xi,rev(xi)),c(pdf.upd,rep(0,1000)),col=rgb(1,0,0,0.5),border=rgb(1,0,0,0.5))  
    if(i==1) legend("topright",c("FishLife","Update"),pch=22,pt.cex = 2,pt.bg=c("grey",rgb(1,0,0,0.5)),bty="n")
  }
  dev.off()
  
  # Summarize results
  sd.fl = sqrt(diag(covs[parms,parms]))
  mean.fl= mu[parms]
  mu.fl = ifelse(parms!="logitboundh",exp(mean.fl),from_logith(mean.fl))
  lcl.fl = ifelse(parms!="logitboundh",exp(mean.fl-1.96*sd.fl),from_logith(mean.fl-1.96*sd.fl))
  ucl.fl = ifelse(parms!="logitboundh",exp(mean.fl+1.96*sd.fl),from_logith(mean.fl+1.96*sd.fl))
  sd.fl = sqrt(diag(covs.upd[parms,parms]))
  mean.upd= mu.upd[parms]
  mu.upd = ifelse(parms!="logitboundh",exp(mean.upd),from_logith(mean.upd))
  lcl.upd = ifelse(parms!="logitboundh",exp(mean.upd-1.96*sd.fl),from_logith(mean.upd-1.96*sd.upd))
  ucl.upd = ifelse(parms!="logitboundh",exp(mean.upd+1.96*sd.fl),from_logith(mean.upd+1.96*sd.upd))
  out = data.frame(trait=pars,mu.fl,lcl.fl,ucl.fl,mu.upd,lcl.upd,ucl.upd)
  out[,-1] = round(out[,-1],4)
  
  return(list(traits=out,vcm=updfl))
  
} 
#-------------------------------------------------------------------------------------------
# End of function 
#-------------------------------------------------------------------------------------------

#><>><>><>><>
# EXAMPLES
#><>><>><>><>
# White stumpnose (Endemic Sparid, South Africa)
Genus = "Rhabdosargus"
Species = "globiceps"

fl_mcupd(Genus=Genus,Species = Species,Loo=c(50,0.2),tm=c(3,0.2),h=c(0.4,0.9))

# School (Soupfin) shark
Genus = "Galeorhinus"
Species = "galeus"
h=c(0.21,0.6)
fl_mcupd(Genus=Genus,Species = Species,Loo=c(150,0.2),h=h)$traits
