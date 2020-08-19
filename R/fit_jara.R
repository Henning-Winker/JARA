#' fit_jara()
#'
#' Fits JARA model in JAGS and produce output object as list()
#' @param jarainput List of input variables as output by build_jara()
#' @param credibility Confidence Credibility Value with default of 0.95 
#' MCMC settings
#' @param ni number of iterations
#' @param nt thinning interval of saved iterations
#' @param nb burn-in
#' @param nc number of mcmc chains
#' @param save.jara saves jara list as .rdata to output.dir
#' @param save.all saves the all posteriors as .rdata to output.dir (big file)
#' @param saves the all posteriors as .rdata to output.dir (big file)
#' @param save.csvs writes results into csv to output.dir
#' @param output.dir path to save plot. default is getwd()
#' @param peels = NULL, # retrospective 'peel' option
#' @param save.jarafile saves jara model and convergence stats as .txt file (default = TRUE)
#' @param quickmcmc option "test run" jara with short (fast) mcmc chains 
#' @return A result list containing estimates of JARA model input, settings and results
#' @export
#' @author Henning Winker and Richard Sherley 
#' @examples
#' data(jaradata)
#' inp = jrdat$Yellowspot_skate
#' buildjr <- build_jara(I = inp$I,se=inp$SE,GL=12)
#' fitjr <- fit_jara(buildjr)
#' jrplot_fits(fitjr)
#' jrplot_trjfit(fitjr)
#' jrplot_poptrj(fitjr)
#' jrplot_changes(fitjr)
#' jrplot_state(fitjr)
#' jrplot_iucn(fitjr)
 

fit_jara = function(jarainput,credibility=0.95,
                     # MCMC settings
                     ni = 9000, # Number of iterations
                     nt = 2, # Steps saved
                     nb = 2000, # Burn-in
                     nc = 3, # number of chains
                     save.jara = FALSE,
                     save.all = FALSE,
                     save.csvs = FALSE,
                     save.jarafile = TRUE,
                     peels = NULL, 
                    output.dir = getwd(),
                    quickmcmc = FALSE
                    ){
  #write jara model
  jara2jags(jarainput)
  if(quickmcmc==TRUE){
  ni <- 7000 # Number of iterations
  nt <- 1 # Steps saved
  nb <- 1000 # Burn-in
  nc <- 2
  }
  
  
  # mcmc saved
  nsaved = (ni-nb)/nt*nc
  
  cat(paste0("\n","><> Running  JARA as ",jarainput$settings$model.type," model for ",jarainput$settings$assessment," ",jarainput$settings$scenario," <><","\n"))
  

  # jara model data
  jd = jarainput$jagsdata
  n.indices = ncol(jd$y)
  nT = nrow(jd$y)
  years = jarainput$data$yr
  year = jarainput$data$pyr
  n.years= length(years)
  GL = jarainput$settings$GL 
  GL1 = round(GL,0) # rounded for r.recent
  GL3 = round(3*GL,0) # 3 x GL rounded for year steps
  
  
  if(is.null(peels)) peels = 0
  peel.i = NULL
  if(peels > 0){
    peel.i=(length(years)-peels+1) : length(years)
    jd$y[peel.i,]  = NA
  }
  
  # Check if all indices have observations
  checks = exp(jd$y)
  checks[is.na(checks)] = 0
  jd$y[1,which(apply(checks,2,mean)==0)] = mean(jd$y,na.rm=T)
  
  # Initial starting values (new Eq)
  if(jarainput$settings$model.type=="census"){
    inits <- function(){list(mean.r = rnorm(n.indices,0,0.5),isigma2.est=runif(1,20,100), itau2=runif(jarainput$settings$nvar,80,200), logN.est =matrix(log(rbind(jarainput$settings$Ninit,matrix(rep(NA,(nT-1)*n.indices),(nT-1),n.indices))),nT,n.indices))  }
  } else {
    inits <- function(){list(isigma2.est=runif(1,20,100), itau2=runif(jarainput$settings$nvar,80,200), mean.r = rnorm(0,0.2),iq = 1/jarainput$settings$q.init)}
  } 
  
  out = output.dir
  if(file.exists(out)==FALSE) stop("\n","\n","><> output.dir does not exist <><","\n","\n")

  # jara model building conditions
  params = jarainput$settings$params


  ptm <- proc.time()

  jara.mod <- R2jags::jags(jd, inits,params,paste0(tempdir(),"/JARA.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)  # adapt is burn-in

  proc.time() - ptm
  save.time = proc.time() - ptm

  # unpack
  settings = c(jarainput$data,jarainput$jagsdata,jarainput$settings)
  abundance = settings$model.type
  parameters = settings$params
  qs=settings$qs
  indices = names(settings$I)[-1]
  
  
  if(n.years<GL1){
    WARN = paste0("\n","WARNING!!!","\n","-----------------------------------------------------\n",
                  "\n","The projections are based on less than one Generation Time!",
                  "\n","\n","The results are unreliable and MUST be interpreted with caution","\n","\n")
  } else {WARN=""}
  
  Chains1 = paste0("\n",paste0("><> Ran  ",settings$assessment,"  with ",nc," mcmc chains"),"\n")
  Chains2 =  paste0("\n",paste0("><> Each with ",ni, " iterations, a burn-in of ",nb,", a thinning rate of ",nt),"\n")
  Chains3 = paste0("\n",paste0("><> A total of ",nsaved," MCMC iterations were saved","\n"))
  
  
  RunTime =paste0("\n",paste0("><> Run  ",settings$assessment," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n"))
  cat(Chains1)
  cat(Chains2)
  cat(Chains3)
  cat(RunTime)
  
  
  #-----------------------------------------------------------
  # <><<><<><<><<><<><<>< Outputs ><>><>><>><>><>><>><>><>><>
  #-----------------------------------------------------------
  posteriors=jara.mod$BUGSoutput$sims.list
  if(abundance=="census"){sel.par=c(1:2)} else {sel.par=c(1:2,6)}
  par.dat= data.frame(posteriors[parameters[sel.par]])
  geweke = coda::geweke.diag(data.frame(par.dat))
  pvalues <- 2*pnorm(-abs(geweke$z))
  pvalues
  heidle = coda::heidel.diag(data.frame(par.dat))
  
  # Capture Results
  results = round(data.frame(cbind(apply(par.dat,2,median),t(HDInterval::hdi(par.dat,credMass=credibility)))),4)

  
  pars = data.frame(median = results[,1],lci=results[,2],uci=results[,3],Geweke.p=round(pvalues,3),Heidelberger.p=round(heidle[,3],3))
  if(abundance=="relative"){par.names = c("mu.r","sigma.proc",paste0("q.",indices[qs]))} else {
    par.names = c(paste0("r.",indices[qs]),"sigma.proc")  
  }
  row.names(pars) = par.names
  
  r.prj = data.frame(t(rbind(apply(posteriors$r.proj,2,median),HDInterval::hdi(posteriors$r.proj,credMass=credibility))))
  colnames(r.prj) = c("mu","lci","uci")
  if(abundance=="census"){
    rownames(r.prj) = c(paste0("r.",indices[qs]))
  }
  
  if(save.jarafile==TRUE){
  options(max.print=1000000)
  
    
  sink(paste0(output.dir,"/JARAjags_",settings$assessment,"_",settings$scenario,".txt"))
  cat(WARN)
  cat("-----------------------------------------------------\n")
  cat(RunTime)
  cat("-----------------------------------------------------\n")
  
  cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
  cat(paste0("\n","><> Input ",settings$assessment," ",settings$abundance," <><"))
  cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
  print(jarainput$settings)
  cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
  cat(paste0("\n","><> Parameter Estimates ",settings$assessment," ",settings$abundance," <><"))
  cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
  print(pars)
  cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
  cat(paste0("\n","><> JAGS output ",settings$assessment," ",settings$abundance," <><"))
  cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
  print(jara.mod)
  
  sink()
  options(max.print=1000)
  }
  
  
  #----------------------------
  # get individual trends
  #---------------------------
  end.yr = length(years)
  fitted <- lower <- upper <- lo.ppd <- hi.ppd <- matrix(NA,end.yr,(n.indices))
  
  if(abundance=="census"){
    for(i in 1:n.indices){
      for (t in 1:end.yr){
        fitted[t,i] <- median(posteriors$N.est[,t,i])
        lower[t,i] <- HDInterval::hdi(posteriors$N.est[,t,i],credMass=credibility)[1]
        upper[t,i] <- HDInterval::hdi(posteriors$N.est[,t,i],credMass=credibility)[2]
        lo.ppd[t,i] <-HDInterval::hdi(posteriors$ppd[,t,i],credMass=credibility)[1]
        hi.ppd[t,i] <-HDInterval::hdi(posteriors$ppd[,t,i],credMass=credibility)[2]}}
    
    # log-normal bias correction from Sherley et al. (in press)
    Nbias.correct <- array(0, dim=c((((ni-nb)*nc)/nt),nT,n.indices))
      for(i in 1:n.indices){
        for (t in 1:nT){
          Nbias.correct[,t,i] = exp(log(posteriors$N.est[,t,i])-0.5*var(log(posteriors$N.est[,t,i])))
        }}
    
    Nfit <- Nlow <- Nhigh <- as.numeric()
    
    pop.posterior = NULL
    # get total pop size
    for (t in 1:nT){
      pop.posterior = cbind(pop.posterior,rowSums(Nbias.correct[,t,]))
      Nfit[t] = mean(rowSums(Nbias.correct[,t,]))# o<
      Nlow[t] = HDInterval::hdi(rowSums(Nbias.correct[,t,]),credMass=credibility)[1]# o<
      Nhigh[t] = HDInterval::hdi(rowSums(Nbias.correct[,t,]),credMass=credibility)[2]# o<
    }  
    } else {
    q.adj = apply(posteriors$q,2,median)
    fitted <- lower <- upper <- as.null()
    for (t in 1:end.yr){
      fitted[t] <- median(posteriors$Y.est[,t])
      lower[t] <- HDInterval::hdi(posteriors$Y.est[,t],credMass=credibility)[1]
      upper[t] <- HDInterval::hdi(posteriors$Y.est[,t],credMass=credibility)[2]
    }
    
    lo.ppd <- hi.ppd <- matrix(NA,end.yr,(n.indices))
    for(i in 1:n.indices){
      for (t in 1:end.yr){
        lo.ppd[t,i] <- HDInterval::hdi(posteriors$ppd[,t,i],credMass=credibility)[1]
        hi.ppd[t,i] <- HDInterval::hdi(posteriors$ppd[,t,i],credMass=credibility)[2]}}
    
    
      # Total population
      Nfit <- Nlow <- Nhigh <- as.numeric()
      pop.posterior = NULL
      # get total pop size
      for (t in 1:nT){
      pop.posterior = cbind(pop.posterior,posteriors$Y.est[,t])  
      Nfit[t] = median(posteriors$Y.est[,t])
      Nlow[t] = HDInterval::hdi(posteriors$Y.est[,t],credMass=credibility)[1]
      Nhigh[t] = HDInterval::hdi(posteriors$Y.est[,t],credMass=credibility)[2]}
    }
      
  
  

  
  #------------------------------------------------------------------
  # Goodness of fit Statistics
  #------------------------------------------------------------------
  
  DIC =round(jara.mod$BUGSoutput$DIC,1)
  
  # get residuals
  Resids = NULL
  for(i in 1:n.indices){
    if(abundance=="census") Resids =rbind(Resids,log(settings$I[,i+1])-log(fitted[,i]))   
    if(abundance=="relative") Resids =rbind(Resids,log(settings$I[,qs[i]+1]/q.adj[i])-log(fitted))   
  }
  
  Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
  RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/(Nobs-1)),1)
  # Produce statistice describing the Goodness of the Fit
  GOF = data.frame(Stastistic = c("Nobs","RMSE","DIC"),Value = c(Nobs,RMSE,DIC))
  
  # Save Obs,Fit,Residuals
  diags = NULL
    for(i in 1:n.indices){
      
      Yr = years
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      
      if(abundance=="relative"){
      exp.i = fitted[is.na(settings$I[,qs[i]+1])==F]*q.adj[i]
      expLCI.i = lower[is.na(settings$I[,qs[i]+1])==F]*q.adj[i]
      expUCI.i = upper[is.na(settings$I[,qs[i]+1])==F]*q.adj[i]
      ppdLCI.i = lo.ppd[is.na(settings$I[,qs[i]+1])==F,i]
      ppdUCI.i = hi.ppd[is.na(settings$I[,qs[i]+1])==F,i]
      } else {
      exp.i = fitted[is.na(settings$I[,qs[i]+1])==F,i]
      expLCI.i = lower[is.na(settings$I[,qs[i]+1])==F,i]
      expUCI.i = upper[is.na(settings$I[,qs[i]+1])==F,i]
      ppdLCI.i = lo.ppd[is.na(settings$I[,qs[i]+1])==F,i]
      ppdUCI.i = hi.ppd[is.na(settings$I[,qs[i]+1])==F,i]
      } 
      obs.i = settings$I[,qs[i]+1][is.na(settings$I[,qs[i]+1])==F]
      sigma.obs.i = (apply(posteriors$TOE[,,i],2,quantile,c(0.5)))[is.na(settings$y[,i])==F]
      
      yr.i = Yr[is.na(settings$I[,qs[i]+1])==F]
      diags = rbind(diags,data.frame(assessment=settings$assessment,scenario=settings$scenario,name=names(settings$I)[qs[i]+1],year=yr.i,obs=obs.i,obs.err=sigma.obs.i,hat=exp.i,lci=expLCI.i,uci=expUCI.i,lpp=ppdLCI.i,upp=ppdUCI.i,residual=log(obs.i)-log(exp.i),retro.peels=peels))
    } 
    
  # predicted abundance trajectories
  
  trj=data.frame(assessment=settings$assessment,scenario=settings$scenario,name="global", yr=year,Type=abundance,estimation=ifelse(year %in% years,rep("fit",length(year)) ,rep("prj",length(year))), mu=Nfit,lci=Nlow,uci=Nhigh,lpp=NA,upp=NA)
  
  for(i in 1:n.indices){
    Yr = year
    Yr = min(Yr):max(Yr)
    yr = Yr-min(years)+1
    if(abundance=="relative"){
      exp.i = apply(posteriors$Y.est,2,median)*q.adj[i]
      expLCI.i =HDInterval::hdi(posteriors$Y.est,credMass=credibility)[1,]*q.adj[i]
      expUCI.i = HDInterval::hdi(posteriors$Y.est,credMass=credibility)[2,]*q.adj[i]
      ppdLCI.i = HDInterval::hdi(posteriors$ppd[,,i],credMass=credibility)[1,]
      ppdUCI.i = HDInterval::hdi(posteriors$ppd[,,i],credMass=credibility)[2,]
      
    } else {
      exp.i = apply(posteriors$N.est[,,i],2,median)
      expLCI.i = HDInterval::hdi(posteriors$N.est[,,i],credMass=credibility)[1,]
      expUCI.i = HDInterval::hdi(posteriors$N.est[,,i],credMass=credibility)[2,]
      ppdLCI.i = HDInterval::hdi(posteriors$ppd[,,i],credMass=credibility)[1,]
      ppdUCI.i = HDInterval::hdi(posteriors$ppd[,,i],credMass=credibility)[2,]
      
      } 
      
    trj=rbind(trj,data.frame(assessment=settings$assessment,scenario=settings$scenario,name=names(settings$I)[qs[i]+1], yr=year,Type=abundance,estimation=ifelse(year %in% years,rep("fit",length(year)) ,rep("prj",length(year))), mu=exp.i,lci=expLCI.i,uci=expUCI.i,lpp=ppdLCI.i,upp=ppdUCI.i))
  } 
  
  if(abundance=="census"){
    rs =  as.matrix(apply(posteriors$r.tot[,1:(n.years-1)],1,median))
    if(n.years>GL1) rs = cbind(rs,apply(posteriors$r.tot[,(n.years-round(GL,0)+1):(n.years-1)],1,median))
    if(n.years>2*GL1) rs = cbind(rs,apply(posteriors$r.tot[,(n.years-round(GL*2,0)+1):(n.years-1)],1,median))
    if(n.years>3*GL1) rs = cbind(rs,apply(posteriors$r.tot[,(n.years-round(GL*3,0)+1):(n.years-1)],1,median))
    
    
  } else {
    rs =  as.matrix(apply(posteriors$r[,1:(n.years-1)],1,median)) 
    if(n.years>GL) rs = cbind(rs,apply(posteriors$r[,(n.years-round(GL,0)+1):(n.years-1)],1,median))
    if(n.years>2*GL) rs = cbind(rs,apply(posteriors$r[,(n.years-round(GL*2,0)+1):(n.years-1)],1,median))
    if(n.years>3*GL) rs = cbind(rs,apply(posteriors$r[,(n.years-round(GL*3,0)+1):(n.years-1)],1,median))
  }
  
  change = (apply(posteriors$Ntot[,(settings$mp.assess[2]-1):(settings$mp.assess[2]+1)],1,median)/apply(posteriors$Ntot[,(settings$mp.assess[1]-1):(settings$mp.assess[1]+1)],1,median)-1)*100
  cnam = c("All.yrs","1GL","2GL","3GL")  
  trends =  data.frame(change,rs)
  colnames(trends ) = c("pop.change",paste0("r.",cnam[1:ncol(rs)]))
  
  pop.posterior = data.frame(  pop.posterior)
  colnames(pop.posterior) = Yr
  #jara = list(settings=Settings,data=jags.data,parameters=pars,trends=trends,abundance=abundance.est,perc.risk=data.frame(CR=CR,EN=EN,VU=VU,LC=LC) ,status=status)
  

  #-------------------------------
  # summarize results in jara list
  #-------------------------------
  jara = list()
  jara$assessment  = settings$assessment
  jara$scenario = settings$scenario
  jara$yr = years
  jara$pyr = year
  jara$indices = indices
  jarainput$settings$peel.i = peel.i
  jara$settings = c(jarainput$jagsdata,jarainput$settings)
  jara$inputseries = list(I=jarainput$data$I,se=jarainput$data$se)
  jara$pars=pars
  jara$stats = GOF
  jara$fits = diags  
  jara$residuals = Resids
  jara$trj = trj
  jara$posteriors = trends
  jara$r.prj = r.prj
  jara$pop.posterior = pop.posterior
  
  if(save.jara==TRUE){
    save(jara,file=paste0(output.dir,"/",settings$assessment,"_",settings$scenario,"_jara.rdata"))
  }
  
  # Safe posteriors (Produces large object!)
  if(save.all==TRUE) save(posteriors,file=paste0(output.dir,"/",settings$assessment,"_",settings$scenario,"_posteriors.rdata"))
  
  if(save.csvs==TRUE){
    # Save results
    write.csv(diags,paste0(output.dir,"/Fits_",settings$assessment,"_",settings$scenario,".csv"),row.names = FALSE)
    write.csv(data.frame(pars),paste0(output.dir,"/Estimates_",settings$assessment,"_",settings$scenario,".csv"))
    write.csv(trj,paste0(output.dir,"/Trajectories_",settings$assessment,"_",settings$scenario,".csv"),row.names = FALSE)
    write.csv(jara$stats,paste0(output.dir,"/GoodnessFit_",settings$assessment,"_",settings$scenario,".csv"),row.names = FALSE)
    write.csv(trends,paste0(output.dir,"/Posteriors_",settings$assessment,"_",settings$scenario,".csv"),row.names = FALSE)
    
  }
  
  
  
  return(jara)

  } # END of fit_jara()
  
  