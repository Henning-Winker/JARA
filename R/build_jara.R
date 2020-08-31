#' build_jara function
#'
#' Creates a data list with JABBA input and settings to be passed to fit_jara()
#' @param I abundance indices/counts,  require data.frame(year, I.1,I.2,...,I.N)
#' @param se optional log standard error (CV) time series,requires data.frame(year, se.1,se.2,...,se.N)
#' @param assessment = "species.X",
#' @param scenario = "s1",
#' @param model.type abundance data type c("relative","census")
#' @param GL Generation length (default n.years/3)
#' @param start.year subsetting option for start year
#' @param end.year subsetting option for end year
# Variance settings
#' @param fixed.obsE minimum plausible observation error (fixed)
#' @param variance.weighting option "equal", "model", "fixed" or user assinged vector e.g. c(1,1,2) 
#' @param sigma.proc.fixed option to fix the process error e.g = 0.1 (default FALSE)
#' @param proc.pen advanced user setting to penalize extreme process error deviations
# Porjection settings
#' @param proj.mod specified by theta c("theta10","logistic")
#' @param pk.prior lognormal prior of population to K c(mu,lod.sd,yr) envoked during projections  
#' @param pk.yr reference year for the depletion prior. NULL ratio to max pred. abundance across all years 
#' @param pk.i option to specify specify different depletion levels for colonies in the census model 
#' @param proj.r  # rate of change for c("all","GL1", years), "all" is default
#' @param proj.stoch allows for projections with process error c(TRUE, FALSE), FALSE is default
#' @param proj.yrs.user option to overwrite GL and costomize projections for forecasting  
#' @return List to be used as data input to JARA JAGS model
#' @export
#' @author Henning Winker
#' @examples 
#' data(jaradata)
#' inp = jrdat$Afr_penguin
#' buildjr <- build_jara(I = inp$I,GL=9,model.type="census")
#' jrplot_indices(buildjr)

build_jara <- function(I = NULL, se = NULL,assessment = "Unnamed",
              scenario = "s1",model.type = c("relative","census")[1],
              GL=NULL, 
              start.year = NA,
              end.year = NA,
              fixed.obsE = NULL,
              variance.weighting = "equal",
              sigma.proc.fixed = FALSE,
              proc.pen = NULL,  
              proj.mod = c("theta10","logistic","theta.value")[1],
              pk.prior = c(0.25,0.1),
              pk.yr = NULL,
              pk.i = NULL,
              proj.r = c("all","GL1","year")[1],
              proj.yrs.user = NULL,
              proj.stoch = FALSE){
  
  
  #-------------------------
  # Prepare input data
  #-------------------------
                  
                if(is.null(GL)&is.null(proj.yrs.user)){GL = floor((nrow(I)-3)/3)} else {GL=GL[1]}
                # If proj.yrs.user is provided it overwrites GL 
                if(is.null(proj.yrs.user)==FALSE){GL = floor((nrow(I)+proj.yrs.user-1)/3)} else {GL=GL[1]}
  
                GL1 = round(GL,0) # rounded for r.recent
                GL3 = round(3*GL,0) # 3 x GL rounded for year steps
                
                
                if(is.null(fixed.obsE)){
                fixed.obsE = ifelse(is.null(se),0.15,0.01)
                }
                if(is.null(proc.pen)){
                  proc.pen = c(ifelse(model.type=="census",1,0.5),log(ifelse(model.type=="census",0.1,0.5)),ifelse(model.type=="census",log(5),log(2)))
                }
                
                cat(paste0("\n","><> Prepare input data <><","\n","\n"))
                # Do checks 
                # check for zeros 
                find0 = which(as.matrix(I[,-1])==0,arr.ind = T)
                
                if(nrow(find0)>0){
                  warning("\n","There are zeros in the times series!","\n", 
                          "These will replaced by constants C = 0.1*mean(Index[i]) to prevent log(0) errors","\n",
                          "Please carefully check residual patterns and the origin of these zeros (true/false)","\n"
                          )
                  find0[,2] = find0[,2]+1 
                  muI = apply(as.matrix(I[,unique(find0[,2])]),2,mean,na.rm=T)
                  I[find0] = muI*0.1  
                }
                dat = data.frame(I)
                
                
                indices = names(dat)[2:ncol(dat)]
                n.indices = max(length(indices),1)
                if(is.na(start.year)) start.year = dat[1,1]   
                if(is.na(end.year)) end.year = dat[nrow(dat),1]   
                dat = dat[dat[,1]>=as.numeric(start.year),]  
                dat = dat[dat[,1]<=as.numeric(end.year),]
                
                # check for NAs in first years 
                Itest = dat
                Itest[is.na(Itest)] = 0
                na.rows = as.numeric(apply(as.matrix(Itest[,-1]),1,sum))
                na.st = which(na.rows>0)[1]
                if(na.st>1 ){
                  start.year = Itest[na.st,1]
                  
                  warning("\n","Only NAs in first year!","\n", 
                          "The first ",na.st-1," years are excluded from the dataset to prevent errors","\n",
                          "New start.year = ",start.year,"\n")
                  dat = dat[dat[,1]>=as.numeric(start.year),]  
                  
                  
                }
                
                
                if(is.null(se)==TRUE){ SE.I=FALSE} else {SE.I=TRUE}
                
                if(SE.I==TRUE){
                  se = se[se[,1]>=as.numeric(start.year),]  
                  se = se[se[,1]<=as.numeric(end.year),]  
                }
                
                sigma.obs.est =TRUE
                if(variance.weighting[1]=="equal"){
                  sets.var = rep(1,ncol(dat)-1)
                } else if(variance.weighting[1]=="model"){
                  sets.var = 1:(ncol(dat)-1)  
                } else if(variance.weighting[1]=="fixed"){ 
                  sets.var = rep(1,ncol(dat)-1)
                  sigma.obs.est =FALSE} else {
                  if(length(variance.weighting)!=(ncol(dat)-1)){
                    stop("variance assingment mismatch with availble indices")
                  }
                  sets.var = variance.weighting[1:(ncol(dat)-1)]}
                nvar = length(unique(sets.var))
                
                years=dat[,1]
                styr = min(years)
                endyr = max(years)
                n.years = length(years)
                # process error settings
                if(sigma.proc.fixed==FALSE){ 
                  #------------------------------------------
                  sigma.proc = TRUE
                  
                  igamma = c(0.001,0.001) #specify inv-gamma parameters
                  
                  # Process error check
                  gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
                  # check mean process error + CV
                  mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
                  
                  # check CV
                  round(c(mu.proc,CV.proc),3)
                  quantile(sqrt(gamma.check),c(0.1,0.9))
                }else{
                  sigma.proc = as.numeric(sigma.proc.fixed) #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
                }
                
                
                
                #------------------------------------
                # Index count, biomass or relative
                #------------------------------------
                conv.I = as.numeric(as.matrix(dat[,-1]))
                I_y=matrix(conv.I,nrow=n.years,ncol=n.indices)
                
                
                if(SE.I==FALSE){
                  se = dat  
                  conv.se = as.numeric(as.matrix(dat[,-1]))
                  se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
                } else{
                  conv.se = as.numeric(as.matrix(se[,-1]))
                  se2 = matrix(ifelse(is.na(conv.se),rep(0.01,length(conv.se)),conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
                }
                
                #------------------------------------
                # Projection Settings
                #------------------------------------
                
                if(proj.mod %in% c("theta10","logistic")){
                  theta = ifelse(proj.mod=="theta10",10,1)
                } else {
                  theta = as.numeric(proj.mod)
                }
                
                # Carrying capacity settings
                if(is.null(pk.yr)) pk.yr= years
                pk.y = which(years%in%pk.yr) 
                
                if(model.type=="census"){
                  pk.cv = rep(as.numeric(pk.prior[2]),n.indices)
                  
                  if(is.null(pk.i)==TRUE | length(pk.i)==1){ 
                  pk.mu = rep(pk.prior[1],n.indices)
                  } else {
                    pk.mu = pk.i
                  }
                  
                  if(n.indices!=length(pk.i))  cat("ERROR: pk.i value(s) must match the number of subpopulations") 
                } else {
                  pk.mu = pk.prior[1]
                  pk.cv = pk.prior[2]
                }
                
                # Adding projections if needed
                
                # prediction years 
                pyears = max(GL3+2-n.years,1) # Fixed
                
                # set up year vectors
                year <- years[1]:(years[length(years)] + pyears)
                yr = years[1]:(years[length(years)])
                
                nT = length(year)
                mp.assess = c(nT-GL3,nT-1)     #start and end midpoints for population trend assessment
                
                #define the r type used for projections
                prjr.type = "robs"  # DEFAULT 
                
                if(proj.r =="all"){
                  prjr = 1:n.years
                } else if(proj.r=="GL1" & GL1 < n.years) {
                  prjr = (n.years-GL1+1):n.years  
                } else if(proj.r=="GL1" & GL1 >= n.years){
                  proj.r = 1:n.years  
                } else {
                  prjr.st = round(as.numeric(proj.r),0)  
                  prjr = (years[n.years]-prjr.st+1):n.years
                } 
                
                # creat prediction matrix for observations
                I_y= rbind(I_y,matrix(NA,pyears,(n.indices))) # now always
                # No Zero allowed in log modell
                I_y[I_y==0]=1
                se2= rbind(se2,matrix(0.001,pyears,(n.indices))) # now always
                
                #---------------------
                # Index color palette
                #---------------------
                cols = as.character(c('#e6194b', "#3cb44b", "#ffe119",
                                              "#0082c8","#f58231", "#911eb4",
                                              "#46f0f0", "#f032e6", "#d2f53c",
                                              "#fabebe", "#008080","#e6beff", "#aa6e28",rainbow(10)))
                
                cat("\n","><> Assuming a Generation length of GL =",round(GL,2),"years <><","\n","\n")
                
                Ninit = NULL
                # get starting values Ninit
                for(i in 1:n.indices){
                  Ninit[i] = na.omit(I_y[,i])[1]  
                }
                
                if(model.type=="census"){
                
                  cat("\n","><> Setting up JARA for census data <><","\n","\n")  
                # Bundle data for jags
                  jags.data <- list(y = log(I_y+10^-20),SE2=se2, T = length(year),nI=n.indices,Ninit=Ninit,sets.var=sets.var,nvar=nvar,pk.mu=pk.mu,pk.cv=pk.cv,pk.y=pk.y,EY = n.years,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),igamma=igamma,penSig=0,prjr=prjr,proc.pen=proc.pen,theta=theta)
                  # order of indices
                  qs = 1:n.indices
                  # Parameters monitored
                  parameters <- c("mean.r","sigma","logNtot", "N.est","Ntot","r.tot","r.proj","TOE","ppd","K")
                }  
                                
                if(model.type=="relative"){
                  cat("\n","><> Setting up JARA for relative abundance indices <><","\n","\n")  
                  
                  #find first time-series with first CPUE
                  q1.y = c(1:length(year))[is.na(apply(I_y,1,mean,na.rm=TRUE))==FALSE][1] #first year with CPUE
                  q1.I =  which.max(I_y[q1.y,])
                  # Order of indices
                  qs = c(q1.I,c(1:(ncol(I_y)))[-q1.I])
                  
                  qI_y = as.matrix(I_y[,qs])
                  qse2 = as.matrix(se2[,qs])
                  
                  q.init = 1
                  if(n.indices>1) for(i in 2:n.indices){q.init[i] = mean(qI_y[,i],na.rm=TRUE)/mean(qI_y[,1],na.rm=TRUE)}
                  
                  # Bundle data
                  jags.data <- list(y = log(qI_y),SE2=qse2, logY1 = log(qI_y[1,1]), T = length(year),EY = n.years,nI=n.indices,sets.var=sets.var,nvar=nvar,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),igamma=igamma,penSig=0,pk.mu=pk.mu,pk.cv=pk.cv,pk.y=pk.y,prjr=prjr,proc.pen=proc.pen,theta=theta)
                  
                  # Parameters monitored
                  parameters <- c("mean.r", "sigma","r", "Y.est","Ntot","q","r.proj","TOE","ppd","K")
                  
                }
                
                cat("\n","><> Compile JARA input <><","\n","\n")  
                #--------------------------
                # Capture Settings
                #--------------------------
                jarainput = list()
                jarainput$data = list()
                jarainput$jagsdata = list()
                jarainput$settings = list()
                jarainput$data$yr = years
                jarainput$data$pyr = year
                jarainput$data$I = dat
                jarainput$data$se = se
                jarainput$jagsdata = jags.data
                jarainput$settings$model.type = model.type 
                jarainput$settings$assessment = assessment
                jarainput$settings$scenario = scenario
                jarainput$settings$params = parameters
                jarainput$settings$GL = GL
                jarainput$settings$SE.I = SE.I
                jarainput$settings$sigma.proc = sigma.proc
                jarainput$settings$sigma.obs.est= sigma.obs.est
                jarainput$settings$fixed.obsE = fixed.obsE
                jarainput$settings$variance.weighting = variance.weighting
                jarainput$settings$sets.var = sets.var
                jarainput$settings$nvar = nvar
                jarainput$settings$sigma.proc.fixed = sigma.proc.fixed  
                jarainput$settings$proc.pen = proc.pen
                # Projection stuff
                jarainput$settings$proj.mod =  proj.mod
                jarainput$settings$prjr.type = prjr.type
                jarainput$settings$proj.r = proj.r
                jarainput$settings$proj.stoch =proj.stoch
                jarainput$settings$pk.prior = pk.prior
                jarainput$settings$pk.yr = pk.yr
                jarainput$settings$pk.i = pk.i
                jarainput$settings$pk.mu =pk.mu
                jarainput$settings$pk.cv =pk.cv
                jarainput$settings$proj.yrs.user = proj.yrs.user
                jarainput$settings$Ninit = Ninit
                jarainput$settings$qs = qs
                
                if(model.type!="census") {jarainput$settings$q.init = q.init} else {
                  jarainput$settings$q.init = "Not available for census model"
                }
                jarainput$settings$mp.assess = mp.assess
                jarainput$settings$cols = cols
                #-------------------------------------------------------------------
                # write JAGS MODEL
                
                #jabba2jags(jbinput)
                
                return(jarainput)
                
                
} # end of build_jara()
  
            