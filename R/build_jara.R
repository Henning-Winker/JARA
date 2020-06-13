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
#' @param sigma.obs.est Estimates additional variance 
#' @param fixed.obsE minimum plausible process error (fixed)
#' @param sigma.proc.fixed option to fix the process error (default FALSE)
#' @param proc.pen advanced user setting to penalize extreme process error deviations
# Porjection settings
#' @param Klim penalty to restrict extrem values during projections c(TRUE, FALSE)
#' @param K.manual option to specify a carrying capacity for each census time series
#' @param prjr.type  # rate of change for c("all","GL1", years), "all" is default
#' @param proj.stoch allows for projections with process error c(TRUE, FALSE), FALSE is default
#' @return List to be used as data input to JARA JAGS model
#' @export
#' @author Henning Winker  
build_jara <- function(
              I = NULL, #abundance indices/counts,  require data.frame(year, I.1,I.2,...,I.N)
              se = NULL,
              assessment = "Unnamed",
              scenario = "s1",
              model.type = c("relative","census")[1],
              GL=NULL, # Generation length
              start.year = NA,
              end.year = NA, 
              sigma.obs.est = TRUE, 
              fixed.obsE = NULL,
              sigma.proc.fixed = FALSE,
              proc.pen = NULL,  
              Klim = FALSE,
              K.manual = FALSE,
              prjr.type = c("all","GL1","year")[1], 
              proj.stoch = FALSE){
  
  
  #-------------------------
  # Prepare input data
  #-------------------------
  
                if(is.null(GL)){GL = floor((nrow(I)-3)/3)} else {GL=GL[1]}
                
                GL1 = round(GL,0) # rounded for r.recent
                GL3 = round(3*GL,0) # 3 x GL rounded for year steps
                
                if(is.null(fixed.obsE)){
                fixed.obsE = ifelse(is.null(se),0.15,0.01)
                }
                if(is.null(proc.pen)){
                  proc.pen = c(ifelse(model.type=="census",1,0.5),log(ifelse(model.type=="census",0.1,0.5)),ifelse(model.type=="census",log(5),log(2)))
                }
                
                cat(paste0("\n","><> Prepare input data <><","\n","\n"))
                dat = I
                indices = names(dat)[2:ncol(dat)]
                n.indices = max(length(indices),1)
                if(is.na(start.year)) start.year = dat[1,1]   
                if(is.na(end.year)) end.year = dat[nrow(dat),1]   
                dat = subset(dat,dat[,1]>=as.numeric(start.year))  
                dat = subset(dat,dat[,1]<=as.numeric(end.year))  
                if(is.null(se)==TRUE){ SE.I=FALSE} else {SE.I=TRUE}
                
                if(SE.I==TRUE){
                    se = subset(se,se[,1]>=as.numeric(start.year))  
                  se = subset(se,se[,1]<=as.numeric(end.year))  
                }
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
                
                # Carrying capacity settings
                if(is.null(K.manual)[1]==FALSE){
                  if(K.manual[1]==FALSE) K.manual=NULL 
                }
                if(Klim==FALSE) Ks = 5
                if(Klim==TRUE & is.null(K.manual)==TRUE) Ks = 1.25
                if(Klim==TRUE & is.null(K.manual)==FALSE){
                  if(is.numeric(K.manual)==FALSE){
                    cat("ERROR: K.manual value(s) must be numeric")
                  }
                  Ks = K.manual  
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
                if(model.type=="census"){
                  if(is.null(K.manual)==TRUE | length(Ks)==1){ Ks = rep(Ks,n.indices)} 
                  
                  if(n.indices!=length(Ks))  cat("ERROR: K.manual value(s) must match the number of subpopulations") 
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
                if(prjr.type=="all" | prjr.type=="mean"){
                  prjr = 1:n.years
                } else if(prjr.type=="GL1" & GL1 >=n.years) {
                  prjr = (n.years-GL1+1):n.years  
                } else if(prjr.type=="GL1" & GL1 < n.years){
                  prjr = 1:n.years  
                } else {
                  prjr.st = round(as.numeric(prjr.type),0)  
                  prjr = (n.years-prjr.st+1):n.years
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
                  jags.data <- list(y = log(I_y+10^-20),SE2=se2, T = length(year),nI=n.indices,Ninit=Ninit,Ks=Ks,EY = n.years,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),igamma=igamma,penSig=0,prjr=prjr,proc.pen=proc.pen)
                  # order of indices
                  qs = 1:n.indices
                  # Parameters monitored
                  parameters <- c("mean.r","sigma","logNtot", "N.est","Ntot","r.tot","r.proj","TOE")
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
                  jags.data <- list(y = log(qI_y),SE2=qse2, logY1 = log(qI_y[1,1]), T = length(year),EY = n.years,nI=n.indices,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),igamma=igamma,penSig=0, Ks = Ks,prjr=prjr,proc.pen=proc.pen)
                  
                  # Parameters monitored
                  parameters <- c("mean.r", "sigma","r", "Y.est","Ntot","q","r.proj","TOE")
                  
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
                jarainput$data$I = I
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
                jarainput$settings$prjr.type = prjr.type
                jarainput$settings$proj.stoch =proj.stoch
                jarainput$settings$Ninit = Ninit
                jarainput$settings$qs = qs
                jarainput$settings$proc.pen = proc.pen
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
  
            