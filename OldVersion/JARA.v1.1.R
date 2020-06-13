

##############################################################################
# JARA: Just Another Redlisting Assessment
# Bayesian state-space redlisting tool 
# Model File: Run from Prime File (Settings)
# Developed by: Henning Winker* & Richard Sherley**
# * Department of Agriculture, Forestry and Fisheries (DAFF), ZAF \
# Email: henning.winker@gmail.com
# **University of Exeter, UK (richard.sherley@gmail.com)
# Email: richard.sherley@gmail.com
##############################################################################

cat(paste0("\n","><>><>><>><>><>><>><>><>><>><><>><>><>><>><>><>"))
cat(paste0("\n","><> Set up model ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Define objects to make sure they exist if not included in Prime file
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><
if(exists("posteriors")) rm(posteriors)
if(exists("jara")) rm(jara)
if(exists("start.year")==FALSE) start.year = NA   
if(exists("end.year")==FALSE) end.year = NA   
if(is.na(start.year)) start.year = dat[1,1]   
if(is.na(end.year)) end.year = dat[nrow(dat),1]   
if(exists("SE.I")==FALSE) SE.I = FALSE   
if(exists("K")==FALSE) K=FALSE   
if(exists("A1")==FALSE) A1=FALSE 
if(exists("plot.width")==FALSE) plot.width=5   
if(exists("plot.cex")==FALSE) plot.cex=1   
if(exists("sigma.est")==FALSE) sigma.est=TRUE   
if(exists("fixed.obsE")==FALSE) fixed.obsE=0.1   
if(exists("prjr.type")==FALSE) prjr.type="all"   
if(exists("proj.stoch")==FALSE) proj.stoch=FALSE   
if(exists("proc.pen")==FALSE){
  if(abundance=="census"){
    proc.pen = c(1,log(0.1),log(5))   
  }else{  
    proc.pen = c(0.5,log(0.5),log(2))   
  }}
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
# Process Error
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
#Estimate set sigma.proc == True
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
  sigma.proc = 0.1 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
}
if(is.null(K.manual)[1]==FALSE){
if(K.manual[1]==FALSE) K.manual=NULL 
}
# JABBA color palette
jabba.colors = as.character(c('#e6194b', "#3cb44b", "#ffe119",
                              "#0082c8","#f58231", "#911eb4",
                              "#46f0f0", "#f032e6", "#d2f53c",
                              "#fabebe", "#008080","#e6beff", "#aa6e28",rainbow(10)))



######################################################
# Prep data and Setup time horizon for trend analysis
######################################################
cat(paste0("\n","><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> Run Model ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>","\n","\n"))
# setwd(paste(File))
output.dir = paste0(File,"/",assessment,"/output",run)
dir.create(output.dir,showWarnings = FALSE)
# For only using data corresponding to 3 x generation  
# Set GL3.dat =TRUE and GL3.dat =FALSE otherwise 

GL1 = round(GL,0) # rounded for r.recent
GL3 = round(3*GL,0) # 3 x GL rounded for year steps

GL3.dat = FALSE
if(GL3.dat==TRUE) dat = dat[(nrow(dat)-GL3):nrow(dat),]

cat(paste0("\n","><> Prepare input data <><","\n","\n"))
indices = names(dat)[2:ncol(dat)]
n.indices = max(length(indices),1)

dat = subset(dat,dat[,1]>=as.numeric(start.year))  
dat = subset(dat,dat[,1]<=as.numeric(end.year))  
if(SE.I==TRUE){
se = subset(se,se[,1]>=as.numeric(start.year))  
se = subset(se,se[,1]<=as.numeric(end.year))  
}
years=dat[,1]
styr = min(years)
endyr = max(years)
n.years = length(years)

# Index count, biomass or relative
conv.I = as.numeric(as.matrix(dat[,-1]))
I_y=matrix(conv.I,nrow=n.years,ncol=n.indices)


if(SE.I==FALSE){
  se = dat  
  conv.se = as.numeric(as.matrix(dat[,-1]))
  se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
} else{
  conv.se = as.numeric(as.matrix(se[,-1]))
  #conv.se = sqrt(conv.se^2+fixed.obsE^2) 
  se2 = matrix(ifelse(is.na(conv.se),rep(0.01,length(conv.se)),conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
}

#------------------------------------------------------
# All input data plot + CIs
#------------------------------------------------------

Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
png(file = paste0(output.dir,"/AbundanceData_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)
options(warn=-1)
Ylim = c(0, max(exp(log(dat[,-1])+1.96*sqrt(se2))  ,na.rm =T))
plot(dat[,1],dat[,1],type="n",xlim=c(min(dat[,1]-1),max(dat[,1]+1)),ylab=ifelse(abundance=="census","Counts","Abunance Index"),xlab="Year",ylim=Ylim, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
for(j in 1:ncol(I_y)){
  plotCI(x = dat[,1]+runif(1,-0.2,0.2),y = dat[,j+1],type="b",ui=ifelse(is.na(dat[,j+1]),NA,exp(log(dat[,j+1])+1.96*sqrt(se2[,j]))) ,li=ifelse(is.na(dat[,j+1]),NA,exp(log(dat[,j+1])-1.96*sqrt(se2[,j]))),gap=0.1,pch=21,cex=1.2,col=grey(0.4,0.7),pt.bg=jabba.colors[j],add=T)
}
axis(1,at=seq(min(dat[,1]),max(dat[,1])+5,ceiling(length(dat[,1])/8)),tick=seq(min(dat[,1]),max(dat[,1]),ceiling(length(dat[,1])/8)),cex.axis=0.9)

posl = c(max(dat[1:3,-1],na.rm=T),max(dat[(length(years)):length(years),-1],na.rm=T))
legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(names(dat)[2:(n.indices+1)]), lty = c(1, rep(n.indices)), lwd = c(rep(-1,n.indices)),pch=c(rep(21,n.indices)), pt.bg = c(jabba.colors[1:n.indices]), bty = "n", cex = 0.9,y.intersp = 0.8)
dev.off()
options(warn=0)


# Carrying capacity settings
if(Klim==FALSE) Ks = 5
if(Klim==TRUE & is.null(K.manual)==TRUE) Ks = 1.25
if(Klim==TRUE & is.null(K.manual)==FALSE){
  if(is.numeric(K.manual)==FALSE){
    cat("ERROR: K.manual value(s) must be numeric")
   }
  Ks = K.manual  
} 


if(abundance=="census"){
  if(is.null(K.manual)==TRUE | length(Ks)==1){ Ks = rep(Ks,n.indices)} 
  
  if(n.indices!=length(Ks))  cat("ERROR: K.manual value(s) must match the number of subpopulations") 
  }

#Start Year
styr = years[1]
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
se2= rbind(se2,matrix(0.1,pyears,(n.indices))) # now always




#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Setup model for population counts
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


Ninit = NULL
# get starting values Ninit
for(i in 1:n.indices){
Ninit[i] = na.omit(I_y[,i])[1]  
}

if(abundance=="census"){
# Bundle data for jags
jags.data <- list(y = log(I_y+10^-20),SE2=se2, T = length(year),nI=n.indices,Ninit=Ninit,Ks=Ks,EY = n.years,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),igamma=igamma,penSig=0,prjr=prjr,proc.pen=proc.pen)

# Initial values for jags
inits <- function(){list(mean.r = rnorm(n.indices,0,0.5),isigma2.est=runif(1,20,100), itau2=runif(1,80,200), logN.est =matrix(log(rbind(Ninit,matrix(rep(NA,(nT-1)*n.indices),(nT-1),n.indices))),nT,n.indices))  }

qs = 1:n.indices

# Parameters monitored
parameters <- c("mean.r","sigma","logNtot", "N.est","Ntot","r.tot","r.proj")

sink(paste0(output.dir,"/jara.jags"))
cat("
    model {
    # Priors and constraints
    for(i in 1:nI)
    {
    mean.r[i] ~ dnorm(0, 0.001) 
    logN.est[1,i] ~ dnorm(log(Ninit[i]),pow(0.5,-2))   # Prior for initial population size with CV =100%
    }
    
    
    
    ")

if(sigma.proc==TRUE){
  cat("
      # Process variance
      isigma2 <- isigma2.est 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
      ",append=TRUE)  
}else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)
           
           ",append=TRUE)}

if(sigma.est==TRUE){
  cat("
      # Obsevation variance
      # Observation error
      itau2~ dgamma(0.001,0.001)
      tau2 <- 1/itau2
      
      
      for(i in 1:nI)
      {
      for(t in 1:T)
      {
      var.obs[t,i] <- SE2[t,i]+tau2
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i])
      
      }}
      ",append=TRUE)  
}else{ cat(" 
      # Obsevation variance
           # Observation error
           itau2~ dgamma(2,2)
           tau2 <- 1/itau2
           
           
           for(i in 1:nI)
           {
           for(t in 1:T)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2
           
           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)     
           TOE[t,i] <- sqrt(var.obs[t,i])
           
           }}
           
           ",append=TRUE)}

# Run rest of code  
cat("  
    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])
    pen.sigma <- ifelse(sigma>proc.pen[1],log(sigma)-log(proc.pen[1]),0) 
    penSig  ~ dnorm(pen.sigma,pow(0.2,-2))
    
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    for (i in 1:nI){
    logN.est[t+1,i] <- logN.est[t,i]+r[t,i]
    }}
    # Set last year r to r.mean
    for (i in 1:nI){
    for(t in 1:(EY-1)){
    rdev[t,i] ~ dnorm(0, isigma2) #T(proc.pen[2],proc.pen[3])
    r[t,i] <- mean.r[i]+ rdev[t,i]-0.5*sigma2    #dnorm(mean.r[i], isigma2)
    }}
    
               ",append=TRUE)

if(prjr.type=="mean"){
  cat("  
    for (i in 1:nI){
    r.proj[i] <- mean.r[i]
    } 
      
     ",append=TRUE)}else{
   cat(" 
   for (i in 1:nI){
   r.proj[i] <- mean(mean.r[i]+ rdev[prjr,i]-0.5*sigma2)} 
   ",append=TRUE)  
}


if(proj.stoch==FALSE){
    cat("  
    for (i in 1:nI){
    for(t in EY:(T-1)){
    rdev[t,i] ~ dnorm(0, pow(0.01,-2))
    r[t,i] <- ifelse(logN.est[t,i]>max(logN.est[1:(EY-1),i])+log(Ks[i]),0,r.proj[i]) 
      
    }}
    ",append=TRUE)}else{
      cat(" 
    for (i in 1:nI){
    for(t in EY:(T-1)){
    rdev[t,i] ~ dnorm(0, isigma2)T(proc.pen[2],proc.pen[3])
    r[t,i] <- ifelse(logN.est[t,i]>max(logN.est[1:(EY-1),i])+log(Ks[i]),-0.01,r.proj[i]+ rdev[t,i]-0.5*sigma2) 
    }}
    ",append=TRUE)  
    }
     
    cat("  
    #for(t in 1:(T)){
    #for (i in 1:nI){
    #devK[t,i]  <- ifelse(logN.est[t,i]>log(Ks[i]),logN.est[t,i]-log(Ks[i]),0) # penalty if N > K 
    #penK[t,i] ~ dnorm(devK[t,i],pow(0.1,-2))
    #}}
    
    
    lambda[1] <- 1
    # Observation process
    for (t in 1:T) {
    for(i in 1:nI){
    y[t,i] ~ dnorm(logN.est[t,i], ivar.obs[t,i])
    }}
    
    # Population sizes on real scale
    for (t in 1:T) {
    for(i in 1:nI){
    N.est[t,i] <- exp(logN.est[t,i])
    }
    Ntot[t] <- sum(N.est[t,])
    logNtot[t] <- log(Ntot[t])
    }
    for(t in 2:T)
    {
    r.tot[t] <- logNtot[t]-logNtot[t-1]
    }
    } 
    ",fill = TRUE)
sink()

}

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Setup model for relative abundance
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


if(abundance=="relative"){
   
  #find first time-series with first CPUE
  q1.y = c(1:length(year))[is.na(apply(I_y,1,mean,na.rm=TRUE))==FALSE][1] #first year with CPUE
  q1.I =  which.max(I_y[q1.y,])
  qs = c(q1.I,c(1:(ncol(I_y)))[-q1.I])
  
  qI_y = as.matrix(I_y[,qs])
  qse2 = as.matrix(se2[,qs])
  
  q.init = 1
  if(n.indices>1) for(i in 2:n.indices){q.init[i] = mean(qI_y[,i],na.rm=TRUE)/mean(qI_y[,1],na.rm=TRUE)}
  
  # Bundle data
  jags.data <- list(y = log(qI_y),SE2=qse2, logY1 = log(qI_y[1,1]), T = length(year),EY = n.years,nI=n.indices,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),igamma=igamma,penSig=0, Ks = Ks,prjr=prjr,proc.pen=proc.pen)
  
  # Initial values
  inits <- function(){list(isigma2.est=runif(1,20,100), itau2=runif(1,80,200), mean.r = rnorm(0,0.2),iq = 1/q.init)}
  
  # Parameters monitored
  parameters <- c("mean.r", "sigma","r", "Y.est","Ntot","q","r.proj")
  
  sink(paste0(output.dir,"/jara.jags"))
  cat("
    model {
      
      # Prior specifications  
      eps <- 0.0000000000001 # small constant    
      
      iq[1] ~ dgamma(1000,1000)
      q[1] <-  pow(iq[1],-1)
      logq[1] <- log(1)
      for(i in 2:nI){
      iq[i] ~ dgamma(0.001,0.001)
      q[i] <- pow(iq[i],-1)
      logq[i] <-  log(q[i])
      }
      
      
      ")
  
  if(sigma.proc==TRUE){
    cat("
        # Process variance
        isigma2 <- isigma2.est 
        sigma2 <- pow(isigma2,-1)
        sigma <- sqrt(sigma2)
        fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
        ",append=TRUE)  
  }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
             sigma2 <- pow(isigma2,-1)
             sigma <- sqrt(sigma2)
             
             ",append=TRUE)}
  
  if(sigma.est==TRUE){
    cat("
        # Obsevation variance
        # Observation error
        itau2~ dgamma(0.001,0.001)
        tau2 <- 1/itau2
        
        
        for(i in 1:nI)
        {
        for(t in 1:T)
        {
        var.obs[t,i] <- SE2[t,i]+tau2
        ivar.obs[t,i] <- 1/var.obs[t,i]
        # note total observation error (TOE)     
        TOE[t,i] <- sqrt(var.obs[t,i])
        
        }}
        ",append=TRUE)  
  }else{ cat(" 
      # Obsevation variance
             # Observation error
             itau2~ dgamma(2,2)
             tau2 <- 1/itau2
             
             
             for(i in 1:nI)
             {
             for(t in 1:T)
             {
             var.obs[t,i] <- SE2[t,i] # drop tau2
             fake.tau[t,i] <- tau2
             
             ivar.obs[t,i] <- 1/var.obs[t,i]
             # note total observation error (TOE)     
             TOE[t,i] <- sqrt(var.obs[t,i])
             
             }}
             
             ",append=TRUE)}
  
  # Run rest of code  
  cat("  
      # Process variance prior
      isigma2.est ~ dgamma(igamma[1],igamma[2])
      pen.sigma <- ifelse(sigma>proc.pen[1],log(sigma)-log(proc.pen[1]),0) 
      penSig  ~ dnorm(pen.sigma,pow(0.2,-2))
      
      # Priors and constraints
      logY.est[1] ~ dnorm(logY1, 1)       # Prior for initial population size
      
      mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
      
      # Likelihood
      # State process
      for (t in 1:(EY-1)){
      rdev[t] ~ dnorm(0, isigma2)T(proc.pen[2],proc.pen[3])
      r[t] <- mean.r+rdev[t]-0.5*sigma2   
      logY.est[t+1] <- logY.est[t] + r[t] 
      }

      ",append=TRUE)
  
    if(prjr.type=="mean"){
    cat("  
      r.proj <- mean.r  
      prjr.dummy <- prjr   
    ",append=TRUE)}else{
    cat(" 
        r.proj <- mean(mean.r+rdev[prjr]-0.5*sigma2) 
        ",append=TRUE)  
    }
  
  
     if(proj.stoch==FALSE){
      cat("
      for (t in EY:(T)){
      rdev[t] ~ dnorm(0, pow(0.01,-2))
      r[t] <- ifelse(logY.est[t]>(max(logY.est[1:(EY-1)])+log(Ks)),0,r.proj) 
      logY.est[t+1] <- logY.est[t]+r[t]   
      }
      ",append=TRUE)} else {
      cat("
      for (t in EY:(T)){
            rdev[t] ~ dnorm(0, pow(isigma2,-2))T(0.5*proc.pen[2],0.5*proc.pen[3])
            r[t] <- ifelse(logY.est[t]>(max(logY.est[1:(EY-1)])+log(Ks)),0,r.proj+rdev[t]-0.5*sigma2) 
            logY.est[t+1] <- logY.est[t] + r[t] 
            #devK[t]  <- ifelse(logY.est[t+1]>log(Ks),logY.est[t+1]-log(Ks),0) # penalty if Y > K 
      }
      ",append=TRUE)}  
        
  cat("

      #for(t in (EY):T){ # Apply penalty for projections only
      #penK[t] ~ dnorm(devK[t],pow(0.1,-2))
      #}
      
      # Observation process
      for (t in 1:EY) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      }}
      for (t in (EY+1):T) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i],pow(0.001,-2)) # does not effect projections
      }}
      

      # Population sizes on real scale
      for (t in 1:T) {
      Y.est[t] <- exp(logY.est[t])
      Ntot[t] <- Y.est[t]
      }
      
  } 
      ",fill = TRUE)
  sink()
  
} 
  
 
  
# Call JAGS from R (BRT 3 min)# Call JAGS from R (BRT 3 min)
ptm <- proc.time()
jara.mod <- jags(jags.data, inits, parameters, paste0(output.dir,"/jara.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

proc.time() - ptm
save.time = proc.time() - ptm



RunTime =paste0("\n",paste0("><> Run  ",assessment," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n"))
cat(RunTime)
cat(paste0("\n","><> Produce results output of ",assessment," Run ",run," <><","\n"))


#----------------
# Set up output
#---------------
end.yr=n.years

# remove scientific numbers
options(scipen=999)

if(n.years<GL1){
WARN = paste0("\n","WARNING!!!","\n","-----------------------------------------------------\n",
"\n","The projections are based on less than one Generation Time!",
"\n","\n","The results are unreliable and MUST be interpreted with caution","\n","\n")
} else {WARN=""}



#--------------------------
# Capture Settings
#--------------------------
Settings = NULL
Settings$assessment = assessment
Settings$abundance = abundance
Settings$generation.length=GL
Settings$K.penalty=K
index.se=SE.I
Settings$sigma.obs.est=sigma.est
Settings$sigma.obs.add= fixed.obsE
Settings$sigma.proc.fixed = sigma.proc.fixed 
Settings$prjr.type = prjr.type 
Settings$proj.stoch = proj.stoch
Settings$start.year=start.year
Settings$end.year=end.year
capture.output( Settings, file=paste0(output.dir,"/Settings_",assessment,run,".txt"))
#-------------------------------------------------------------------
# Capture results
posteriors=jara.mod$BUGSoutput$sims.list
if(abundance=="census"){sel.par=c(1:2)} else {sel.par=c(1:2,6)}
par.dat= data.frame(posteriors[parameters[sel.par]])
geweke = geweke.diag(data.frame(par.dat))
pvalues <- 2*pnorm(-abs(geweke$z))
pvalues
heidle = heidel.diag(data.frame(par.dat))

# Capture Results
results = round(t(cbind(apply(par.dat,2,quantile,c(0.025,0.5,0.975)))),3)
pars = data.frame(Median = results[,2],LCI=results[,1],UCI=results[,3],Geweke.p=round(pvalues,3),Heidelberger.p=round(heidle[,3],3))
if(abundance=="relative"){par.names = c("mu.r","sigma.proc",paste0("q.",indices[qs]))} else {
  par.names = c(paste0("r.",indices[qs]),"sigma.proc")  
}
row.names(pars) = par.names
write.csv(pars,paste0(output.dir,"/Parameters_",assessment,".csv"))


options(max.print=1000000)

sink(paste0(output.dir,"/JARAjags_",assessment,".txt"))
cat(WARN)
cat("-----------------------------------------------------\n")
cat(RunTime)
cat("-----------------------------------------------------\n")

cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> Input ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
print(Settings)
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> Parameter Estimates ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
print(pars)
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> JAGS output ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
print(jara.mod)

sink()
options(max.print=1000)


# get individual trends
fitted <- lower <- upper <- mat.or.vec(end.yr,(n.indices))

if(abundance=="census"){
for(i in 1:n.indices){
for (t in 1:end.yr){
  fitted[t,i] <- median(posteriors$N.est[,t,i])
  lower[t,i] <- quantile(posteriors$N.est[,t,i], 0.025)
  upper[t,i] <- quantile(posteriors$N.est[,t,i], 0.975)}}
} else {
fitted <- lower <- upper <- as.null()
  for (t in 1:end.yr){
    fitted[t] <- median(posteriors$Y.est[,t])
    lower[t] <- quantile(posteriors$Y.est[,t], 0.025)
    upper[t] <- quantile(posteriors$Y.est[,t], 0.975)}
}


Nfit <- Nlow <- Nhigh <- as.numeric()
# get total pop size
for (t in 1:nT){
Nfit[t] =  median(posteriors$Ntot[,t])
Nlow[t] = quantile(posteriors$Ntot[,t],0.025)
Nhigh[t] = quantile(posteriors$Ntot[,t],0.975)
}

#Abundance FITS
Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
png(file = paste0(output.dir,"/Fits_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)


m1 <- 0
m2 <- max(c(fitted, upper*1.1), na.rm = TRUE)


if(abundance=="census"){
plot(0, 0, ylim = c(m1, m2), xlim = c(min(years-1),max(years+1)), ylab = "Population Numbers", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- jabba.colors
  
for(i in 1:n.indices) polygon(x = c(years,rev(years)), y = c(lower[,i],upper[end.yr:1,i]), col = gray(runif(1,0.5,0.9),0.5), border = "gray90")

for(i in 1:n.indices)
{
lines(years,fitted[,i], type = "l",col=col_line[i], lwd=1)
points(years,dat[,i+1], bg = col_line[i],pch=21) 
}
posl = c(max(fitted[1,]),max(fitted[length(years),]))
} else {  
plot(0, 0, ylim = c(m1, m2), xlim =  c(min(years-1),max(years+1)), ylab = "Abudance Index", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- jabba.colors
q.adj = apply(posteriors$q,2,median)

polygon(x = c(years,rev(years)), y = c(lower,upper[end.yr:1]), col = "gray", border = "gray90")


for(i in 1:n.indices)
{
points(year,I_y[,qs[i]]/q.adj[i], bg = col_line[i],col=col_line[i],pch=21,type="b") 
}
lines(years,fitted, type = "l",col=1, lwd=2)
posl = c(max(fitted[1]),max(fitted[length(years)]))
}
 
legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(names(dat)[2:(n.indices+1)]), lty = c(1, rep(n.indices)), lwd = c(rep(-1,n.indices)),pch=c(rep(21,n.indices)), pt.bg = c(jabba.colors[1:n.indices]), bty = "n", cex = 0.9,y.intersp = 0.8)
axis(1,at=seq(min(years)-1,max(years)+5,ceiling(n.years/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
dev.off()

#------------------------------------------------------------------
# Goodness of fit Statistics
#------------------------------------------------------------------


DIC =round(jara.mod$BUGSoutput$DIC,1)

# get residuals
Resids = NULL
for(i in 1:n.indices){
  if(abundance=="census") Resids =rbind(Resids,log(dat[,i+1])-log(fitted[,i]))   
  if(abundance=="relative") Resids =rbind(Resids,log(dat[,qs[i]+1]/q.adj[i])-log(fitted))   
  }

Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/(Nobs-1)),1)
# Produce statistice describing the Goodness of the Fit

GOF = data.frame(Stastistic = c("Nobs","RMSE","DIC"),Value = c(Nobs,RMSE,DIC))
write.csv(GOF,paste0(output.dir,"/GoodnessFit_",assessment,".csv"))

# JABBA-residual plot
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=plot.cex)
png(file = paste0(output.dir,"/Residuals_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)
plot(years,years,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),ylab="log Residuals",xlab="Year")
boxplot(Resids,add=TRUE,at=c(years),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)
for(i in 1:n.indices){
  for(t in 1:n.years){
  lines(rep((years+positions[i])[t],2),c(0,Resids[i,t]),col=jabba.colors[i])}
  points(years+positions[i],Resids[i,],col=1,pch=21,bg=jabba.colors[i])}
  mean.res = apply(Resids,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~years),data.frame(Yr=years))
lines(years,smooth.res,lwd=2)
legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
legend("bottomright",paste(runs), legend = c("Loess",paste(indices[qs])), lty = c(1, rep(-1,n.indices)), lwd = c(2, rep(-1,n.indices)),pch=c(-1,rep(21,n.indices)), pt.bg = c(1,col_line), bty = "n", cex = 0.9,y.intersp = 0.8)
dev.off()


#log Index FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/logFits_",assessment,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  if(abundance=="census"){fit = fitted[,i]} else {fit = fitted*q.adj[[i]]}
  
  I.i = I_y[is.na(I_y[,qs[i]])==F,qs[i]]
  yr.i = Yr[is.na(I_y[,qs[i]])==F]
  se.i = sqrt(se2[is.na(I_y[,qs[i]])==F,(i)])
  
  ylim = log(c(min(fit[yr[is.na(I_y[,qs[i]])==F]]*0.8,exp(log(I.i)-1.96*se.i)), max(fit[yr[is.na(I_y[,i])==F]]*1.5,exp(log(I.i)+1.96*se.i))))
  
  # Plot Observed vs predicted CPUE
  plot(years,log(dat[,qs[i]+1]),ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  #polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  
  
  lines(years,log(fit),lwd=2,col=4)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,log(I.i),ui=log(exp(log(I.i)+1.96*se.i)),li=log(exp(log(I.i)-1.96*se.i)),add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,log(I.i),pch=21,xaxt="n",yaxt="n",bg="white")}
  legend('topright',paste(indices[qs[i]]),bty="n",y.intersp = -0.2,cex=1)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Log Abundance"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()







Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
png(file = paste0(output.dir,"/PopTrend_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)

# Total N
m1 <- 0
m2 <- max(Nhigh, na.rm = TRUE)

plot(0, 0, ylim = c(m1, m2), xlim = c(min(year-1),max(year+1)), ylab = paste("Population",ifelse(abundance=="census","size","trend")), xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n",xaxs="i",yaxs="i")
axis(1,at=seq(min(year),max(year)+5,ceiling(length(year)/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)

polygon(x = c(year,rev(year)), y = c(Nlow,rev(Nhigh)), col = gray(0.6,0.3), border = "gray90")
#add trend
polygon(x = c(1:end.yr,end.yr:1), y = c(Nlow[1:end.yr],Nhigh[end.yr:1]), col = gray(0.7,0.3),border = "gray90")
lines(year[end.yr:nT],Nfit[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
lines(years,Nfit[1:end.yr], type = "l",col=1,lwd=2)
if(n.years-3*GL-1>0) lines(rep(year[n.years]-3*GL,2),c(0,m2*.93),lty=2,col=2)
lines(rep(year[n.years]-1*GL,2),c(0,m2*.93),lty=2,col=4,lwd=2)
lines(rep(year[n.years]-2*GL,2),c(0,m2*.93),lty=2,col=3,lwd=2)
lines(rep(year[n.years],2),c(0,m2)*.93,lty=2,col=1,lwd=2)
#lines(rep(year[mp.assess[2]],2),c(0,m2)*.93,lty=5,col=1)
#lines(rep(year[mp.assess[1]],2),c(0,m2)*.93,lty=2)
if(n.years-3*GL-1>0) text(year[n.years]-3*GL,m2*.96,"-3GL",lwd=2)
text(year[n.years]-2*GL,m2*.96,"-2GL")
text(year[n.years]-GL,m2*.96,"-1xGL")
text(year[n.years],m2*.96,paste0(year[n.years]))
dev.off()

#### IUCN plot
Par = list(mfrow=c(1,1),mar = c(0.5, 1, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.3, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/IUCNplot_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)
change = (apply(posteriors$Ntot[,(mp.assess[2]-1):(mp.assess[2]+1)],1,median)/apply(posteriors$Ntot[,(mp.assess[1]-1):(mp.assess[1]+1)],1,median)-1)*100
den = density(change,adjust=2)
x1 = den$x
y1 = den$y

mu.change = round(median(change),1)
sign=""
if(mu.change>0) sign="+"
CR = round(sum(ifelse(change< ifelse(A1,-90,-80),1,0))/length(change)*100,1)
EN = round(sum(ifelse(change> ifelse(A1,-90,-80) & change< ifelse(A1,-70,-50),1,0))/length(change)*100,1)
VU = round(sum(ifelse(change> ifelse(A1,-70,-50) & change< ifelse(A1,-50,-30),1,0))/length(change)*100,1)
LC = round(sum(ifelse(change> -30,1,0))/length(change)*100,1)
Decline = round(sum(ifelse(change> 0,1,0))/length(change)*100,1)
plot(x1,y1,type="n",xlim=c(-100,min(max(30,quantile(change,.99)),1000)),ylim=c(0,max(y1*1.05)),ylab="Density",xlab="",cex.main=0.9,frame=T,xaxt="n",yaxt="n",xaxs="i",yaxs="i")

maxy = max(y1*1.08)
x2 = c(ifelse(A1,-50,-30),1500); y2 = c(0,5)
polygon(c(x2,rev(x2)),c(rep(maxy,2),rev(rep(0,2))),col="green",border="green")
x3 = c(ifelse(A1,-70,-50),x2[1])
polygon(c(x3,rev(x3)),c(rep(maxy,2),rev(rep(0,2))),col="yellow",border="yellow")
x4 = c(ifelse(A1,-90,-80),x3[1])
polygon(c(x4,rev(x4)),c(rep(maxy,2),rep(0,2)),col="orange",border="orange")
x5 = c(-100,x4[1])
polygon(c(x5,rev(x5)),c(rep(maxy,2),rep(0,2)),col="red",border="red")

polygon(c(x1,rev(x1)),c(y1,rep(0,length(y1))),col="grey")
axis(1,at=seq(-100,max(x1,30)+50,ifelse(max(x1,30)>150,50,25)),tick=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)))
mtext(paste("Density"), side=2, outer=T, at=0.55,line=-1.2,cex=1)
mtext(paste("Change (%)"), side=1, outer=T, at=0.5,line=-1.5,cex=1)
legend("right",c(paste0("CR (",CR,"%)"),paste0("EN (",EN,"%)"),
                 paste0("VU (",VU,"%)"),paste0("LC (",LC,"%)")),col=1,pt.bg=c("red","orange","yellow","green"),pt.cex=1.4,pch=22,bg="white",cex=1.1)
text(ifelse(mean(change)< -80,-80,mean(change)),max(y1*1.03),paste0("Change = ",sign,mu.change,"%"),bg="white",cex=1.2)
dev.off()
# End IUCN Plot

Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
png(file = paste0(output.dir,"/AnnualRate_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)
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

lamdas = (exp(rs)-1)*100

lymax=rymax = lxrange = rxrange =NULL # maximum and range for plotting
for(i in 1:ncol(rs)){
den = density(lamdas[,i],adjust=2)
assign(paste0("xl",i),den$x)
assign(paste0("yl",i),den$y)
lymax=c(lymax,max(den$y))
lxrange = c(lxrange,range(den$x))
den = density(rs[,i],adjust=2)
assign(paste0("xr",i),den$x)
assign(paste0("yr",i),den$y)
rymax=c(rymax,max(den$y))
rxrange = c(rxrange,range(den$x))
}
cnam = c("All.yrs","1GL","2GL","3GL")  
 
jcol = c(grey(0.5,0.6),rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
plot(0,0,type="n",ylab="Density",xlab="Annual Rate of Change(%)",xaxt="n",cex.main=0.9,ylim=c(0,1.1*max(lymax)),xlim=quantile(lamdas,c(0.001,0.999)),xaxs="i",yaxs="i") 
for(i in 1:ncol(rs)){
  x = get(paste0("xl",i))
  y = get(paste0("yl",i))
  polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
  mu.lamda = round(median(lamdas[,i]),10)
  lines(rep(mu.lamda,2),c(0,max(y)),col=c(1,4,3,2)[i],lwd=1,lty=1)
  
}
axis(1,at=seq(floor(min(x)),ceiling(max(x)),1),tick=seq(min(year),max(year),5),cex.axis=0.9)
    


abline(v=0,lty=2)
legend("topright", paste0(cnam[1:ncol(rs)]," = ",ifelse(round(apply(lamdas,2,median),2)>0,"+",""),round(apply(lamdas,2,median),2),"%"),pch=15,col=c(jcol),bty="n")
dev.off()

Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
png(file = paste0(output.dir,"/logr_",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)
plot(0,0,type="n",ylab="Density",xlab="log rate r",cex.main=0.9,ylim=c(0,1.1*max(rymax)),xlim=quantile(rs,c(0.001,0.999)),xaxs="i",yaxs="i")
for(i in 1:ncol(rs)){
  x = get(paste0("xr",i))
  y = get(paste0("yr",i))
  polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
  mu.r=round(median(rs[,i]),10)
  lines(rep(mu.r,2),c(0,max(y)),col=c(1,4,3,2)[i],lwd=1,lty=1)
  
}
abline(v=0,lty=2)

legend("topright", paste0(cnam[1:ncol(rs)]," = ",ifelse(round(apply(rs,2,median),3)>0,"+",""),round(apply(rs,2,median),3)),pch=15,col=c(jcol),bty="n")
dev.off()

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Save key posterior results
trends =  data.frame(change,rs)
colnames(trends ) = c("pop.change",paste0("r.",cnam[1:ncol(rs)]))
abundance.est=data.frame(yr=year,Type=abundance,estimation=ifelse(year %in% years,rep("fit",length(year)) ,rep("prj",length(year))), total=Nfit,total.lcl=Nlow,total.ucl=Nhigh)

if(abundance=="census"){
  r.i = posteriors$mean.r
  colnames(r.i) = paste0("r.",names(dat[,-1]))
  abund.i = mat.or.vec(nT,n.indices)  
  for(i in 1:n.indices){
  for (t in 1:nT){
        abund.i[t,i] <- median(posteriors$N.est[,t,i])}}
  
  
  colnames(abund.i) = names(dat[,-1])
  abundance.est = data.frame(abundance.est,abund.i)
  }


categories = c("CR","EN","VU","LC")  
percentages = c(CR,EN,VU,LC)
status= ifelse(which(percentages==max(percentages))==4 & max(percentages)<50,"NT",categories[which(percentages==max(percentages))])
jara = list(settings=Settings,data=jags.data,parameters=pars,trends=trends,abundance=abundance.est,perc.risk=data.frame(CR=CR,EN=EN,VU=VU,LC=LC) ,status=status)
save(jara,file=paste0(output.dir,"/jara.",assessment,".rdata"))

# Save estimated population trajectory
write.csv(abundance.est,paste0(output.dir,"/Pop.trajectory.csv"),row.names = FALSE)
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
rm(proc.pen)
cat(paste0("\n","><> Complete: ",assessment," - Run ",run," - Type: ",abundance," <><"))
if(n.years<GL1){
cat(paste0("\n","WARNING!!!"))
cat(paste0("\n","The projections are based on less than one Generation Time!"))
cat(paste0("\n","The results are unreliable and MUST be interpreted with caution"))
}



