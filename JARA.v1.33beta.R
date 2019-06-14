

##############################################################################
# JARA: Just Another Redlist tool
# Bayesian state-space approach for trend analysis 
# Developed by: Henning Winker (henning.winker@gmail.com)
##############################################################################


cat(paste0("\n","><>><>><>><>><>><>><>><>><>><><>><>><>><>><>><>"))
cat(paste0("\n","><> Set up model ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Define objects to make sure they exist if not included in Prime file
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
if(exists("start.year")==FALSE) start.year = NA   
if(exists("SE.I")==FALSE) SE.I = FALSE   
if(exists("K")==FALSE) K=FALSE   
if(exists("sigma.est")==FALSE) sigma.est=TRUE   
if(exists("fixed.obsE")==FALSE) fixed.obsE=0.1   
if(exists("proj.eval")==FALSE) proj.eval="all"   
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

if(K.manual==FALSE) K.manual=NULL 

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
# Set GT3.dat =TRUE and GT3.dat =FALSE otherwise 

GT1 = round(GT,0) # rounded for r.recent
GT3 = round(3*GT,0) # 3 x GT rounded for year steps

GT3.dat = FALSE
if(GT3.dat==TRUE) dat = dat[(nrow(dat)-GT3):nrow(dat),]

cat(paste0("\n","><> Prepare input data <><","\n","\n"))
indices = names(dat)[2:ncol(dat)]
n.indices = max(length(indices),1)

if(is.na(start.year)==FALSE){
dat = subset(Year>=as.numeric(start.year))  
if(SE.I==TRUE){
se = subset(se,Year>=as.numeric(start.year))  
}  
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
  se2 = matrix(ifelse(is.na(conv.se),rep(0.3,length(conv.se)),conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
}

# Carrying capacity settings
if(Klim==FALSE) Ks = 5
if(Klim==TRUE & is.null(K.manual)==TRUE) Ks = 1.25
if(Klim==TRUE & is.null(K.manual)==FALSE) Ks = K.manual
if(abundance=="census" & is.null(K.manual)==TRUE) Ks = rep(Ks,n.indices)
#Start Year
styr = years[1]
# prediction years 
pyears = ifelse(GT3>(n.years+2),GT3-n.years+2,1) # +2 required for avergaging start + end with y = -1 and +1 
if(GT3==n.years |GT3==n.years+1) pyears = 2

# set up year vectors
year <- years[1]:(years[length(years)] + pyears)
yr = years[1]:(years[length(years)])

nT = length(year)
mp.assess = c(nT-GT3,nT-1)     #start and end midpoints for population trend assessment

#define the r type used for projections
if(prjr.type=="all" | prjr.type=="mean"){
prjr = 1:n.years
} else if(prjr.type=="GT1" & GT1 >=n.years) {
prjr = (n.years-GT1+1):n.years  
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
    r[t,i] <- ifelse(logN.est[t,i]>max(logN.est[1:(EY-1),i])+log(Ks[i]),0,r.proj[i]+ rdev[t,i]-0.5*sigma2) 
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

if(n.years<GT1){
WARN = paste0("\n","WARNING!!!","\n","-----------------------------------------------------\n",
"\n","The projections are based on less than one Generation Time!",
"\n","\n","The results are unreliable and MUST be interpreted with caution","\n","\n")
} else {WARN=""}

input = t(data.frame(assessment,abundance,generation.time=GT,K.penalty=K,index.se=SE.I
                   ,sigma.obs.est=sigma.est,sigma.obs.add= fixed.obsE,sigma.proc.fixed,prjr.type,proj.stoch,start.year=years[1]))
colnames(input) = "Settings"

options(max.print=1000000)


sink(paste0(output.dir,"/JARAout_",assessment,".txt"))
cat(WARN)
cat("-----------------------------------------------------\n")
print(jara.mod)
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> Input ",assessment," ",abundance," <><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>><>><>><>><>><>><>","\n","\n"))
print(input)
cat("-----------------------------------------------------\n")
cat(RunTime)
sink()

options(max.print=1000)


posteriors=jara.mod$BUGSoutput$sims.list


par.dat= data.frame(posteriors[parameters[c(1:2)]])
geweke = geweke.diag(data.frame(par.dat))
pvalues <- 2*pnorm(-abs(geweke$z))
pvalues

heidle = heidel.diag(data.frame(par.dat))

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

dir.create(paste0(getwd(),"/",assessment),showWarnings = FALSE)

#Abundance FITS
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(0.7, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(output.dir,"/Fits_",assessment,".png"), width = 4.5, height = 4, 
    res = 200, units = "in")
par(Par)


m1 <- 0
m2 <- max(c(fitted, upper*1.1), na.rm = TRUE)


if(abundance=="census"){
plot(0, 0, ylim = c(m1, m2), xlim = range(years), ylab = "Population size", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE)
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- jabba.colors
  
for(i in 1:n.indices) polygon(x = c(years,rev(years)), y = c(lower[,i],upper[end.yr:1,i]), col = gray(runif(1,0.5,0.9),0.5), border = "gray90")

for(i in 1:n.indices)
{
lines(years,fitted[,i], type = "l",col=col_line[i], lwd=1)
points(years,dat[,i+1], bg = col_line[i],pch=21) 
}
} else {  
plot(0, 0, ylim = c(m1, m2), xlim = range(years), ylab = "Abudance Index", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE)
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
}
legend('topleft', legend = c("Fit",paste(indices[qs])), lty = c(1, rep(-1,n.indices)), lwd = c(2, rep(-1,n.indices)),pch=c(-1,rep(21,n.indices)), pt.bg = c(1,col_line), bty = "n", cex = 1)

dev.off()

Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(0.7, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(output.dir,"/PopTrend_",assessment,".png"), width = 4.5, height = 4, 
    res = 200, units = "in")
par(Par)

# Total N
m1 <- 0
m2 <- max(Nhigh, na.rm = TRUE)
if(abundance=="census"){
  plot(0, 0, ylim = c(m1, m2), xlim = range(year), ylab = "Total population size", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n")
  axis(1,at=seq(min(year),max(year),5),tick=seq(min(year),max(year),5))
  
  polygon(x = c(year,rev(year)), y = c(Nlow,rev(Nhigh)), col = gray(0.6,0.3), border = "gray90")
  #add trend
  polygon(x = c(1:end.yr,end.yr:1), y = c(Nlow[1:end.yr],Nhigh[end.yr:1]), col = gray(0.7,0.3),border = "gray90")
  lines(year[end.yr:nT],Nfit[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
  lines(years,Nfit[1:end.yr], type = "l",col=1,lwd=2)
  
  
  lines(rep(year[max(mp.assess[1],1)],2),c(0,m2),lty=2)
  lines(rep(year[mp.assess[2]],2),c(0,m2),lty=2)
} else {  
  plot(0, 0, ylim = c(m1, m2), xlim = range(year), ylab = "Abudance Index", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n")
  axis(1,at=seq(min(year),max(year),5),tick=seq(min(year),max(year),5))
  
  polygon(x = c(year,rev(year)), y = c(Nlow,rev(Nhigh)), col = gray(0.6,0.3), border = "gray90")
  #add trend
  polygon(x = c(1:end.yr,end.yr:1), y = c(Nlow[1:end.yr],Nhigh[end.yr:1]), col = gray(0.7,0.3),border = "gray90")
  lines(year[end.yr:nT],Nfit[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
  lines(years,Nfit[1:end.yr], type = "l",col=1,lwd=2)
  
  
  lines(rep(year[max(mp.assess[1],1)],2),c(0,m2),lty=2)
  lines(rep(year[mp.assess[2]],2),c(0,m2),lty=2)
}
dev.off()

#### IUCN plot
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(0.7, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(output.dir,"/IUCNplot_",assessment,".png"), width = 4.5, height = 4, 
    res = 200, units = "in")
par(Par)
change = (apply(posteriors$Ntot[,(mp.assess[2]-1):(mp.assess[2]+1)],1,median)/apply(posteriors$Ntot[,(mp.assess[1]-1):(mp.assess[1]+1)],1,median)-1)*100
den = density(change,adjust=2)
x1 = den$x
y1 = den$y

mu.change = round(median(change),1)
sign=""
if(mu.change>0) sign="+"
CR = round(sum(ifelse(change< -80,1,0))/length(change)*100,1)
EN = round(sum(ifelse(change> -80 & change< -50,1,0))/length(change)*100,1)
VU = round(sum(ifelse(change> -50 & change< -30,1,0))/length(change)*100,1)
NT = round(sum(ifelse(change> -30 & change< -20,1,0))/length(change)*100,1)
LC = round(sum(ifelse(change> -20,1,0))/length(change)*100,1)
Decline = round(sum(ifelse(change> 0,1,0))/length(change)*100,1)
plot(x1,y1,type="n",xlim=c(-100,min(max(30,quantile(change,.99)),1000)),ylim=c(0,max(y1*1.05)),ylab="Density",xlab="Change (%)",cex.main=0.9,frame=F,xaxt="n")
axis(1,at=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)),tick=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)))

maxy = max(y1*1.08)
x2 = c(-20,1500); y2 = c(0,5)
polygon(c(x2,rev(x2)),c(rep(maxy,2),rev(rep(0,2))),col="green",border=0)
x3 = c(-30,-20); y2 = c(0,5)
polygon(c(x3,rev(x3)),c(rep(maxy,2),rev(rep(0,2))),col="greenyellow",border=0)
x4 = c(-50,-30)
polygon(c(x4,rev(x4)),c(rep(maxy,2),rev(rep(0,2))),col="yellow",border=0)
x5 = c(-80,-50)
polygon(c(x5,rev(x5)),c(rep(maxy,2),rep(0,2)),col="orange",border=0)
x6 = c(-100,-80)
polygon(c(x6,rev(x6)),c(rep(maxy,2),rep(0,2)),col="red",border=0)

polygon(c(x1,rev(x1)),c(y1,rep(0,length(y1))),col="grey")
#lines(x1,y1,col=1,lwd=1)
abline(v=mu.change,lty=2)
legend("right",c(paste0("CR (",CR,"%)"),paste0("EN (",EN,"%)"),
                 paste0("VU (",VU,"%)"),paste0("NT (",NT,"%)"),
                 paste0("LC (",LC,"%)")),col=1,pt.bg=c("red","orange","yellow","greenyellow","green"),pt.cex=1.4,pch=22,bg="white",cex=0.9)

text(ifelse(median(change)< -80,-80,median(change)),max(y1*1.02),paste0("Change = ",sign,mu.change,"%"),bg="white",cex=1)

dev.off()


Par = list(mfrow=c(1,2),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(0.7, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(output.dir,"/AnnualRates_",assessment,".png"), width = 8, height = 4, 
    res = 200, units = "in")
par(Par)
if(abundance=="census"){
rs =  as.matrix(apply(posteriors$r.tot[,1:(n.years-1)],1,median))
if(n.years>GT1) rs = cbind(rs,apply(posteriors$r.tot[,(n.years-round(GT,0)+1):(n.years-1)],1,median))
if(n.years>2*GT1) rs = cbind(rs,apply(posteriors$r.tot[,(n.years-round(GT*2,0)+1):(n.years-1)],1,median))
if(n.years>3*GT1) rs = cbind(rs,apply(posteriors$r.tot[,(n.years-round(GT*3,0)+1):(n.years-1)],1,median))


} else {
rs =  as.matrix(apply(posteriors$r[,1:(n.years-1)],1,median)) 
if(n.years>GT) rs = cbind(rs,apply(posteriors$r[,(n.years-round(GT,0)+1):(n.years-1)],1,median))
if(n.years>2*GT) rs = cbind(rs,apply(posteriors$r[,(n.years-round(GT*2,0)+1):(n.years-1)],1,median))
if(n.years>3*GT) rs = cbind(rs,apply(posteriors$r[,(n.years-round(GT*3,0)+1):(n.years-1)],1,median))
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
cnam = c("All.yrs","1xGT","2xGT","3xGT")  
 
jcol = c(grey(0.5,0.6),rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
plot(0,0,type="n",ylab="Density",xlab="Annual rate of change (%)",cex.main=0.9,ylim=c(0,1.1*max(lymax)),xlim=quantile(lamdas,c(0.001,0.999))) 
axis(2,labels = FALSE)
for(i in 1:ncol(rs)){
  x = get(paste0("xl",i))
  y = get(paste0("yl",i))
  polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
  mu.lamda = round(median(lamdas[,i]),10)
  lines(rep(mu.lamda,2),c(0,max(y)),col=c(1,4,3,2)[i],lwd=1,lty=1)
  
}
abline(v=0,lty=2)
legend("topright", paste0(cnam[1:ncol(rs)]," = ",ifelse(round(apply(lamdas,2,median),2)>0,"+",""),round(apply(lamdas,2,median),2),"%"),pch=15,col=c(jcol),bty="n")

plot(0,0,type="n",ylab="Density",xlab="log rate r",cex.main=0.9,ylim=c(0,1.1*max(rymax)),xlim=quantile(rs,c(0.001,0.999)))
axis(2,labels = FALSE)
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

#---------------------------------------------------
# Observed Population Trajectory without Projections
#---------------------------------------------------
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(0.7, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(output.dir,"/PopObsYrs_",assessment,".png"), width = 4.5, height = 4, 
    res = 200, units = "in")
par(Par)

# Total N
m1 <- 0
m2 <- max(Nhigh[1:n.years], na.rm = TRUE)
plot(0, 0, ylim = c(m1, m2), xlim = range(years), ylab = "Total population size", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n")
axis(1,at=seq(min(year),max(year),5),tick=seq(min(year),max(year),5))
polygon(x = c(years,rev(years)), y = c(Nlow[1:n.years],rev(Nhigh[1:n.years])), col = gray(0.5,0.5), border = "gray90")
#add trend
lines(years,Nfit[1:end.yr], type = "l",col=1,lwd=2)
dev.off()

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Save key posterior results
out =  data.frame(change) 
out =  data.frame(out,rs)
colnames(out) = c("Pop.change",paste0("r.",cnam[1:ncol(rs)]))
if(abundance=="census") out=data.frame(out,r.mean=posteriors$mean.r) 
save(out,file=paste0(output.dir,"/jara.",assessment,".rdata"))
# Save estimated population trajectory
abundance.out=data.frame(Year=year, Abudance=Nfit,LCI=Nlow,UCI=Nhigh,Type=abundance,Estimation=ifelse(year %in% years,rep("Obs",length(year)) ,rep("Prj",length(year))))
write.csv(abundance.out,paste0(output.dir,"/Pop.trajectory.csv"),row.names = FALSE)
out1 = data.frame(pos=posteriors$Ntot)
save(out1,file=paste0(output.dir,"/jara.",assessment,"Ntot.rdata")) #O<
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
rm(proc.pen)
cat(paste0("\n","><> Complete: ",assessment," - Run ",run," - Type: ",abundance," <><"))
if(n.years<GT1){
cat(paste0("\n","WARNING!!!"))
cat(paste0("\n","The projections are based on less than one Generation Time!"))
cat(paste0("\n","The results are unreliable and MUST be interpreted with caution"))
}
    