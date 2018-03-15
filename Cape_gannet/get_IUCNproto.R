

##############################################################################
# Bayesian state-space redlisting tool 
# Developed by: Henning Winker
##############################################################################

rm(list=ls())
cat("\014") # clear console# set working directory
gc()
#Just in case
library(gplots)
library(coda)
library(jagsUI) ### THIS IS NEW
library("fitdistrplus")

#library(rjags)
#library(R2jags)

#Just in case
#detach("package:jagsUI", unload=TRUE)
detach("package:R2jags", unload=TRUE)
options(warn=-1) # surpresses plotting warnings

# specify dataset here
SP = 3 # select number of data_set in order 

data_set = c("Afr_penguin","Mountain_zebra","Red_steenbras","Roman")[SP]

# Determine type of abundance data (relative abundance "relative" or population counts "census")
abundance = c("census","census",rep("relative",2))[SP]


# Specify if a csv file with standard errors is provided
SE = c(FALSE,FALSE,rep(TRUE,2))[SP]

# Set up working directory
File = "C:/Work/Research/Redlisttool/get_IUCNstatus"
setwd(File)

# Load dataset
dat = read.csv(paste0("Data/",data_set,".csv"))
datSE = read.csv(paste0("Data/",data_set,ifelse(SE==TRUE,"SE",""),".csv"))


###########################
# Set generation time
###########################

#><> Enter new GT for new data_sets
# Afr Penguin, Mountan Zebra, Red Steenbras 
GT = c(9,12,16,11)[SP]


####################################################
# Set carrying capacities (optional for census data)
####################################################

K = c(TRUE,TRUE,FALSE,FALSE)[SP] # Only relevant if 3GT > number of years
# specify Ks manual here (expert opinion see Mountain Zebra)
K.manual = NULL 
if(SP==2) K.manual = c(170,1200,40,140,1000,100,300,140,150) # Only used for MZ
######################################################
# Prep data and Setup time horizon for trend analysis
######################################################

# For only using data corresponding to 3 x generation  
# Set GT3.dat =TRUE and GT3.dat =FALSE otherwise 

GT1 = round(GT,0) # rounded for r.recent
GT3 = round(3*GT,0) # 3 x GT rounded for year steps

GT3.dat = FALSE
if(GT3.dat==TRUE) dat = dat[(nrow(dat)-GT3):nrow(dat),]

counts = dat[,2:ncol(dat)]
if(abundance=="relative"){# normalize
meancol = apply(counts,2,mean,na.rm=TRUE)
counts = counts / meancol[col(counts)]} 

years = dat[,1]
gr.names=names(counts)
n.years = length(years)
n.groups = ncol(counts)
conv.c = as.numeric(as.matrix(counts))
counts=matrix(conv.c,nrow=n.years,ncol=n.groups)

se = datSE[,2:ncol(dat)]
conv.se = as.numeric(as.matrix(se))
se2 = matrix(ifelse(is.na(conv.se),0.2,conv.se)^2,n.years,n.groups)
if(SE==FALSE) se2 = ifelse(se2>0,0,se2)


if(K==FALSE) Ks = apply(counts,2,max,na.rm=TRUE)*10 
if(K==TRUE & is.null(K.manual)==TRUE) Ks = apply(counts,2,max,na.rm=TRUE)*1.5
if(K==TRUE & is.null(K.manual)==FALSE) Ks = K.manual

#Start Year
st.yr = years[1]
# prediction years 
pyears = ifelse(GT3>(n.years+2),GT3-n.years+2,1) # +2 required for avergaging start + end with y = -1 and +1 

# set up year vectors
year <- years[1]:(years[length(years)] + pyears)
yr = years[1]:(years[length(years)])

nT = length(year)
mp.assess = c(nT-GT3,nT-1)     #start and end midpoints for population trend assessment

# you could use the the assessment period to cut the count data to match 3xGT

# creat prediction matrix for observations
if(pyears>0) counts= rbind(counts,matrix(NA,pyears,n.groups)) # now always
if(pyears>0) se2= rbind(se2,matrix(0.1,pyears,n.groups)) # now always


Ninit = NULL
# get starting values Ninit
for(i in 1:n.groups){
Ninit[i] = na.omit(counts[,i])[1]  
}
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Setup model for population counts
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


if(abundance=="census"){
# Bundle data for jags
jags.data <- list(y = log(counts), T = length(year),nG=n.groups,Ninit=Ninit,Ks=Ks,EY = n.years,penK = matrix(0,nT,n.groups))

# Initial values for jags
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(n.groups), sigma.obs = runif(n.groups, 0, 0.3), logN.est =matrix(log(rbind(Ninit,matrix(rep(NA,(nT-1)*n.groups),(nT-1),n.groups))),nT,n.groups))  }

# Parameters monitored
parameters <- c("mean.r","sigma2.obs", "sigma2.proc","r","logNtot", "N.est","Ntot","lambda")

sink("ssm.jags")
cat("
    model {
    # Priors and constraints
    for(i in 1:nG)
    {
    logN.est[1,i] ~ dnorm(log(Ninit[i]), 1)   # Prior for initial population size with CV =100%
    sigma.obs[i] ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs[i] <- pow(sigma.obs[i], 2)
    tau.obs[i] <- pow(sigma.obs[i], -2)
    mean.r[i] ~ dnorm(0, 0.001)             # Prior for mean growth rate

    }
    
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    for (i in 1:nG){
    r[t,i] ~ dnorm(mean.r[i], tau.proc)
    logN.est[t+1,i] <- ifelse(logN.est[t,i]>log(Ks[i]),log(Ks[i]),logN.est[t,i])+r[t,i] # Ks conditioned predictions
    devK[t,i]  <- ifelse(logN.est[t,i]>log(Ks[i]),logN.est[t,i]-log(Ks[i]),0) # penalty if N > K 
    }}

    
    for(t in 1:EY){
    for (i in 1:nG){
    penK[t,i] ~ dnorm(devK[t,i],1/0.2^2)
    }}

   
    lambda[1] <- 1
   # Observation process
    for (t in 1:T) {
    for(i in 1:nG){
    y[t,i] ~ dnorm(logN.est[t,i], tau.obs[i])
    }}
    
    # Population sizes on real scale
    for (t in 1:T) {
    for(i in 1:nG){
    N.est[t,i] <- exp(logN.est[t,i])
    }
    Ntot[t] <- sum(N.est[t,])
    logNtot[t] <- log(Ntot[t])
    }
    for(t in 2:T)
    {
    lambda[t] <- Ntot[t]/Ntot[t-1]
    }
   } 
   ",fill = TRUE)
sink()

}

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Setup model for relative abundance
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


if(abundance=="relative"){
  # Bundle data for jags
  jags.data <- list(y = log(counts),se2=se2, T = length(year),nG=n.groups,Ninit=Ninit)
  
  # Initial values for jags
  inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.add = runif(n.groups, 0, 0.3), logN.est =matrix(log(rbind(Ninit,matrix(rep(NA,(nT-1)*n.groups),(nT-1),n.groups))),nT,n.groups)[,1])  }
  
  # Parameters monitored
  parameters <- c("mean.r","sigma2.obs", "sigma2.proc","r", "N.est","Ntot","q")
  
  sink("ssm.jags")
  cat("
  model {
       # Priors and constraints
       for(i in 1:nG){
      sigma.add[i] ~ dunif(0, 1)              # Prior for sd of observation process
      sigma2.add[i] <- pow(sigma.add[i], 2)
      }

      for(t in 1:T){
      for(i in 1:nG){
      sigma2.obs[t,i] <- sigma2.add[i]+se2[t,i]
      tau.obs[t,i] <- pow(sigma2.obs[t,i],-1)
      }}
      
      iq[1] <- 1
      q[1] <- 1 
      logq[1] <- log(1)
      for(i in 2:nG){
      iq[i] ~ dgamma(0.001,0.001)
      q[i] <- pow(iq[i],-1)
      logq[i] <-  log(q[i])
      }
      
      
      # Priors and constraints
      logN.est[1] ~ dnorm(log(Ninit[1]), 1)       # Prior for initial population size
      mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
      sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
      sigma2.proc <- pow(sigma.proc, 2)
      tau.proc <- pow(sigma.proc, -2)
      
      # Likelihood
      # State process
      for (t in 1:(T-1)){
      r[t] ~ dnorm(mean.r, tau.proc)
      logN.est[t+1] <- logN.est[t] + r[t] }
      
    
      # Observation process
      for (t in 1:T) {
      for(i in 1:nG){
      y[t,i] ~ dnorm(logN.est[t]+logq[i], tau.obs[t,i])
      }}
      
      # Population sizes on real scale
      for (t in 1:T) {
      N.est[t] <- exp(logN.est[t])
      Ntot[t] <- N.est[t] 
      }
      
      } 
      ",fill = TRUE)
  sink()
  
} 
  
  ni <- 20000
  nt <- 5
  nb <- 5000
  nc <- 3
  
# Call JAGS from R (BRT 3 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,verbose=TRUE,parallel=TRUE)
  
# remove scientific numbers
options(scipen=999)

# Summarize posteriors
print(ssm, digits = 3)

end.yr=n.years

# get individual trends
fitted <- lower <- upper <- mat.or.vec(end.yr,n.groups)

if(abundance=="census"){
for(i in 1:n.groups){
for (t in 1:end.yr){
  fitted[t,i] <- mean(ssm$sims.list$N.est[,t,i])
  lower[t,i] <- quantile(ssm$sims.list$N.est[,t,i], 0.025)
  upper[t,i] <- quantile(ssm$sims.list$N.est[,t,i], 0.975)}}
} else {
fitted <- lower <- upper <- as.null()
  for (t in 1:end.yr){
    fitted[t] <- mean(ssm$sims.list$N.est[,t])
    lower[t] <- quantile(ssm$sims.list$N.est[,t], 0.025)
    upper[t] <- quantile(ssm$sims.list$N.est[,t], 0.975)}
}


Nfit <- Nlow <- Nhigh <- as.numeric()
# get total pop size
for (t in 1:nT){
Nfit[t] =  mean(ssm$sims.list$Ntot[,t])
Nlow[t] = quantile(ssm$sims.list$Ntot[,t],0.025)
Nhigh[t] = quantile(ssm$sims.list$Ntot[,t],0.975)
}

dir.create(paste0(getwd(),"/",data_set),showWarnings = FALSE)

#Abundance FITS
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(1, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(getwd(),"/",data_set,"/Fits_",data_set,".png"), width = 5, height = 4, 
    res = 200, units = "in")
par(Par)

m1 <- 0
m2 <- max(c(fitted, counts, upper), na.rm = TRUE)



if(abundance=="census"){
plot(0, 0, ylim = c(m1, m2), xlim = range(years), ylab = "Population size", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE)
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- rainbow(n.groups)
  
for(i in 1:n.groups) polygon(x = c(years,rev(years)), y = c(lower[,i],upper[end.yr:1,i]), col = gray(runif(1,0.5,0.9),0.5), border = "gray90")
for(i in 1:n.groups)
{
points(year,counts[,i], bg = col_line[i],pch=21) 
lines(years,fitted[,i], type = "l",col=col_line[i], lwd=2)}
} else {  
plot(0, 0, ylim = c(m1, m2), xlim = range(years), ylab = "Abudance Index", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE)
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- rainbow(n.groups)
q = apply(ssm$sims.list$q,2,mean)
polygon(x = c(years,rev(years)), y = c(lower,upper[end.yr:1]), col = gray(runif(1,0.5,0.9),0.5), border = "gray90")
lines(years,fitted, type = "l",col=1, lwd=2)
for(i in 1:n.groups)
{
#if(SE==FALSE) points(year,counts[,i]*q[i], bg = col_line[i],pch=21) 
#if(SE==TRUE) plotCI(year,counts[,i]*q[i],uiw=se[,i], pt.bg = col_line[i],pch=21,add=T,gap=0.)  
points(year,counts[,i]*q[i], bg = col_line[i],pch=21) 
}}
legend('topleft', legend = c("Counts",paste(gr.names[1:n.groups])), lty = c(-1, rep(1,n.groups)), lwd = c(-1, rep(2,n.groups)),pch=c(21,rep(-1,n.groups)), col = c(1,col_line), bty = "n", cex = 1)

dev.off()

Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(1, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(getwd(),"/",data_set,"/PopTrend_",data_set,".png"), width = 5, height = 4, 
    res = 200, units = "in")
par(Par)

# Total N
m1 <- 0
m2 <- max(Nhigh, na.rm = TRUE)

plot(0, 0, ylim = c(m1, m2), xlim = range(year), ylab = "Total population size", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n")
axis(1,at=seq(min(year),max(year),5),tick=seq(min(year),max(year),5))

polygon(x = c(year,rev(year)), y = c(Nlow,rev(Nhigh)), col = gray(0.6,0.3), border = "gray90")
#add trend
polygon(x = c(1:end.yr,end.yr:1), y = c(Nlow[1:end.yr],Nhigh[end.yr:1]), col = gray(0.7,0.3),border = "gray90")
lines(year[end.yr:nT],Nfit[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
lines(years,Nfit[1:end.yr], type = "l",col=1,lwd=2)
dev.off()

#### IUCN plot
Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(1, 0.7, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
png(file = paste0(getwd(),"/",data_set,"/IUCNplot_",data_set,".png"), width = 4, height = 4, 
    res = 200, units = "in")
par(Par)
change = (apply(ssm$sims.list$Ntot[,(mp.assess[2]-1):(mp.assess[2]+1)],1,mean)/apply(ssm$sims.list$Ntot[,(mp.assess[1]-1):(mp.assess[1]+1)],1,mean)-1)*100
den = density(change)
x1 = den$x
y1 = den$y

mu.change = round(mean(change),1)
sign=""
if(mu.change>0) sign="+"
CR = round(sum(ifelse(change< -80,1,0))/length(change)*100,1)
EN = round(sum(ifelse(change> -80 & change< -50,1,0))/length(change)*100,1)
VU = round(sum(ifelse(change> -50 & change< -30,1,0))/length(change)*100,1)
LC = round(sum(ifelse(change> -30,1,0))/length(change),0)
plot(x1,y1,type="n",xlim=c(-100,min(max(30,quantile(change,.99)),1000)),ylim=c(0,max(y1*1.05)),ylab="Density",xlab="Change (%)",cex.main=0.9,frame=F,xaxt="n")
axis(1,at=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)),tick=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)))

maxy = max(y1*1.08)
x2 = c(-30,1500); y2 = c(0,5)
polygon(c(x2,rev(x2)),c(rep(maxy,2),rev(rep(0,2))),col=rgb(0.,1,0.,0.1),border=0)
x2 = c(-50,-30)
polygon(c(x2,rev(x2)),c(rep(maxy,2),rev(rep(0,2))),col="yellow",border=0)
x3 = c(-80,-50)
polygon(c(x3,rev(x3)),c(rep(maxy,2),rep(0,2)),col="orange",border=0)
x4 = c(-100,-80)
polygon(c(x4,rev(x4)),c(rep(maxy,2),rep(0,2)),col="red",border=0)

polygon(c(x1,rev(x1)),c(y1,rep(0,length(y1))),col="grey")

legend("right",c(paste0("CR (",CR,"%)"),paste0("EN (",EN,"%)"),
                 paste0("VU (",VU,"%)")),col=1,pt.bg=c("red","orange","yellow"),pt.cex=1.4,pch=22,bg="white",cex=0.9)

text(ifelse(mean(change)< -80,-80,mean(change)),max(y1*1.02),paste0("Change = ",sign,mu.change,"%"),bg="white",cex=1)
dev.off()

