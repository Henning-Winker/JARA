

##############################################################################
# JARA: Just Another Redlisting Assessment
# Bayesian state-space redlisting tool 
# Retrospective Analysis Setup
# Developed by: Henning Winker* & Richard Sherley**
# * Department of Agriculture, Forestry and Fisheries (DAFF), ZAF \
# Email: henning.winker@gmail.com
# **University of Exeter, UK (richard.sherley@gmail.com)
# Email: richard.sherley@gmail.com
##############################################################################

rm(list=ls())
cat("\014") # clear console# set working directory
gc()

# required packages
library(gplots)
library(coda)
library(rjags)
library(R2jags)
library("fitdistrplus")


#----------------------------------------------------------------
# Setup working directories and output folder labels 
#-----------------------------------------------------------------
# Set Working directory file, where assessments are stored 
File = "C:/Work/Research/Github/JARA/"
# Set working directory for JABBA R source code
JARA.file = "C:/Work/Research/GitHub/JARA"
# JABBA version
version = "v1.1"

# Read assessment species information file
sp.assess = read.csv(paste0(File,"/jara.assessments.csv"))

print(data.frame(Assessment.list=sp.assess[,1]))

# Select assessment species
sp = 1  # Here African Penguin
spsel= sp.assess[sp,]

# RETROSPECTIVES LOOP
nback = 0:11 # number of years back
for(step in 1:length(nback)){

#-----------------------------------------------------------------
# Set up JARA
#-----------------------------------------------------------------
  assessment = spsel$assessment # assessment name
  run = spsel$run 
  abundance = spsel$abundance # c("census","relative")
  GL = spsel$generation.length # numeric generation length (years)
  Klim = spsel$Klim # c(TRUE,FALSE)
  K.manual =  spsel$K.manual # c(TRUE,FALSE), if TRUE provide vector of K for each pop
  if(spsel[,1] =="Mountain_zebra") K.manual = c(1.55,1.01,0.85,1.75,1.18,1.02,1.27,0.88,0.95)
  SE.I = spsel$index.se #  c(TRUE,FALSE), if TRUE provide assessment_se.csv file
  sigma.est = spsel$sigma.obs.est # c(TRUE,FALSE)
  fixed.obsE = spsel$sigma.obs.add # numeric, fixed component of observation error 
  sigma.proc.fixed = spsel$sigma.proc.fixed # numeric or TRUE otherwise
  prjr.type = (spsel$project.r) # c("all","GL1", years) # "all" is default
  proj.stoch = spsel$stochastic.projection # c(FALSE, TRUE), FALSE is default
  start.year = spsel$start.year # numeric, NA takes all years
  end.year = spsel$end.year # numeric, NA takes all years
  A1 = spsel$A1
  plot.width = 5 # default is 5 inches
  #--------------------------------------------------
# Read csv files
#--------------------------------------------------
# Load dataset
dat = read.csv(paste0(File,"/",assessment,"/",assessment,".csv"))#[,-3]

if(SE.I ==TRUE){
  se = read.csv(paste0(File,"/",assessment,"/",assessment,"_se.csv"))
}

#---------------------------------------------------------------------
# option to manually exclude abundance series or change starting year
#---------------------------------------------------------------------
# Do Retrospective removal of indices
end.year = max(dat[,1])
if(nback[step]>0) end.year = end.year-nback[step] 
run = end.year
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Execute model and produce output
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
if(abundance=="census"){
# MCMC settings
ni <- 6000 # Number of iterations
nt <- 1 # Steps saved
nb <- 1000 # Burn-in
} else {
  ni <- 12000 # Number of iterations
  nt <- 2 # Steps saved
  nb <- 1000 # Burn-in
}

nc <- 2 # number of chains
nsaved = (ni-nb)/nt*nc

source(paste0(JARA.file,"/JARA.",version,".R")) 
}

#-----------------------------------------------------
# Compile Retrospective runs output
#-----------------------------------------------------
# read original data
dat = read.csv(paste0(File,"/",assessment,"/",assessment,".csv"))#[,-3]
jara.retro = list()
runs = NULL
for(step in 1:(length(nback))){
end.year = max(dat[,1])
if(nback[step]>0) end.year = end.year-nback[step] 
run = end.year
runs = c(runs,run)
get.jara = paste0(File,"/",assessment,"/output",run,"/jara.",assessment,".rdata")
load(get.jara,verbose=TRUE)
if(step==1){
jara.retro$change = data.frame(run=jara$trends$pop.change) 
jara.retro$abundance = data.frame(run,jara$abundance) 
jara.retro$perc.risk = data.frame(run,jara$perc.risk) 
jara.retro$status = data.frame(run,status=jara$status) 
}else{
jara.retro$change = cbind(jara.retro$change,data.frame(run=jara$trends$pop.change)) 
jara.retro$abundance = rbind(jara.retro$abundance, data.frame(run,jara$abundance)) 
jara.retro$perc.risk = rbind(jara.retro$perc.risk,data.frame(run,jara$perc.risk)) 
jara.retro$status = rbind(jara.retro$status,data.frame(run,status=jara$status)) 
}}
colnames(jara.retro$change) = runs

save(jara.retro,file =paste0(File,"/",assessment,"/jara.retro.",assessment,".rdata"))

#-------------------------------------------------------------------
# Produce Retro-IUCN Plot
#------------------------------------------------------------------
# Load rdata (to avoid rerunning)

plot.width = 5
load(paste0(File,"/",assessment,"/jara.retro.",assessment,".rdata"),verbose=T)
par.save = par
# Extract residual for by scenario and index 
d = jara.retro
names(d)

Par = list(mfrow=c(1,1),mar = c(3, 3, 1.2, 0.1),oma=c(0, 0, 0, 4), mgp =c(2.,0.5,0), tck = -0.02,cex=1)
png(file = paste0(File,"/",assessment,"/RetroPosteriors.",assessment,".png"), width = plot.width, height = 4, 
    res = 200, units = "in")
par(Par)
ylim=c(-100,min(max(30,quantile(unlist(d$change),0.99)),1000))
xlim=c(0.5,length(runs)+0.5)
  plot(0,0,type="n",xlab="Assessment year",xlim=xlim,ylim=ylim,axes=F,xaxs = "i",yaxs="i",ylab="Change(%)")  
  axis(1,at=seq(0,length(runs)+1,1),labels=rev(c(min(runs-1),runs,max(runs+1))),cex.axis=0.8,mgp=c(2,0.5,0))
  xall = unlist(d$change)
  axis(2,at=seq(-100,max(xall,30)+50,ifelse(max(xall,30)>150,50,25)),tick=seq(-100,max(x1,30)+50,ifelse(max(xall,30)>150,50,25)),cex.axis=0.8,mgp=c(2,0.5,0))
  for(j in 1:length(runs)){
  change = d$change[,rev(runs) ==runs[j]]
  den = density(change,adjust=2)
  y1 = den$x
  x1 = den$y/max(den$y)+j-0.5
  
  lc = c(ifelse(A1,-50,-30))
  polygon(c(x1[y1>=lc],rep(min(x1),length(y1[y1>=lc]))),c(y1[y1>=lc],rev(y1[y1>=lc])),col="green",border="green")
  vu = c(ifelse(A1,-70,-50))
  polygon(c(x1[y1<lc & y1>=vu],rep(min(x1),length(y1[y1<lc & y1>=vu]))),c(y1[y1<lc & y1>=vu],rev(y1[y1<lc & y1>=vu])),col="yellow",border="yellow")
  en =ifelse(A1,-90,-80)
  polygon(c(x1[y1<vu & y1>=en],rep(min(x1),length(y1[y1<vu & y1>=en]))),c(y1[y1<vu & y1>=en],rev(y1[y1<vu & y1>=en])),col="orange",border="orange")
  polygon(c(x1[y1<en],rep(min(x1),length(y1[y1<en]))),c(y1[y1<en],rev(y1[y1<en])),col="red",border="red")
  polygon(c(x1,rep(min(x1),length(x1))),c(y1,rev(y1)))
  }
  text(1:length(runs),max(ylim)*0.95,rev(d$status[,2]),cex=0.8)
  legend(par('usr')[2]*1.01, quantile(par('usr')[3:4],0.6), bty='n', xpd=NA,
         c("LC","VU","EN","CR"),pch=15,col=c("green","yellow","orange","red"),pt.cex=2,cex=0.9)
  
  
  dev.off()  
  par=par.save
  
  
  cols=rainbow(length(runs))
  # Retro Trends 
  Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=1)
  png(file = paste0(File,"/",assessment,"/RetroTrends.",assessment,".png"), width = 7, height = 4, 
      res = 200, units = "in")
  par(Par)
  trend = d$abundance[runs[1]==d$abundance$run,]
  yr = trend$yr
  ymax = max(d$abundance$total,trend$total.ucl)
  
  plot(0, 0, ylim = c(0, ymax), xlim = c(min(yr-1),max(yr+1)), ylab = "Total population size", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n",xaxs="i",yaxs="i")
  axis(1,at=seq(min(yr),max(yr)+5,5),tick=seq(min(yr),max(yr),5))
  polygon(c(yr,rev(yr)),c(trend$total.lcl,rev(trend$total.ucl)),col="grey",border="grey")
  for(j in 2:length(runs)){
  ts = d$abundance[runs[j]==d$abundance$run,]
  lines(ts$yr,ts$total,col=cols[j-1],lwd=2,lty=1)
  }
  lines(yr,trend$total,lwd=2)
  legend(ifelse(trend$total[1]<trend$total[length(yr)],"topleft","topright"),paste(runs),col=c(1,cols),bty="n",cex=0.7,lwd=c(2))
  dev.off()
  
  