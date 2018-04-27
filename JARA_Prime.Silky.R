

##############################################################################
# JARA: Just Another Redlisting Assessment
# Bayesian state-space redlisting tool 
# Developed by: Henning Winker
##############################################################################

#rm(list=ls())
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
# Set working directory for JARA R source code
JARA.file = "C:/Work/Research/GitHub/JARA"
# JARA version
version = "v1.2beta"
# select number of data_set in assessments 
SP = 1

# save output
dir.create(paste0(File,"Silky_Global"))
yrs = 1992:2017
silky = silky_se = matrix(NA,length(yrs),3)

for(SP in 1:3){
# Set Assessment file: assement folder within File that includes .csv input files
assessment = c("Silky_Pac","Silky_Atl","Silky_Global")[SP]
# add specifier for assessment (File names of outputs)


#-----------------------------------------------------------------
# Determine whether census data or relative abudance indices
#-----------------------------------------------------------------

# Determine type of abundance data (relative abundance (e.g CPUE) or population counts "census")
abundance = c("relative")

#-------------------------- 
# Set generation time GT
#--------------------------

#><> Enter new GT for new data_sest
GT = c(7.5,5,15)[SP]


#----------------------------------------------------
# Set carrying capacities (optional for census data)
#----------------------------------------------------

K = c(FALSE) # Only relevant if 3GT > number of years
# specify Ks manual here (expert opinNULLion see Mountain Zebra)

K.manual = NULL

#--------------------------------------------------
# Read csv files
#--------------------------------------------------

# Use SEs from csv file for abudance indices (TRUE/FALSE)
SE.I = c(FALSE,FALSE,TRUE)[SP]

# Load dataset
dat = read.csv(paste0(File,assessment,"/",assessment,".csv"))


if(SE.I ==TRUE){
  se = read.csv(paste0(File,"/",assessment,"/",assessment,"_se.csv"))
}



#--------------------------------------------------
# option to exclude abundance series 
#--------------------------------------------------
if(SP==1) dat = dat[,1:3]
if(SP==3){
  dat = dat[,-1]
  se = se[,-1]
}


#--------------------------------------------------
# Variance settings
#--------------------------------------------------

# Observation Error

#To Estimate additional observation variance set sigma.add = TRUE
sigma.est = TRUE

# As option for data-weighing
# minimum fixed observation error for each variance set (optional choose 1 value for both)
fixed.obsE = 0.25 # Important if SE.I is not availble

# Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
# Process Error
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
#Estimate set sigma.proc == True
sigma.proc = TRUE
# Determines if process error deviation are estimated for all years (TRUE)  
# or only from the point the first abundance index becomes available (FALSE)
proc.dev.all = FALSE 
#------------------------------------------
if(sigma.proc == TRUE){
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


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Execute model and produce output
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# MCMC settings
ni <- 8000 # Number of iterations
nt <- 1 # Steps saved
nb <- 500 # Burn-in
nc <- 2 # number of chains
nsaved = (ni-nb)/nt*nc

source(paste0(JARA.file,"/JARA.",version,".R")) 

if(SP<3){
silky[,1] = silky_se[,1]= yrs
silky[yrs %in% year,SP+1] = apply(posteriors$Ntot,2,median)[1:length(year)]  
silky_se[yrs %in% year,SP+1] = apply(log(posteriors$Ntot),2,sd)[1:length(year)]  
}
if(SP==2){
write.csv(silky,paste0(File,"Silky_Global/","Silky_Global.csv"))  
write.csv(silky_se,paste0(File,"Silky_Global/","Silky_Global_se.csv"))  
}

}


