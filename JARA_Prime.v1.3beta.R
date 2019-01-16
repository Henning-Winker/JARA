

##############################################################################
# JARA: Just Another Redlisting Assessment
# Bayesian state-space redlisting tool 
# Developed by: Henning Winker
##############################################################################

rm(list=ls())
cat("\014") # clear console# set working directory
#gc()

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
version = "v1.3beta"

# Read assessment species information file
sp.assess = read.csv(paste0(JARA.file,"/jara.assessments.csv"))

print(data.frame(Assessment.list=sp.assess[,1]))

# Select assessment species
for(sp in 3:3){
spsel= sp.assess[sp,]

#-----------------------------------------------------------------
# Set up JARA
#-----------------------------------------------------------------
assessment = spsel$assessment # assessment name
run = spsel$run 
abundance = spsel$abundance # c("census","relative")
GT = spsel$generation.time # numeric generation time (years)
Klim = spsel$Klim # c(TRUE,FALSE)
K.manual =  spsel$K.manual # c(TRUE,FALSE), if TRUE provide vector of K for each pop
SE.I = spsel$index.se #  c(TRUE,FALSE), if TRUE provide assessment_se.csv file
sigma.est = spsel$sigma.obs.est # c(TRUE,FALSE)
fixed.obsE = spsel$sigma.obs.add # numeric, fixed component of observation error 
sigma.proc.fixed = spsel$sigma.proc.fixed # numeric or TRUE otherwise
prjr.type = paste(spsel$project.r) # c("mean","all","GT1", years) # "all" is default
proj.stoch = spsel$stochastic.projection # c(FALSE, TRUE), FALSE is default
start.year = spsel$start.year # numeric, NA takes all years
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


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Execute model and produce output
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
if(abundance=="census"){
# MCMC settings
ni <- 1100 # Number of iterations
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




