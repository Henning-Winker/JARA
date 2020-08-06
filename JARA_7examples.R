
##############################################################################
# JARA: Just Another Redlisting Assessment
# Bayesian state-space redlisting tool 
# Prime File (Settings)  
# Developed by: Henning Winker* Nathan Pacoureau & Richard Sherley**
# * JRC - European Commission, Ispra, Italy
# Email: henning.winker@gmail.com
# **University of Exeter, UK (richard.sherley@gmail.com)
# Email: richard.sherley@gmail.com
# 
##############################################################################

rm(list=ls())
cat("\014") # clear console# set working directory

library(JARA)
data("jaradata")
dat = jrdat
# Set Working directory file, where assessments are stored 
File = "C:/Work/Research/Github/JARAruns"
# Read assessment species information file
sp.assess = jara.examples
# Check species
print(data.frame(Assessment.list=sp.assess[,1]))

# Select assessment species by row numbers
specs = 1:nrow(sp.assess) # select all species
specs=1

# Loop over selected species
for(sp in specs){
spsel= sp.assess[sp,]

#-----------------------------------------------------------------
# # Extract info for build_jara() for each example
#-----------------------------------------------------------------
ts = dat[[paste0(spsel$assessment)]] # get time series
assessment = spsel$assessment # assessment name
run = spsel$run 
abundance = spsel$abundance # c("census","relative")
GL = spsel$generation.length # numeric generation length (years)
pk.i = NULL
if(spsel[,1] =="Mountain_zebra") pk.i = 1/c(1.55,1.01,0.85,1.75,1.18,1.02,1.27,0.88,0.95)
fixed.obsE = spsel$sigma.obs.add # numeric, fixed component of observation error 
# Organise folders 
output = file.path(File,assessment)
dir.create(output,showWarnings = F)

#--------------------------------------------------
# Build JARA model
#--------------------------------------------------

jara.input = build_jara(I=ts$I,se=ts$SE,model.type = abundance,
                        assessment = assessment,scenario = "noK",
                        GL=GL,fixed.obsE=fixed.obsE,
                        pk.prior = c(0.0001,0.1),
                        pk.i = NULL        # Carrying capacitiy example for Mountain Zebra
                        )
# Check input
jrplot_indices(jara.input,as.png = T,output.dir = output)

# Run JARA with csv file output
fit = fit_jara(jarainput =jara.input,save.csvs = T,output.dir = output)
# plot all routine output
jara_plots(fit,output.dir = output)
}
#-----------------
# The End




