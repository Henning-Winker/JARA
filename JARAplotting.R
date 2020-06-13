#------------------------------------
# Examples for easy plotting in JARA
# Henning Winker (henning.wiker@gmail.com)
# JRC - European Commission
#------------------------------------


library(JARA)
data("jaradata")

# Set working directory
setwd("C:/Work/Research/GitHub/JARA")
# create subfolder
output.dir = "JARAplotting"
dir.create(file.path(getwd(),output.dir),showWarnings = F)

# Example Smoothhound Shark
shsh.input = build_jara(I=dat$SmoothhoundShark$I,se=dat$SmoothhoundShark$SE,model.type = "relative",assessment = "Smoothhound",GL=13.7)
# check trawl survey indices
jrplot_indices(shsh.input)
# Now save as png to output dir
jrplot_indices(shsh.input,as.png=T,output.dir = output.dir)
# fit jata
shsh.fit = fit_jara(shsh.input,quickmcmc=F) # quick test run
# now plot overall fit 
jrplot_trjfit(shsh.fit)
# you can save this also as jpg
dev.print(jpeg,paste0(output.dir,"/shsh.trj.jpg"), width = 5, height = 4.5, res = 300, units = "in")
# or use the inbuilt save option
jrplot_trjfit(shsh.fit,as.png = T,output.dir = output.dir)
# check out the plotting options
?jrplot_trjfit
# another useful function is jrpar() to setup the graphic parameters 
jrpar(mfrow=c(1,1),labs=T,plot.cex = 0.8)
jrplot_trjfit(shsh.fit)
dev.print(jpeg,paste0(output.dir,"/shsh.trj.jpg"), width = 5, height = 4.5, res = 300, units = "in")

# Lets look at Striped Marlin (relative)
# create JARA input
mls.input = build_jara(I=dat$StripedMarlin_IO_CPUE$I,se=dat$StripedMarlin_IO_CPUE$SE,model.type = "relative",assessment = "MLS",GL=5.5)
# Check input
jrplot_indices(mls.input)
# fit JARA model
mls.fit = fit_jara(mls.input,quickmcmc=T) # quick test run
jrplot_trjfit(mls.fit)
# Test Colour option
jrplot_trjfit(mls.fit,cols=rainbow(6))

# now if we want to show two examples in the same plot we use the option add=T
jrpar(mfrow = c(1,2),plot.cex = 0.8)
jrplot_indices(shsh.input,add=T)
mtext(c("a)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
jrplot_indices(mls.input,add=T)
mtext(c("b)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
# save to file
dev.print(png,paste0(output.dir,"/Indices2x1.png"), width = 8, height = 4, res = 300, units = "in")
# note the legend may distorted if using R studio to avoid this use of 
#e.g. windows(width=8,height=4) or quartz()

# We can do the same with fits
jrpar(mfrow = c(1,2),plot.cex = 0.8)
jrplot_trjfit(shsh.fit,add=T)
mtext(c("a)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
jrplot_trjfit(mls.fit,add=T)
mtext(c("b)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
dev.print(png,paste0(output.dir,"/trj2x1.png"), width = 8, height = 4, res = 300, units = "in")

# Save Fits for Readme file illustration
jrplot_fits(mls.fit,as.png = T,output.dir = output.dir)
jrplot_logfits(mls.fit,as.png = T,output.dir = output.dir)

# Now take a look at the census data for Penguin
# African Penguin
# create JARA input
ap.input = build_jara(I=dat$Afr_penguin$I,assessment="AfrPenguin",model.type = "census",GL=9.9,fixed.obsE = 0.15)
# Check input indices
jrplot_indices(ap.input)
# fit JARA model
ap.fit = fit_jara(ap.input,quickmcmc=F)

# and check the fits and population trajectory
jrpar(mfrow = c(1,2),plot.cex = 0.7)
jrplot_trjfit(ap.fit,add=T)
mtext(c("a)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
jrplot_poptrj(ap.fit,add=T)
mtext(c("b)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
dev.print(png,paste0(output.dir,"/AP2x1.png"), width = 8, height = 4, res = 300, units = "in")

# we can quickly compare to a relative abudance example
jrpar(mfrow = c(1,2),plot.cex = 0.7)
jrplot_trjfit(shsh.fit,add=T)
mtext(c("a)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
jrplot_poptrj(shsh.fit,add=T)
mtext(c("b)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
dev.print(png,paste0(output.dir,"/SHSH_trend.png"), width = 8, height = 4, res = 300, units = "in")

# Probably the most important JARA result is
# the change over time 
jrpar(mfrow = c(1,2),plot.cex = 0.7)
jrplot_changes(ap.fit,add=T)
mtext(c("a)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
jrplot_iucn(ap.fit,add=T)
mtext(c("b)"),side=3,outer=F,cex=1.1,line=0.5,adj=-0.1)
dev.print(png,paste0(output.dir,"/APstatus2x1.png"), width = 8, height = 4, res = 300, units = "in")

#-----------------------------------------------------
# An alternative option is to read in the JARA plots from a completed assessment folder
# Examples of joining JARA outputs into multiplots
#-----------------------------------------------------

# All plot will be save in the respective species assessment folders
# For example run Yellowspotted Skate 
# create a folder
assessment = "YellowspottedSkate"
newoutput.dir = assessment 
dir.create("YellowspottedSkate",showWarnings = F)
ysk=build_jara(I=dat$Yellowspot_skate$I,se=dat$Yellowspot_skate$SE,model.type = "relative",assessment = assessment,GL=12)
# fit
fit.ysk = fit_jara(ysk)

# Plot all output in newoutput.dir
jara_plots(fit.ysk,as.png = T,output.dir = newoutput.dir)

#----------------------------------------------------------------        
# Produce IUCN Supplement Information "style" Plot
#----------------------------------------------------------------
DIMs=c(5,4) # Plot dimension
sp=8  # Select species  
Plot = c("IUCN_SoupfinShark")
DIMs=c(2*5,2*4)/1.5
# Choos plot types
plots = c("TrjFits_","PopTrend_","AnnualRate_","IUCNplot_")
j=1
runs  = 1 # run name
# layout the plots into a matrix  columns, by row
library(png)
par(mar=rep(0,4),omi= c(0, 0, 0, 0),cex=1.3) # no margins
layout(matrix(1:4, ncol=2, byrow=TRUE))

for(j in 1:length(plots)){
  run = runs[1]

  # plot path
  get_plot = file.path(getwd(),assessment,paste0(plots[j],assessment,"_s1.png"))
  # load image
  img <- readPNG(paste0(get_plot))
  
  # do the plotting
  plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
  rasterImage(img,0,0,1,1)
  legend("topleft",paste0("(",letters[j],")"),bty="n",cex=1.2,x.intersp = -0.55,y.intersp=-0.2)
}

# Write plot
dev.print(png,paste0(output.dir,"/",assessment,"_status.png"), width = DIMs[1], height = DIMs[2], res = 300, units = "in")

