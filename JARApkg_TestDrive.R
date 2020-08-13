
#setwd("C:/Work/Research/GitHub/JARA")
#library(devtools)
#load_all()
#document()
#check()
#data("jaradata")


#======================================
# INSTALL JARA
library(devtools)
install_github("Henning-Winker/JARA")
#======================================

# CHECK IT OUT!
library(JARA)
data("jaradata")

# TEST DRIVE

# Striped Marlin (relative)
# create JARA input
mls.input = build_jara(I=jrdat$StripedMarlin_IO_CPUE$I,se=jrdat$StripedMarlin_IO_CPUE$SE,model.type = "relative",assessment = "MLS",GL=5.5)
# Check input
jrplot_indices(mls.input)

mls.input = build_jara(I=jrdat$StripedMarlin_IO_CPUE$I,se=jrdat$StripedMarlin_IO_CPUE$SE,model.type = "relative",assessment = "MLS",GL=5.5)
jarainput = mls.input = build_jara(I=jrdat$StripedMarlin_IO_CPUE$I,se=jrdat$StripedMarlin_IO_CPUE$SE,model.type = "relative",assessment = "MLS",GL=5.5)


# fit JARA model
mls.fit = fit_jara(mls.input,quickmcmc=T) # quick test run
# re-run with 3 chains (default)
mls.fit = fit_jara(mls.input,quickmcmc=F)

jrplot_trjfit(mls.fit)
# Test Colour option
jrplot_trjfit(mls.fit,cols=rainbow(6))
# Fits
jrplot_fits(mls.fit)
jrplot_logfits(mls.fit)

jrplot_poptrj(mls.fit)
jrplot_r(mls.fit)
mls.threat = jrplot_iucn(mls.fit)
mls.threat$perc.risk
mls.threat$status

jrpar(mfrow=c(2,2))
jrplot_trjfit(mls.fit,add=T)
jrplot_poptrj(mls.fit,add=T)
jrplot_changes(mls.fit,add=T)
jrplot_iucn(mls.fit,add=T)


# African Penguin
# create JARA input
ap.input = build_jara(I=jrdat$Afr_penguin$I,assessment="AfrPenguin",model.type = "census",GL=9.9,fixed.obsE = 0.15)
# Check input indices
jrplot_indices(ap.input)
# fit JARA model
ap.fit = fit_jara(ap.input,quickmcmc=T)

jrpar(mfrow=c(2,1))
jrplot_indices(ap.input,add=T)
jrplot_trjfit(ap.fit,add=T)

jrpar(mfrow=c(2,2))
jrplot_trjfit(ap.fit,add=T)
jrplot_poptrj(ap.fit,add=T)
jrplot_changes(ap.fit,add=T)
jrplot_iucn(ap.fit,add=T)

jrplot_fits(ap.fit)
# Choose indices to look at Robben and Dassen
jrplot_trjfit(ap.fit,indices=ap.fit$indices[1:2])
jrplot_fits(ap.fit,indices=ap.fit$indices[1:2],single.plots = F,add=F)
# Manually pick Namibian Islands
jrplot_fits(ap.fit,indices=c("Mercury","Ichaboe","Halifax","Possession"))

# For Rich
jrpar(mfrow=c(4,3),lab=F)
jrplot_fits(ap.fit,indices=ap.fit$indices,single.plots = T,add=T)
# only islands
jrpar(mfrow=c(3,3),labs=F)
jrplot_fits(ap.fit,indices=ap.fit$indices[-c(3,4)],single.plots = T,add=T)
jrplot_logfits(ap.fit,indices=ap.fit$indices[-c(3,4)],single.plots = T,add=T)
ap.threat = jrplot_iucn(ap.fit)
ap.threat$perc.risk
ap.threat$status

# Show both assessments 
jrpar(mfrow=c(2,3))
jrplot_poptrj(ap.fit,add=T)
jrplot_changes(ap.fit,add=T)
legend("topleft","Afr. Penguin",cex=1.,bty="n",y.intersp = -.2,x.intersp = -0.5)
ap.threat = jrplot_iucn(ap.fit,add=T)
jrplot_poptrj(mls.fit,add=T)
jrplot_changes(mls.fit,add=T)
legend("topleft","Striped Marlin",cex=1.,bty="n",y.intersp = -.2,x.intersp = -0.5)
mls.threat = jrplot_iucn(mls.fit,add=T)

