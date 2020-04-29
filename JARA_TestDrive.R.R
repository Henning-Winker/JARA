

library(devtools)
load_all()
document()
check()
data("jaradata")


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
mls.input = build_jara(I=dat$StripedMarlin_IO_CPUE$I,dat$StripedMarlin_IO_CPUE$se,model.type = "relative",assessment = "MLS",GL=5.5)
# Check input
jrplot(mls.input)

# fit JARA model
mls.fit = fit_jara(mls.input,quickmcmc=T) # quick test run
# re-run with 3 chains (default)
mls.fit = fit_jara(mls.input,quickmcmc=F)

jrplot_trjfit(mls.fit)
jrplot_fits(mls.fit)
jrplot_logfits(mls.fit)

jrplot_poptrj(mls.fit)
mls.threat = jrplot_iucn(mls.fit)
mls.threat$perc.risk
mls.threat$status

# African Penguin
# create JARA input
ap.input = build_jara(I=dat$Afr_penguin$I,assessment="AfrPenguin",model.type = "census",GL=9.9)
# fit JARA model
ap.fit = fit_jara(ap.input,quickmcmc=T)

jrplot_trjfit(ap.fit)
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

# 
jrpar(mfrow=c(2,2))
jrplot_poptrj(ap.fit,add=T)
ap.threat = jrplot_iucn(ap.fit,add=T)
legend("topright","Afr. Penguin",cex=1.2,bty="n",y.intersp = -.2)
jrplot_poptrj(mls.fit,add=T)
ap.threat = jrplot_iucn(mls.fit,add=T)
legend("topright","Striped Marlin",cex=1.2,bty="n",y.intersp = -.2)

