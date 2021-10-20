
# reinstall to update
#devtools::install_github("henning-winker/JARA")

library(JARA)
setwd("C:/Work/Research/GitHub/JARA_timeblock")
load(file="Sphyrna_tiburo.rdata",verbose=T)

# build JARA normally with large GL to see effect of projections
jb = build_jara(index,se,fixed.obsE = 0.01,GL=GL+5)
# build JARA with time block
tb = build_jara(index,se,fixed.obsE = 0.01,GL=GL+5,timeblock = 1994)
# FIT both
fjb = fit_jara(jb)
ftb = fit_jara(tb)
# The fits are only slightly effected by accounting for addition variation
# However, the post time-block r is now used for projections
jrpar(mfrow=c(2,2))
jrplot_trjfit(fjb,add=T)
jrplot_trjfit(ftb,add=T)
jrplot_poptrj(fjb,add=T)
jrplot_poptrj(ftb,add=T)
mtext(c("No Time-Block","Time-Block 1994"),3,at=c(0.27,0.78),outer=T,cex=0.9)
# check how proj r differs now
fjb$r.prj
ftb$r.prj # mean.r estimated for 1994+

# What you after is the effect size
jrplot_timeblock(fjb) # no time-block - no plot
jrplot_timeblock(ftb) 

# the posterior is here
ftb$timeblock

# You can also check r
jrplot_timeblock(ftb,type="r") # no time-block - no plot

# get r post impact
ftb$pars[1,1]
ftb$r.prj[1]




