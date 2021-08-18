

library(JARA)
setwd("C:/Work/Research/GitHub/JARA_timeblock")
load(file="Sphyrna_tiburo.rdata",verbose=T)

# build JARA normally with large GL to see effect of projections
jb = build_jara(I,se,fixed.obsE = 0.01,GL=GL+5)
# build JARA with time block
tb = build_jara(I,se,fixed.obsE = 0.01,GL=GL+5,timeblock = 1994)
# Fit both
fjb = fit_jara(jb)
ftb = fit_jara(tb)
# The fits are only slightly effected by accounting for addition variation
jrplot_poptrj(fjb)
jrplot_poptrj(ftb)
# However, the post time-block r is now used for projections
jrplot_poptrj(fjb)
jrplot_poptrj(ftb)
# What you after is the effect size
jrplot_timeblock(fjb) # no time-block - no plot
jrplot_timeblock(ftb) # no time-block - no plot



fitntb = fit_jara(notb,quickmcmc = T)
fittb = fit_jara(tb,quickmcmc = T)
jara = fittb
jrplot_poptrj(fitntb) 
jrplot_poptrj(fittb) 

