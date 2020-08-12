#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# library(JARA)
# Code to produce figures for JARA Paper 
#
# github.com/henning-winker/JARA 
#
# Henning Winker (henning.winker@gmail.com) 
# Joint Research Centre of the European Commission
# Work Email: Henning.Winker@ec.europa.eu 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

setwd("C:/Work/Research/MS_JARA2020/MSplots")
# Install JARA
#devtools::install_github("henning-winker/jara")

# load pakage
library(JARA)
# load data
# jara.examples # inputs
# jrdat # Indices

specs = data.frame(Species=jara.examples$assessment)
specs
# Run Yellowspotted skate 
i = 8
sp = jara.examples$assessment[i]
pars =  jara.examples[i,]
dat = jrdat$Yellowspot_skate
# Build JARA input 
jr1 = build_jara(I=dat$I,se=dat$SE,
                 assessment=sp,
                 scenario = "MS",
                 GL=pars$generation.length)
# Fit JARA
fit1 = fit_jara(jr1,save.jara = T)   

# Run Whitespotted Smoothhound 
i = 7
sp = jara.examples$assessment[i]
pars =  jara.examples[i,]
dat = jrdat$SmoothhoundShark
# Build JARA input 
jr2 = build_jara(I=dat$I,se=dat$SE,
                 assessment=sp,
                 scenario = "MS",
                 GL=pars$generation.length)
# Fit JARA
fit2 = fit_jara(jr2,save.jara = T)   

#*************************
# Figure 1: Indices
#*************************
l = 1.2
plname = "Fig1"
pwidth = 8
pheight = 7
res=300
adj=-.15
windows(width=pwidth,height=pheight) # Sorry Rich
jrpar(mfrow=c(2,2),labs=T,plot.cex=1)
jrplot_indices(jr1,add=T)
mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_indices(jr2,add=T)
mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_trjfit(fit1,add=T)
mtext(text=paste0(letters[3],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_trjfit(fit2,add=T)
mtext(text=paste0(letters[4],")"), xpd=NA, side=3, adj=adj, cex=1)
mtext(c("Yellowspotted Skate","Whitespotted Houndshark"), side=3, outer=T,line= -0.2,cex=1.,c(0.27,0.77))

dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")

#*************************
# Figure 2: Decision plot
#*************************
plname = "Fig2"
pwidth = 8
pheight = 8.5
res=300
windows(width=pwidth,height=pheight) # Sorry Rich
jrpar(mfrow=c(3,2),labs=T,plot.cex=1)
jrplot_poptrj(fit1,add=T)
mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_poptrj(fit2,add=T)
mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_changes(fit1,add=T)
mtext(text=paste0(letters[3],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_changes(fit2,add=T)
mtext(text=paste0(letters[4],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_iucn(fit1,add=T)
mtext(text=paste0(letters[5],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_iucn(fit2,add=T)
mtext(text=paste0(letters[6],")"), xpd=NA, side=3, adj=adj, cex=1)
mtext(c("Yellowspotted Skate","Whitespotted Houndshark"), side=3, outer=T,line= -0.2,cex=1.,c(0.27,0.77))

dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")

#--------------------------
# Do hindcasting
#-------------------------
runhc = FALSE # need to run one time
if(runhc==T){
hc1 = jara_hindcast(jr1,speedup = F,save.jarafile = F)
hc2 = jara_hindcast(jr2,speedup = F,save.jarafile = F)
save(hc1,file="hc1.rdata")
save(hc2,file="hc2.rdata")} else{
  load("hc1.rdata",verbose=T)
  load("hc2.rdata",verbose=T)
}


#*************************
# Figure 3: Retrospectives
#*************************
plname = "Fig3"
pwidth = 8
pheight = 8.5
res=300
adj=-.15
windows(width=pwidth,height=pheight) # Sorry Rich
jrpar(mfrow=c(3,2),labs=T,plot.cex=1)
jrplot_state(fit1,add=T)
mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_state(fit2,add=T)
mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_retrobias(hc1,add=T)
mtext(text=paste0(letters[3],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_retrobias(hc2,add=T,legend.loc = "topleft")
mtext(text=paste0(letters[4],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_retroiucn(hc1,add=T,add.legend = F)
mtext(text=paste0(letters[5],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_retroiucn(hc2,add=T,add.legend = F)
mtext(text=paste0(letters[6],")"), xpd=NA, side=3, adj=adj, cex=1)
mtext(c("Yellowspotted Skate","Whitespotted Houndshark"), side=3, outer=T,line= -0.2,cex=1.,c(0.27,0.77))

dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")


#******************************************
# Figure 4: Blue Marlin CPUE vs Assessment
#******************************************
specs
# Run CPUEs 
i = 5
sp = jara.examples$assessment[i]
pars =  jara.examples[i,]
dat = jrdat$Bluemarlin_Atl_CPUE
# Build JARA input 
jr3 = build_jara(I=dat$I,se=dat$SE,
                 assessment=sp,
                 scenario = "MS",
                 GL=pars$generation.length)
# Fit JARA
fit3 = fit_jara(jr3,save.jara = T)   

# Run Assessment output 
i = 6
sp = jara.examples$assessment[i]
pars =  jara.examples[i,]
dat = jrdat$Bluemarlin_ICCAT
# Build JARA input 
jr4 = build_jara(I=dat$I,se=dat$SE,
                 assessment=sp,
                 scenario = "MS",
                 GL=pars$generation.length)
# Fit JARA
fit4= fit_jara(jr4,save.jara = T)   

# Make Figure
plname = "Fig4"
pwidth = 8
pheight = 8.5
res=300
windows(width=pwidth,height=pheight) # Sorry Rich
jrpar(mfrow=c(3,2),labs=T,plot.cex=1)
jrplot_trjfit(fit3,add=T)
mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_trjfit(fit4,add=T)
mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_poptrj(fit3,add=T)
mtext(text=paste0(letters[3],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_poptrj(fit4,add=T)
mtext(text=paste0(letters[4],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_iucn(fit3,add=T)
mtext(text=paste0(letters[5],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_iucn(fit4,add=T)
mtext(text=paste0(letters[6],")"), xpd=NA, side=3, adj=adj, cex=1)
mtext(c("CPUE Indices","Stock Assessment"), side=3, outer=T,line= -0.2,cex=1.,c(0.27,0.77))

dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")

#******************************************
# Figure 5: Cape Gannet
#******************************************
specs
# Run CPUEs 
i = 2
sp = jara.examples$assessment[i]
pars =  jara.examples[i,]
dat = jrdat$Cape_gannet
# Build JARA input 
jr5 = build_jara(I=dat$I,se=dat$SE,
                 assessment=sp,
                 model.type = "census",
                 scenario = "MS",
                 GL=pars$generation.length)
# Fit JARA
fit5 = fit_jara(jr5,save.jara = T)   

# do user projections
jr5p = build_jara(I=dat$I,se=dat$SE,
                 assessment=sp,
                 model.type = "census",
                 scenario = "MS",proj.yrs.user = 20,
                 pk.prior = c(0.8,0.1))
# Assumes that those colonies that currently increasing have 
# a maximum growth potential + 20% relative their historical maximum
# i.e. Algoa Bay and Lamberts Bay and Malgas
fit5$pars
# Interesting story line
# Fit JARA
fit5p = fit_jara(jr5p,save.jara = T)   
# check
jrplot_poptrj(fit5)
jrplot_poptrj(fit5p)

# do retrospectives (10 years back)
runhc = FALSE # need to run one time
if(runhc==T){
  hc5 = jara_hindcast(jr5,peels = 0:10,speedup = F,save.jarafile = F)
  save(hc5,file="hc5.rdata")
  } else{
    load("hc5.rdata",verbose=T)
  }

# Make Figure
plname = "Fig5"
pwidth = 8
pheight = 8.5
res=300
windows(width=pwidth,height=pheight) # Sorry Rich
jrpar(mfrow=c(3,2),labs=T,plot.cex=1)
jrplot_trjfit(fit5,add=T)
mtext(text=paste0(letters[1],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_changes(fit5,add=T)
mtext(text=paste0(letters[2],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_poptrj(fit5,add=T)
mtext(text=paste0(letters[3],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_iucn(fit5,add=T)
mtext(text=paste0(letters[4],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_poptrj(fit5p,add=T)
mtext(text=paste0(letters[5],")"), xpd=NA, side=3, adj=adj, cex=1)
jrplot_retroiucn(hc5,add=T,add.legend = F)
mtext(text=paste0(letters[6],")"), xpd=NA, side=3, adj=adj, cex=1)
mtext(c("CPUE Indices","Stock Assessment"), side=3, outer=T,line= -0.2,cex=1.,c(0.27,0.77))

dev.print(tiff,paste0(getwd(),"/",plname,"_hires.tiff"), width = pwidth, height = pheight, res = res, units = "in")
dev.print(jpeg,paste0(getwd(),"/",plname,".jpg"), width = pwidth, height = pheight, res = res, units = "in")

# End (for now)
