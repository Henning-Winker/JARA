
#------------------------------
# Plot process error deviation
#------------------------------

proc.dev = apply(posteriors$r.tot[,-1],2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/ProcDev_",assessment,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = c(min(-0.22,proc.dev),max(0.22,proc.dev))#range(proc.dev)*1.1
cord.x <- c(years,rev(years))
cord.y <- c(proc.dev[1,],rev(proc.dev[3,]))
# Process Error
plot(years,proc.dev[2,],ylab="Process Error Deviates r[t]",xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,proc.dev[2,],lwd=2)
lines(years,rep(0,length(years)),lty=5)
dev.off()