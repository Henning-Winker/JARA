

#' jrplot_pars()
#'
#' Set the par() to options suitable for JARA multi plots   
#' @param mfrow determines plot frame set up
#' @param plot.cex cex graphic option
#' @export
jrpar <- function(mfrow=c(1,1),plot.cex=1,mai=c(0.35,0.15,0,.15),labs=TRUE){
if(labs)  mai=c(0.45,0.45,0.15,.15)
  par(list(mfrow=mfrow,mai = mai, mgp =c(2.,0.5,0),omi = c(0.2,0.25,0.2,0) + 0.1, tck = -0.02,cex=0.8))
}


#' jrplot_iucn
#'
#' IUCN posterior plot for A1 or A2   
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param criteria A1 or A2 for decline
#' @param add if TRUE par is not called to enable manual multiplots
#' @param plot.cex cex graphic option
#' @param iucn.cols to use iucn color recommendation or a brighter version if FALSE
#' @param criteria option to choose between IUCN A1 or A2 thresholds (A2 is default)  
#' @param Plot if FALSE then only threat status stats are returned 
#' @return IUCN classification 
#' @export
jrplot_iucn <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,criteria=c("A2","A1")[1],iucn.cols=TRUE,add=FALSE,Plot=TRUE){
  
  cat(paste0("\n","><> jrplot_iucn() - return % threat classification <><","\n"))
  change= jara$posteriors$pop.change
  den = stats::density(change,adjust=2)
  x1 = den$x
  y1 = den$y
  A1 = ifelse(criteria=="A1",TRUE,FALSE)
  mu.change = round(median(change),1)
  sign=""
  if(mu.change>0) sign="+"
  CR = round(sum(ifelse(change< ifelse(A1,-90,-80),1,0))/length(change)*100,1)
  EN = round(sum(ifelse(change> ifelse(A1,-90,-80) & change< ifelse(A1,-70,-50),1,0))/length(change)*100,1)
  VU = round(sum(ifelse(change> ifelse(A1,-70,-50) & change< ifelse(A1,-50,-30),1,0))/length(change)*100,1)
  LC = round(sum(ifelse(change> -30,1,0))/length(change)*100,1)
  Decline = round(sum(ifelse(change< 0,1,0))/length(change)*100,1)
  if(iucn.cols==T){
    cols = c("#60C659","lightgreen","#F9E814","#FC7F3F","#D81E05")[c(1,3:5)] # green to red
  } else {
    cols=c("green","lightgreen","yellow","orange","red")[c(1,3:5)]  
  }
  
  if(Plot==TRUE){
    Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
    if(as.png==TRUE){png(file = paste0(output.dir,"/IUCNplot_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    if(add==FALSE) par(Par)
    
    plot(x1,y1,type="n",xlim=c(-100,min(max(30,quantile(change,.99)),1000)),ylim=c(0,max(y1*1.1)),ylab="",xlab="",cex.main=0.9,frame=T,xaxt="n",yaxt="n",xaxs="i",yaxs="i")
    maxy = max(y1*1.11)
    x2 = c(ifelse(A1,-50,-30),1500); y2 = c(0,5)
    polygon(c(x2,rev(x2)),c(rep(maxy,2),rev(rep(0,2))),col=cols[1],border=cols[1])
    x3 = c(ifelse(A1,-70,-50),x2[1])
    polygon(c(x3,rev(x3)),c(rep(maxy,2),rev(rep(0,2))),col=cols[2],border=cols[2])
    x4 = c(ifelse(A1,-90,-80),x3[1])
    polygon(c(x4,rev(x4)),c(rep(maxy,2),rep(0,2)),col=cols[3],border=cols[3])
    x5 = c(-100,x4[1])
    polygon(c(x5,rev(x5)),c(rep(maxy,2),rep(0,2)),col=cols[4],border=cols[4])
    
    polygon(c(x1,rev(x1)),c(y1,rep(0,length(y1))),col="grey")
    axis(1,at=seq(-100,max(x1,30)+50,ifelse(max(x1,30)>150,50,25)),tick=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)))
    mtext(paste("Density"), side=2, outer=F,line=1.9,cex=1)
    mtext(paste("Change (%)"), side=1, outer=F,line=1.9,cex=1)
    legend("right",c(paste0("CR (",CR,"%)"),paste0("EN (",EN,"%)"),
                     paste0("VU (",VU,"%)"),paste0("LC (",LC,"%)")),col=1,pt.bg=c("red","orange","yellow","green"),pt.cex=1.4,pch=22,bg="white",cex=1.1)
    text(ifelse(mean(change)< -80,-80,mean(change)),max(y1*1.03),paste0("Change = ",sign,mu.change,"%"),bg="white",cex=1.2)
    
    if(as.png==TRUE) dev.off()
  } # End IUCN Plot = TRUE
  categories = c("Pr.Decl","change3GL","CR","EN","VU","LC")  
  percentages = c(CR,EN,VU,LC)
  status= ifelse(which(percentages==max(percentages))==4 & max(percentages)<50,"NT",categories[3:6][which(percentages==max(percentages))])
  perc.risk = data.frame(perc.risk=c(Decline,mu.change,percentages))
  rownames(perc.risk)=categories
  out = list()
  out$perc.risk = perc.risk
  out$status = status
  
  return(out)
}


#' jrplot_trjfit()
#'
#' Plots the estimated and predicted population trajectors on the same scale    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @param indices names of indices to plot (default = "all")
#' @export
jrplot_trjfit <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE,indices="all"){
  
  cat(paste0("\n","><> jrplot_trjfit() - fit to trajectories <><","\n"))
  
  
  Nt = jara$trj[jara$trj$name=="global" & jara$trj$estimation=="fit",]
  year = Nt$yr
  end.yr = which(year==max(Nt$yr[Nt$estimation=="fit"])) 
  nT = length(year)
  years= year[1:end.yr]
  n.years = length(years)
  abundance = jara$settings$model.type
  GL = jara$settings$GL
  years = jara$yr
  if(indices[1]=="all"){
    indices = unique(jara$fits$name)
    n.indices = jara$settings$nI
  } else {
    if(length(indices[indices%in%unique(jara$fits$name)])<1) stop("non-existent index name provided")
    indices = indices[indices%in%unique(jara$fits$name)]
    n.indices = length(indices)
  }
  
  
  
  
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(as.png==TRUE){png(file = paste0(output.dir,"/TrjFits_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  
  # Total N
  m1 <- 0
  m2 <- max(jara$trj[jara$trj$name!="global",]$uci, na.rm = TRUE)
  
  if(abundance=="census"){
    plot(0, 0, ylim = c(m1, m2), xlim = c(min(years-1),max(years+1)), ylab = "Population Numbers", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
    cs = sample(seq(80,90,1))
    cols=paste0("gray",cs)
    col_line <- jara$settings$cols
    for(i in 1:n.indices){ 
      Nfit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$estimation=="fit",]
      polygon(x = c(years,rev(years)), y = c(Nfit$lci,rev(Nfit$uci)), col = gray(runif(1,0.5,0.9),0.5), border = "gray90")
    }
    for(i in 1:n.indices){
      Nfit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$estimation=="fit",]
      lines(years,Nfit$mu, type = "l",col=col_line[i], lwd=1)
      fits= jara$fits[jara$fits$name%in%indices[i],]
      points(fits$year,fits$obs, bg = col_line[i],pch=21) 
    }
    posl = c(max(Nt[Nt$yr==min(years),"uci"]),max(Nt[Nt$yr==max(years),"uci"]))
  } else {  
    plot(0, 0, ylim = c(m1, m2), xlim =  c(min(years-1),max(years+1)), ylab = "Abudance Index", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
    cs = sample(seq(80,90,1))
    cols=paste0("gray",cs)
    col_line <- jara$settings$cols
    q.adj = jara$pars$median[-c(1:2)]
    
    polygon(x = c(years,rev(years)), y = c(Nt$lci,rev(Nt$uci)), col = "gray", border = "gray90")
    
    
    for(i in 1:n.indices)
    {
      fits= jara$fits[jara$fits$name==indices[i],]
      points(fits$year,fits$obs/q.adj[i], bg = col_line[i],col=col_line[i],lty=2,pch=21,type="o") 
     }
    lines(years,Nt$mu, type = "l",col=1, lwd=2)
  }
  posl = c(max(Nt[Nt$yr==min(years),"uci"]),max(Nt[Nt$yr==max(years),"uci"]))


   legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(indices), lty = c(1, rep(n.indices)), lwd = c(rep(-1,n.indices)),pch=c(rep(21,n.indices)), pt.bg = c(col_line[1:n.indices]), bty = "n", cex = 0.9,y.intersp = 0.8)
   axis(1,at=seq(min(years)-1,max(years)+5,ceiling(n.years/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
   if(as.png==TRUE) dev.off()

} # end of jrplot_trjfit()

#' jrplot_poptrj()
#'
#' Plots the estimated population trajector relative to GL   
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_poptrj <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE){
  
  cat(paste0("\n","><> jrplot_poptrj() - plots trends against GL blocks <><","\n"))
  
  Nt = jara$trj[jara$trj$name=="global",]
  year = Nt$yr
  end.yr = which(year==max(Nt$yr[Nt$estimation=="fit"])) 
  nT = length(year)
  years= year[1:end.yr]
  n.years = length(years)
  abundance = jara$settings$model.type
  GL = jara$settings$GL
  
    Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
    if(as.png==TRUE){png(file = paste0(output.dir,"/PopTrend_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    if(add==FALSE) par(Par)
  
  # Total N
  m1 <- 0
  m2 <- max(Nt$uci, na.rm = TRUE)
  
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(year-1),max(year+1)), ylab = paste("Population",ifelse(abundance=="census","size","trend")), xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n",xaxs="i",yaxs="i")
  axis(1,at=seq(min(year),max(year)+5,ceiling(length(year)/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
  
  polygon(x = c(year,rev(year)), y = c(Nt$lci,rev(Nt$uci)), col = gray(0.6,0.3), border = "gray90")
  #add trend
  polygon(x = c(years,rev(years)), y = c(Nt$lci[1:end.yr],Nt$uci[end.yr:1]), col = gray(0.7,0.5),border = "gray90")
  lines(year[end.yr:nT],Nt$mu[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
  lines(years,Nt$mu[1:end.yr], type = "l",col=1,lwd=2)
  if(n.years-3*GL-1>0) lines(rep(year[n.years]-3*GL,2),c(0,m2*.93),lty=2,col=2)
  lines(rep(year[n.years]-1*GL,2),c(0,m2*.93),lty=2,col=4,lwd=2)
  lines(rep(year[n.years]-2*GL,2),c(0,m2*.93),lty=2,col=3,lwd=2)
  lines(rep(year[n.years],2),c(0,m2)*.93,lty=2,col=1,lwd=2)
  if(n.years-3*GL-1>0) text(year[n.years]-3*GL,m2*.96,"-3GL",lwd=2)
  text(year[n.years]-2*GL,m2*.96,"-2GL")
  text(year[n.years]-GL,m2*.96,"-1xGL")
  text(year[n.years],m2*.96,paste0(year[n.years]))
  
  if(as.png==TRUE) dev.off()
} # end of jrplot_poptrj
  

#' jrplot_fits
#'
#' Plots observed and fitted indices with expexted CIs (dark grey) 
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param width plot width
#' @param height plot hight
#' @param plot.cex graphic option
#' @param indices names of indices to plot (default = "all")
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_fits <- function(jara, output.dir=getwd(),as.png=FALSE,single.plots=FALSE,width=NULL,height=NULL,indices="all",add=FALSE){
  
    cat(paste0("\n","><> jrplot_fits() - fits to abudance indices <><","\n"))
    years = jara$yr
    N = length(years)
    if(indices[1]=="all"){
      indices = unique(jara$fits$name)
      n.indices = jara$settings$nI
    } else {
      if(length(indices[indices%in%unique(jara$fits$name)])<1) stop("non-existent index name provided")
      indices = indices[indices%in%unique(jara$fits$name)]
      n.indices = length(indices)
    }
    series = 1:jara$settings$nI
    CPUE = jara$settings$y
    check.yrs = abs(apply(jara$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(i in 1:n.indices){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/Fits",jara$assessment,"_",jara$scenario,"_",indices[i],".png"), width = width, height = height,
                             res = 200, units = "in")}
        if(add==FALSE){
        if(as.png==TRUE | i==1) par(Par)
        }  
        # set observed vs predicted CPUE
        Yr = jara$yr
        Yr = min(Yr):max(Yr)
        yr = Yr-min(years)+1
        
        fit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$yr%in%years,c("mu","lci","uci")]
        mufit = mean(fit[,1])
        fit = fit/mufit
        fit.hat = fit[,1]/mufit
        
        cpue.i = jara$fits[jara$fits$name%in%indices[i],]$obs
        yr.i =  jara$fits[jara$fits$name%in%indices[i],]$year
        se.i =  jara$fits[jara$fits$name%in%indices[i],]$obs.err
        
        ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
        
        cord.x <- c(Yr,rev(Yr))
        cord.y <- c(fit[yr,2],rev(fit[yr,3]))
        
        # Plot Observed vs predicted CPUE
        plot(years,fit[,1],ylab="",xlab="",ylim=ylim,xlim=range(jara$yr),type='n',xaxt="n",yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)
        polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
        
        lines(Yr,fit[yr,1],lwd=2,col=1)
        if(jara$settings$SE.I  ==TRUE | max(jara$settings$SE2)>0.005){ 
          iv = c(-0.25,0.25)
          for(t in 1:length(yr.i)){
            lines(rep(yr.i[t],2),c(exp(log(cpue.i[t])-1.96*se.i[t])/mufit,exp(log(cpue.i[t])+1.96*se.i[t])/mufit))  
            lines(yr.i[t]+iv,rep(exp(log(cpue.i[t])-1.96*se.i[t])/mufit,2))  
            lines(yr.i[t]+iv,rep(exp(log(cpue.i[t])+1.96*se.i[t])/mufit,2))
          }  
        }
        points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")
        
        legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
        mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=0.7,cex=1)
        mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
        if(as.png==TRUE) dev.off()
      }
      } else {
        
        if(is.null(width)) width = 7
        if(is.null(height)) height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0)
        Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/Fits_",jara$assessment,"_",jara$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                             res = 200, units = "in")}
        par(Par)
        
        for(i in 1:n.indices){
          # set observed vs predicted CPUE
          Yr = jara$yr
          Yr = min(Yr):max(Yr)
          yr = Yr-min(years)+1
          
          fit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$yr%in%years,c("mu","lci","uci")]
          mufit = mean(fit[,1])
          fit = fit/mufit
          fit.hat = fit[,1]/mufit
          
          cpue.i = jara$fits[jara$fits$name%in%indices[i],]$obs
          yr.i =  jara$fits[jara$fits$name%in%indices[i],]$year
          se.i =  jara$fits[jara$fits$name%in%indices[i],]$obs.err
          
          ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
          
          cord.x <- c(Yr,rev(Yr))
          cord.y <- c(fit[yr,2],rev(fit[yr,3]))
          # Plot Observed vs predicted CPUE
          plot(years,fit[,1],ylab="",xlab="",ylim=ylim,xlim=range(jara$yr),type='n',xaxt="n",yaxt="n")
          axis(1,labels=TRUE,cex=0.8)
          axis(2,labels=TRUE,cex=0.8)
          polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
          lines(Yr,fit[yr,1],lwd=2,col=1)
          if(jara$settings$SE.I  ==TRUE | max(jara$settings$SE2)>0.005){ 
            iv = c(-0.25,0.25)
            for(t in 1:length(yr.i)){
            lines(rep(yr.i[t],2),c(exp(log(cpue.i[t])-1.96*se.i[t])/mufit,exp(log(cpue.i[t])+1.96*se.i[t])/mufit))  
            lines(yr.i[t]+iv,rep(exp(log(cpue.i[t])-1.96*se.i[t])/mufit,2))  
            lines(yr.i[t]+iv,rep(exp(log(cpue.i[t])+1.96*se.i[t])/mufit,2))
            }  
          }
            points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")
          
        legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
        }
        mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=0.75,cex=1)
        mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
        if(as.png==TRUE){dev.off()}
      }
      
  } # End of CPUE plot function

#' log(index) fits
#'
#' Plot of fitted CPUE indices on log-scale (r4ss-style)
#'
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param width plot width
#' @param height plot hight
#' @param indices names of indices to plot (default = "all")
#' @param add if TRUE par is not called to enable manual multiplots
#' @export

jrplot_logfits <- function(jara, output.dir=getwd(),as.png=FALSE,single.plots=FALSE,width=NULL,height=NULL,indices="all",add=FALSE){
    cat(paste0("\n","><> jrplot_logfits()  <><","\n"))
    
  years = jara$yr
  N = length(years)
  if(indices[1]=="all"){
    indices = unique(jara$fits$name)
    n.indices = jara$settings$nI
  } else {
    if(length(indices[indices%in%unique(jara$fits$name)])<1) stop("non-existent index name provided")
    indices = indices[indices%in%unique(jara$fits$name)]
    n.indices = length(indices)
  }
  series = 1:jara$settings$nI
  CPUE = jara$settings$y
  check.yrs = abs(apply(jara$residuals,2,sum,na.rm=TRUE))
  cpue.yrs = years[check.yrs>0]
  
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(i in 1:n.indices){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/logFits",jara$assessment,"_",jara$scenario,"_",indices[i],".png"), width = width, height = height,
                             res = 200, units = "in")}
        if(add==FALSE){
        if(as.png==TRUE | i==1) par(Par)
        }
        Yr = jara$yr
        Yr = min(Yr):max(Yr)
        yr = Yr-min(years)+1
        cpue.i = jara$fits[jara$fits$name%in%indices[i],]$obs
        hat.i = jara$fits[jara$fits$name%in%indices[i],]$hat
        yr.i =  jara$fits[jara$fits$name%in%indices[i],]$year
        se.i =  jara$fits[jara$fits$name%in%indices[i],]$obs.err
        
        
        ylim = (c(min((log(cpue.i)-1.96*se.i)), max((log(cpue.i)+1.96*se.i))))
        
        # Plot Observed vs predicted CPUE
        plot(yr.i,log(cpue.i),ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)
        
        lines(yr.i,log(hat.i),lwd=2,col=4)
        if(jara$settings$SE.I  ==TRUE | max(jara$settings$SE2)>0.005){ 
          iv = c(-0.25,0.25)
          for(t in 1:length(yr.i)){
            lines(rep(yr.i[t],2),c((log(cpue.i[t])-1.96*se.i[t]),(log(cpue.i[t])+1.96*se.i[t])))  
            lines(yr.i[t]+iv,rep((log(cpue.i[t])-1.96*se.i[t]),2))  
            lines(yr.i[t]+iv,rep((log(cpue.i[t])+1.96*se.i[t]),2))
          }}  
          points(yr.i,log(cpue.i),pch=21,xaxt="n",yaxt="n",bg="white")
          
        legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=1)
        mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=0.7,cex=1)
        mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
        if(as.png==TRUE){dev.off()}
      }
    } else { # single.plots = F
      if(is.null(width)) width = 7
      if(is.null(height)) height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0)
      Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
      if(as.png==TRUE){png(file = paste0(output.dir,"/logFits_",jara$assessment,"_",jara$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                           res = 200, units = "in")}
      par(Par)
      for(i in 1:n.indices){
        Yr = jara$yr
        Yr = min(Yr):max(Yr)
        yr = Yr-min(years)+1
        cpue.i = jara$fits[jara$fits$name%in%indices[i],]$obs
        hat.i = jara$fits[jara$fits$name%in%indices[i],]$hat
        yr.i =  jara$fits[jara$fits$name%in%indices[i],]$year
        se.i =  jara$fits[jara$fits$name%in%indices[i],]$obs.err
        yr.i = Yr[is.na(CPUE[,i])==F]
        
        
        ylim = (c(min((log(cpue.i)-1.96*se.i)), max((log(cpue.i)+1.96*se.i))))
        
        # Plot Observed vs predicted CPUE
        plot(yr.i,log(cpue.i),ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)
        
        lines(yr.i,log(hat.i),lwd=2,col=4)
        if(jara$settings$SE.I  ==TRUE | max(jara$settings$SE2)>0.005){ 
          iv = c(-0.25,0.25)
          for(t in 1:length(yr.i)){
            lines(rep(yr.i[t],2),c((log(cpue.i[t])-1.96*se.i[t]),(log(cpue.i[t])+1.96*se.i[t])))  
            lines(yr.i[t]+iv,rep((log(cpue.i[t])-1.96*se.i[t]),2))  
            lines(yr.i[t]+iv,rep((log(cpue.i[t])+1.96*se.i[t]),2))
          }}  
        points(yr.i,log(cpue.i),pch=21,xaxt="n",yaxt="n",bg="white")
        
        legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=1)
      }
      mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=0.75,cex=1)
      mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE){dev.off()}
      
    }

} # End of logfit
