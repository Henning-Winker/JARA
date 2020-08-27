#' jrplot_pars()
#'
#' Set the par() to options suitable for JARA multi plots   
#' @param mfrow determines plot frame set up
#' @param mai graphical par for plot margins
#' @param omi graphical par for outer plot margins
#' @param plot.cex cex graphic option
#' @export
jrpar <- function(mfrow=c(1,1),plot.cex=1,mai=NULL,omi=c(0,0,0,0),labs=TRUE){
  if(is.null(mai)){
  if(labs){
    omi = omi+0.05
    mai=c(0.5,0.45,0.15,.15)} else {
      omi = omi+0.2
      mai=c(0.35,0.2,0,.15)}
    }
  par(list(mfrow=mfrow,mai = mai, mgp =c(2.,0.5,0),omi = omi + 0.1, tck = -0.02,cex=0.8))
}

#' jrplot_indices
#'
#' Plot mean rates of change (%) over 1, 2 anf 3 Generation lengths    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param criteria A1 or A2 for decline
#' @param add if TRUE par is not called to enable manual multiplots
#' @param add.legend option to add legend
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param cols option to choose own colour palette
#' @export
jrplot_indices <- function(jarainput, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,xlab="Year",ylab="Abundance index",add=FALSE,add.legend=TRUE,cols=NULL){

years =   jarainput$data$yr  
abundance =jarainput$settings$model.type
dat = jarainput$data$I
se = as.matrix(sqrt(jarainput$jagsdata$SE2)[1:length(years),])
y = as.matrix(jarainput$jagsdata$y[1:length(years),])
nI = ncol(y) 
if(is.null(cols)) cols = jarainput$settings$cols
Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
if(as.png==TRUE){png(file = paste0(output.dir,"/AnnualRate_",jarainput$settings$assessment,"_",jarainput$settings$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
if(add==FALSE) par(Par)

Ylim = c(0, max(exp(jarainput$jagsdata$y+1.96*sqrt(jarainput$jagsdata$SE2)*1.05)  ,na.rm =T))
plot(years,years,type="n",xlim=c(min(years-1),max(years+2)),ylab=ylab,xlab=xlab,ylim=Ylim, frame = TRUE,xaxs="i",yaxs="i",xaxt="n")
iv = c(-0.25,0.25)
for(j in 1:nI){
  ds =runif(1,-0.2,0.2)
  for(t in 1:length(years)){
    lines(rep(years[t]+ds,2),c(exp(y[t,j]-1.96*se[t,j]),exp(y[t,j]+1.96*se[t,j])))
    lines(years[t]+ds+iv,rep(exp(y[t,j]-1.96*se[t,j]),2))  
    lines(years[t]+ds+iv,rep(exp(y[t,j]+1.96*se[t,j]),2))
    }  
    lines(years+ds,dat[,j+1],type="b",pch=21,cex=1.2,col=grey(0.4,0.7),bg=cols[j])
}

axis(1,at=seq(min(dat[,1]),max(dat[,1]),ceiling(length(dat[,1])/8)),tick=seq(min(dat[,1]),max(dat[,1]),ceiling(length(dat[,1])/8)),cex.axis=0.9)

posl = c(max(dat[1:3,-1],na.rm=T),max(dat[(length(years)):length(years),-1],na.rm=T))
if(add.legend) legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(names(dat)[2:(nI+1)]), lty = c(1, rep(nI)), lwd = c(rep(-1,nI)),pch=c(rep(21,nI)), pt.bg = c(cols[1:nI]), bty = "n", cex = 0.9,y.intersp = 0.8)
if(as.png==TRUE) dev.off()
} # End of index plot

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
#' @param legend.cex lengend size cex graphic option
#' @param iucn.cols to use iucn color recommendation or a brighter version if FALSE
#' @param criteria option to choose between IUCN A1 or A2 thresholds (A2 is default)  
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param Plot if FALSE then only threat status stats are returned 
#' @return IUCN classification 
#' @author Henning Winker, Richard Sherley and Nathan Pacoureau
#' @export
jrplot_iucn <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,legend.cex=0.9,criteria=c("A2","A1")[1],iucn.cols=TRUE,ylab="Density",xlab="Change (%)",add=FALSE,Plot=TRUE){
  
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
  
    plot(x1,y1,type="n",xlim=c(-100,min(max(30,quantile(change,.99)),1000)),ylim=c(0,max(y1*1.1)),ylab=ylab,xlab=xlab,cex.main=0.9,frame=TRUE,xaxt="n",yaxt="n",xaxs="i",yaxs="i")
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
    legend("right",c(paste0("CR (",CR,"%)"),paste0("EN (",EN,"%)"),
                     paste0("VU (",VU,"%)"),paste0("LC (",LC,"%)")),col=1,pt.bg=c("red","orange","yellow","green"),pt.cex=1.2,pch=22,bg="white",cex=legend.cex,y.intersp = 0.8,x.intersp = 0.8)
    text(ifelse(mean(change)< -80,-80,mean(change)),max(y1*1.05),paste0("Change = ",sign,mu.change,"%"),bg="white",cex=legend.cex+0.1)
    
    if(as.png==TRUE) dev.off()
  } # End IUCN Plot = TRUE
  categories = c("Pr.Decl","change3GL","CR","EN","VU","LC")  
  percentages = c(CR,EN,VU,LC)
  status= ifelse(which(percentages==max(percentages))==4 & max(percentages)<50,"NT",categories[3:6][which(percentages==max(percentages))])
  perc.risk = data.frame(perc.risk=c(Decline,mu.change,percentages))
  rownames(perc.risk)=categories
  out = list()
  out$perc.risk = t(perc.risk)
  out$status = status
  
  return(out)
}



#' jrplot_retroiucn
#'
#' Retrospective IUCN posterior plot for A1 or A2   
#' @param hc output list from jara_hindcast()
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param criteria A1 or A2 for decline
#' @param add if TRUE par is not called to enable manual multiplots
#' @param plot.cex cex graphic option
#' @param iucn.cols to use iucn color recommendation or a brighter version if FALSE
#' @param criteria option to choose between IUCN A1 or A2 thresholds (A2 is default)  
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param Plot if FALSE then only threat status stats are returned 
#' @param add.legend option add legend
#' @return Retrospectice IUCN classification 
#' @author Henning Winker, Richard Sherley and Nathan Pacoureau
#' @export
jrplot_retroiucn <- function(hc, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,criteria=c("A2","A1")[1],iucn.cols=TRUE,ylab="Change (%)",xlab="Year",add=FALSE,Plot=TRUE,add.legend=TRUE){
  
  cat(paste0("\n","><> jrplot_retroiucn() - returns retrospective analysis of threat classification <><","\n"))
  A1 = ifelse(criteria=="A1",TRUE,FALSE)
  # settings
  runs = hc$peels
  d = hc$posteriors
  ymax1 = quantile(hc$posteriors$pop.change[hc$posteriors$level==0],0.995)
  ymax2 = quantile(hc$posteriors$pop.change[hc$posteriors$level>0],0.985)
  ymax = max(ymax1,ymax2)
  ylim=c(-100,min(max(30,ymax),1000))
  xlim=c(0.5,length(hc$peels)+0.49)
  xall = hc$posteriors$pop.change
  if(iucn.cols==T){
    cols = c("#60C659","lightgreen","#F9E814","#FC7F3F","#D81E05")[c(1,3:5)] # green to red
  } else {
    cols=c("green","lightgreen","yellow","orange","red")[c(1,3:5)]  
  }

  if(Plot==TRUE){
    Par = list(mfrow=c(1,1),mar = c(3, 3, 0, 0), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.6),mex=0.8, tck = -0.02,cex=plot.cex)
    if(as.png==TRUE){png(file = paste0(output.dir,"/IUCNretro_",hc$assessment,"_",hc$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    if(add==FALSE) par(Par)
    
    plot(0,0,type="n",xlab=xlab,xlim=xlim,ylim=ylim,axes=F,xaxs = "i",yaxs="i",ylab=ylab)  
    axis(1,at=seq(0,length(hc$peels),1),labels=max(hc$yr)-rev(c(hc$peels,max(hc$peels+1))),cex.axis=0.8,mgp=c(2,0.5,0))
    axis(2,at=seq(-100,max(xall,30)+50,ifelse(max(xall,30)>150,50,25)),tick=seq(-100,max(xall,30)+50,ifelse(max(xall,30)>150,50,25)),cex.axis=0.8,mgp=c(2,0.5,0))
    box()
    out = NULL
    for(j in 1:length(runs)){
      #change = log(d[rev(runs) ==runs[j],]$pop.change+100)
      change = d[d$level ==rev(runs)[j],]$pop.change
      change=ifelse(change>ymax,max(ymax+50,45),change)
      den = stats::density(change,adjust=1)
      
      #y1 = exp(den$x)-100
      y1 = den$x
      x1 = den$y/max(den$y)+j-0.5
      # get status
      categories = c("Pr.Decl","change3GL","CR","EN","VU","LC")  
      CR = round(sum(ifelse(change< ifelse(A1,-90,-80),1,0))/length(change)*100,1)
      EN = round(sum(ifelse(change> ifelse(A1,-90,-80) & change< ifelse(A1,-70,-50),1,0))/length(change)*100,1)
      VU = round(sum(ifelse(change> ifelse(A1,-70,-50) & change< ifelse(A1,-50,-30),1,0))/length(change)*100,1)
      LC = round(sum(ifelse(change> -30,1,0))/length(change)*100,1)
      Change3xGL = round(median(change),3)# round(sum(ifelse(change< 0,1,0))/length(change)*100,1)
      percentages = c(CR,EN,VU,LC)
      status= ifelse(which(percentages==max(percentages))==4 & max(percentages)<50,"NT",categories[3:6][which(percentages==max(percentages))])
      out = rbind(out,data.frame(Year=max(hc$yr)-rev(runs)[j], Change3xGL,CR,EN,VU,LC,status))
      
      lc = c(ifelse(A1,-50,-30))
      polygon(c(x1[y1>=lc],rep(min(x1),length(y1[y1>=lc]))),c(y1[y1>=lc],rev(y1[y1>=lc])),col=cols[1] ,border=cols[1])
      vu = c(ifelse(A1,-70,-50))
      polygon(c(x1[y1<=lc & y1>=vu],rep(min(x1),length(y1[y1<=lc & y1>=vu]))),c(y1[y1<lc & y1>=vu],rev(y1[y1<lc & y1>=vu])),col=cols[2],border=cols[2])
      en =ifelse(A1,-90,-80)
      polygon(c(x1[y1<vu & y1>=en],rep(min(x1),length(y1[y1<vu & y1>=en]))),c(y1[y1<vu & y1>=en],rev(y1[y1<vu & y1>=en])),col=cols[3],border=cols[3])
      polygon(c(x1[y1<en],rep(min(x1),length(y1[y1<en]))),c(y1[y1<en],rev(y1[y1<en])),col=cols[4],border=cols[4])
      polygon(c(x1,rep(min(x1),length(x1))),c(y1,rev(y1)))
    }
    abline(h=0,lty=2)
    
    text(1:length(runs),par('usr')[4],(out$status),cex=0.8,pos=1,offset = 0.2)
    if(add.legend)legend(par('usr')[2]*1.01, quantile(par('usr')[3:4],0.6), bty='n', xpd=NA,
           c("LC","VU","EN","CR"),pch=15,col=c(cols),pt.cex=2,cex=0.9)
    
    
    if(as.png==TRUE) dev.off()
  } # End IUCN Plot = TRUE
  
  return(out)
}

#' jrplot_retrobias() to plot retrospective pattern
#'
#' Plots retrospective pattern of population trend
#' @param hc output from jara_hindast()
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param add if TRUE par is not called to enable manual multiplots
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex plotting parameter
#' @param Xlim  allows to "zoom-in" requires speficiation Xlim=c(first.yr,last.yr)
#' @param cols option to add colour palette 
#' @param legend.loc location of legend
#' @param add.legend option to add legend
#' @param show.rho shows rho statistic in plot
#' @return Mohn's rho statistic for several quantaties
#' @export
jrplot_retrobias <- function(hc,output.dir=getwd(),as.png=FALSE,width=5,height=4.5,add=FALSE,xlab="Year",ylab="Abundance",plot.cex=1,Xlim=NULL,cols=NULL,legend.loc="topright",add.legend=TRUE,show.rho = TRUE){
  
  cat(paste0("\n","><> jrplot_retrobias() - retrospective analysis <><","\n"))
  if(is.null(cols)) cols=hc$settings$cols
  Nt = hc$trj[hc$trj$name=="global",]
  Nt0 = Nt[Nt$level==0,]
  year = hc$y
  end.yr = which(year==max(Nt$yr[Nt$estimation=="fit"])) 
  nyrs = length(year)
  suby = 1:nyrs
  peels = hc$peels
  if(is.null(Xlim)) Xlim = c(min(year)-0.99,max(year)+0.99)
  
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(as.png==TRUE){png(file = paste0(output.dir,"/Retrobias_",hc$assessment,"_",hc$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  # Total N
  m1 <- 0
  m2 <- max(c(Nt0$uci,Nt$mu), na.rm = TRUE)*1.1
  
  plot(0, 0, ylim = c(m1, m2), xlim = Xlim, ylab = ylab, xlab = xlab, col = "black", type = "n", lwd = 2, frame = TRUE,xaxt="n",xaxs="i",yaxs="i")
  axis(1,at=seq(min(year),max(year),ceiling(length(year)/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
  
  polygon(x = c(year,rev(year)), y = c(Nt0$lci[suby],rev(Nt0$uci[suby])), col = gray(0.4,0.5), border = gray(0.4,0.5))
  #add trend
  lines(year,Nt$mu[suby], type = "l",col=1,lwd=2)
  diags = data.frame(rho=NULL,hcrho=NULL)
  for(i in 2:length(hc$peels)){
  Nti  = Nt[Nt$level==peels[i],]
  sub0 = 1:(end.yr-peels[i])
  sub1 = 1:(end.yr-peels[i]+1)
  sub2 = (end.yr-peels[i]+1)
  lines(year[sub0],Nti$mu[sub0], type = "l",col=cols[i-1],lwd=2)
  lines(year[sub1],Nti$mu[sub1], type = "l",col=cols[i-1],lwd=2,lty=2)
  points(year[sub2],Nti$mu[sub2],bg=cols[i-1],pch=21,cex=0.8)
  rho = (Nti$mu[sub2-1]-Nt0$mu[sub2-1])/Nt0$mu[sub2-1]
  hcrho = (Nti$mu[sub2]-Nt0$mu[sub2])/Nt0$mu[sub2]
  diags = rbind(diags,data.frame(rho=rho,hcrho=hcrho))
  
  }
  lines(year,Nt$mu[suby], type = "l",col=1,lwd=1)
  
  if(add.legend) legend(legend.loc,paste(year[nyrs-peels]),col=c(1,cols),bty="n",cex=0.7,pt.cex=0.7,lwd=c(2,rep(1,length(peels))))
  diags = rbind(diags,data.frame(rho=mean(diags$rho),hcrho=mean(diags$hcrho)))
  mrho = round(diags[nrow(diags),1],2)
  mhcrho = round(diags[nrow(diags),2],2)
  
  if(show.rho) 
    legend("top",legend=paste0("Bias = ",mrho," (",mhcrho,")"),bty="n",y.intersp = -0.2,cex=0.9)
  
  row.names(diags) <- c(year[nyrs-peels[-1]],"mu")   

  if(as.png==TRUE) dev.off()

  return(diags)
} # end of Retrospective Plot





#' jrplot_changes
#'
#' Plot mean rates of change (%) over 1, 2 anf 3 Generation lengths    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param legend.cex legend sizr graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
#' @author Henning Winker and Richard Sherley

jrplot_changes <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,ylab ="Density",xlab="Annual rate of change (%)",plot.cex=1,legend.cex=0.9,add=FALSE){
  
  cat(paste0("\n","><> jrplot_change() - %change over  1, 2, 3 x GL <><","\n"))
  
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(as.png==TRUE){png(file = paste0(output.dir,"/AnnualRate_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  rs = jara$posteriors[,-1]
  
  lamdas = (exp(rs)-1)*100
  
  lymax=rymax = lxrange = rxrange =NULL # maximum and range for plotting
  for(i in 1:ncol(rs)){
    den = stats::density(lamdas[,i],adjust=2)
    assign(paste0("xl",i),den$x)
    assign(paste0("yl",i),den$y)
    lymax=c(lymax,max(den$y))
    lxrange = c(lxrange,range(den$x))
    den = stats::density(rs[,i],adjust=2)
    assign(paste0("xr",i),den$x)
    assign(paste0("yr",i),den$y)
    rymax=c(rymax,max(den$y))
    rxrange = c(rxrange,range(den$x))
  }
  cnam = c("All.yrs","1G","2G","3G")  
  
  jcol = c(grey(0.5,0.6),rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
  plot(0,0,type="n",ylab=ylab,xlab=xlab,xaxt="n",cex.main=0.9,ylim=c(0,1.1*max(lymax)),xlim=quantile(as.matrix(lamdas),c(0.001,0.999)),xaxs="i",yaxs="i") 
  for(i in 1:ncol(rs)){
    x = get(paste0("xl",i))
    y = get(paste0("yl",i))
    polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
    mu.lamda = round(median(lamdas[,i]),10)
    lines(rep(mu.lamda,2),c(0,max(y)),col=c(1,4,3,2)[i],lwd=1,lty=1)
    
  }
  axis(1,at=seq(floor(min(lamdas)),ceiling(max(lamdas)),1),tick=seq(min(x),max(x),5),cex.axis=0.9)
  abline(v=0,lty=2)
  legend("topright", paste0(cnam[1:ncol(rs)]," = ",ifelse(round(apply(lamdas,2,median),2)>0,"+",""),round(apply(lamdas,2,median),2),"%"),pch=15,col=c(jcol),bty="n",cex=legend.cex,y.intersp = 0.8,x.intersp = 0.8)
  if(as.png==TRUE) dev.off()
} # End rate of change plot  


#' jrplot_state
#'
#' Plots current or projected population change relative to the first year    
#' @param jara output list from fit_jara
#' @param type final year type = c("current","projected","both")
#' @param ref.yr year(s) used as reference; default is avg. first 3 years
#' @param extinction threshold for extiction classification default or ref.yr
#' @param credibility Confidence interval credibility calue with default of 0.95 
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param legend.cex legend size graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @return Relative state 
#' @export
#' @author Henning Winker 
jrplot_state <- function(jara, type=NULL,ref.yr=NULL,
                         extinction=0.01,credibility=0.95, output.dir=getwd(),
                         as.png=FALSE,width=5,height=4.5,ylab = "Density",xlab="Relative state",plot.cex=1,legend.cex=0.9,add=FALSE){
  
  cat(paste0("\n","><> jrplot_state() - %change relative to reference year <><","\n"))
  
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(as.png==TRUE){png(file = paste0(output.dir,"/State_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  pdyn = jara$pop.posterior
  yrs = 1:ncol(pdyn)
  nyrs = length(yrs)
  yr = as.numeric(names(pdyn))
  if(is.null(ref.yr)) ref.yr = yr[1:3]
  end.yr = max(jara$yr) 
  prj.yr = max(jara$pyr)
  pop.ref = apply(pdyn[,which(yr%in%ref.yr)],1,mean)
  states =  cbind(pdyn[,which(yr%in%end.yr)]/pop.ref,pdyn[,which(yr%in%prj.yr)]/pop.ref)
  
  if(is.null(type)){
    type=ifelse(prj.yr-end.yr<3,"current","both") 
  }
  
  
  lymax=rymax = lxrange = lxmax =NULL # maximum and range for plotting
  for(i in 1:2){
    if(i == 1 & type =="current" | type== "both" |i == 2 & type =="projected"){
    den = stats::density(states[,i],adjust=2)
    assign(paste0("xl",i),den$x)
    assign(paste0("yl",i),den$y)
    lymax=c(lymax,max(den$y))
    lxmax = c(lxmax,quantile(states[,i],0.99))
    }}
  
  lxrange = ifelse(lxrange<0,0,lxrange)
  jcol = c(grey(0.4,0.6),rgb(1,0,0,0.6))
  plot(0,0,type="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex.main=0.9,ylim=c(0,1.22*max(lymax)),xlim=c(0,max(lxmax,1.1)),xaxs="i",yaxs="i",frame=FALSE) 
  for(i in 2:1){
    if(i == 1 & type =="current" | type== "both" |i == 2 & type =="projected"){
    x = get(paste0("xl",i))
    y = get(paste0("yl",i))
    polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
    mu = round(median(states[,i]),10)
    lines(rep(mu,2),c(0,max(lymax*c(1.05,1.0)[i])),col=c(1,2)[i],lwd=1,lty=c(1))
    text(max(mu,0.05),max(lymax*c(1.11,1.05)[i]),c(end.yr,prj.yr)[i],cex=0.9)
  }}
  axis(1,at=seq(0,ceiling(max(states)),0.2),cex.axis=0.9)
  axis(2,cex.axis=0.9)
  
  lines(rep(1,2),c(0,max(lymax*c(1.13,1.0)[1])),col=1,lwd=1,lty=2)
  text(1,max(lymax*c(1.18)),min((ref.yr)),cex=0.9)
  
  cnam = c(paste0("Cur = ",round(median(states[,1]),2)),paste0("Proj = ",round(median(states[,2]),2)))
  if(type =="current") type.id = 1 
  if(type =="projected") type.id = 2 
  if(type =="both") type.id = 1:2 
  legend("right",cnam[type.id],pch=15,col=c(jcol),box.col = "white",cex=legend.cex,y.intersp = 0.8,x.intersp = 0.8)
  mu =apply(states,2,quantile,c(0.5))
  quants = rbind(mu,HDInterval::hdi(states,credMass=credibility))
  box()
  state = NULL
  state$state = data.frame(State=cnam,year=c(end.yr,prj.yr),median=quants[1,],lci=quants[2,],uci=quants[3,])
  if(type=="current"){state$prob.pextinct = "Requires projection horizon"} else {
    state$prob.extinct = sum(ifelse(states[,2]<extinction,1,0))/nrow(states)
  }
  
  if(as.png==TRUE) dev.off()
  return(state)
} # End rate of change plot  



#' jrplot_r
#'
#' Plot mean rates r = log(lamda) over 1, 2 anf 3 Generation lengths    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_r <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,xlab="r",ylab="Density",plot.cex=1,add=FALSE){
  
  cat(paste0("\n","><> jrplot_r() - r = log(lamda) over  1, 2, 3 x G <><","\n"))
  
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(as.png==TRUE){png(file = paste0(output.dir,"/rchange_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                       res = 200, units = "in")}
  if(add==FALSE) par(Par)
  
  rs = jara$posteriors[,-1]
  
  lamdas = rs
  
  lymax=rymax = lxrange = rxrange =NULL # maximum and range for plotting
  for(i in 1:ncol(rs)){
    den = stats::density(lamdas[,i],adjust=2)
    assign(paste0("xl",i),den$x)
    assign(paste0("yl",i),den$y)
    lymax=c(lymax,max(den$y))
    lxrange = c(lxrange,range(den$x))
    den = stats::density(rs[,i],adjust=2)
    assign(paste0("xr",i),den$x)
    assign(paste0("yr",i),den$y)
    rymax=c(rymax,max(den$y))
    rxrange = c(rxrange,range(den$x))
  }
  cnam = c("All.yrs","1G","2G","3G")  
  
  jcol = c(grey(0.5,0.6),rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
  plot(0,0,type="n",ylab=ylab,xlab=xlab,xaxt="n",cex.main=0.9,ylim=c(0,1.1*max(lymax)),xlim=quantile(as.matrix(lamdas),c(0.001,0.999)),xaxs="i",yaxs="i") 
  for(i in 1:ncol(rs)){
    x = get(paste0("xl",i))
    y = get(paste0("yl",i))
    polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
    mu.lamda = round(median(lamdas[,i]),10)
    lines(rep(mu.lamda,2),c(0,max(y)),col=c(1,4,3,2)[i],lwd=1,lty=1)
    
  }
  axis(1,at=seq(floor(min(x)),ceiling(max(x)),0.1),tick=seq(min(x),max(x),0.1),cex.axis=0.9)
  abline(v=0,lty=2)
  legend("topright", paste0(cnam[1:ncol(rs)]," = ",ifelse(round(apply(lamdas,2,median),2)>0,"+",""),round(apply(lamdas,2,median),3)),pch=15,col=c(jcol),bty="n")
  if(as.png==TRUE) dev.off()
} # End rate of change plot  


#' jrplot_trjfit()
#'
#' Plots the estimated and predicted population trajectors on the same scale    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @param indices names of indices to plot (default = "all")
#' @param cols option to choose own colour palette
#' @export
jrplot_trjfit <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,ylab="Abundance",xlab="Year",plot.cex=1,add=FALSE,indices="all",cols=NULL){
  
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
    plot(0, 0, ylim = c(m1, m2), xlim = c(min(years-1),max(years+2)), ylab = ylab, xlab = xlab, col = "black", type = "n", lwd = 2, frame = TRUE,xaxs="i",yaxs="i",xaxt="n")
    cs = sample(seq(80,90,1))
    
    colci=paste0("gray",cs)
    if(is.null(cols)) cols <- jara$settings$cols
    for(i in 1:n.indices){ 
      Nfit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$estimation=="fit",]
      polygon(x = c(years,rev(years)), y = c(Nfit$lci,rev(Nfit$uci)), col = gray(runif(1,0.5,0.9),0.5), border = "gray90")
    }
    for(i in 1:n.indices){
      Nfit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$estimation=="fit",]
      lines(years,Nfit$mu, type = "l",col=cols[i], lwd=1)
      fits= jara$fits[jara$fits$name%in%indices[i],]
      points(fits$year,fits$obs, bg = cols[i],pch=21) 
    }
    posl = c(max(Nt[Nt$yr==min(years),"uci"]),max(Nt[Nt$yr==max(years),"uci"]))
  } else {  
    plot(0, 0, ylim = c(m1, m2), xlim =  c(min(years-1),max(years+2)), ylab = ylab, xlab = xlab, col = "black", type = "n", lwd = 2, frame = TRUE,xaxs="i",yaxs="i",xaxt="n")
    cs = sample(seq(80,90,1))
    colci=paste0("gray",cs)
    if(is.null(cols)) cols <- jara$settings$cols
    q.adj = jara$pars$median[-c(1:2)]
    
    polygon(x = c(years,rev(years)), y = c(Nt$lci,rev(Nt$uci)), col = "gray", border = "gray90")
    
    
    for(i in 1:n.indices)
    {
      fits= jara$fits[jara$fits$name==indices[i],]
      points(fits$year,fits$obs/q.adj[i], bg = cols[i],col=cols[i],lty=2,pch=21,type="p") 
     }
    lines(years,Nt$mu, type = "l",col=1, lwd=2)
  }
  posl = c(max(Nt[Nt$yr==min(years),"uci"]),max(Nt[Nt$yr==max(years),"uci"]))


   legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(indices), lty = c(1, rep(n.indices)), lwd = c(rep(-1,n.indices)),pch=c(rep(21,n.indices)), pt.bg = c(cols[1:n.indices]), bty = "n", cex = 0.9,y.intersp = 0.8)
   axis(1,at=seq(min(years),max(years),ceiling(n.years/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
   if(as.png==TRUE) dev.off()

} # end of jrplot_trjfit()

#' jrplot_poptrj()
#'
#' Plots the estimated population trajector relative to GL   
#' @param jara output list from fit_jara
#' @param plotGL TRUE/FALSE indicates Generation Length
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_poptrj <- function(jara,plotGL =NULL, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,xlab="Year",ylab="Abundance",plot.cex=1,add=FALSE){
  
  cat(paste0("\n","><> jrplot_poptrj() - plots trends against GL blocks <><","\n"))
  
  Nt = jara$trj[jara$trj$name=="global",]
  year = Nt$yr
  end.yr = which(year==max(Nt$yr[Nt$estimation=="fit"])) 
  nT = length(year)
  years= year[1:end.yr]
  n.years = length(years)
  abundance = jara$settings$model.type
  GL = jara$settings$GL
  if(is.null(plotGL) & is.null(jara$settings$proj.yrs.user)){
    plotGL = TRUE
  } else if(is.null(plotGL) & is.null(jara$settings$proj.yrs.user)==FALSE){
    plotGL = FALSE
  }
  
    Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
    if(as.png==TRUE){png(file = paste0(output.dir,"/PopTrend_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    if(add==FALSE) par(Par)
  
  # Total N
  m1 <- 0
  m2 <- max(Nt$uci, na.rm = TRUE)
  
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(year-1),max(year+2)), ylab =ylab, xlab = xlab, col = "black", type = "n", lwd = 2, frame = TRUE,xaxt="n",xaxs="i",yaxs="i")
  axis(1,at=seq(min(year),max(year),ceiling(length(year)/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
  
  polygon(x = c(year,rev(year)), y = c(Nt$lci,rev(Nt$uci)), col = gray(0.6,0.3), border = "gray90")
  #add trend
  polygon(x = c(years,rev(years)), y = c(Nt$lci[1:end.yr],Nt$uci[end.yr:1]), col = gray(0.7,0.5),border = "gray90")
  lines(year[end.yr:nT],Nt$mu[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
  lines(years,Nt$mu[1:end.yr], type = "l",col=1,lwd=2)
  if(n.years-3*GL-1>0) lines(rep(year[n.years]-3*GL,2),c(0,m2*.93),lty=2,col=2)
  if(plotGL){
  lines(rep(year[n.years]-1*GL,2),c(0,m2*.93),lty=2,col=4,lwd=2)
  lines(rep(year[n.years]-2*GL,2),c(0,m2*.93),lty=2,col=3,lwd=2)
  if(n.years-3*GL-1>0) text(year[n.years]-3*GL,m2*.96,"3G",lwd=2)
  text(year[n.years]-2*GL,m2*.96,"2G")
  text(year[n.years]-GL,m2*.96,"1G")
  } else {
    cat("\n","><> Not showing GLs as assessment period is specified by user","\n")
  }
  lines(rep(year[n.years],2),c(0,m2)*.93,lty=2,col=1,lwd=2)
  text(year[n.years],m2*.96,paste0(year[n.years]))
  
  if(as.png==TRUE) dev.off()
} # end of jrplot_poptrj
  

#' jrplot_fits
#'
#' Plots observed and fitted indices with expexted CIs (dark grey) 
#' @param jara output list from fit_jara
#' @param ppd show posterior predictive distribution
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param width plot width
#' @param height plot hight
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param plot.cex graphic option
#' @param indices names of indices to plot (default = "all")
#' @param index.label show index name in plot
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_fits <- function(jara,ppd=TRUE, output.dir=getwd(),as.png=FALSE,single.plots=FALSE,width=NULL,height=NULL,ylab="Normalized index",xlab="Year",indices="all",index.label=TRUE,add=FALSE){
  
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
        cord.lpp <- c(fit[yr,"lpp"],rev(fit[yr,"upp"]))
        
        # Plot Observed vs predicted CPUE
        plot(years,fit[,1],ylab="",xlab="",ylim=ylim,xlim=range(jara$yr),type='n',xaxt="n",yaxt="n")
        axis(1,labels=TRUE,cex=0.8)
        axis(2,labels=TRUE,cex=0.8)
        if(ppd) polygon(cord.x,cord.lpp,col=grey(0.5,0.5),border=0,lty=2)
        polygon(cord.x,cord.y,col=grey(ifelse(ppd,0.4,0.6),0.8),border=grey(0.4,0.8),lty=2)
        
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
        
        if(index.label==TRUE)legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
        mtext(paste(xlab), side=1, outer=TRUE, at=0.5,line=0.7,cex=1)
        mtext(paste(ylab), side=2, outer=TRUE, at=0.5,line=1,cex=1)
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
          
          fit = jara$trj[jara$trj$name%in%indices[i] & jara$trj$yr%in%years,c("mu","lci","uci","lpp","upp")]
          mufit = mean(fit[,1])
          fit = fit/mufit
          fit.hat = fit[,1]/mufit
          
          cpue.i = jara$fits[jara$fits$name%in%indices[i],]$obs
          yr.i =  jara$fits[jara$fits$name%in%indices[i],]$year
          se.i =  jara$fits[jara$fits$name%in%indices[i],]$obs.err
          
          ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
          
          cord.x <- c(Yr,rev(Yr))
          cord.y <- c(fit[yr,2],rev(fit[yr,3]))
          cord.lpp <- c(fit[yr,"lpp"],rev(fit[yr,"upp"]))
          
          # Plot Observed vs predicted CPUE
          plot(years,fit[,1],ylab="",xlab="",ylim=ylim,xlim=range(jara$yr),type='n',xaxt="n",yaxt="n")
          axis(1,labels=TRUE,cex=0.8)
          axis(2,labels=TRUE,cex=0.8)
          if(ppd) polygon(cord.x,cord.lpp,col=grey(0.5,0.5),border=0,lty=2)
          
          polygon(cord.x,cord.y,col=grey(ifelse(ppd,0.4,0.6),0.8),border=grey(0.4,0.8),lty=2)
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
          
            if(index.label==TRUE) legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
        }
        mtext(paste(xlab), side=1, outer=TRUE, at=0.5,line=0.75,cex=1)
        mtext(paste(ylab), side=2, outer=TRUE, at=0.5,line=1,cex=1)
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
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param indices names of indices to plot (default = "all")
#' @param add if TRUE par is not called to enable manual multiplots
#' @export

jrplot_logfits <- function(jara, output.dir=getwd(),as.png=FALSE,single.plots=FALSE,width=NULL,height=NULL,ylab="Log(index)",xlab="Year",indices="all",add=FALSE){
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
        mtext(paste(xlab), side=1, outer=TRUE, at=0.5,line=0.7,cex=1)
        mtext(paste(ylab), side=2, outer=TRUE, at=0.5,line=1,cex=1)
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
      mtext(paste(xlab), side=1, outer=TRUE, at=0.5,line=0.75,cex=1)
      mtext(paste(ylab), side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE){dev.off()}
      
    }

} # End of logfit

#' joint residual plot
#'
#' plots residuals for all indices as boxplot with a loess showing systematic trends
#'
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param width plot width
#' @param height plot hight
#' @param cols option to add colour palette 
#' @export
jrplot_residuals <- function(jara,output.dir=getwd(),as.png = FALSE,add=FALSE,ylab="log Residuals",xlab="Year", width = 5, height = 3.5,cols=NULL){
  
    cat(paste0("\n","><> jrplot_residuals() - Joint residual plot  <><","\n"))
    if(is.null(cols)) cols = jara$settings$cols 
    years = jara$yr
    check.yrs = abs(apply(jara$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    Resids = jara$residuals
    Yr = jara$yr
    n.years = length(Yr)
    n.indices = jara$settings$nI
    indices = unique(jara$diags$name)
    series = 1:jara$settings$nI
    
    # Joint-residual plot
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    if(as.png==TRUE){png(file = paste0(output.dir,"/Residuals_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                         res = 200, units = "in")}
    if(add==FALSE) par(Par)
    
    
    plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab=ylab,xlab=xlab)
    boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
    abline(h=0,lty=2)
    positions=runif(n.indices,-0.2,0.2)
    
    for(i in 1:n.indices){
      for(t in 1:n.years){
        lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=cols[i])}
      points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=cols[i])}
    mean.res = apply(Resids,2,mean,na.rm =TRUE)[Yr%in%cpue.yrs]
    smooth.res = predict(stats::loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
    lines(cpue.yrs,smooth.res,lwd=2)
    # get degree of freedom
    Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
    
    RMSE = round(jara$stats[2,2],2)
    
    legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
    legend('bottomright',c(paste(jara$indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jara$settings$col[series],1),lwd=c(rep(-1,n.indices),2))
    if(as.png==TRUE){dev.off()}
} # End of functions


#' JARA runs test plots
#'
#' Residual diagnostics with runs test p-value and 3xsigma limits
#' @param jara output list from fit_jara
#' @param index option to plot specific indices (numeric & in order)
#' @param alternative hypothesis undermixing "less" or both "two-sided"
#' @param output.dir directory to save plots
#' @param add if true par() is surpressed within the plot function
#' @param as.png save as png file of TRUE
#' @param single.plots if TRUE plot invidual fits else make multiplot
#' @param ylab option to change y-axis label
#' @param xlab option to change x-axis label
#' @param width plot width
#' @param height plot hight
#' @export
jrplot_runstest <- function(jara,indices="all",alternative="less", output.dir=getwd(),add=FALSE,as.png=FALSE,single.plots=FALSE,ylab="Residuals",xlab="Year",width=NULL,height=NULL){
  
    cat(paste0("\n","><> jrplot_runstest()   <><","\n"))
    
  
  years = jara$yr
  N = length(years)
  if(indices[1]=="all"){
    indices = unique(jara$fits$name)
    n.indices = jara$settings$nI
    index = 1:n.indices
  } else {
    if(length(indices[indices%in%unique(jara$fits$name)])<1) stop("non-existent index name provided")
    indices = indices[indices%in%unique(jara$fits$name)]
    n.indices = length(indices)
    index = which(unique(jara$fits$name)%in%indices)
  }
  series = 1:jara$settings$nI
  CPUE = jara$settings$y
  check.yrs = abs(apply(jara$residuals,2,sum,na.rm=TRUE))
  cpue.yrs = years[check.yrs>0]
  
  Resids = jara$residuals

    
    if(single.plots==TRUE){
      if(is.null(width)) width = 5
      if(is.null(height)) height = 3.5
      for(i in 1:n.indices){
        Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
        if(as.png==TRUE){png(file = paste0(output.dir,"/ResRunsTests_",jara$assessment,"_",jara$scenario,"_",indices[i],".png"), width = width, height = height,
                             res = 200, units = "in")}
        
        if(add==FALSE){
          if(as.png==TRUE | i==1) par(Par)
        }
        
        
        resid = (Resids[index[i],is.na(Resids[index[i],])==F])
        res.yr = years[is.na(Resids[index[i],])==F]
        runstest = jr_runs(x=as.numeric(resid),type="resid",altenative = alternative)
        # CPUE Residuals with runs test
        plot(res.yr,rep(0,length(res.yr)),type="n",ylim=c(min(-1,runstest$sig3lim[1]*1.25),max(1,runstest$sig3lim[2]*1.25)),lty=1,lwd=1.3,xlab="Year",ylab="Residuals")
        abline(h=0,lty=2)
        lims = runstest$sig3lim
        cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
        rect(min(years-1),lims[1],max(years+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
        for(j in 1:length(resid)){
          lines(c(res.yr[j],res.yr[j]),c(0,resid[j]))
        }
        points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
        legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
        #mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
        #mtext(expression(log(cpue[obs])-log(cpue[pred])), side=2, outer=TRUE, at=0.5,line=1,cex=1)
        if(as.png==TRUE){dev.off()}
      } # end of loop
    } else { # single.plot = FALSE
      if(is.null(width)) width = 7
      if(is.null(height)) height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0)
      Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
      if(as.png==TRUE){png(file = paste0(output.dir,"/ResRunsTests_",jara$assessment,"_",jara$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0),
                           res = 200, units = "in")}
      if(add==FALSE) par(Par)
      for(i in 1:n.indices){
        resid = (Resids[index[i],is.na(Resids[index[i],])==F])
        res.yr = years[is.na(Resids[index[i],])==F]
        runstest = jr_runs(x=as.numeric(resid),type="resid")
        # CPUE Residuals with runs test
        plot(res.yr,rep(0,length(res.yr)),type="n",ylim=c(min(-1,runstest$sig3lim[1]*1.25),max(1,runstest$sig3lim[2]*1.25)),lty=1,lwd=1.3,xlab="Year",ylab="Residuals")
        abline(h=0,lty=2)
        lims = runstest$sig3lim
        cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
        rect(min(years-1),lims[1],max(years+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
        for(j in 1:length(resid)){
          lines(c(res.yr[j],res.yr[j]),c(0,resid[j]))
        }
        points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
        legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
        
      }
      mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
      mtext("Residuals", side=2, outer=TRUE, at=0.5,line=1,cex=1)
      if(as.png==TRUE){dev.off()}
    }
    
    

  
} # end of runstest plot function


#' wrapper jara_plots function
#'
#' plots all routine JARA plots to output.dir if as.png=TRUE (default)
#'
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @export
jara_plots = function(jara,output.dir = getwd(),as.png=TRUE,statusplot ="kobe"){
  
  jrplot_fits(jara,as.png=as.png,output.dir=output.dir) 
  jrplot_logfits(jara,as.png=as.png,output.dir=output.dir)
  jrplot_trjfit(jara,as.png=as.png,output.dir=output.dir)
  jrplot_poptrj(jara,as.png=as.png,output.dir=output.dir)
  jrplot_r(jara,as.png=as.png,output.dir=output.dir)
  jrplot_changes(jara,as.png=as.png,output.dir=output.dir)
  jrplot_state(jara,as.png=as.png,output.dir=output.dir)
  jrplot_iucn(jara,as.png=as.png,output.dir=output.dir)
}



