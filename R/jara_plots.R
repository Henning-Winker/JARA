#' jrpar()
#'
#' Set the par() to options suitable for ss3diags multi plots   
#' @param mfrow determines plot frame set up
#' @param plot.cex cex graphic option
#' @param mai graphical par for plot margins
#' @param labs if TRUE margins are narrow 
#' @export
jrpar <- function(mfrow=c(1,1),plot.cex=1,mai=c(0.5,0.5,0.1,.1),omi = c(0.,0.,0.,0)+ 0.05,labs=FALSE){
  if(labs==F){
    mai=c(0.25,0.25,0.15,.15)
    omi = c(0.3,0.35,0.2,0.2)}
  par(list(mfrow=mfrow,mai = mai, mgp =c(2.,0.5,0),omi =omi, tck = -0.02,cex=plot.cex))
}



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
#' @param plot.cex cex graphic option
#' @param cols option to choose own colour palette
#' @export
jrplot_indices <- function(jarainput, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE,cols=NULL){

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

Ylim = c(0, max(exp(jarainput$jagsdata$y+1.96*sqrt(jarainput$jagsdata$SE2))  ,na.rm =T))
plot(years,years,type="n",xlim=c(min(years-1),max(years+1)),ylab=ifelse(abundance=="census","Counts","Abunance Index"),xlab="Year",ylim=Ylim, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
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

axis(1,at=seq(min(dat[,1]),max(dat[,1])+5,ceiling(length(dat[,1])/8)),tick=seq(min(dat[,1]),max(dat[,1]),ceiling(length(dat[,1])/8)),cex.axis=0.9)

posl = c(max(dat[1:3,-1],na.rm=T),max(dat[(length(years)):length(years),-1],na.rm=T))
legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(names(dat)[2:(nI+1)]), lty = c(1, rep(nI)), lwd = c(rep(-1,nI)),pch=c(rep(21,nI)), pt.bg = c(cols[1:nI]), bty = "n", cex = 0.9,y.intersp = 0.8)
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
#' @param iucn.cols to use iucn color recommendation or a brighter version if FALSE
#' @param criteria option to choose between IUCN A1 or A2 thresholds (A2 is default)  
#' @param Plot if FALSE then only threat status stats are returned 
#' @return IUCN classification 
#' @author Henning Winker, Richard Sherley and Nathan Pacoureau
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


#' jrplot_changes
#'
#' Plot mean rates of change (%) over 1, 2 anf 3 Generation lengths    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
#' @author Henning Winker and Richard Sherley

jrplot_changes <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE){
  
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
  cnam = c("All.yrs","1GL","2GL","3GL")  
  
  jcol = c(grey(0.5,0.6),rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
  plot(0,0,type="n",ylab="Density",xlab="Annual Rate of Change(%)",xaxt="n",cex.main=0.9,ylim=c(0,1.1*max(lymax)),xlim=quantile(as.matrix(lamdas),c(0.001,0.999)),xaxs="i",yaxs="i") 
  for(i in 1:ncol(rs)){
    x = get(paste0("xl",i))
    y = get(paste0("yl",i))
    polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
    mu.lamda = round(median(lamdas[,i]),10)
    lines(rep(mu.lamda,2),c(0,max(y)),col=c(1,4,3,2)[i],lwd=1,lty=1)
    
  }
  axis(1,at=seq(floor(min(lamdas)),ceiling(max(lamdas)),1),tick=seq(min(x),max(x),5),cex.axis=0.9)
  abline(v=0,lty=2)
  legend("topright", paste0(cnam[1:ncol(rs)]," = ",ifelse(round(apply(lamdas,2,median),2)>0,"+",""),round(apply(lamdas,2,median),2),"%"),pch=15,col=c(jcol),bty="n")
  if(as.png==TRUE) dev.off()
} # End rate of change plot  


#' jrplot_state
#'
#' Plots current or projected population change relative to the first year    
#' @param jara output list from fit_jara
#' @param type final year type = c("current","projected","both")
#' @param ref.yr year(s) used as reference; default is avg. first 3 years
#' @param extinction threshold for extiction classification default < 1% or ref.yr
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @return Relative state 
#' @export
#' @author Henning Winker and Richard Sherley

jrplot_state <- function(jara, type=NULL,ref.yr=NULL,extinction=0.01, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE){
  
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
  plot(1,1)
  
  lxrange = ifelse(lxrange<0,0,lxrange)
  jcol = c(grey(0.4,0.6),rgb(1,0,0,0.6))
  plot(0,0,type="n",ylab="Density",xlab=paste("Relative State"),xaxt="n",cex.main=0.9,ylim=c(0,1.15*max(lymax)),xlim=c(0,max(lxmax)),xaxs="i",yaxs="i") 
  for(i in 2:1){
    x = get(paste0("xl",i))
    y = get(paste0("yl",i))
    polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=jcol[i],border=0)
    mu = round(median(states[,i]),10)
    lines(rep(mu,2),c(0,max(lymax*c(1.05,1.0)[i])),col=c(1,2)[i],lwd=1,lty=c(1))
    text(mu,max(lymax*c(1.1,1.05)[i]),c(end.yr,prj.yr)[i],cex=0.9)
  }
  axis(1,at=seq(0,ceiling(max(states)),0.2),cex.axis=0.9)
  lines(rep(1,2),c(0,max(lymax*c(1.05,1.0)[1])),col=1,lwd=1,lty=2)
  text(1,max(lymax*c(1.1)),min((ref.yr)),cex=0.9)
  
  #cnam = c(bquote("Cur"[.(end.yr) ~ "/" ~ .(min(ref.yr))] ~ "=" ~ .(round(median(states[,1]),2))),
  #         bquote("Prj"[.(end.yr) ~ "/" ~ .(min(ref.yr))]~ "=" ~ .(round(median(states[,2]),2))))
  cnam = c(paste0("Cur = ",round(median(states[,1]),2)),paste0("Proj = ",round(median(states[,2]),2)))
  if(type =="current") type.id = 1 
  if(type =="projected") type.id = 2 
  if(type =="both") type.id = 1:2 
  legend("topright",cnam[type.id],pch=15,col=c(jcol),bty="n")
  quants =apply(states,2,quantile,c(0.5,0.05,0.98))
  state = NULL
  state$state = data.frame(State=cnam,year=c(end.yr,prj.yr),median=quants[1,],lci=quants[2,],uci=quants[3,])
  if(type=="current"){state$prob.pextinct = "Requires projection horizon"} else {
    state$prob.extinct = sum(ifelse(states[,2]<extinction,1,0))/nrow(states)
  }
  
  if(as.png==TRUE) dev.off()
} # End rate of change plot  



#' jrplot_r
#'
#' Plot mean rates r = log(lamda) over 1, 2 anf 3 Generation lengths    
#' @param jara output list from fit_jara
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param width plot width
#' @param height plot hight
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_r <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE){
  
  cat(paste0("\n","><> jrplot_r() - r = log(lamda) over  1, 2, 3 x GL <><","\n"))
  
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
  cnam = c("All.yrs","1GL","2GL","3GL")  
  
  jcol = c(grey(0.5,0.6),rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
  plot(0,0,type="n",ylab="Density",xlab="r",xaxt="n",cex.main=0.9,ylim=c(0,1.1*max(lymax)),xlim=quantile(as.matrix(lamdas),c(0.001,0.999)),xaxs="i",yaxs="i") 
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
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @param indices names of indices to plot (default = "all")
#' @param cols option to choose own colour palette
#' @export
jrplot_trjfit <- function(jara, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE,indices="all",cols=NULL){
  
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
    plot(0, 0, ylim = c(m1, m2), xlim =  c(min(years-1),max(years+1)), ylab = "Abudance Index", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
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
   axis(1,at=seq(min(years)-1,max(years)+5,ceiling(n.years/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
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
#' @param plot.cex cex graphic option
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_poptrj <- function(jara,plotGL =NULL, output.dir=getwd(),as.png=FALSE,width=5,height=4.5,plot.cex=1,add=FALSE){
  
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
  
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(year-1),max(year+1)), ylab = paste("Population",ifelse(abundance=="census","size","trend")), xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxt="n",xaxs="i",yaxs="i")
  axis(1,at=seq(min(year),max(year)+5,ceiling(length(year)/8)),tick=seq(min(year),max(year),5),cex.axis=0.9)
  
  polygon(x = c(year,rev(year)), y = c(Nt$lci,rev(Nt$uci)), col = gray(0.6,0.3), border = "gray90")
  #add trend
  polygon(x = c(years,rev(years)), y = c(Nt$lci[1:end.yr],Nt$uci[end.yr:1]), col = gray(0.7,0.5),border = "gray90")
  lines(year[end.yr:nT],Nt$mu[end.yr:nT], type = "l",col=2, lwd=2,lty=5)
  lines(years,Nt$mu[1:end.yr], type = "l",col=1,lwd=2)
  if(n.years-3*GL-1>0) lines(rep(year[n.years]-3*GL,2),c(0,m2*.93),lty=2,col=2)
  if(plotGL){
  lines(rep(year[n.years]-1*GL,2),c(0,m2*.93),lty=2,col=4,lwd=2)
  lines(rep(year[n.years]-2*GL,2),c(0,m2*.93),lty=2,col=3,lwd=2)
  if(n.years-3*GL-1>0) text(year[n.years]-3*GL,m2*.96,"-3GL",lwd=2)
  text(year[n.years]-2*GL,m2*.96,"-2GL")
  text(year[n.years]-GL,m2*.96,"-1xGL")
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
#' @param plot.cex graphic option
#' @param indices names of indices to plot (default = "all")
#' @param index.label show index name in plot
#' @param add if TRUE par is not called to enable manual multiplots
#' @export
jrplot_fits <- function(jara,ppd=TRUE, output.dir=getwd(),as.png=FALSE,single.plots=FALSE,width=NULL,height=NULL,indices="all",index.label=TRUE,add=FALSE){
  
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


#' wrapper jara_plots function
#'
#' plots all routine JARA plots to output.dir if as.png=TRUE (default)
#'
#' @param jara output list from fit_jabba
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
  jrplot_iucn(jara,as.png=as.png,output.dir=output.dir)
}



