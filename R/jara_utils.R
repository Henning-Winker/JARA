#' rqpois()
#'
#' Random generator for quasi-poisson variables    
#' @param n number of random observations
#' @param lambda vector of (non-negative) expected means
#' @param phi over-dispersion (variance) parameter 
#' @return random quasi-poisson observations
#' @export
rqpois = function(n, lambda, phi) {
  mu = lambda
  k = mu/phi/(1-1/phi)
  r = stats::rnbinom(n, mu = mu, size = k)
  return(r)
}

#' assign_iucn()
#'
#' Function to assign IUCN status based on decline and criteria     
#' @param change measured change over 3 generations 
#' @param criteria A1 or A2 for decline
#' @return IUCN status
#' @export
assign_iucn <- function(change,criteria=c("A2","A1")[1]){
  A = ifelse(criteria=="A1",1,2)
  status = ifelse(change< -80,"CR",ifelse(change>= c(-90,-80)[A] & change< c(-70,-50)[A],"EN",ifelse(change>= c(-70,-50)[A] & change< c(-50,-30)[A],"VU","LC"))) 
  return(status)
}  

#' iucn_estimators()
#'
#' Current estimators based on two points (p2) and log-linear regression (reg)      
#' @param I obserseved abundance index data.frame(yr,y)  
#' @param GL generation length
#' @param criteria A1 or A2 for decline
#' @return change over 3 x GL for methods p2 and reg
#' @export

iucn_estimators = function(I=NULL,GL=10,criteria=c("A2","A1")[1]){
  d. = data.frame(I)[,1:2]
  colnames(d.) = c("yr","y")
  p2obs = d.[nrow(d.),2]/d.[1,2]
  fit.reg = lm(log(y)~yr,d.)
  pr.reg = exp(predict(fit.reg))
  regobs = as.numeric(pr.reg[nrow(d.)]/pr.reg[1])
  p2change = (as.numeric(p2obs^(3*GL/length(d.$yr)))-1)*100
  regchange = (as.numeric(regobs^(3*GL/length(d.$yr)))-1)*100
  
  status = data.frame()
  out = list()
  out$reg.r = data.frame(est=coef(fit.reg)[2],se=summary(fit.reg)$coef[2,2])
  out$reg.pr = data.frame(yr=d.$yr,yobs=d.$y,yhat=pr.reg)
  out$change = data.frame(p2=p2change,reg=regchange)
  out$status = data.frame(p2=assign_iucn(p2change,criteria),reg=assign_iucn(regchange,criteria))
  return(out)
}  
  
#' iucn_frame
#'
#' Empty iucn color plot   
#' @param xylim determines upper x and y lims  
#' @param criteria A1 or A2 for decline
#' @param add if TRUE par is not called to enable manual multiplots
#' @param plot.cex cex graphic option
#' @param legend.cex lengend size cex graphic option
#' @param iucn.cols to use iucn color recommendation or a brighter version if FALSE
#' @param criteria option to choose between IUCN A1 or A2 thresholds (A2 is default)  
#' @author Henning Winker, Richard Sherley and Nathan Pacoureau
#' @export
iucn_frame <- function(xylim=c(100,1),plot.cex=1,legend.cex=0.9,criteria=c("A2","A1")[1],iucn.cols=TRUE,add=FALSE){
  
  A1 = ifelse(criteria=="A1",TRUE,FALSE)
  
  if(iucn.cols==T){
    cols = c("#60C659","lightgreen","#F9E814","#FC7F3F","#D81E05")[c(1,3:5)] # green to red
  } else {
    cols=c("green","lightgreen","yellow","orange","red")[c(1,3:5)]  
  }
  
  Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,0.5,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
  if(add==FALSE) par(Par)
  x1 = c(-100,xylim[1])
  plot(c(-100,xylim[1]),c(0,xylim[2]),type="n",ylab="Density",xlab="Change (%)",cex.main=0.9,frame=T,xaxt="n",yaxt="n",xaxs="i",yaxs="i")
  x2 = c(ifelse(A1,-50,-30),1500); y2 = c(0,5)
  polygon(c(x2,rev(x2)),c(rep(xylim[2],2),rev(rep(0,2))),col=cols[1],border=cols[1])
  x3 = c(ifelse(A1,-70,-50),x2[1])
  polygon(c(x3,rev(x3)),c(rep(xylim[2],2),rev(rep(0,2))),col=cols[2],border=cols[2])
  x4 = c(ifelse(A1,-90,-80),x3[1])
  polygon(c(x4,rev(x4)),c(rep(xylim[2],2),rep(0,2)),col=cols[3],border=cols[3])
  x5 = c(-100,x4[1])
  polygon(c(x5,rev(x5)),c(rep(xylim[2],2),rep(0,2)),col=cols[4],border=cols[4])
  
  axis(1,at=seq(-100,max(x1,30)+50,ifelse(max(x1,30)>150,50,25)),tick=seq(-100,max(x1,30),ifelse(max(x1,30)>150,50,25)))
  
} # End IUCN frame 



#' Function to do runs.test and 3 x sigma limits
#'
#' runs test is conducted with library(randtests)
#' @param x residuals from CPUE fits
#' @param type only c("resid","observations")
#' @param mixing c("less","greater","two.sided"). Default less is checking for postive autocorrelation only    
#' @return runs p value and 3 x sigma limits
#' @export
#' @author Henning Winker (JRC-EC) and Laurance Kell (Sea++)
jr_runs <- function(x,type=NULL,mixing="less") {
  if(is.null(type)) type="resid"
  if(type=="resid"){
    mu = 0}else{mu = mean(x, na.rm = TRUE)}
  alternative=c("two.sided","left.sided")[which(c("two.sided", "less")%in%mixing)]
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){
    # Make the runs test non-parametric
    runstest = randtests::runs.test(x,threshold = 0,alternative = alternative)
    if(is.na(runstest$p.value)) p.value =0.001
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }
  
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}


#' mixed.trend()
#'
#' Function to compile combined index
#' @param jara fit_jara output
#' @param run name qualifier of data.frame
#' @param refyr sets index relative to a reference year
#' @param type c("mu","pr"), note the probability is only available of "mixed.trends" 
#' @return data.frame
#' @export
#' @author Henning Winker (JRC-EC) 
mixed.trend <- function(jara,run="joint",refyr=FALSE,type=c("mu","pr")[1]){
  df =jara$pop.posterior
  if(type=="pr") df =jara$prop.posterior
  if(jara$settings$mixed.scale=="geomean")
      Obs = exp(aggregate(log(obs)~year,jara$fit,mean)$`log(obs)`)
  if(jara$settings$mixed.scale=="mean")
      Obs = (aggregate(obs~year,jara$fit,mean)$obs)
  
  joint = data.frame(Year=jara$yr,Obs=Obs,
                     n=aggregate(obs~year,jara$fits,length)$obs,
                     t(apply(df,2,quantile,c(0.025,0.25,0.5,0.75,0.975)))[1:length(jara$yr),])
  colnames(joint) = c("Year","Obs","n","5%","25%","50%","75%","95%")
  if(refyr){
    joint$Obs=joint$Obs/joint$`50%`[1]
    joint[,-c(1:2)] = joint[,-c(1:2)]/joint$`50%`[1]
  }
  joint$run = run
  return(joint)
}

#' dfidx()
#'
#' Function to compile fits into obs and fit data.frames
#' @param jara fit_jara output
#' @param run name qualifier of data.frame
#' @return data.frame
#' @export
#' @author Henning Winker (JRC-EC) 
dfidx = function(jara,run="obs"){
  dat=list()
  dat$i = jara$trj
  dat$i = dat$i[dat$i$name!="global" & dat$i$estimation=="fit",]
  dat$i$run=run
  dat$o = obs=jara$fits
  dat$o$run=run
  return(dat)
}


#' iucn 
#'
#' IUCN computes IUCN status for A1 or A2   
#' @param jara output list from fit_jara
#' @param criteria option to choose between IUCN A1 or A2 thresholds (A2 is default)  
#' @return IUCN classification 
#' @author Henning Winker, Richard Sherley and Nathan Pacoureau
#' @export
iucn <- function(jara, criteria=c("A2","A1")[1]){
  
  change= jara$posteriors$pop.change
  A1 = ifelse(criteria=="A1",TRUE,FALSE)
  mu.change = round(median(change),1)
  sign=""
  if(mu.change>0) sign="+"
  CR = round(sum(ifelse(change< ifelse(A1,-90,-80),1,0))/length(change)*100,1)
  EN = round(sum(ifelse(change> ifelse(A1,-90,-80) & change< ifelse(A1,-70,-50),1,0))/length(change)*100,1)
  VU = round(sum(ifelse(change> ifelse(A1,-70,-50) & change< ifelse(A1,-50,-30),1,0))/length(change)*100,1)
  LC = round(sum(ifelse(change> -30,1,0))/length(change)*100,1)
  Decline = round(sum(ifelse(change< 0,1,0))/length(change)*100,1)
  
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

#' jrroc
#'
#' ROC curve Function
#' @param x observed measure
#' @param y true or reference measure
#' @param rfp reference point threshold
#' @return roc list
#' @export 
jrroc <- function(x,y,rfp=NULL){
  if(is.null(rfp)) rfp = 0
  Y <- y[order(x, decreasing=F)]
  curve=data.frame(x=x[order(x, decreasing=F)],
                   fpr=cumsum(!Y)/sum(!Y), # false positives
                   tpr=cumsum(Y)/sum(Y)) # true positives
  pt = curve[which(abs(curve$x-rfp)==min(abs(curve$x-rfp)))[1],2:3]
  return(list(curve=curve[,2:3],pt=pt,xo=curve[,1]))
}
