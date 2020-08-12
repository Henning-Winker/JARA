#' JARA hindcasting function
#'
#' Wrapper to coduct histcasts for retrospective analysis and cross-validation
#' @param jarainput List of input variables as output by build_jabba()
#' MCMC settings
#' @param ni number of iterations
#' @param nt thinning interval of saved iterations
#' @param nb burn-in
#' @param nc number of mcmc chains
#' @param peels sequence of retro spective peels default 0:5
#' @param save.jara saves jara list as .rdata to output.dir
#' @param save.all saves the all posteriors as .rdata to output.dir (big file)
#' @param saves the all posteriors as .rdata to output.dir (big file)
#' @param save.csvs writes results into csv to output.dir
#' @param output.dir path to save plot. default is getwd()
#' @param save.jarafile saves jara model and convergence stats as .txt file (default = TRUE)
#' @return A result list containing estimates of JARA model input, settings and results
#' @param save.hc Save hindcast list output as .rdata
#' @param plotall if TRUE makes jara_plots() for each run    
#' @param speedup Reduces MCMC after setting runs 2+ inits to first "full" reference run
#' @param rm.yr if TRUE remove last year, if FALSE replace obs with NA
#' @return hc containing estimates of key joint results from all hindcast run 
#' @export
jara_hindcast = function(jarainput,
                          # MCMC settings
                          ni = 9000, # Number of iterations
                          nt = 2, # Steps saved
                          nb = 2000, # Burn-in
                          nc = 3, # number of chains
                          peels = 0:5, # retro peel option
                          save.jara = TRUE,
                          save.all = FALSE,
                          save.csvs  = FALSE,
                          output.dir = getwd(),
                          save.jarafile = TRUE,
                          save.hc = TRUE, 
                          plotall = FALSE, 
                          speedup = TRUE,
                          rm.yr = TRUE){
  
  
  # hindcast object define object  
  hc = list(assessment= jarainput$settings$assessment,scenario = jarainput$settings$scenario, yr=jarainput$data$yr,pyr = jarainput$data$pyr,peels=peels,trj = NULL,fits=NULL,posteiors=NULL,settings=NULL)
  hc$settings$cols = jarainput$settings$cols
  
  Scenario = jarainput$settings$scenario
  for(i in 1:length(peels)){
    if(rm.yr==TRUE){
      index = jarainput$data$I
      subset = 1:(nrow(index)-peels[i])
      if(jarainput$settings$SE.I){
        se= jarainput$data$se
      } else {
        se=NULL
      }
      jrin = build_jara(
      I = jarainput$data$I, 
      se = se,
      assessment = jarainput$settings$assessment,
      scenario = peels[i],
      model.type = jarainput$settings$model.type,
      GL=jarainput$settings$GL, # Generation length
      start.year = NA,
      end.year = max(jarainput$data$yr)-peels[i], 
      sigma.obs.est = jarainput$settings$sigma.obs.est, 
      fixed.obsE = jarainput$settings$fixed.obsE,
      sigma.proc.fixed = jarainput$settings$sigma.proc.fixed,
      proc.pen = jarainput$settings$proc.pen,  
      proj.mod = jarainput$settings$proj.mod,
      pk.prior = jarainput$settings$pk.prior,
      pk.yr = jarainput$settings$pk.yr,
      pk.i = jarainput$settings$pk.i,
      proj.r = jarainput$settings$proj.r,
      proj.yrs.user = jarainput$settings$proj.yrs.user,
      proj.stoch = jarainput$settings$proj.stoch)
      retros = 0
    } else {
      jrin = jarainput
      jrin$settings$scenario = peels[i]
      retros = peels[i]
    }
    
    if(i == 1 | speedup==FALSE){
      mci = ni
      mct = nt
      mcb = nb
      mcc = nc
          } else {
      mci = 4500
      mct = 1
      mcb = 1000
      mcc = 3
      init.values=TRUE
    } 
    fithc = fit_jara(jrin,save.jara=save.jara,save.jarafile = save.jarafile,
                       output.dir=output.dir,
                      save.all = save.all,
                     save.csvs  = save.csvs,
                      ni = mci, # Number of iterations
                      nt = mct, # Steps saved
                      nb = mcb, # Burn-in
                      nc = mcc, # number of chains
                      peels = retros) # retro peel option
    
    hc$trj = rbind(hc$trj,data.frame(factor=fithc$fits[1,1],level=fithc$fits[1,2],fithc$trj[,-c(1,2)])) 
    hc$fits = rbind(hc$fits,fithc$fits)
    hc$posteriors= rbind(hc$posteriors ,data.frame(factor=fithc$fits[1,1],level=fithc$fits[1,2],pop.change=fithc$posteriors[,1]))
    
    if(plotall==TRUE){
      jara_plots(fithc,output.dir = output.dir)
    }
    
  } # end of loop
  if(save.hc==TRUE){
    save(hc,file=paste0(output.dir,"/hc_",Scenario,".rdata"))  
  }
  return(hc)
} #end of hindcast function 
