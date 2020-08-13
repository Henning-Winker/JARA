#' jara2jags()
#'
#' Writes JAGS code for JARA into a temporary directory
#' @param jarainput JARA input object from build_jara()
#' @export
#' @author Henning Winker
jara2jags = function(jarainput){

  if(jarainput$settings$model.type=="census"){
  # JAGS MODEL
  sink(paste0(tempdir(),"/jara.jags"))
    cat("
    model {
    # Priors and constraints
    for(i in 1:nI)
    {
    mean.r[i] ~ dnorm(0, 0.001) 
    logN.est[1,i] ~ dnorm(log(Ninit[i]),pow(0.5,-2))   # Prior for initial population size with CV =100%
    }
    
    ")
    
    if(jarainput$settings$sigma.proc==TRUE){
      cat("
      # Process variance
      isigma2 <- isigma2.est 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
      ",append=TRUE)  
    }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)
           
           ",append=TRUE)}
    
    if(jarainput$settings$sigma.obs.est==TRUE){
      cat("
      # Obsevation variance
      for(j in 1:nvar)
      {
      itau2[j]~ dgamma(0.001,0.001)
      tau2[j] <- 1/itau2[j]
      }
      
      for(i in 1:nI)
      {
      for(t in 1:T)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]]
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i])
      }
      }
      ",append=TRUE)  
    }else{ cat(" 
      # Obsevation variance
           # Observation error
           itau2~ dgamma(2,2)
           tau2 <- 1/itau2
           note1 <- nvar
           note2 <- sets.var
           
           for(i in 1:nI)
           {
           for(t in 1:T)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2
           
           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)     
           TOE[t,i] <- sqrt(var.obs[t,i])
           
           }}
           
           ",append=TRUE)}
    
    # Run rest of code  
    cat("  
    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])
    pen.sigma <- ifelse(sigma>proc.pen[1],log(sigma)-log(proc.pen[1]),0) 
    penSig  ~ dnorm(pen.sigma,pow(0.2,-2))
    
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    for (i in 1:nI){
    logN.est[t+1,i] <- logN.est[t,i]+r[t,i]
    }}
    # Set last year r to r.mean
    for (i in 1:nI){
    for(t in 1:(EY-1)){
    rdev[t,i] ~ dnorm(0, isigma2) #T(proc.pen[2],proc.pen[3])
    # Theta-Logistic
    r[t,i] <- mean.r[i]+ rdev[t,i]- 0.5*sigma2
    }}
    
    
    # Carrying Capacity for Projections only
    for(i in 1:nI){
    logpk[i] ~ dnorm(log(pk.mu[i]),pow(pk.cv[i],-2))
    logK[i] <- log(1)-logpk[i]+max(logN.est[pk.y,i])
    K[i] <- exp(logK[i])
    }
    
               ",append=TRUE)
    
    if(jarainput$settings$prjr.type=="mean"){
      cat("  
    for (i in 1:nI){
    r.proj[i] <- mean.r[i]
    } 
      
     ",append=TRUE)}else{
       cat(" 
   for (i in 1:nI){
   r.proj[i] <- mean(mean.r[i]+ rdev[prjr,i]-0.5*sigma2)} 
   ",append=TRUE)  
     }
    
    
    if(jarainput$settings$proj.stoch==FALSE){
      cat("  
    for (i in 1:nI){
    for(t in EY:(T-1)){
    rdev[t,i] ~ dnorm(0, pow(0.01,-2))
    r[t,i] <- r.proj[i]*(1-pow(exp(logN.est[t,i]-logK[i]),theta))
      
    }}
    ",append=TRUE)}else{
      cat(" 
    for (i in 1:nI){
    for(t in EY:(T-1)){
    rdev[t,i] ~ dnorm(0, isigma2)T(proc.pen[2],proc.pen[3])
    r[t,i] <- r.proj[i]*(1-pow(exp(logN.est[t,i]-logK[i]),theta))+ rdev[t,i]-0.5*sigma2) 
    }}
    ",append=TRUE)  
    }
    
    cat("  
    
    lambda[1] <- 1
    # Observation process
    for (t in 1:T) {
    for(i in 1:nI){
    y[t,i] ~ dnorm(logN.est[t,i], ivar.obs[t,i])
    ppd[t,i] ~ dlnorm(logN.est[t,i], ivar.obs[t,i])
    }}
    
    # Population sizes on real scale
    for (t in 1:T) {
    for(i in 1:nI){
    N.est[t,i] <- exp(logN.est[t,i])
    }
    Ntot[t] <- sum(N.est[t,])
    logNtot[t] <- log(Ntot[t])
    }
    for(t in 2:T)
    {
    r.tot[t] <- logNtot[t]-logNtot[t-1]
    }
    } 
    ",fill = TRUE)
    sink()
  } else { #--------------------------Relative Abudance Model-------------------------------#
    
    sink(paste0(tempdir(),"/jara.jags"))
    
    cat("
    model {
      
      # Prior specifications  
      eps <- 0.0000000000001 # small constant    
      
      iq[1] ~ dgamma(1000,1000)
      q[1] <-  pow(iq[1],-1)
      logq[1] <- log(1)
      for(i in 2:nI){
      iq[i] ~ dgamma(0.001,0.001)
      q[i] <- pow(iq[i],-1)
      logq[i] <-  log(q[i])
      }
      
      ")
    
    if(jarainput$settings$sigma.proc==TRUE){
      cat("
        # Process variance
        isigma2 <- isigma2.est 
        sigma2 <- pow(isigma2,-1)
        sigma <- sqrt(sigma2)
        fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
        ",append=TRUE)  
    }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
             sigma2 <- pow(isigma2,-1)
             sigma <- sqrt(sigma2)
             
             ",append=TRUE)}
    
    if(jarainput$settings$sigma.obs.est==TRUE){
      cat("
        # Obsevation variance
      for(j in 1:nvar)
      {
      itau2[j]~ dgamma(0.001,0.001)
      tau2[j] <- 1/itau2[j]
      }
      
      for(i in 1:nI)
      {
      for(t in 1:T)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]]
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i])
      }
      }
        ",append=TRUE)  
    }else{ cat(" 
      # Obsevation variance
             # Observation error
             itau2~ dgamma(2,2)
             tau2 <- 1/itau2
             note1 <- nvar
             note2 <- sets.var
             
             for(i in 1:nI)
             {
             for(t in 1:T)
             {
             var.obs[t,i] <- SE2[t,i] # drop tau2
             fake.tau[t,i] <- tau2
             
             ivar.obs[t,i] <- 1/var.obs[t,i]
             # note total observation error (TOE)     
             TOE[t,i] <- sqrt(var.obs[t,i])
             
             }}
             
             ",append=TRUE)}
    
    # Run rest of code  
    cat("  
      # Process variance prior
      isigma2.est ~ dgamma(igamma[1],igamma[2])
      pen.sigma <- ifelse(sigma>proc.pen[1],log(sigma)-log(proc.pen[1]),0) 
      penSig  ~ dnorm(pen.sigma,pow(0.2,-2))
      
      # Priors and constraints
      logY.est[1] ~ dnorm(logY1, 1)       # Prior for initial population size
      
      mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
      
      # Likelihood
      # State process
      for (t in 1:(EY-1)){
      rdev[t] ~ dnorm(0, isigma2)T(proc.pen[2],proc.pen[3])
      r[t] <- mean.r+rdev[t]-0.5*sigma2   
      logY.est[t+1] <- logY.est[t] + r[t] 
      }
      
      # Carrying Capacity for Projections only
      logpk ~ dnorm(log(pk.mu),pow(pk.cv,-2))
      logK <- log(1)-logpk+max(logY.est[pk.y])
      K <- exp(logK)
      #logK <- log(pow(pk,-1))+max(I[pk.y])
      #K <- exp(logK)
      ",append=TRUE)
    
    if(jarainput$settings$prjr.type=="mean"){
      cat("  
      r.proj <- mean.r  
      prjr.dummy <- prjr   
    ",append=TRUE)}else{
      cat(" 
        r.proj <- mean(mean.r+rdev[prjr]-0.5*sigma2) 
        ",append=TRUE)  
    }
    
    
    if(jarainput$settings$proj.stoch==FALSE){
      cat("
      for (t in EY:(T)){
      rdev[t] ~ dnorm(0, pow(0.01,-2))
      r[t] <- r.proj*(1-pow(exp(logY.est[t]-logK),theta)) 
      logY.est[t+1] <- logY.est[t]+r[t]   
      }
      ",append=TRUE)} else {
        cat("
      for (t in EY:(T)){
            rdev[t] ~ dnorm(0, pow(isigma2,-2))T(0.5*proc.pen[2],0.5*proc.pen[3])
            r[t] <- r.proj*(1-pow(exp(logY.est[t]-logK),theta))+rdev[t]-0.5*sigma2 
            logY.est[t+1] <- logY.est[t] + r[t] 
          
      }
      ",append=TRUE)}  
    
    cat("

      #for(t in (EY):T){ # Apply penalty for projections only
      #penK[t] ~ dnorm(devK[t],pow(0.1,-2))
      #}
      
      # Observation process
      for (t in 1:EY) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      ppd[t,i] ~ dlnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      }}
      for (t in (EY+1):T) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i],pow(0.001,-2)) # does not effect projections
      ppd[t,i] ~ dlnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      }}
      

      # Population sizes on real scale
      for (t in 1:T) {
      Y.est[t] <- exp(logY.est[t])
      Ntot[t] <- Y.est[t]
      }
      
     } 
      ",fill = TRUE)
    sink()
  } # end model.type = "relative  
    
} # end of jara2jags()  


