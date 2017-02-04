################################
# Population based Monte Carlo #
################################

source("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana\\simple_oscillator.R")

N = 5 #number of chains
T_N = seq(0.2,1, length = 5) #temperature ladder

log_prior <- function(th){
  out = dmvnorm(th, mean = c(2,1), sigma = diag(1,nrow=2), log=T)
  #out = log(dgamma(th[1],16,8)) + log(dgamma(th[2],8,8))
  return(out)
}

log_likelihood <- function (th,y_obs,y0){
  n = length(y_obs)
  times = seq(0,30, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode23")[,3]
  
  out = 0
  
  #define the likelihood
  for(i in (1:n))
  {
    out <- out + dnorm(y_obs[i],mean = y_mod[i], sd = 0.5,log=T)
  }
  return(out)
}

log_target <- function(th,y_obs,y0,t_n){
  out = t_n*log_likelihood(th,y_obs,y0) + log_prior(th)
  return(out)
}

population_MCMC <- function(niter, burnin, th0, T_N ,Sig, y0, p_m,log_target)
{
  # define the vector that contains the output MCMC sample
  th <- NULL
  #alfa = numeric(niter)
  #alfa[1] = 1
  
  nacp = 0 # number of accepted moves
  # Start from th0
  for(i in 2:(niter))
  {
    for(j in 1:length(T_N)){
      p0 = runif(1,0,1)
      #### Choose between crossover and local move
      if(p0<p_m){
        #this is the local change
        delta = as.vector(rmvnorm(1, mean = th0[j], sig = Sig))
        lacp <- log_target(th = delta, y_obs = y_obs, y0 = y0, t_n=T_N[j])
        lacp <- lacp - log_target(th = th0[j], y_obs = y_obs, y0 = y0, t_n=T_N[j])
        
        lgu <- log(runif(1))  
        if(lgu < lacp)
        {
          th0[j] <- delta
          nacp = nacp + 1
        }
      }else{
        #this is the crossover operation
        
      }
    }
    
    
    ## propose delta: the mean is given by the previous value 
    ## N.B: delta = c(th,alfa_i), we initially use a discrete uniform for alfa_i
    
    th1 = as.vector(rmvnorm(1, mean = th0[1:2], sig = Sig))
    alfa1 =   sample(x=c(0.2,0.4,0.6,0.8,1), size=1, replace=TRUE, prob=c(1,0,0,0,0))              
    
    delta <- c(th1,alfa1)
    
    # First consider the accept/reject ratio
    # numerator
    lacp <- log_target(th = delta[1:2], y_obs = y_obs, y0 = y0, alfa = delta[3])
    # denominator
    lacp <- lacp - log_target(th = th0[1:2], y_obs = y_obs, y0 = y0, alfa = th0[3])
    
    #lacp <- min(0, lacp)  
    # Note: The proposal is symmetrical and therefore does not appear in the value of acceptance / rejection!
    
    lgu <- log(runif(1))  # lgu is negative
    ## if u < acp accept the move
    if(lgu < lacp)
    {
      th0 <- delta
      nacp = nacp + 1
    }
    
    #alfa[i] = th0[3]
    
    if(i>burnin & th0[3] == 1)
    {
      th=rbind(th,th0)
      #cat("# accepted move at iteration =", i, "\n")
    }
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter,"alfa=",th0[3], "th = ", th0[1:2], "\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}
