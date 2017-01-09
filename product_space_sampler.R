#########################
# Product Space Sampler #
#########################

# Define the taget density: p(th|y,alfa) = prod (p(th|y,alfa_i))
# where p(th|y,alfa_i) is proportional to p(y|th)^alfa_i*p(th)

L = 11
  
alfa = numeric(L)

alfa = seq(0,1,length=L)

# A function defining the prior

prior <- function(th){
  return(dgamma(th[1],1,1)*dgamma(th[2],1,1))
}

# A function defining the likelihood

likelihood <- function (th,y_obs,y0){
  n = dim(y_obs)[1]
  m = dim(y_obs)[2]
  times = seq(0,30, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode45")[,2:3]
  
  out = 1
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out * dmvnorm(y_obs[i,],mean = y_mod[i,], sigma = diag(0.5,2),log=F)
  }
  return(out)
}

# A function defining the target density

density <- function(th,y_obs,y0,alfa,likelihood,prior){
    return(likelihood(th,y_obs,y0)^alfa * prior(th) )
}





############

target_density_sampler <- function(niter, burnin, th0, Sig,y0,alfa,likelihood, prior)
{
  # define the vector that contains the output MCMC sample
  th <- matrix(nrow= ceiling((niter-burnin)/thin), ncol = 2)
  L = length(alfa)
  
  nacp = 0 # number of accepted moves
  # Start from th0
  for(i in 1:(niter))
  {
    ## propose delta: the mean is given by the previous value 
    ## N.B: delta = c(th,alfa_i), we initially use a discrete uniform for alfa_i
    delta <- c(as.vector(rmvnorm(1, mean = th0, sig = Sig)),sample(alfa,1))
    
    # First consider the accept/reject ratio
    # numerator
    acp <- density(th = delta[1:2], y_obs = y_obs, y0 = y0, alfa = delta[3], likelihood, prior)
    # denominator
    acp <- acp / density(th = th0[1:2], y_obs = y_obs, y0 = y0, alfa = th0[3], likelihood, prior)

    lacp <- min(0, lacp)  
    # Note: The proposal is symmetrical and therefore does not appear in the value of acceptance / rejection!
    
    lgu <- log(runif(1))  # lgu is negative
    ## if u < acp accept the move
    if(lgu < lacp)
    {
      th0 <- delta
      nacp = nacp + 1
    }
    
    
    if(i>burnin & th0[3] == 1)
    {
      # MODIFY HERE: insert as last element th0[1:2]
      th[(i-burnin)/thin,] = th0
    }
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter,"th=",tail(th,1) ,"\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}
}
