#########################
# Product Space Sampler #
#########################
source("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana\\simple_oscillator.R")
# Define the taget density: p(th|y,alfa) = prod (p(th|y,alfa_i))
# where p(th|y,alfa_i) is proportional to p(y|th)^alfa_i*p(th)

L = 11
  
alfa = numeric(L)

alfa = seq(0,1,length=L)

# A function defining the prior

prior <- function(th){
  return(dgamma(th[1],16,8)*dgamma(th[2],8,8))
}

log_prior <- function(th){
  out = log(dgamma(th[1],16,8)) + log(dgamma(th[2],8,8))
  return(out)
}

# A function defining the likelihood

likelihood <- function (th,y_obs,y0){
  n = dim(y_obs)[1]
  m = dim(y_obs)[2]
  times = seq(0,30, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode45")[,3]
  
  out = 1
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out * dnorm(y_obs[i,],mean = y_mod[i,], sigma = 0.5 ,log=F)
  }
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

# A function defining the target density

target <- function(th,y_obs,y0,alfa,likelihood,prior){
    return(likelihood(th,y_obs,y0)^alfa * prior(th) )
}

log_target <- function(th,y_obs,y0,alfa,log_likelihood,log_prior){
  out = alfa*log_likelihood(th,y_obs,y0) + log_prior(th)
  return(out)
}
###########
# Plot the log_target for different choices of alfa

grid_k3 = seq(1.3,6,by=0.1)
grid_k4 = seq(0.3,5,by=0.1)

alfa = 0;

plot_grid = matrix(nrow=length(grid_k3), ncol=length(grid_k4))
for(i in (1:length(grid_k3))){
  for(j in (1:length(grid_k4)))
    plot_grid[i,j] = log_target(th=c(grid_k3[i],grid_k4[j]), y_obs = y_obs, y0=y0, alfa=alfa,log_likelihood = log_likelihood,log_prior = log_prior)
}

x11()
persp(grid_k3, grid_k4, plot_grid, theta = 30)

alfa = 0.4

plot_grid = matrix(nrow=length(grid_k3), ncol=length(grid_k4))
for(i in (1:length(grid_k3))){
  for(j in (1:length(grid_k4)))
    plot_grid[i,j] = log_target(th=c(grid_k3[i],grid_k4[j]), y_obs = y_obs, y0=y0, alfa=alfa,log_likelihood = log_likelihood,log_prior = log_prior)
}

x11()
persp(grid_k3, grid_k4, plot_grid, theta = 30)

alfa = 0.75

plot_grid = matrix(nrow=length(grid_k3), ncol=length(grid_k4))
for(i in (1:length(grid_k3))){
  for(j in (1:length(grid_k4)))
    plot_grid[i,j] = log_target(th=c(grid_k3[i],grid_k4[j]), y_obs = y_obs, y0=y0, alfa=alfa,log_likelihood = log_likelihood,log_prior = log_prior)
}

x11()
persp(grid_k3, grid_k4, plot_grid, theta = 30)

alfa = 1

plot_grid = matrix(nrow=length(grid_k3), ncol=length(grid_k4))
for(i in (1:length(grid_k3))){
  for(j in (1:length(grid_k4)))
    plot_grid[i,j] = log_target(th=c(grid_k3[i],grid_k4[j]), y_obs = y_obs, y0=y0, alfa=alfa,log_likelihood = log_likelihood,log_prior = log_prior)
}

x11()
persp(grid_k3, grid_k4, plot_grid, theta = 30)


############

target_density_sampler <- function(niter, burnin, th0, Sig,y0,alfa,log_target,log_likelihood, log_prior)
{
  # define the vector that contains the output MCMC sample
  th <- NULL
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
    lacp <- log_target(th = delta[1:2], y_obs = y_obs, y0 = y0, alfa = delta[3], log_likelihood, log_prior)
    # denominator
    lacp <- lacp - log_target(th = th0[1:2], y_obs = y_obs, y0 = y0, alfa = th0[3], log_likelihood, log_prior)

    #lacp <- min(0, lacp)  
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
      th=rbind(th,th0)
    }
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter,"th=",tail(th,1) ,"\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}

Sig = matrix(data = c(1.788331e-05, -0.0000462595, -4.625950e-05, 0.0001702542),nrow=2,ncol=2)
niter=16000
burnin=10000

th0 = c(1.5,1)
alfa = 0.5
th.post <- target_density_sampler(niter = niter, burnin = burnin, th0 = th0, Sig = Sig,y0 = y0,alfa = alfa,log_likelihood, log_prior)
dim(th.post)

th.post.mc <- mcmc(th.post, start = 0, end = niter, thin = 1)

plot(th.post.mc)

geweke.plot(th.post.mc, frac1 = 0.1, frac2 = 0.5, nbins = 20 )