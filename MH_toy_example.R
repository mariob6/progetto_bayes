#############################################
## Metropolis Hastings on the simple model ##
#############################################
library(mvtnorm)
library(deSolve)

# Model: dx/dt = k1/(36 + k2* y) - k3
#        dy/dt = k4*x - k5

#Define a function for the evaluation of the system

circ_oscillator <- function(time,state,params){
 k1 = params[1]
 k2 = params[2]
 k3 = params[3]
 k4 = params[4]
 k5 = params[5]
 
 dx = k1/(36 + k2* state[2]) - k3
 dy = k4*state[1] - k5
 return(list(c(dx,dy)))
}
# Set the "true" parameters

k1 = 72
k2 = 1
k3 = 2
k4 = 1
k5 = 1

t_end = 30

t_obs=seq(0,t_end, by=0.5)

#Set the "true" solutions evaluated in t_obs using ODE45

y0 = c(7,-10)
y_true = ode(y0,t_obs,circ_oscillator,c(k1,k2,k3,k4,k5),method = "ode45")[,2:3]


N = dim(y_true)[1]

noise = rmvnorm(N-1,mean=c(0,0),sigma=diag(0.5,2))
y_obs = y_true
y_obs[2:N,] = y_true[2:N,] + noise
y_obs = as.matrix(y_obs)


windows()
par(mfrow = c(1,2))
plot(t_obs,y_obs[,1])
lines(t_obs,y_true[,1],col="red")
plot(t_obs,y_obs[,2])
lines(t_obs,y_true[,2],col="red")

## We assume k1,k2,k5 known and assign a Gamma prior with parameters 1,1 to k3,k4

posterior <- function (th,y_obs,y0){
  n = dim(y_obs)[1]
  m = dim(y_obs)[2]
  times = seq(0,30, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode45")[,2:3]
  
  out = 0
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out + dmvnorm(y_obs[i,],mean = y_mod[i,], sigma = diag(0.5,2),log=T)
  }
  #Multiply times the prior
  out = out + log(dunif(th[1],0.1,5)) + log(dgamma(th[2],0.1,5))
  return(out)
}

#Plot the posterior on a grid

grid.k3 = seq(0.1,5,length=50)
grid.k4 = seq(0.1,5,length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,y0)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}

max(post_grid)
min(post_grid)

post_grid = post_grid - max(post_grid)

post_grid = exp(post_grid)

x11()
# Contour fun: Create a contour plot, or add contour lines to an existing plot
contour(grid.k3, grid.k4, post_grid, lwd = 3, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")

x11()
persp(grid.k3, grid.k4, post_grid, theta = 30)



######

# Toy Example, measuring only y

######


t_obs=seq(0,t_end, by=0.5)

#Set the "true" solutions evaluated in t_obs using ODE45

starting_point = c(7,-10)
y_true = ode(starting_point,t_obs,circ_oscillator,c(k1,k2,k3,k4,k5),method = "ode45")[,3]


N = length(y_true)

noise = rnorm(N-1,mean=0,sd=1)
y_obs = y_true
y_obs[2:N] = y_true[2:N] + noise



windows()
plot(t_obs,y_obs)
lines(t_obs,y_true,col="red")

## We assume k1,k2,k5 known and assign a Gamma prior with parameters 1,1 to k3,k4

posterior <- function (th,y_obs,y0){
  n = length(y_obs)
  times = seq(0,60, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode45")[,3]
  
  out = 0
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out + dnorm(y_obs[i],mean = y_mod[i], sd = 0.5,log=T)
  }

  #Multiply times the prior
  out = out + log(dgamma(th[1],1,1)) + log(dgamma(th[2],1,1))
  return(out)
}

#Plot the posterior on a grid


grid.k3 = seq(0.01,4,length=50)
grid.k4 = seq(0.01,4,length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,starting_point)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}

max(post_grid)
min(post_grid)

post_grid = post_grid - max(post_grid)

post_grid = exp(post_grid)


x11()
# Contour fun: Create a contour plot, or add contour lines to an existing plot
contour(grid.k3, grid.k4, post_grid, lwd = 3, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")

x11()
persp(grid.k3, grid.k4, post_grid, theta = 30)


##########################
## Metropolis Hastings  ##
##########################

# Define the proposal as gaussian with suitable mean and variance

?optim
ottimo <- optim(par = c(mean(y_obs), log(sd(y_obs))), 
                fn = posterior, y_obs = y_obs,
                y0 = y0, 
                hessian = T, control = list(fnscale=-1))

map <- ottimo$par 
map
# ...but also an "estimate" of the input Hessian matrix around the maximum
# Through this estimate I define the matrix of variance and covariance of the proposal
Sig <- - solve(ottimo$hessian)  
# Estimation of the covariance
Sig

gamma = 1000

Sig = gamma * Sig

metro_hastings <- function(niter, burnin, thin, th0, Sig,y0)
{
  
  # define the vector that contains the output MCMC sample
  th <- matrix(nrow= ceiling((niter-burnin)/thin), ncol = 2)
  
  nacp = 0 # number of accepted moves
  # Start from th0
  for(i in 1:(niter))
  {
    ## propose delta: the mean is given by the previous value
    delta <- as.vector(rmvnorm(1, mean = th0, sig = Sig))
    
    # First consider the log of accept/reject ratio
    # numerator
    lacp <- posterior(th = delta, y_obs = y_obs, y0 = y0)
    # denominator
    lacp <- lacp - posterior(th = th0, y_obs = y_obs, y0 = y0)

    if(is.nan(lacp))
      lacp = 0
    lacp <- min(0, lacp)  
    # Note: The proposal is symmetrical and therefore does not appear in the value of acceptance / rejection!
    
    lgu <- log(runif(1))  # lgu is negative
    ## if u < acp accept the move
    if(lgu < lacp)
    {
      th0 <- delta
      nacp = nacp + 1
    }
  
    
    if(i>burnin & (i-burnin)%%thin==0)
    {
      th[(i-burnin)/thin,] = th0
    }
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter,"th=",tail(th,1) ,"\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}

niter = 60000
burnin = 10000
thin = 10
# At the end the chain will contain
# (niter-burnin)/thin = 5000 observations

#th0 = map
th0 = c(0.5,1)
th.post <- metro_hastings(niter = niter, burnin = burnin, thin = thin, th0 = th0, Sig = Sig, y0 = starting_point)

dim(th.post)

