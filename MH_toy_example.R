#############################################
## Metropolis Hastings on the simple model ##
#############################################
library(mvtnorm)

# Model: dx/dt = k1/(36 + k2* y) - k3
#        dy/dt = k4*x - k5

# Set the "true" parameters

k1 = 72
k2 = 1
k3 = 2
k4 = 1
k5 = 1

t_obs=seq(0,60, by=0.5)

#Set the "true" solutions evaluated in t_obs using ODE45

y_true = read.table('toy_example_sol.csv',sep=",")

N = dim(yy)[1]

noise = rmvnorm(N,mean=c(0,0),sigma=diag(0.5,2))
y_obs = y_true + noise

windows()
par(mfrow = c(1,2))
plot(t_obs,y_obs[,1])
lines(t_obs,y_true[,1],col="red")
plot(t_obs,y_obs[,2])
lines(t_obs,y_true[,2],col="red")

## We assume k1,k2,k5 known and assign a Gamma prior with parameters 1,1 to k3,k4

posterior <- function (th,y){
  n = dim(y)[1]
  m = dim(y)[2]
  
  out = 0
  #define the likelihood
  
  for(i in 1:n)
  {
    out <- out + lambda - log(1 + exp(-2*lambda)*(y[i]-mu)^2)
  }
  
}





