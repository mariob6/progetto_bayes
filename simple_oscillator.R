#### Define the simple oscillator model

library(mvtnorm)
library(deSolve)
library(coda)
#library(rgl)

set.seed(09012017)
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

y0 = c(7,-10)
y_true = ode(y0,t_obs,circ_oscillator,c(k1,k2,k3,k4,k5),method = "ode45")[,3]


N = length(y_true)

noise = rnorm(N-1,mean=0,sd=0.45)
y_obs = y_true
y_obs[2:N] = y_true[2:N] + noise