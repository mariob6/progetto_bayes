#############################################
## Metropolis Hastings on the simple model ##
#############################################

setwd("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana")

library(mvtnorm)
library(deSolve)
library(coda)

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

## We assume k1,k2,k5 known and assign a bivariate normal with mean (2,1) and sigma = diag(0.1,2)

posterior <- function (th,y_obs,y0){
  n = dim(y_obs)[1]
  m = dim(y_obs)[2]
  times = seq(0,30, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode45")[,2:3]
  
  diff = max(abs(y_mod - y_obs))
  
  #cat("diff" = diff, "\n")
  
  out = 0
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out + dmvnorm(y_obs[i,],mean = y_mod[i,], sigma = diag(0.5,2),log=T)
  }
  #Multiply times the prior
  out = out + dmvnorm(th, mean = c(2,1), sigma = diag(0.1,nrow=2), log=T)
  return(out) + 132
}

#Plot the posterior on a grid

grid.k3 = seq(1.85,2.25,length=50)
grid.k4 = seq(0.85,1.25,length=50)

post_grid = matrix(nrow=50,ncol=50)

plot(y_obs[,2],col="red",type="l")

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,y0)
    #th = c(grid.k3[k],grid.k4[j])
    #lines(ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode45")[,3],col="blue")
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}


x11()
# Contour fun: Create a contour plot, or add contour lines to an existing plot
contour(grid.k3, grid.k4, post_grid, lwd = 3, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")

x11()
persp(grid.k3, grid.k4, post_grid, theta = 30)

max(post_grid)

#### Different Grid

grid.k3 = seq(1.5,2.5,length=50)
grid.k4 = seq(0.5,1.5,length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,y0)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}

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

## We assume k1,k2,k5,sigma(=0.5) known and assign a Gamma prior with parameters 1,1 to k3,k4

posterior <- function (th,y_obs,y0){
  n = length(y_obs)
  times = seq(0,60, by=0.5)
  
  y_mod = ode(y0,times,circ_oscillator,c(72,1,th[1],th[2],1),method = "ode23")[,3]
  
  out = 0
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out + dnorm(y_obs[i],mean = y_mod[i], sd = 0.5,log=T)
  }

  #Multiply times the prior
  out = out + dmvnorm(th, mean = c(2,1), sigma = diag(0.1,nrow=2), log=T)
  return(out)
}

#Plot the posterior on a grid


grid.k3 = seq(1,4,length=50)
grid.k4 = seq(1,4,length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,starting_point)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}


x11()
# Contour fun: Create a contour plot, or add contour lines to an existing plot
contour(grid.k3, grid.k4, post_grid, lwd = 3, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")

x11()
persp(grid.k3, grid.k4, post_grid, theta = 30)


##########################
## Metropolis Hastings  ##
##########################

# Define the proposal as gaussian with suitable variance

ottimo <- optim(par = c(mean(y_obs), log(sd(y_obs))), 
                fn = posterior, y_obs = y_obs,
                y0 = y0, 
                hessian = T, control = list(fnscale=-1))


Sig <- - solve(ottimo$hessian)  

Sig

gamma = 10

Sig = gamma * Sig

#Sig = diag(c(0.5,0.5))

metro_hastings <- function(niter, burnin, thin, th0, Sig,y0)
{
  
  # define the vector that contains the output MCMC sample
  th <- matrix(nrow= ceiling(burnin/500 + (niter-burnin)/thin), ncol = 2)
  
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

    lacp <- min(0, lacp)  
    # Note: The proposal is symmetrical and therefore does not appear in the value of acceptance / rejection!
    
    lgu <- log(runif(1))  # lgu is negative
    ## if u < acp accept the move
    if(lgu < lacp)
    {
      th0 <- delta
      nacp = nacp + 1
    }
    
    if(i<burnin & (i-burnin)%%500==0){
      th[i/500,] = th0
    }

    if(i>burnin & (i-burnin)%%thin==0)
    {
      th[burnin/500 + (i-burnin)/thin,] = th0
    }
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter,"alfa=",th0[3], "th = ", th0, "\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}

niter = 16000
burnin = 0
thin = 1



th0 = c(2.5,2)
th.post <- metro_hastings(niter = niter, burnin = burnin, thin = thin, th0 = th0, Sig = Sig, y0 = starting_point)

dim(th.post)

th.post.mc <- mcmc(th.post, start = burnin + 1, end = niter, thin = thin)

plot(th.post.mc)

geweke.plot(th.post.mc, frac1 = 0.1, frac2 = 0.5, nbins = 20 )

#th.post = th
plot(th.post)

grid.k3 = seq(min(th.post[,1]),max(th.post[,1]),length=50)
grid.k4 = seq(min(th.post[,2]),max(th.post[,2]),length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,starting_point)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}


contour(grid.k3, grid.k4, post_grid, lwd = 1, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")
points(th.post)


###### Start 10 MC from different initial points #####
th0_vec = matrix(nrow = 10, ncol = 2)
th0_vec[,1] = runif(10,0.5,5)
th0_vec[,2] = runif(10,0.5,5)

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

finalMatrix <- foreach(i=1:10, .packages = c("mvtnorm","deSolve"), .combine=cbind) %dopar% {
  tempMatrix =  metro_hastings(niter = niter, burnin = burnin, thin = thin, th0 = th0_vec[i,], Sig = Sig, y0 = starting_point)
  tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}

#finalMatrix<-NULL
#for(i in 1:10){
#  tempMatrix =  metro_hastings(niter = niter, burnin = burnin, thin = thin, th0 = th0_vec[i,], Sig = Sig, y0 = starting_point)
#  finalMatrix = cbind(finalMatrix,tempMatrix)
#}                                        
stopCluster(cl)
write.csv(finalMatrix, file="output_10_chains_par.csv")



finalMatrix=read.csv(file="output_10_chains_par.csv")
finalMatrix=finalMatrix[,2:21]

colmax = apply(finalMatrix,2,max)
colmin = apply(finalMatrix,2,min)

range_k3 = c(min(colmin[1],colmin[3],colmin[5],colmin[7],colmin[9],colmin[11],colmin[13],colmin[15],colmin[17],colmin[19]),
             max(colmax[1],colmax[3],colmax[5],colmax[7],colmax[9],colmax[11],colmax[13],colmax[15],colmax[17],colmax[19]))


range_k4 = range_k3;

grid.k3 = seq(range_k3[1],range_k3[2],length=50)
grid.k4 = seq(range_k4[1],range_k4[2],length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,starting_point)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}

graphics.off()

contour(grid.k3, grid.k4, post_grid, lwd = 1, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")

for(i in seq(1,19,by=2)){
  points(finalMatrix[,i],finalMatrix[,i+1])
}







###################
# Try the MH sampler on a simple case: y' = 3*0.219 t^2 + 2* 0.5287*t_obs-0.805 
#                                      y  = k1*t^3 + k2*t^2 - k3*t + k4

t0 = -2
t_end = 2
N = 20
t_obs = seq(t0,t_end,length=N)

#Set the true parameters

k1 = 0.219
k2 = 0.5287
k3 = -0.805
k4 = 0

y_true = k1*t_obs^3 + k2*t_obs^2 + k3*t_obs + k4
y_obs= k1*t_obs^3 + k2*t_obs^2 + k3*t_obs +k4 + rnorm(N,0,0.3)


windows()
plot(t_obs,y_obs)
lines(t_obs,y_true,col="red")


## We assume k1,k4,sigma(=0.5) known and assign a normal prior with parameters 0, 5 to k2,k3 

posterior <- function (th,y_obs,y0){
  n = length(y_obs)
  times = seq(-2,2, length = n)
  
  y_mod = 0.219*times^3 + th[1]*times^2 + th[2]*times + 0 
  
  out = 0
  
  #define the likelihood
  for(i in 1:n)
  {
    out <- out + dnorm(y_obs[i],mean = y_mod[i], sd = 0.5,log=T)
  }
  
  #Multiply times the prior
  out = out + dnorm(th[1], mean = 0.5, sd = 1, log = T) + dnorm(th[2], mean = -0.5, sd = 5, log = T)
  return(out + 10.87)
}

grid.k2 = seq(-2,2,length=50)
grid.k3 = seq(-2,2,length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k2[k],grid.k3[j]),y_obs,starting_point)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}


x11()
# Contour fun: Create a contour plot, or add contour lines to an existing plot
contour(grid.k2, grid.k3, post_grid, lwd = 3, labcex = 1.1, col = "blue", main = "Posterior density", 
        xlab = "k2", ylab = "k3")

x11()
persp(grid.k2, grid.k3, post_grid, theta = 30)


# Define the proposal as gaussian with suitable variance

ottimo <- optim(par = c(mean(y_obs), log(sd(y_obs))), 
                fn = posterior, y_obs = y_obs,
                y0 = y0, 
                hessian = T, control = list(fnscale=-1))


Sig <- - solve(ottimo$hessian)  

Sig

gamma = 100

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
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter ,"\n")
    
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
th0 = c(0.5,-0.7)
th.post <- metro_hastings(niter = niter, burnin = burnin, thin = thin, th0 = th0, Sig = Sig, y0 = starting_point)

dim(th.post)


th.post.mc <- mcmc(th.post, start = burnin + 1, end = niter, thin = thin)
plot(th.post.mc)
geweke.plot(th.post.mc, frac1 = 0.1, frac2 = 0.5, nbins = 20 )