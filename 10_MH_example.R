##############################################################
## Launching 10 MH chains in parallel on Girolami's example ##
##############################################################
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

starting_point = c(7,-10)
y0=starting_point
y_true = ode(starting_point,t_obs,circ_oscillator,c(k1,k2,k3,k4,k5),method = "ode45")[,3]


N = length(y_true)

noise = rnorm(N-1,mean=0,sd=1)
y_obs = y_true
y_obs[2:N] = y_true[2:N] + noise



windows()
plot(t_obs,y_obs)
lines(t_obs,y_true,col="red")

## We assume k1,k2,k5,sigma(=0.5) known and assign a Gamma prior with parameters 16,8 to k3
##                                                                               8,8 to k4

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
  out = out + log(dgamma(th[1],shape = 8, rate = 4)) + log(dgamma(th[2],shape = 8, rate = 4))
  return(out)
}

#Plot the posterior on a grid


grid.k3 = seq(1.6,5.5,by=0.1)
grid.k4 = seq(0.1,5,by=0.1)

post_grid = matrix(nrow=length(grid.k3),ncol=length(grid.k4))

for(k in 1:length(grid.k3)){
  for(j in 1:length(grid.k4)){
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

ottimo <- optim(par = c(1,1), 
                fn = posterior, y_obs = y_obs,
                y0 = y0, 
                hessian = T, control = list(fnscale=-1))


Sig <- - solve(ottimo$hessian)  

Sig

gamma = 1

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


niter = 16000
burnin = 10000
thin = 10


###### Start 10 MC from different initial points #####
th0_vec = matrix(nrow = 10, ncol = 2)
th0_vec[,1] = runif(10,1.6,5)
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



#finalMatrix=read.csv(file="output_10_chains_par.csv")
#finalMatrix=finalMatrix[,2:21]

colmax = apply(finalMatrix,2,max)
colmin = apply(finalMatrix,2,min)

range_k3 = c(min(colmin[1],colmin[3],colmin[5],colmin[7],colmin[9],colmin[11],colmin[13],colmin[15],colmin[17],colmin[19]),
             max(colmax[1],colmax[3],colmax[5],colmax[7],colmax[9],colmax[11],colmax[13],colmax[15],colmax[17],colmax[19]))


range_k4 = c(min(colmin[2],colmin[4],colmin[6],colmin[8],colmin[10],colmin[12],colmin[14],colmin[16],colmin[18],colmin[20]),
             max(colmax[2],colmax[4],colmax[6],colmax[8],colmax[10],colmax[12],colmax[15],colmax[16],colmax[18],colmax[20]))


grid.k3 = seq(-2,2,length=50)
grid.k4 = seq(-2,2,length=50)

post_grid = matrix(nrow=50,ncol=50)

for(k in 1:50){
  for(j in 1:50){
    post_grid[k,j] = posterior(c(grid.k3[k],grid.k4[j]),y_obs,starting_point)
  }
  if(k%%10 == 0)
    cat("k=",k, "\n")
}

graphics.off()

contour(grid.k3, grid.k4, post_grid, lwd = 1, labcex = 1.1,zlim=c(-2e05,0), col = "blue", main = "Posterior density", 
        xlab = "k3", ylab = "k4")

for(i in seq(1,17,by=2)){
  points(finalMatrix[,i],finalMatrix[,(i+1)])
}

this.post <- finalMatrix[,15:16]
this.post.mcmc <- mcmc(this.post,start=burnin+1,end=niter,thin=thin)
plot(this.post.mcmc)



log_prior <- function(th){
  #out = dmvnorm(th, mean = c(2,1), sigma = diag(0.1,nrow=2), log=T)
  out = log(dgamma(th[1],1,1)) + log(dgamma(th[2],1,1))
  return(out)
}

# A function defining the likelihood

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

log_target <- function(th,y_obs,y0,alfa){
  out = alfa*log_likelihood(th,y_obs,y0) + log_prior(th)
  return(out)
}



grid_k3 = seq(1.3,5.1,by=0.1)
grid_k4 = seq(0.1,4,by=0.1)

alfa = 1

plot_grid = matrix(nrow=length(grid_k3), ncol=length(grid_k4))
for(i in (1:length(grid_k3))){
  for(j in (1:length(grid_k4)))
    plot_grid[i,j] = log_target(th=c(grid_k3[i],grid_k4[j]), y_obs = y_obs, y0=y0, alfa=alfa)
}

nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(plot_grid, nbcol)
p<-persp3d(grid_k3, grid_k4, plot_grid, theta = 50,phi=25, col=color[zcol],
           ticktype="detailed", zlim = c(-2e05,1),xlab = "k3", ylab= "k4", zlab="",axes=TRUE)
rgl.snapshot("posterior.png")


image(grid_k3, grid_k4, plot_grid)
contour(grid_k3, grid_k4, plot_grid, lwd = 1, labcex = 1.1,xlim=c(1.2,5.1),zlim=c(-2e05,0), col = "black", 
        xlab = "k3", ylab = "k4",drawlabels=FALSE,add=TRUE)


for(i in seq(1,17,by=2)){
  points(finalMatrix[,i],finalMatrix[,(i+1)],cex=1.5)
}
for(i in seq(1,17,by=2)){
  lines(finalMatrix[,i],finalMatrix[,(i+1)])
}

