################################
# Population based Monte Carlo #
################################
setwd("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana")
source("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana\\simple_oscillator.R")
#source("/home/mario/Scrivania/progetto_bayes/simple_oscillator.R")
#require(compiler)
enableJIT(3)
N = 50 #number of chains
T_N = numeric(N) #temperature ladder

for(n in 1:(N-1)){
  T_N[n]=n^5/(1.3*N)^5
}

T_N[N]=1

log_prior <- function(th){
  #out = dmvnorm(th, mean = c(2,1), sigma = diag(1,nrow=2), log=T)
  out = log(dgamma(th[1],1,1)) + log(dgamma(th[2],1,1))
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

population_MCMC <- function(niter, burnin,thin ,th0, T_N ,Sig, y0, p_m,log_target, parallel)
{ 
  # th0 will be updated at each step, th will contail the output of interest (that is, when T_N = 1)
  th <- matrix(nrow= ceiling((niter-burnin)/thin), ncol=2)
  N = length(T_N)
  nacp = 0 # number of accepted moves

  for(i in 1:(niter))
  {
      for(j in 1:N){
        p0 = runif(1,0,1)
        #### Choose between crossover and local move
        if(p0 <= p_m){
          #this is the local change
          delta = as.vector(rmvnorm(1, mean = th0[j,], sig = Sig))
          lacp <- log_target(th = delta, y_obs = y_obs, y0 = y0, t_n=T_N[j])
          lacp <- lacp - log_target(th = th0[j,], y_obs = y_obs, y0 = y0, t_n=T_N[j])
          #cat(lacp,"\n")
          lgu <- log(runif(1))  
          if(!(is.na(lacp)) & lgu < lacp)
          {
            th0[j,] <- delta
            nacp = nacp + 1
          }
        }else{
          #this is the crossover operation
          lm = sample(N,2)
          l = lm[1]
          m = lm[2]
          c=sample(2,1)
          #crossover
          th1 = th0[l,]
          th1[c:2] = th0[m, c:2]
          th2 = th0[m,]
          th2[c:2] = th0[l,c:2]
          
          #acceptance-rejection
            lacp = log_likelihood(th=th1,y_obs=y_obs,y0=y0)*T_N[l] + log_likelihood(th=th2,y_obs=y_obs,y0=y0)*T_N[m]
            lacp = lacp -  (log_likelihood(th=th0[l,],y_obs=y_obs,y0=y0)*T_N[l] + log_likelihood(th=th0[m,],y_obs=y_obs,y0=y0)*T_N[m])
            lgu <- log(runif(1))  
            if(lgu < lacp)
            {
              th0[l,] <- th1
              th0[m,] <- th2
              nacp = nacp + 1
            }
        }
      }
    
    # Try to exchange theta_l and theta_m where m = l+1 or m= l-1 if l=! 1 and l=! length(T_N)
    for(k in 1:N){
      l = sample(x=(1:N),size=1,prob = rep(1/N,N))
      if(l>1 && l< length(T_N)){
        u = runif(1)
        m = l + 1*(u<=0.5) - 1*(u>0.5)
      }else if(l==1){
        m=2
      }else{m=N-1}
      
      lacp = log_likelihood(th=th0[m,],y_obs=y_obs,y0=y0)*T_N[l] + log_likelihood(th=th0[l,],y_obs=y_obs,y0=y0)*T_N[m]
      lacp = lacp -  (log_likelihood(th=th0[m,],y_obs=y_obs,y0=y0)*T_N[m] + log_likelihood(th=th0[l,],y_obs=y_obs,y0=y0)*T_N[l])
      lgu <- log(runif(1))  
      if(lgu < lacp)
      {
        th_aux = th0[l,]
        th0[l,] <- th0[m,]
        th0[m,] <- th_aux
        nacp = nacp + 1
      }
    }
    if(i>burnin & (i-burnin)%%thin==0){
      th[(i-burnin)/thin,] = th0[N,]
    }
    
    if(i%%1==0) cat("*** Iteration number ", i,"/", niter, "\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}

parallel = FALSE


niter = 2000
burnin = 0
thin = 1 
Sig = matrix(data = c(0.05, 0.01, 0.01, 0.05),nrow=2,ncol=2)

th0 = matrix( c(runif(N,0,5),runif(N,0,5)),ncol=2, byrow=T)
th0[N,] = c(3.5,3.5)

th.post <- population_MCMC(niter = niter, burnin=burnin, thin = thin ,th0=th0, T_N=T_N ,Sig=Sig, y0=y0, p_m=0.75,log_target=log_target, parallel = parallel)
dim(th.post)
write.table(th.post, file = "output_pop_MCMC_1602v3.txt",row.names = F)
#th.post<-read.table(file="output_pop_MCMC1402v2.txt",header=T)

# Plotting the markov chain in the state space

grid_k3 = seq(1,4.5,by=0.1)
grid_k4 = seq(0.5,3.5,by=0.1)
t_n = 1

plot_grid = matrix(nrow=length(grid_k3), ncol=length(grid_k4))
for(i in (1:length(grid_k3))){
  for(j in (1:length(grid_k4)))
    plot_grid[i,j] = log_target(th=c(grid_k3[i],grid_k4[j]), y_obs = y_obs, y0=y0, t_n = t_n)
}

persp(grid_k3,grid_k4,plot_grid,zlim=c(-1e05,1))
points(th.post,pch=16)
contour(grid_k3,grid_k4,plot_grid,zlim=c(-2e05,0))


th.post.mc <- mcmc(th.post[1000:2000 ,], start = 1000+1, end = niter, thin = thin)

x11()
plot(th.post.mc)

