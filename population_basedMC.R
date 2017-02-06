################################
# Population based Monte Carlo #
################################

source("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana\\simple_oscillator.R")
#source("/home/mario/Scrivania/progetto_bayes/simple_oscillator.R")

if(parallel){
  library(foreach)
  library(doParallel)
  cores=detectCores()
  if (cores == 1){
    cat("cannot perform parallel on 1 cpu \n")
    parallel = FALSE
  }
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  if (cores==2){
    parallel = FALSE
  }
  
}

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

population_MCMC <- function(niter, burnin,thin ,th0, T_N ,Sig, y0, p_m,log_target, parallel)
{ 
  # th0 will be updated at each step, th will contail the output of interest (that is, when T_N = 1)
  th <- matrix(nrow= ceiling((niter-burnin)/thin), ncol=2)

  nacp = 0 # number of accepted moves

  for(i in 1:(niter))
  {
      for(j in 1:length(T_N)){
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
          
        }
      }
    
    # Try to exchange theta_l and theta_m where m = l+1 or m= l-1 if l=! 1 and l=! length(T_N)
    
    l = sample(x=(1:length(T_N)),size=1,prob = rep(1/length(T_N),length(T_N)))
    if(l>1 && l< length(T_N)){
      u = runif(1)
      m = l + 1*(u<=0.5) - 1*(u>0.5)
    }else if(l==1){
      m=2
    }else{m=length(T_N)-1}
    
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
    
    if(i>burnin & (i-burnin)%%thin==0){
      th[(i-burnin)/thin,] = th0[length(T_N),]
    }
    
    if(i%%1000==0) cat("*** Iteration number ", i,"/", niter, "\n")
  }
  cat("Acceptance rate =", nacp/niter, "\n")
  return(th)
}

parallel = FALSE


niter = 100000
burnin = 1000
thin = 10 
Sig = matrix(data = c(0.05, 0, 0, 0.05),nrow=2,ncol=2)

th0 = matrix( rep(c(1,0.5),length(T_N)),ncol=2, byrow=T)

th.post <- population_MCMC(niter = niter, burnin=burnin, thin = thin ,th0=th0, T_N=T_N ,Sig=Sig, y0=y0, p_m=1,log_target=log_target, parallel = parallel)
dim(th.post)
write.table(th.post, file = "output_pop_MCMC22.txt",row.names = F)
th.post.mc <- mcmc(th.post, start = burnin+ 1, end = niter, thin = thin)

x11()
plot(th.post.mc)
