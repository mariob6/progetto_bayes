#######################################################
# Linear Example (Section 4 of the paper by Girolami) #

setwd("C:\\Users\\mario\\Desktop\\UNIVERSITA'\\Progetti\\bayesiana")

library("MCMCpack")
library("MASS")
set.seed(21122016)
N=25
t_obs=runif(25,-2,2)
t_obs=sort(t_obs)
y_obs=0.219*t_obs^3 + 0.5287*t_obs^2-0.805*t_obs + rnorm(N,0,1)

t_true=seq(-2.02,2.02,by=0.1)
y_true=0.219*t_true^3 + 0.5287*t_true^2-0.805*t_true
windows()
plot(t_obs,y_obs,pch=16)
lines(t_true,y_true,col="red")

K=10;

# Set Uniform prior over the models

prior_models=rep(1/K,K)

# Bayesian Model-Averaged mean prediction 
# E[y*|t*,y_obs,t_obs] = mu* = sum_over_all_models( phi^T * mu_k* posterior(M_k|y_obs,t_obs))
# phi^T = kx1 vector with the k basis of model M_k evaluated at t*
# mu_k = Sigma_k^(-1)* PHI^T *y_obs
#         PHI = Nxk matrix with basis responses at points t_obs
#         Sigma_k = I_k + PHI^T*PHI

# Posterior(M_k|y_obs,t_obs) = (p(y_obs|t_obs,M_k)*prior_models(k))/sum_k(p(y_obs|t_obs,M_k)*prior_models(k))

# p(y_obs|t_obs,M_k) = integrated likelihood

integrated_likelihood <- function(N,Sigma,PHI,y_obs)
{
  out = gamma(1+N/2)/sqrt(pi^N*det(Sigma)) * (1 + 0.5* t(y_obs)%*%(diag(N)-PHI%*%solve(Sigma)%*%t(PHI))%*%y_obs)^(-1-N/2)
  return(out)
}


# Setting the parameters for all the models, this implementation might be highly inefficient when
# the number of models is high

phi = matrix(nrow=N,ncol=1)
phi[,1] = t_obs;
PHI_1 = phi;
Sigma_1 = 1 + t(phi)%*%phi
mu_1 = solve(Sigma_1)%*%t(PHI_1)%*%y_obs


for(k in 2:K){
  #for different basis functions, define the following line accordingly
    phi = cbind(phi,t_obs^k)
  #create PHI_k
    assign(paste("PHI",k,sep="_"),phi)
  #create Sigma_k
    Sigma = diag(k) + t(phi)%*%phi
    assign(paste("Sigma",k,sep="_"), Sigma)
}

#setting the normalizing constant for the posterior

posterior_denominator=0;

for(k in 1:K){
  posterior_denominator = posterior_denominator + integrated_likelihood(N,get(paste("Sigma",k,sep="_")),get(paste("PHI",k,sep="_")),y_obs)
}

posterior_models=numeric(K);

for(k in 1:K){
  posterior_models[k] = integrated_likelihood(N,get(paste("Sigma",k,sep="_")),get(paste("PHI",k,sep="_")),y_obs)/posterior_denominator
}


# Bayesian Model Averaged Mean Predictior
# make predictions, evaluate mu_star

mu_star=0
L = length(t_true)
y_pred = numeric(L)

for(i in 1:L){
  t_star=t_true[i]
  for(k in 1:K){
    phi = numeric(k)
    for(j in 1:k){
      phi[j] = t_star^j
    }
    mu = solve(get(paste("Sigma",k,sep="_")))%*%t(get(paste("PHI",k,sep="_")))%*%y_obs
    y_pred[i] = y_pred[i] + t(phi)%*%mu%*%posterior_models[k]
  }
}

windows()
plot(t_obs,y_obs)
lines(t_true,y_true,col="red")
lines(t_true,y_pred,col="dark green")
windows()
barplot(posterior_models)


##################################
## Caso lineare, base di Fourier #
##################################

## se y_obs e y_true sono combinazioni lineari delle basi di fourier, nessun problema.
## se invece y_obs e y_true sono scelti come nel caso precedente, il risultato non è entusiasmante
## - aumentare grado del modello?


set.seed(21122016)
N=20
t_obs=runif(N,0,1)
t_obs=sort(t_obs)
#y_obs=0.219*t_obs^3 + 0.5287*t_obs^2-0.805*t_obs + rnorm(N,0,0.1)

y_obs = sin(2*pi*t_obs) + 4*cos(4*pi*t_obs) + rnorm(N,0,0.1)

t_true=seq(-.02,1.02,by=0.01)
#y_true=0.219*t_true^3 + 0.5287*t_true^2-0.805*t_true
y_true = sin(2*pi*t_true) + 4*cos(4*pi*t_true)
#windows()
plot(t_obs,y_obs,pch=16)
lines(t_true,y_true,col="red")

K=19;

fourier_basis <-function(index,seno){
  #params: -index: integer cointaining the index of the fourier basis to be evaluated
  #        -seno: boolean = 1 if we want to evaluate sin, 0 for cos
  function(x) {
    sqrt(2)*sin(2*pi*index*x)*seno + sqrt(2)*cos(2*pi*index*x)*(1-seno)
  }
 
}

prior_models=rep(1/K,K)

integrated_likelihood <- function(N,Sigma,PHI,y_obs)
{
  out = gamma(1+N/2)/sqrt(pi^N*det(Sigma)) * (1 + 0.5* t(y_obs)%*%(diag(N)-PHI%*%solve(Sigma)%*%t(PHI))%*%y_obs)^(-1-N/2)
  return(out)
}


phi = matrix(nrow=N,ncol=K)

phi[,1] = fourier_basis(0,0)(t_obs);
phi[,2] = fourier_basis(1,0)(t_obs);
phi[,3] = fourier_basis(1,1)(t_obs);
phi[,4] = fourier_basis(2,0)(t_obs);
phi[,5] = fourier_basis(2,1)(t_obs);
phi[,6] = fourier_basis(3,0)(t_obs);
phi[,7] = fourier_basis(3,1)(t_obs);
phi[,8] = fourier_basis(4,0)(t_obs);
phi[,9] = fourier_basis(4,1)(t_obs);
phi[,10] = fourier_basis(5,0)(t_obs);
phi[,11] = fourier_basis(5,1)(t_obs);
phi[,12] = fourier_basis(6,0)(t_obs);
phi[,13] = fourier_basis(6,1)(t_obs);
phi[,14] = fourier_basis(7,0)(t_obs);
phi[,15] = fourier_basis(7,1)(t_obs);
phi[,16] = fourier_basis(8,0)(t_obs);
phi[,17] = fourier_basis(8,1)(t_obs);
phi[,18] = fourier_basis(9,0)(t_obs);
phi[,19] = fourier_basis(9,1)(t_obs);
phi[,20] = fourier_basis(10,0)(t_obs);

indici <- matrix(nrow=N,ncol=2)
indici[1,] = c(0,0);
indici[2,] = c(1,0);
indici[3,] = c(1,1);
indici[4,] = c(2,0);
indici[5,] = c(2,1);
indici[6,] = c(3,0);
indici[7,] = c(3,1);
indici[8,] = c(4,0);
indici[9,] = c(4,1);
indici[10,] = c(5,0);
indici[11,] = c(5,1);
indici[12,] = c(6,0);
indici[13,] = c(6,1);
indici[14,] = c(7,0);
indici[15,] = c(7,1);
indici[16,] = c(8,0);
indici[17,] = c(8,1);
indici[18,] = c(9,0);
indici[19,] = c(9,1);
indici[20,] = c(10,0);


for(k in 1:K){
  #create PHI_k
  assign(paste("PHI",k,sep="_"),phi[,1:k])
  #create Sigma_k
  Sigma = diag(k) + t(phi[,1:k])%*%phi[,1:k]
  assign(paste("Sigma",k,sep="_"), Sigma)
}

#setting the normalizing constant for the posterior

posterior_denominator=0;

for(k in 1:K){
  posterior_denominator = posterior_denominator + integrated_likelihood(N,get(paste("Sigma",k,sep="_")),get(paste("PHI",k,sep="_")),y_obs)
}


posterior_models=numeric(K);

for(k in 1:K){
  posterior_models[k] = integrated_likelihood(N,get(paste("Sigma",k,sep="_")),get(paste("PHI",k,sep="_")),y_obs)/posterior_denominator
}


# Bayesian Model Averaged Mean Predictior

L = length(t_true)
y_pred = numeric(L)

for(i in 1:L){
  t_star=t_true[i]
  for(k in 1:K){
    phi = numeric(k)
    for(j in 1:k) {
     index = indici[j,1];
     seno = indici[j,2];
     phi[j] = fourier_basis(index,seno)(t_star)
    }
    mu = solve(get(paste("Sigma",k,sep="_")))%*%t(get(paste("PHI",k,sep="_")))%*%y_obs
    y_pred[i] = y_pred[i] + t(phi)%*%mu%*%posterior_models[k]
  }
}

#windows()
plot(t_obs,y_obs)
lines(t_true,y_true,col="red")
lines(t_true,y_pred,col="dark green")
windows()
barplot(posterior_models)


