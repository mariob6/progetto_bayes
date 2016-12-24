library("MCMCpack")
library("MASS")
set.seed(21122016)
N=25
t_obs=runif(25,-2,2)
t_obs=sort(t_obs)
y_obs=0.219*t_obs^3 + 0.5287*t_obs^2-0.805*t_obs + rnorm(N,0,0.5)

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
  #create mu_k, useful for the predictive
    mu = solve(Sigma)%*%t(phi)%*%y_obs
    assign(paste("mu",k,sep="_"),mu)
    MU[k]=mu
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

mu_ <- function(t_star,k,mu_k){
  out = 0;
  phi = numeric(k)
  for(j in 1:k){
    phi[j] = t_star^j
  }
 # for(j in 1:k){
    out = out + t(phi)%*%mu_k*posterior_models[j]
  #}
  return(out)
}

#make predictions, evaluate mu_star

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
plot(t_true,y_pred)
windows()
barplot(posterior_models)

