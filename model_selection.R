###################
# Model Selection #
###################

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

polinomail_model_prediction <- function(K,t_obs,y_obs){
  #param K is the order of the model
  prior_models=rep(1/K,K)
  
  N = length(y_obs)
  
  integrated_likelihood <- function(N,Sigma,PHI,y_obs)
  {
    out = gamma(1+N/2)/sqrt(pi^N*det(Sigma)) * (1 + 0.5* t(y_obs)%*%(diag(N)-PHI%*%solve(Sigma)%*%t(PHI))%*%y_obs)^(-1-N/2)
    return(out)
  }
  
  phi = matrix(nrow=N,ncol=1)
  phi[,1] = t_obs;
  PHI_1 = phi;
  Sigma_1 = 1 + t(phi)%*%phi
  mu_1 = solve(Sigma_1)%*%t(PHI_1)%*%y_obs
  
if(K>1){
  for(k in 2:K){
    #for different basis functions, define the following line accordingly
    phi = cbind(phi,t_obs^k)
    #create PHI_k
    assign(paste("PHI",k,sep="_"),phi)
    #create Sigma_k
    Sigma = diag(k) + t(phi)%*%phi
    assign(paste("Sigma",k,sep="_"), Sigma)
  }
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
      phi_ = numeric(k)
      for(j in 1:k){
        phi_[j] = t_star^j
      }
      mu = solve(get(paste("Sigma",k,sep="_")))%*%t(get(paste("PHI",k,sep="_")))%*%y_obs
      y_pred[i] =t(phi_)%*%mu
    }
  }
  
  windows()
  plot(t_obs,y_obs,main = paste(K))
  lines(t_true,y_true,col="red")
  lines(t_true,y_pred,col="dark green")
 # windows()
 # barplot(posterior_models)
  
}


for(k in 2:10){
  k
  polinomail_model_prediction(k,t_obs,y_obs)
}

