library("MCMCpack")
set.seed(21122016)
N=25
t_obs=runif(25,-2,2)
t_obs=sort(t)
y_obs=0.219*t_obs^3 + 0.5287*t_obs^2-0.805*t_obs + rnorm(N,0,0.5)

t_true=seq(-2,2,by=0.1)
y_true=0.219*t_true^3 + 0.5287*t_true^2-0.805*t_true
plot(t_obs,y_obs,pch=16)
lines(t_true,y_true,col="red")

rnorminvgamma<-function(n,mu,alpha,gamma){
  if(length(mu)!=n)
    fprintf("error:dimensions must agree")
  sigma <- rinvgamma(alpha, gamma)
  x <- rnorm(mu, sigma*diag(n))
  data.frame(x = x)
}

predictive <- function(k,)