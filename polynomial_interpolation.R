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

K=10


# Store in each row k, the residuals sum of squares and adjusted R^2 

results=matrix(nrow=10,ncol=2)

#fit <- lm(y_obs  ~ t_obs)
#results[1,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2))
results[2,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3))
results[3,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3) + I(t_obs^4))
results[4,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3) + I(t_obs^4) + I(t_obs^5))
results[5,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3)+ I(t_obs^4) + I(t_obs^5)+ I(t_obs^6))
results[6,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3)+ I(t_obs^4) + I(t_obs^5)+ I(t_obs^6) + I(t_obs^7))
results[7,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3)+ I(t_obs^4) + I(t_obs^5)+ I(t_obs^6) + I(t_obs^7)+ I(t_obs^8))
results[8,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3)+ I(t_obs^4) + I(t_obs^5)+ I(t_obs^6) + I(t_obs^7)+ I(t_obs^8) + I(t_obs^9))
results[9,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3)+ I(t_obs^4) + I(t_obs^5)+ I(t_obs^6) + I(t_obs^7)+ I(t_obs^8) + I(t_obs^9)+ I(t_obs^10))
results[10,]=c(fit$residuals%*%fit$residuals,summary(fit)[[9]])

results = results[2:10,]

par(mfrow=c(1,2))
plot(results[,1])
plot(results[,2])


## based on R^2 & Rss I would choose model of order 7
## I use it to interpolate

fit <- lm(y_obs  ~ t_obs + I(t_obs^2) + I(t_obs^3)+ I(t_obs^4) + I(t_obs^5)+ I(t_obs^6) + I(t_obs^7))
coeff =fit$coefficients

L = length(t_true)
y_pred = numeric(L)
for(i in 1:L){
  t_star=t_true[i]
  tt = c(1,t_star,t_star^2,t_star^3,t_star^4,t_star^5,t_star^6,t_star^7)
  y_pred[i] = coeff%*%tt
}

plot(t_obs,y_obs)
lines(t_true,y_true,col="red")
lines(t_true,y_pred,col="dark green")