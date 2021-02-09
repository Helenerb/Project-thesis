#   ----  first attempt at inference of Lee-Carter models for mortality 
#  predictions with inlabru  ----   

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)

N = 10000

x = sample(1:100, N, replace = TRUE)   # 1000 samples of ages 1-100
t = sample(1:100, N, replace = TRUE)   # 1000 samples of years 1901-2000: represented as 1-100

obs = data.frame(x,t)

#   model parameters for underlying models:

tau.iid = 1/0.1**2   #  Standard deviation of 0.1 for the iid effects (beta)
tau.rw = 1/0.2**2    #  Standard deviation of 0.2 for the random walk
tau.epsilon = 1/0.05**2   #  Standard deviation of 0.05 for the noise 

kappa = -cos(((1:100 - 20)* pi)/80)
kappa = kappa - mean(kappa)

alpha = -2.0   #  Intercept - perhaps make this iid in the future? 

phi = 0.025  #  Drift of random walk

#  sample synthetic data:
beta = rnorm(100, 0, sqrt(1/tau.iid))  # should it not depend on t??
beta = 1/100 + beta - mean(beta)   # sum to 1

beta.x = beta[x]

phi.t = phi*t

kappa.t = kappa[t]

m.epsilon = matrix(rnorm(100*100, 0, sqrt(1/tau.epsilon)), 100, 100)
# how to structure samples from these?  

eta = alpha + beta.x*phi.t + beta.x*kappa.t + m.epsilon[x,t]  # is this the correct representation of eta?
# no, you now have all combinations of the 1000 samples of x and t each. You want the grid representation. 

# add everything to the dataframe instead of specifying beta.x and kappa.t etc separately. 
# expand the obs dataframe:

#obs = cbind(obs, epsilon = rnorm(1000, 0, sqrt(1/tau.epsilon)))  # add epsilon
obs = cbind(obs, beta = beta[as.vector(obs$x)])
obs = cbind(obs, kappa = kappa[obs$t])
obs = cbind(obs, alpha = rep(alpha, N))
obs = cbind(obs, phi.t = phi*obs$t)

epsilon.map <- function(x,t){
  return(m.epsilon[x, t])
}

obs = cbind(obs, epsilon = apply(obs, 1, function(row) epsilon.map(row['x'], row['t'])))

eta.func <- function(alpha, beta, phit, kappa){
  return(alpha + beta*phit + beta*kappa)
}

obs = cbind(obs, eta = apply(obs, 1, function(row) eta.func(row['alpha'], row['beta'], row['phi.t'], row['kappa'])))
obs = cbind(obs, at.risk = rep(1000, N))  # add constant at risk:

# sample actual observations:
y = rpois(N, obs$at.risk*exp(obs$eta))
obs = cbind(obs, y = y)

# add extra t to the observations for the sake of inlabru:
obs = cbind(obs, t1 = obs$t)

# plot observations:
ggplot(data = obs, aes(x=x, y=t, fill = y)) + geom_tile()


#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A = matrix(1, nrow = 100, ncol = 1)  #  not sure if you did this correctly
e = rep(1,1)

comp = ~ Intercept + 
  phi(t, model = "linear") + 
  beta(x, model = "iid", extraconst = list(A = A, e = e)) + 
  kappa(t1, model = "rw1", values = 1:100, constr = TRUE) + 
  epsilon(tx, model = "iid")
  
  

