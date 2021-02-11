# reduced version of Lee-Carter model

#   ----  first attempt at inference of Lee-Carter models for mortality 
#  predictions with inlabru  ----   
library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)

N = 10000

nx = 10
nt = 10

at.risk = 1000

x = sample(1:nx, N, replace = TRUE)   # 1000 samples of ages 1-100
t = sample(1:nt, N, replace = TRUE)   # 1000 samples of years 1901-2000: represented as 1-100

obs = data.frame(x,t)

#   model parameters for underlying models:

tau.iid = 1/0.1**2   #  Standard deviation of 0.1 for the iid effects (beta)
# change tau.epsilon from 1/0.1**2 to 1/0.01**2 to see if it was causing too much noise 
tau.epsilon = 1/0.01**2   #  Standard deviation of 0.01 for the noise 

kappa = cos(((1:nt - 3)* pi)/6)
kappa = kappa - mean(kappa)

alpha = -2.0   #  Intercept - perhaps make this iid in the future? 
#alpha = 0 # try what happens if I change the intercept

phi = 0.025  #  Drift of random walk

#  sample synthetic data:
beta = rnorm(nx, 0, sqrt(1/tau.iid))  # should it not depend on t??
beta = 1/nx + beta - mean(beta)   # sum to 1

m.epsilon = matrix(rnorm(nx*nt, 0, sqrt(1/tau.epsilon)), nx, nt)
# how to structure samples from these?  

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

eta.func <- function(alpha, beta, phit, kappa, epsilon){
  return(alpha + beta*phit + beta*kappa + epsilon)
}

obs = cbind(obs, eta = apply(obs, 1,
                             function(row) eta.func(row['alpha'], row['beta'],
                                                    row['phi.t'], row['kappa'],
                                                    row['epsilon'])))
obs = cbind(obs, at.risk = rep(at.risk, N))  # add constant at risk

# sample actual observations:
y = rpois(N, obs$at.risk*exp(obs$eta))
obs = cbind(obs, y = y)

# add extra t to the observations for the sake of inlabru:
obs = cbind(obs, t1 = obs$t)

# add xt index for epsilon effect in inlabru:
xt.func <- function(x,t){
  return(x*100 + t)
}
obs = cbind(obs, xt = apply(obs, 1, function(row) xt.func(row['x'], row['t'])))

# plot observations:
ggplot(data = obs, aes(x=t, y=x, fill = y)) + geom_tile()


#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)  #  not sure if you did this correctly
e.vec = 1

pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.2,0.8)))

comp = ~ Intercept + 
  phi(t, model = "linear") + 
  beta(x, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior) + 
  epsilon(xt, model = "iid")  # change from iid to rw1

form.1 = y ~ Intercept + beta*phi + beta*kappa + epsilon
likelihood.1 = like(formula = form.1, family = "poisson", data = obs)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res = bru(components = comp,
          likelihood.1, 
          options = list(verbose = F,
                         num.threads = "1:1",
                         control.compute = c.c,
                         control.inla = list(int.strategy = "eb"))) 

gg.beta.true = ggplot(data = obs, aes(x = x, y = beta)) + geom_point(color = "hotpink") + ggtitle("True beta")
gg.kappa.true = ggplot(data = obs, aes(x = t, y = kappa)) + geom_line(color = "hotpink") + ggtitle("True kappa")
gg.phi.true = ggplot(data = obs, aes(x = t, y = phi)) + geom_line(color = "hotpink") + ggtitle("True phi")
gg.epsilon.true = ggplot(data = obs, aes(x = x, y = t, fill = epsilon)) + geom_tile() + ggtitle("True epsilon")
ggplot(data = obs, aes(x = epsilon)) + geom_density()

gg.beta = ggplot(data = cbind(res$summary.random$beta, beta.true = beta[res$summary.random$beta$ID]), aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_point(aes(y = mean), color = "lightskyblue") + 
  geom_point(aes(y = beta.true), color = "dodgerblue1")
gg.beta  

data.kappa = cbind(res$summary.random$kappa, kappa.true = kappa[res$summary.random$kappa$ID])
gg.kappa = ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_line(aes(y = mean), color = "lightskyblue") + 
  geom_line(aes(y = kappa.true), color = "dodgerblue1")
gg.kappa

gg.epsilon.true
data.epsilon = res$summary.random$epsilon
data.epsilon = cbind(data.epsilon, x = apply(data.epsilon,1, function(row) row['ID']%/%100 ))
data.epsilon = cbind(data.epsilon, t = apply(data.epsilon, 1, function(row) row['ID']%%100))
gg.epsilon = ggplot(data = data.epsilon, aes(x = x, y = t, fill = mean)) + geom_tile() + ggtitle("Simulated epsilon")
gg.epsilon

# density plots of true and simulated epsilon:
data.epsilon.density = rbind(data.frame(epsilon = data.epsilon$mean, sim = "T"), data.frame(epsilon = obs$epsilon, sim = "F"))
gg.epsilon.density = ggplot(data = data.epsilon.density, aes(x = epsilon, color = sim)) + geom_density()
gg.epsilon.density

#  results of hyperparameters:
cat("Precision for beta: ")
cat("True value: ", tau.iid)
cat("Simulated value: ", res$summary.hyperpar$mean[1])

cat("Precision for kappa: \n", "\"True\" value: ", tau.rw,"\n Simulated value: ",
    res$summary.hyperpar$mean[2])

cat("Precision for epsiloin: \n", "True value: ", tau.epsilon,"\n Simulated value: ",
    res$summary.hyperpar$mean[3])

# density plot of true eta and predicted eta:
eta.sim = res$summary.linear.predictor$mean
eta.sim = eta.sim - log(at.risk)
data.eta.density = rbind(data.frame(eta = obs$eta, sim = "F"),
                         data.frame(eta = eta.sim, sim = "T"))
gg.eta.density = ggplot(data = data.eta.density, aes(x = eta, color = sim)) + geom_density()
gg.eta.density

cat("Intercept: ", res$summary.fixed$mean[1] - log(at.risk))
