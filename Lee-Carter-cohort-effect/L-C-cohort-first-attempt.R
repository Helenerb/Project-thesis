# R-script for first attempt at inference with inlabru of 
# Lee-Carter models with a cohort effect. 
# The model is inspired by Wiśniowski, A., Smith, P.W.F., Bijak, J. et al., in 
# Bayesian Population Forecasting: Extending the Lee-Carter Method.

# in this file we will follow the previous underlying model as closely as possible
# to be better able to see whether or not the additional cohort effect 
# affected the preformance of the inlabru inference. 

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)

N = 10000

nx = 20
nt = 20

t_min_x.min = 1-nx
t_min_x.max = nt-1
n.t_min_x = t_min_x.max + abs(t_min_x.min) + 1

at.risk = 1000

x = sample(1:nx, N, replace = TRUE)   # 1000 samples of ages 1-100
t = sample(1:nt, N, replace = TRUE)   # 1000 samples of years 1901-2000: represented as 1-100

t_min_x = t-x + abs(t_min_x.min) + 1

obs = data.frame(x,t, t_min_x)

#   underlying model parameters the same as before:
tau.iid = 1/0.1**2   # = 100
tau.epsilon = 1/0.01**2   # =10000

kappa = 2*cos((1:nx)*pi/20)
kappa = kappa - mean(kappa)

alpha = cos(((1:nx - 3)* pi)/6)
alpha = alpha - mean(alpha)

phi = 0.025  #  Drift of random walk

beta = rnorm(nx, 0, sqrt(1/tau.iid))
beta = 1/nx + beta - mean(beta)

m.epsilon = matrix(rnorm(nx*nt, 0, sqrt(1/tau.epsilon)), nx, nt)

# introduce kohort effect:
gamma = (t_min_x.min:t_min_x.max)**3/(3*10**3)
gamma = gamma - mean(gamma)

obs = cbind(obs, beta = beta[as.vector(obs$x)])
obs = cbind(obs, kappa = kappa[obs$t])
obs = cbind(obs, alpha = alpha[obs$x])
obs = cbind(obs, phi.t = phi*obs$t)
obs = cbind(obs, gamma = gamma[obs$t_min_x])

epsilon.map <- function(x,t){
  return(m.epsilon[x, t])
}

obs = cbind(obs, epsilon = apply(obs, 1, function(row) epsilon.map(row['x'], row['t'])))

eta.func <- function(alpha, beta, phit, kappa, epsilon, gamma){
  return(alpha + beta*phit + beta*kappa + gamma + epsilon)
}

obs = cbind(obs, eta = apply(obs, 1,
                             function(row) eta.func(row['alpha'], row['beta'],
                                                    row['phi.t'], row['kappa'],
                                                    row['epsilon'], row['gamma'])))
# add constant at risk
obs = cbind(obs, at.risk = rep(at.risk, N))  

# sample observations:
y = rpois(N, obs$at.risk*exp(obs$eta))
obs = cbind(obs, y = y)

# add extra t to the observations for the sake of inlabru:
obs = cbind(obs, t1 = obs$t)
obs = cbind(obs, x1 = obs$x)

# add xt index for epsilon effect in inlabru:
xt.func <- function(x,t){
  return(x*100 + t)
}
obs = cbind(obs, xt = apply(obs, 1, function(row) xt.func(row['x'], row['t'])))

# Ready to start inference:

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)  #  not sure if you did this correctly
e.vec = 1

pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.2,0.8)))

comp = ~ -1 + 
  alpha(x, model = "rw1", constr = TRUE, hyper = pc.prior) + 
  phi(t, model = "linear") + 
  beta(x1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior) + 
  gamma(t_min_x, model = "rw1", values = 1:n.t_min_x, constr = TRUE) + 
  epsilon(xt, model = "iid")

form.1 = y ~ -1 + alpha + beta*phi + beta*kappa + gamma + epsilon
likelihood.1 = like(formula = form.1, family = "poisson", data = obs, E = at.risk)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)
res = bru(components = comp,
          likelihood.1, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c))

res = bru_rerun(res)

gg.compare <- function(data, title){
  gg = ggplot(data = data, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
    geom_line(aes(y = mean), color = "lightskyblue") + 
    geom_line(aes(y = true_val), color = "dodgerblue1") + 
    ggtitle(title)
  gg
  return(gg)
}

data.beta = cbind(res$summary.random$beta, true_val = beta[res$summary.random$beta$ID])
gg.beta <- gg.compare(data=data.beta, title="Beta"); gg.beta

data.kappa = cbind(res$summary.random$kappa, true_val = kappa[res$summary.random$kappa$ID])
gg.kappa <- gg.compare(data = data.kappa, title = "Kappa"); gg.kappa

data.alpha = cbind(res$summary.random$alpha, true_val = alpha[res$summary.random$alpha$ID])
gg.alpha <- gg.compare(data = data.alpha, title = "Alpha"); gg.alpha

data.eta = cbind(res$summary.linear.predictor, true_val = obs$eta)

eta.sim = res$summary.linear.predictor$mean
data.eta.density = rbind(data.frame(eta = obs$eta, sim = "F"),
                         data.frame(eta = eta.sim, sim = "T"))
gg.eta.density = ggplot(data = data.eta.density, aes(x = eta, color = sim)) + geom_density() + ggtitle("Eta density")
gg.eta.density






