# Lee-Carter model with alpha modelled as effect of x
# In this script, kappa is modelled as realisations from a random walk.  

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)

seed = 226
set.seed(seed)

N = 5000
general.title = paste("N = ", N, "seed = ", seed)

nx = 10
nt = 10

at.risk = 1000

x = sample(1:nx, N, replace = TRUE)   # 1000 samples of ages 1-100
t = sample(1:nt, N, replace = TRUE)   # 1000 samples of years 1901-2000: represented as 1-100

obs = data.frame(x,t)

#   model parameters for underlying models:

tau.iid = 1/0.1**2   #  Standard deviation of 0.1 for the iid effects (beta)
tau.epsilon = 1/0.01**2   #  Standard deviation of 0.01 for the noise 
tau.rw = 1/0.1**2

# model kappa as a random walk:
kappa <- devs <- rnorm(nt, mean = 0, sd = sqrt(1/tau.rw))
for (i in 2:nt){
  kappa[i] = kappa[i-1] + devs[i]
}
kappa = kappa - mean(kappa)

alpha = cos(((1:nx - 3)* pi)/6)
alpha = alpha - mean(alpha)

phi = -0.25 

#  sample synthetic data:
beta = rnorm(nx, 0, sqrt(1/tau.iid))  # should it not depend on t??
beta = 1/nx + beta - mean(beta)   # sum to 1

# note: name all 
obs = obs %>% 
  mutate(beta = beta[as.vector(obs$x)],
         kappa = kappa[as.vector(obs$t)],
         alpha = alpha[as.vector(obs$x)],
         phi.t = phi*obs$t,
         phi = phi,
         at.risk = at.risk,
         epsilon = rnorm(n = N, 0, sqrt(1/tau.epsilon))) %>%
  mutate(eta = alpha + beta*phi.t + beta*kappa + epsilon) %>% # linear predictor
  mutate(y.o = rpois(N, at.risk*exp(eta))) %>%                 # simulate data
  mutate(t1 = t, x1 = x)  %>%                                # add extra t and x to the observations for the sake of inlabru:
  mutate(xt = seq_along(t))

#  plot underlying model:

obs %>% ggplot() + geom_point(aes(x = x, y = beta)) + ggtitle(paste("Underlying beta", general.title))
obs %>% ggplot() + geom_point(aes(x = t, y = kappa)) + ggtitle(paste("Underlying kappa", general.title))
obs %>% ggplot() + geom_line(aes(x = x, y = alpha)) + ggtitle(paste("Underlying alpha", general.title))
obs %>% ggplot() + geom_line(aes(x = t, y = phi.t)) + ggtitle(paste("Underlying phi", general.title))

# plot observations:xs
ggplot(data = obs, aes(x=t, y=x, fill = y.o)) + geom_tile()


#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)  #  not sure if you did this correctly
e.vec = 1

pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.1)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.2, 0.8)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.02, 0.1)))

# note: change names of components, to ensure no mix-up with global variables and 
# variables in the observation.

comp = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", constr = TRUE, hyper = pc.prior.alpha) + 
  phi(t, model = "linear", prec.linear = 0.25) +
  beta(x1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior.kappa) +
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)  # change from iid to rw1

form.1 = y.o ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon

likelihood.1 = like(formula = form.1, family = "poisson", data = obs, E = at.risk)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha.c = alpha, beta.c = beta, kappa.c = kappa, phi.t = phi*(1:nt))

res = bru(components = comp,
          likelihood.1, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c
          )) 

#res = bru_rerun(res)

cat(general.title)
res$summary.fixed
res$summary.hyperpar


data.alpha = cbind(res$summary.random$alpha, alpha.true = alpha[res$summary.random$alpha$ID])
  ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_line(aes(y = mean), color = "lightskyblue") + 
  geom_line(aes(y = alpha.true), color = "dodgerblue1") + 
  ggtitle(paste("Alpha: ", general.title))

data.beta = cbind(res$summary.random$beta, beta.true = beta[res$summary.random$beta$ID])
ggplot(data = data.beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_point(aes(x = ID, y = mean), color = "lightskyblue") + 
  geom_line(aes(x = ID, y = beta.true), color = "dodgerblue1") + 
  ggtitle(paste("Beta: ", general.title))

data.kappa = cbind(res$summary.random$kappa, kappa.true = kappa[res$summary.random$kappa$ID])
  ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_line(aes(y = mean), color = "lightskyblue") + 
  geom_point(aes(y = kappa.true), color = "dodgerblue1") + 
  ggtitle(paste("Kappa: ", general.title))


# density plot of true eta and predicted eta:
data.frame({eta.sim = res$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = obs$eta) %>%
  ggplot() + geom_point(aes(x = eta.sim, y = true.eta)) + 
  ggtitle(paste("Eta: ", general.title))

data.eta.density = rbind(data.frame(eta = obs$eta, sim = "F"), data.frame(eta = eta.sim, sim = "T"))
  ggplot(data = data.eta.density, aes(x = eta, color = sim)) + 
  geom_density() + 
  ggtitle(paste("Eta density", general.title))


