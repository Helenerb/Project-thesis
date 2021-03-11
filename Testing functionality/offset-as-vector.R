# This script will be used to investigate how to enter a vector as an offset:

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)

seed = 324
set.seed(seed)

N = 1000

general.title = paste("N = ", N, "seed = ", seed, " : L-C-alpha-x")

nx = 10
nt = 10

#at.risk = 1000
at.risk.vec = rep(1000, N)

x = sample(1:nx, N, replace = TRUE)   
t = sample(1:nt, N, replace = TRUE)   

obs = data.frame(x,t)

#   model parameters for underlying models:

tau.iid = 1/0.1**2   #  Standard deviation of 0.1 for the iid effects (beta)
tau.epsilon = 1/0.01**2   #  Standard deviation of 0.01 for the noise 

# change kappa to something better suited to model real-life 
#kappa = cos((1:nt)*pi/8)
#kappa = kappa - mean(kappa)
# Note: here, phi is less than sd(kappa) - might cause problems?

# old kappa that caused problems:
#kappa = 2*cos((1:nt)*pi/20)
#kappa = cos((1:nt)*pi/20)
#kappa = sin((1:nt)*pi/20)  # 26.02:1117
kappa = 0.3*cos((1:nt)*pi/5)  # 
#kappa = 0.5*cos((1:nt)*pi/3)  # roughly the same sd as the above
kappa = kappa - mean(kappa)

# change this into an effect of x
alpha = cos(((1:nx - 3)* pi)/6)
alpha = alpha - mean(alpha)

# attempt to make alpha as similar to beta as possible and see if they mix 
# alpha <- devs <- rnorm(nt, mean = 0, sd = sqrt(1/tau.iid))
# for (i in 2:nt){
#   alpha[i] = alpha[i-1] + devs[i]
# }
# alpha = alpha - mean(alpha)

#phi = -0.25  # gave good result
phi = -0.5  # increse to be bigger than sd(kappa)
#phi = 0.025  #  Drift of random walk

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
         at.risk = at.risk.vec,
         epsilon = rnorm(n = N, 0, sqrt(1/tau.epsilon))) %>%
  mutate(eta = alpha + beta*phi.t + beta*kappa + epsilon) %>% # linear predictor
  mutate(y.o = rpois(N, at.risk*exp(eta))) %>%                 # simulate data
  mutate(t1 = t, x1 = x)  %>%                                # add extra t and x to the observations for the sake of inlabru:
  mutate(xt = seq_along(t))

# plot observations:
ggplot(data = obs, aes(x=t, y=x, fill = y.o)) + geom_tile() + ggtitle(paste("Observations: ", general.title))


#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)  #  not sure if you did this correctly
e.vec = 1

#pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.2,0.8)))
pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.3,0.8)))
#pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.1,0.1)))
pc.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1,0.1)))
pc.prior.small <- list(prec = list(prior = "pc.prec", param = c(0.02, 0.1)))

comp = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", constr = TRUE, hyper = pc.alpha) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior) +
  epsilon(xt, model = "iid", hyper = pc.prior.small)

form.1 = y.o ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon

likelihood.1 = like(formula = form.1, family = "poisson", data = obs, E = at.risk.vec)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha = alpha, beta = beta, kappa = kappa, phi.t = phi*(1:nt))

res = bru(components = comp,
          likelihood.1, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c
          )) 

res = bru_rerun(res)

cat(general.title)
res$summary.fixed
res$summary.hyperpar

data.alpha = cbind(res$summary.random$alpha, alpha.true = alpha[res$summary.random$alpha$ID])
ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  geom_point(aes(y = alpha.true, color = "True value")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = "lightskyblue", "True value" = "dodgerblue1") ) +
  ggtitle(paste("Alpha: ", general.title))

data.beta = cbind(res$summary.random$beta, beta.true = beta[res$summary.random$beta$ID])
ggplot(data = data.beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_point(aes(x = ID, y = mean, color = "Estimated")) + 
  geom_point(aes(x = ID, y = beta.true, color = "True value")) +
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = "lightskyblue", "True value" = "dodgerblue1") ) + 
  ggtitle(paste("Beta: ", general.title))

data.kappa = cbind(res$summary.random$kappa, kappa.true = kappa[res$summary.random$kappa$ID])
ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  geom_point(aes(y = kappa.true, color = "True value")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = "lightskyblue", "True value" = "dodgerblue1") ) + 
  ggtitle(paste("Kappa: ", general.title))

data.phi = data.frame(cbind(ID = 1:nt, 
                            mean = res$summary.fixed$mean[2]*1:nt,
                            X0.025quant = res$summary.fixed$`0.025quant`[2]*1:nt,
                            X0.975quant = res$summary.fixed$`0.975quant`[2]*1:nt,
                            phi.true = phi*1:nt))
ggplot(data = data.phi, aes(x = ID)) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), fill = "lightskyblue1") + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  geom_point(aes(y = phi.true, color = "True value")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = "lightskyblue", "True value" = "dodgerblue1") ) + 
  ggtitle(paste("Phi: ", general.title))

# density plot of true eta and predicted eta:
data.frame({eta.sim = res$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = obs$eta) %>%
  ggplot() + geom_point(aes(x = eta.sim, y = true.eta)) + 
  ggtitle(paste("Eta: ", general.title))

data.eta.density = rbind(data.frame(eta = obs$eta, sim = "Simulated"), data.frame(eta = eta.sim, sim = "True value"))
ggplot(data = data.eta.density, aes(x = eta, color = sim)) + 
  geom_density() + 
  ggtitle(paste("Eta density", general.title))




