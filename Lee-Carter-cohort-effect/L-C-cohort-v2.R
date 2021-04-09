# Lee-Carter model with cohort effect
# Based on the Lee-Carter model that has proven to give the best inference results - 
# this should be considered when evaluating the fit of the model. 

# Note: as both the calendar years (t) and the ages (x) are numbered from 1 to N, 
# some of the cohorts (t-x) will be negative. 
# The mapping between the cohort value and its index in gamma is cohort.val - cohort.min + 1

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)

seed = 324
set.seed(seed)

palette.basis <- c('#70A4D4', '#ECC64B', '#A85150', '#607A4D', '#026AA1')
palette.light <- c('#ABC9E6', '#F3DC90', '#C38281', '#86A46F', '#40BBFD')

N = 1000
general.title = paste("N = ", N, "seed = ", seed)

#nx = 10
#nt = 10
nx = 20
nt = 20

n.cohort = (nt - 1) + abs(1-nx) + 1
cohort.min = 1-nx
cohort.max = nt-1

at.risk = 1000

x = sample(1:nx, N, replace = TRUE)   # 1000 samples of ages 1-100
t = sample(1:nt, N, replace = TRUE)   # 1000 samples of years 1901-2000: represented as 1-100
cohort = t-x

obs = data.frame(x,t,cohort)

#   model parameters for underlying models:

tau.iid = 1/0.1**2   #  Precision of iid beta: 100
tau.epsilon = 1/0.01**2   #  Precision of error term: 10000

#kappa = 0.3*cos((1:nt)*pi/5)
#kappa = sin((1:nt)*pi/20)  #26.02:1300
kappa = 0.5*cos((1:nt)*pi/3)   # 25.02:14, conf 3.1

kappa = kappa - mean(kappa)

#alpha = cos(((1:nx - 3)* pi)/6)
alpha = cos(((1:nx)* pi)/8)  # conf 3.1
alpha = alpha - mean(alpha)

#gamma = 0.2*(cohort.min:cohort.max) + sin(cohort.min:cohort.max/2)
#gamma = 0.2*(cohort.min:cohort.max) + sin(cohort.min:cohort.max)
#gamma = 0.2*(cohort.min:cohort.max) + sin(cohort.min:cohort.max/3)
#gamma = 0.5*(0.2*(cohort.min:cohort.max) + sin(cohort.min:cohort.max/3))
gamma = -0.5*(0.1*(cohort.min:cohort.max) + cos((cohort.min:cohort.max - 2)/4))  # conf 3.1
gamma = gamma - mean(gamma)  #center around zero

phi = -0.5   # conf 3.1

#  sample synthetic data:
beta = rnorm(nx, 0, sqrt(1/tau.iid))  # conf 3.1
beta = 1/nx + beta - mean(beta)   # sum to 1

# note: name all 
obs = obs %>% 
  mutate(beta = beta[as.vector(obs$x)],
         kappa = kappa[as.vector(obs$t)],
         alpha = alpha[as.vector(obs$x)],
         gamma = gamma[as.vector(obs$cohort - cohort.min + 1)],
         phi.t = phi*obs$t,
         phi = phi,
         at.risk = at.risk,
         epsilon = rnorm(n = N, 0, sqrt(1/tau.epsilon))) %>%
  mutate(eta = alpha + beta*phi.t + beta*kappa + gamma + epsilon) %>% # linear predictor
  mutate(y.o = rpois(N, at.risk*exp(eta))) %>%                 # simulate data
  mutate(t1 = t, x1 = x)  %>%                                # add extra t and x to the observations for the sake of inlabru:
  mutate(xt = seq_along(t))

#  plot underlying model:

ggplot(data = obs) + geom_point(aes(x = x, y = beta)) + ggtitle(paste("Underlying beta", general.title))
ggplot(data = obs) + geom_point(aes(x = t, y = kappa)) + ggtitle(paste("Underlying kappa", general.title))
ggplot(data = obs) + geom_line(aes(x = x, y = alpha)) + ggtitle(paste("Underlying alpha", general.title))
ggplot(data = obs) + geom_line(aes(x = t, y = phi.t)) + ggtitle(paste("Underlying phi", general.title))
ggplot(data = obs) + geom_line(aes(x = cohort, y = gamma)) + ggtitle(paste("Underlying gamma", general.title))

# plot observations:xs
ggplot(data = obs, aes(x=t, y=x, fill = y.o)) + geom_tile()


#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)  #  not sure if you did this correctly
e.vec = 1

# pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.1)))
# pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.8)))
# pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.02, 0.1)))
# pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.8, 0.8)))

# attempt with less informative priors:
pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.4)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.5)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))
pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.5)))

# note: change names of components, to ensure no mix-up with global variables and 
# variables in the observation.

comp = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", constr = TRUE, hyper = pc.prior.alpha) + 
  phi(t, model = "linear", mean.linear = -0.5, prec.linear = 0.25) +
  beta(x1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior.kappa) +
  gamma(cohort, model = "rw1", values = cohort.min:cohort.max, constr = TRUE, hyper = pc.prior.gamma) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

form.1 = y.o ~ -1 + Int + alpha + beta*phi + beta*kappa + gamma + epsilon

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

res = bru_rerun(res)

cat(general.title)
res$summary.fixed
res$summary.hyperpar

data.alpha = cbind(res$summary.random$alpha, alpha.true = alpha[res$summary.random$alpha$ID])
gg.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(y = alpha.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Alpha", x = "x", y='')

data.beta = cbind(res$summary.random$beta, beta.true = beta[res$summary.random$beta$ID])
gg.beta <- ggplot(data = data.beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(x = ID, y = beta.true, color = "True value")) +
  geom_point(aes(x = ID, y = mean, color = "Estimated")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Beta", x = "x", y='')

data.kappa = cbind(res$summary.random$kappa, kappa.true = kappa[res$summary.random$kappa$ID])
gg.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(y = kappa.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Kappa", x = "t", y='')

data.gamma = cbind(res$summary.random$gamma, gamma.true = gamma[res$summary.random$gamma$ID - cohort.min + 1])
gg.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(y = gamma.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Gamma", x = "t-x", y='')

data.phi = data.frame(cbind(ID = 1:nt, 
                            mean = res$summary.fixed$mean[2]*1:nt,
                            X0.025quant = res$summary.fixed$`0.025quant`[2]*1:nt,
                            X0.975quant = res$summary.fixed$`0.975quant`[2]*1:nt,
                            phi.true = phi*1:nt))
gg.phi <- ggplot(data = data.phi, aes(x = ID)) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), fill = palette.light[1]) + 
  geom_point(aes(y = phi.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "Method",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Phi * t", x = "t", y='')


# density plot of true eta and predicted eta:
data.eta <- data.frame({eta.sim = res$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = obs$eta)
gg.eta <- ggplot(data.eta) + geom_point(aes(x = eta.sim, y = true.eta)) + 
  labs(x="Simulated", y="True", title = "Eta")

(gg.alpha | gg.beta | gg.kappa)/(gg.phi | gg.gamma | gg.eta) + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

data.eta.density = rbind(data.frame(eta = obs$eta, sim = "Simulated"), data.frame(eta = eta.sim, sim = "True value"))
ggplot(data = data.eta.density, aes(x = eta, color = sim)) + 
  geom_density() + 
  ggtitle(paste("Eta density", general.title))






