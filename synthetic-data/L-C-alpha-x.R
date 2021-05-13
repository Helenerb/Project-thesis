# Lee-Carter model with alpha modelled as effect of x
# use same model as in previouys attempts, only add alpha as an effect 
# of x, as this is more similar to the model used in other papers. 

# load workspace image
load("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Workspaces/L-C-alpha-x.RData")

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)

# define palette:
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1',
                   '#5d8060', '#D7B36A', '#826133', '#A85150')

seed = 324
set.seed(seed)

N = 500

general.title = paste("N = ", N, "seed = ", seed, " : L-C-alpha-x")

nx = 10
nt = 10

at.risk = 1000

x = sample(1:nx, N, replace = TRUE)   
t = sample(1:nt, N, replace = TRUE)   

obs = data.frame(x,t)

#   model parameters for underlying models:

#  Standard deviation of 0.1 for the iid effects (beta)
tau.iid = 1/0.1**2  # config 2.1

#  Standard deviation of 0.01 for the noise 
tau.epsilon = 1/0.01**2   # config 2.1

# comment in and out different model choices to test:
#kappa = cos((1:nt)*pi/8)

#kappa = 2*cos((1:nt)*pi/20)
kappa = cos((1:nt)*pi/20)  #  config 2.2
#kappa = sin((1:nt)*pi/20)  # 26.02:1117
#kappa = 0.3*cos((1:nt)*pi/5)  # config 2.1
#kappa = 0.5*cos((1:nt)*pi/3)  # roughly the same sd as the above
kappa = kappa - mean(kappa)

# change this into an effect of x
alpha = cos(((1:nx - 3)* pi)/6)  # config 2.1
alpha = alpha - mean(alpha)

# alpha <- devs <- rnorm(nt, mean = 0, sd = sqrt(1/tau.iid))
# for (i in 2:nt){
#   alpha[i] = alpha[i-1] + devs[i]
# }
# alpha = alpha - mean(alpha)

#phi = -0.25 
phi = -0.5  # config 2.1
#phi = 0.025

#  sample synthetic data:
beta = rnorm(nx, 0, sqrt(1/tau.iid))  # config 2.1
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

#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)
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

likelihood.1 = like(formula = form.1, family = "poisson", data = obs, E = at.risk)

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

#res = bru_rerun(res)

# new plot set-up:

p.Int <- ggplot(data.frame(res$marginals.fixed)) + 
  geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = res$summary.fixed, aes(xintercept = mean[1], color = "Estimated")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of intercept", y = " ", title = "Intercept")

p.Int

data.alpha = cbind(res$summary.random$alpha, alpha.true = alpha[res$summary.random$alpha$ID])
p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(title="Alpha", x = "x", y='')

p.alpha

data.beta = cbind(res$summary.random$beta, beta.true = beta[res$summary.random$beta$ID])
p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "x", y = "beta", title = "Beta")

p.beta

data.kappa = cbind(res$summary.random$kappa, kappa.true = kappa[res$summary.random$kappa$ID])
p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa

p.phi <- ggplot(data.frame(res$marginals.fixed)) + 
  geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = res$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill = "Estimated")) + 
  geom_vline(data = obs, aes(xintercept = phi, color = "True value", fill = "True value")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi

data.eta <- data.frame({eta.sim = res$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = obs$eta)
p.eta <- ggplot(data = data.eta) +
  geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
  labs(x="Estimated eta", y="True value for eta", title = "Eta")
p.eta

# configuration 2.1 --> basic LC model with alpha as an effect of x
# p.alpha.x <- (p.alpha | p.beta | p.kappa)/(p.phi | p.eta) +
#   plot_layout(guides = "collect") & 
#   plot_annotation(title = "Estimated random effects for Lee-Carter model, with synthetic data")
# p.alpha.x
# 
# ggsave('effects-LC-synthetic.png',
#        plot = p.alpha.x,
#        device = "png",
#        path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Figures',
#        height = 5, width = 8, 
#        dpi = "retina"
# )

# configuration 2.2 --> LC model displaying identifiability issues:
p.LC.identifiability <- (p.alpha | p.beta | p.kappa)/(p.phi | p.eta) +
  plot_layout(guides = "collect") &
  plot_annotation(title = "Estimated random effects for Lee-Carter model, with synthetic data")
p.LC.identifiability

ggsave('effects-LC-synthetic-identifiability.png',
       plot = p.LC.identifiability,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Figures',
       height = 5, width = 8,
       dpi = "retina"
)

# old plot set-up:

cat(general.title)
res$summary.fixed
res$summary.hyperpar

data.alpha = cbind(res$summary.random$alpha, alpha.true = alpha[res$summary.random$alpha$ID])
p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(y = alpha.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Alpha", x = "x", y='')

data.beta = cbind(res$summary.random$beta, beta.true = beta[res$summary.random$beta$ID])
p.beta <- ggplot(data = data.beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(x = ID, y = beta.true, color = "True value")) +
  geom_point(aes(x = ID, y = mean, color = "Estimated")) + 
  scale_color_manual(name = "",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title="Beta", x = "x", y = '')

data.kappa = cbind(res$summary.random$kappa, kappa.true = kappa[res$summary.random$kappa$ID])
p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.light[1]) + 
  geom_point(aes(y = kappa.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title = "Kappa", x = "t", y = '')

data.phi = data.frame(cbind(ID = 1:nt, 
                 mean = res$summary.fixed$mean[2]*1:nt,
                 X0.025quant = res$summary.fixed$`0.025quant`[2]*1:nt,
                 X0.975quant = res$summary.fixed$`0.975quant`[2]*1:nt,
                 phi.true = phi*1:nt))
p.phi <- ggplot(data = data.phi, aes(x = ID)) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), fill = palette.light[1]) + 
  geom_point(aes(y = phi.true, color = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated")) + 
  scale_color_manual(name = "",
                     breaks = c("Estimated", "True value"),
                     values = c("Estimated" = palette.basis[1], "True value" = palette.basis[2]) ) +
  labs(title = "Phi", x = "t", y='')

# density plot of true eta and predicted eta:
data.eta <- data.frame({eta.sim = res$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = obs$eta)
p.eta <- ggplot(data = data.eta) + geom_point(aes(x = eta.sim, y = true.eta)) + 
  labs(x="Simulated", y="True", title = "Eta")
  
(p.alpha | p.beta | p.kappa) / (p.phi | p.eta) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

data.eta.density = rbind(data.frame(eta = obs$eta, sim = "Simulated"), data.frame(eta = eta.sim, sim = "True value"))
ggplot(data = data.eta.density, aes(x = eta, color = sim)) + 
  geom_density() + 
  ggtitle(paste("Eta density", general.title))

# save workspace image 
save.image("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Workspaces/L-C-alpha-x.RData")




