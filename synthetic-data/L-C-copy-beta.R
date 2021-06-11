#   Basic, reduced Lee-Carter with the copy-method on beta:
# reduced version of Lee-Carter model

# load workspace:
load("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Workspaces/L-C-copy-beta.RData")

#   ----  first attempt at inference of Lee-Carter models for mortality 
#  predictions with inlabru  ----   
library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)

# define palettes:
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1',
                   '#5d8060', '#D7B36A', '#826133', '#A85150')

seed = 324
set.seed(seed)

N = 500

nx = 10
nt = 10

at.risk = 1000

x = sample(1:nx, N, replace = TRUE)   # 1000 samples of ages 1-100
t = sample(1:nt, N, replace = TRUE)   # 1000 samples of years 1901-2000: represented as 1-100

obs = data.frame(x,t)

#   model parameters for underlying models:

tau.iid = 1/0.1**2   #  Standard deviation of 0.1 for the iid effects (beta)
tau.epsilon = 1/0.01**2   #  Standard deviation of 0.01 for the noise 

kappa = cos(((1:nt - 3)* pi)/6)
kappa = kappa - mean(kappa)

alpha = -2.0   #  model alpha as a constant to make the model as simple as possible. 

phi = 0.025  #  Drift of random walk

#  sample synthetic data:
beta = rnorm(nx, 0, sqrt(1/tau.iid))  # should it not depend on t??
beta = 1/nx + beta - mean(beta)   # sum to 1

obs = obs %>% 
  mutate(beta = beta[as.vector(obs$x)],
         kappa = kappa[obs$t],
         alpha = rep(alpha, N),
         phi.t = phi*obs$t,
         phit = phi,
         at.risk = at.risk,
         epsilon = rnorm(n = N, 0, sqrt(1/tau.epsilon))) %>%
  mutate(eta = alpha + beta*phi.t + beta*kappa + epsilon) %>% # linear predictor
  mutate(y = rpois(N, at.risk*exp(eta))) %>%                 # simulate data
  mutate(t1 = t, x1 = x)  %>%                                # add extra t and x to the observations for the sake of inlabru:
  mutate(xt = seq_along(t))


# plot observations:
ggplot(data = obs, aes(x=t, y=x, fill = y)) + geom_tile()


#   ----  Start defining the inlabru model components  ----   

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = nx)  #  not sure if you did this correctly
e.vec = 1

#pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.2,0.8)))
pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.5)))

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

# components with two beta, beta2 a copy of beta1
comp.copy = ~ Intercept + 
  phi(t, model = "linear") + 
  beta1(x, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior) + 
  beta2(x1, model = "iid", copy="beta1") +
  epsilon(xt, model = "iid")

form.copy = y ~ Intercept + beta1*phi + beta2*kappa + epsilon
likelihood.copy = like(formula = form.copy, family = "poisson", data = obs, E = at.risk)


runtime.copy = system.time({res.copy = bru(components = comp.copy,
          likelihood.copy, 
          options = list(verbose = F,
                         bru_verbose = 1,
                         num.threads = "1:1",
                         control.compute = c.c)) })

runtime.copy.rerun = system.time({res.copy = bru_rerun(res.copy)})
cat("Runtime for the copy method: ")
print(runtime.copy)
cat("\n runtime for rerun: ")
print(runtime.copy.rerun)

comp.single = ~ Intercept + 
  phi(t, model = "linear") + 
  beta(x, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t1, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior) + 
  epsilon(xt, model = "iid")

form.single = y ~ Intercept + beta*phi + beta*kappa + epsilon
likelihood.single = like(formula = form.single, family = "poisson", data = obs, E = at.risk)

runtime.single = system.time({res.single = bru(components = comp.single,
               likelihood.single, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c))}) 

runtime.single.rerun = system.time({res.single = bru_rerun(res.single)})

cat("Runtime for the single-beta method: ")
print(runtime.single)
cat("\n runtime for rerun: ")
print(runtime.single.rerun)

# new and modern plots B)
data.beta <- data.frame(rbind(res.copy$summary.random$beta1, res.copy$summary.random$beta2, res.single$summary.random$beta)) %>%
  mutate(type = rep(c("Copied 1","Copied 2","Single"), each = 10))

gg.beta = ggplot(data.beta) + geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant, color =type), position=position_dodge(width=0.5)) +
  scale_color_manual(name = "",
                     values = c(palette.basis[1], palette.basis[1], palette.basis[2])) + 
  labs(title="Beta", x = "x", y = '')
gg.beta

data.kappa <- data.frame(rbind(res.copy$summary.random$kappa,  res.single$summary.random$kappa)) %>%
  mutate(type = rep(c("Copied","Single"), each = 10))

gg.kappa <- ggplot(data.kappa) + geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant, color =type), position=position_dodge(width=0.5)) +
  scale_color_manual(name = "", values = palette.basis) + 
  labs(title="Kappa", x = "t", y = '')

data.frame(res.single$marginals.fixed) %>%
  rbind(data.frame(res.copy$marginals.fixed)) %>%
  mutate(type = rep(c("Single", "Copied"), each = 75))

p.phi <- ggplot(data.frame(res.single$marginals.fixed) %>%
                  rbind(data.frame(res.copy$marginals.fixed)) %>%
                  mutate(type = rep(c("Single", "Copied"), each = 75)))  + 
  geom_area(aes(x = phi.x, y = phi.y, fill = type), alpha = 0.4) + 
  geom_vline(data = res.single$summary.fixed, aes(xintercept = mean[2], color = "Single", fill = "Single")) + 
  geom_vline(data = res.copy$summary.fixed, aes(xintercept = mean[2], color = "Copied", fill = "Copied")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi

p.alpha <- ggplot(data.frame(res.single$marginals.fixed) %>%
                  rbind(data.frame(res.copy$marginals.fixed)) %>%
                  mutate(type = rep(c("Single", "Copied"), each = 75)))  + 
  geom_area(aes(x = Intercept.x, y = Intercept.y, fill = type), alpha = 0.4) + 
  geom_vline(data = res.single$summary.fixed, aes(xintercept = mean[1], color = "Single", fill = "Single")) + 
  geom_vline(data = res.copy$summary.fixed, aes(xintercept = mean[1], color = "Copied", fill = "Copied")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of alpha", y = " ", title = "Alpha")
p.alpha

data.eta <- data.frame(eta.sim.s = res.single$summary.linear.predictor$mean[1:N]) %>%
  mutate(true.eta = obs$eta) %>%
  mutate(eta.sim.c = res.copy$summary.linear.predictor$mean[1:N])
p.eta <- ggplot(data = data.eta) +
  geom_point(aes(x = eta.sim.s, y = eta.sim.c), color = palette.basis[1]) + 
  labs(x="Estimated eta with single beta", y="Estimated eta with copied beta", title = "Eta")
p.eta

p.copy.beta <- (p.alpha | gg.beta | gg.kappa)/(p.phi | p.eta) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects using the copied and the single beta implementation")
p.copy.beta

ggsave('copy-beta.png',
       plot = p.copy.beta,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plots of distribution of hyperparameter:

p.prec.beta <- ggplot() + 
  geom_area(data = data.frame(res.single$marginals.hyperpar) %>%
              filter(Precision.for.beta.x < 200),
            aes(x = Precision.for.beta.x, y = Precision.for.beta.y, fill = "Single"), alpha = 0.4) + 
  geom_vline(data = res.single$summary.hyperpar, aes(xintercept = mean[1], color = "Single", fill = "Single")) + 
  geom_area(data = data.frame(res.copy$marginals.hyperpar) %>%
              filter(Precision.for.beta1.x < 200),
            aes(x = Precision.for.beta1.x, y = Precision.for.beta1.y, fill = "Copied"), alpha = 0.4) + 
  geom_vline(data = res.copy$summary.hyperpar, aes(xintercept = mean[1], color = "Copied", fill = "Copied")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of precision", y = " ", title = "Precision for beta")
p.prec.beta
# single, copy: 2 values above 1000


p.prec.kappa <- ggplot() + 
  geom_area(data = data.frame(res.single$marginals.hyperpar) %>%
              filter(Precision.for.kappa.x < 30),
            aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, fill = "Single"), alpha = 0.4) + 
  geom_vline(data = res.single$summary.hyperpar, aes(xintercept = mean[2], color = "Single", fill = "Single")) + 
  geom_area(data = data.frame(res.copy$marginals.hyperpar) %>%
              filter(Precision.for.kappa.x < 30),
            aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, fill = "Copied"), alpha = 0.4) + 
  geom_vline(data = res.copy$summary.hyperpar, aes(xintercept = mean[2], color = "Copied", fill = "Copied")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of precision", y = " ", title = "Precision for kappa")
p.prec.kappa
# copy, single : 2 values above 200

p.prec.epsilon <- ggplot() + 
  geom_area(data = data.frame(res.single$marginals.hyperpar) %>%
              filter(Precision.for.epsilon.x < 100000),
            aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, fill = "Single"), alpha = 0.4) + 
  geom_vline(data = res.single$summary.hyperpar, aes(xintercept = mean[3], color = "Single", fill = "Single")) + 
  geom_area(data = data.frame(res.copy$marginals.hyperpar) %>%
              filter(Precision.for.epsilon.x < 100000),
            aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, fill = "Copied"), alpha = 0.4) + 
  geom_vline(data = res.copy$summary.hyperpar, aes(xintercept = mean[3], color = "Copied", fill = "Copied")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of precision", y = " ", title = "Precision for epsilon")
p.prec.epsilon


p.hyperpars <- (p.prec.beta | p.prec.kappa | p.prec.epsilon) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom') &
  plot_annotation(title = "Estimated hyperparameters using the copied and the single beta implementation")
p.hyperpars

ggsave('hyperparameters-LC-copy.png',
       plot = p.hyperpars,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Figures',
       height = 5, width = 8,
       dpi = "retina"
)

# save workspace image
save.image("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/synthetic-data/Workspaces/L-C-copy-beta.RData")
