#   Basic, reduced Lee-Carter with the copy-method on beta:
# reduced version of Lee-Carter model

#   ----  first attempt at inference of Lee-Carter models for mortality 
#  predictions with inlabru  ----   
library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)

# define palettes:
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#A85150', '#026AA1')
palette.light <- c('#ABC9E6', '#F3DC90', '#86A46F', '#C38281', '#40BBFD')

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

pc.prior <- list(prec = list(prior = "pc.prec", param = c(0.2,0.8)))

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

cat("Runtime for the copy method: ")
print(runtime.single)
cat("\n runtime for rerun: ")
print(runtime.single.rerun)

# new and modern plots B)
data.beta <- data.frame(rbind(res.copy$summary.random$beta1, res.copy$summary.random$beta2, res.single$summary.random$beta)) %>%
  mutate(type = rep(c("copy1","copy2","single"), each = 10))

gg.beta = ggplot(data.beta) + geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant, color =type), position=position_dodge(width=0.5)) +
  scale_color_manual(name = "", values = palette.basis) + 
  labs(title="Beta", x = "x", y = '')

data.kappa <- data.frame(rbind(res.copy$summary.random$kappa,  res.single$summary.random$kappa)) %>%
  mutate(type = rep(c("copy1","single"), each = 10))

gg.kappa <- ggplot(data.kappa) + geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant, color =type), position=position_dodge(width=0.5)) +
  scale_color_manual(name = "", values = palette.basis) + 
  labs(title="Kappa", x = "t", y = '')

gg.beta | gg.kappa


# # summary of intercept and phi:
# res.copy$summary.fixed
# res.single$summary.fixed
# 
# # summary of hyperparameters:
# res.copy$summary.hyperpar
# res.single$summary.hyperpar
# 
# data.frame(rbind(res.copy$summary.random$beta1, res.copy$summary.random$beta2, res.single$summary.random$beta)) %>% 
#   mutate(type = rep(c("copy1","copy2","single"), each = 10)) %>%
#   mutate(beta.true = rep(beta, 3)) %>%
#   ggplot() + geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant, color =type), position=position_dodge(width=0.5)) +
#   geom_point(aes(x = ID, y = beta.true))
#   ggtitle("beta")
#   
# data.frame(rbind(res.copy$summary.random$kappa,  res.single$summary.random$kappa)) %>%
#   mutate(type = rep(c("copy1","single"), each = 10)) %>%
#   mutate(kappa.true = rep(kappa, 2)) %>%
#   ggplot() + geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant, color =type), position=position_dodge(width=0.5)) +
#   geom_point(aes(x = ID, y = kappa.true))
#   ggtitle("kappa")
# 
# data.frame(cbind(copy = res.copy$summary.linear.predictor$mean, single = res.single$summary.linear.predictor$mean)) %>%
#   ggplot() + geom_point(aes(x = copy, y = single)) + ggtitle("Linear predictor")
#   
# 
# gg.beta.true = ggplot(data = obs, aes(x = x, y = beta)) + geom_point(color = "hotpink") + ggtitle("True beta"); gg.beta.true
# gg.kappa.true = ggplot(data = obs, aes(x = t, y = kappa)) + geom_line(color = "hotpink") + ggtitle("True kappa")
# gg.phi.true = ggplot(data = obs, aes(x = t, y = phi)) + geom_line(color = "hotpink") + ggtitle("True phi")
# gg.epsilon.true = ggplot(data = obs, aes(x = x, y = t, fill = epsilon)) + geom_tile() + ggtitle("True epsilon")
# ggplot(data = obs, aes(x = epsilon)) + geom_density()
# 
# gg.beta1.c = ggplot(data = cbind(res.copy$summary.random$beta1, beta.true = beta[res.copy$summary.random$beta1$ID]), aes(x = ID)) + 
#   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
#   geom_point(aes(y = mean), color = "lightskyblue") + 
#   geom_point(aes(y = beta.true), color = "dodgerblue1") + 
#   ggtitle("beta1 coopy")
# 
# gg.beta2.c = ggplot(data = cbind(res.copy$summary.random$beta2, beta.true = beta[res.copy$summary.random$beta2$ID]), aes(x = ID)) + 
#   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
#   geom_point(aes(y = mean), color = "lightskyblue") + 
#   geom_point(aes(y = beta.true), color = "dodgerblue1") + 
#   ggtitle("beta2 copy")
# 
# gg.beta.s = ggplot(data = cbind(res.single$summary.random$beta, beta.true = beta[res.single$summary.random$beta$ID]), aes(x = ID)) + 
#   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
#   geom_point(aes(y = mean), color = "lightskyblue") + 
#   geom_point(aes(y = beta.true), color = "dodgerblue1") + 
#   ggtitle("beta single")
# gg.beta.s
# 
# (gg.beta1.c | gg.beta2.c | gg.beta.s)
# 
# data.kappa.c = cbind(res.copy$summary.random$kappa, kappa.true = kappa[res.copy$summary.random$kappa$ID])
# gg.kappa.c = ggplot(data = data.kappa.c, aes(x = ID)) + 
#   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
#   geom_line(aes(y = mean), color = "lightskyblue") + 
#   geom_line(aes(y = kappa.true), color = "dodgerblue1") + 
#   ggtitle("Kappa copy")
# 
# data.kappa.s = cbind(res.single$summary.random$kappa, kappa.true = kappa[res.single$summary.random$kappa$ID])
# gg.kappa.s = ggplot(data = data.kappa.s, aes(x = ID)) + 
#   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "lightskyblue1") + 
#   geom_line(aes(y = mean), color = "lightskyblue") + 
#   geom_line(aes(y = kappa.true), color = "dodgerblue1") + 
#   ggtitle("Kappa single")
# gg.kappa.s
# 
# (gg.kappa.c | gg.kappa.s)
# 
# 
# gg.epsilon.true
# data.epsilon.c = res.copy$summary.random$epsilon
# data.epsilon.c = cbind(data.epsilon.c, x = apply(data.epsilon.c,1, function(row) row['ID']%/%100 ))
# data.epsilon.c = cbind(data.epsilon.c, t = apply(data.epsilon.c, 1, function(row) row['ID']%%100))
# gg.epsilon.c = ggplot(data = data.epsilon.c, aes(x = x, y = t, fill = mean)) + geom_tile() + ggtitle("Simulated epsilon")
# 
# data.epsilon.s = res.single$summary.random$epsilon
# data.epsilon.s = cbind(data.epsilon.s, x = apply(data.epsilon.s,1, function(row) row['ID']%/%100 ))
# data.epsilon.s = cbind(data.epsilon.s, t = apply(data.epsilon.s, 1, function(row) row['ID']%%100))
# gg.epsilon.s = ggplot(data = data.epsilon.s, aes(x = x, y = t, fill = mean)) + geom_tile() + ggtitle("Simulated epsilon")
# gg.epsilon.s
# 
# (gg.epsilon.true | gg.epsilon.c | gg.epsilon.s)
# 
# # density plots of true and simulated epsilon:
# data.epsilon.density.c = rbind(data.frame(epsilon = data.epsilon.c$mean, sim = "Simulated"), data.frame(epsilon = obs$epsilon, sim = "True values"))
# gg.epsilon.density.c = ggplot(data = data.epsilon.density.c, aes(x = epsilon, color = sim)) + 
#   geom_density() + ggtitle("Model with copied beta")
# 
# data.epsilon.density.s = rbind(data.frame(epsilon = data.epsilon.s$mean, sim = "Simulated"), data.frame(epsilon = obs$epsilon, sim = "True values"))
# gg.epsilon.density.s = ggplot(data = data.epsilon.density.s, aes(x = epsilon, color = sim)) + 
#   geom_density() + ggtitle("Model with single beta")
# 
# (gg.epsilon.density.c | gg.epsilon.density.s)
# 
# #  results of hyperparameters:
# cat("Precision for beta: ")
# cat("True value: ", tau.iid)
# cat("Simulated value copy: ", res.copy$summary.hyperpar$mean[1])
# cat("Simulated value single: ", res.single$summary.hyperpar$mean[1])
# 
# cat("Precision for kappa: \n", "\n Simulated value copy: ",
#     res.copy$summary.hyperpar$mean[2], "\n Simulated value single: ", res.single$summary.hyperpar$mean[2])
# 
# cat("Precision for epsilon: \n", "True value: ", tau.epsilon,"\n Simulated value copy: ",
#     res.copy$summary.hyperpar$mean[3], "\n Simulated value single: ", res.single$summary.hyperpar$mean[3])
# 
# # density plot of true eta and predicted eta:
# eta.sim.c = res.copy$summary.linear.predictor$mean
# eta.sim.c = eta.sim.c
# data.eta.density.c = rbind(data.frame(eta = obs$eta, sim = "True values"),
#                          data.frame(eta = eta.sim.c, sim = "Simulated"))
# gg.eta.density.c = ggplot(data = data.eta.density.c, aes(x = eta, color = sim)) + 
#   geom_density() + ggtitle("Model with copied beta")
# gg.eta.density.c
# 
# eta.sim.s = res.single$summary.linear.predictor$mean
# data.eta.density.s = rbind(data.frame(eta = obs$eta, sim = "True values"),
#                            data.frame(eta = eta.sim.s, sim = "Simulated"))
# gg.eta.density.s = ggplot(data = data.eta.density.s, aes(x = eta, color = sim)) + 
#   geom_density() + ggtitle("Model with single beta")
# gg.eta.density.s
# 
# (gg.eta.density.c | gg.eta.density.s)
# 
# cat("Intercept: ", res$summary.fixed$mean[1] - log(at.risk))
