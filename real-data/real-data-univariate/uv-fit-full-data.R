# fit the LCC and LC models to the full range of lung and stomach cancer data,
# to obtain the 

load("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Workspaces/ws_uv-fit-full-data.RData")

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(lubridate)
library(readxl)

population <- read_excel("population-germany.xlsx", sheet=1,
                         col_names = c("age", "male", "female", "total"),
                         skip = 6)
years <- population %>% select(age) %>% 
  filter(!is.na(as.Date(age, format="%d.%m.%Y", optional = TRUE))) 
population <- population %>% slice(1:1848) %>% 
  mutate(year = rep(as.vector(years$age), each=88)) %>%
  filter(year != age) %>% filter(age != "Insgesamt") %>%
  mutate(age.number =  parse_number(age)) %>%
  mutate(age.number = replace(age.number, age == "unter 1 Jahr", 0)) %>%
  mutate(age.int = 5*(age.number%/%5)) %>%
  mutate(age.int = paste(age.int, "-", age.int + 4)) %>%
  mutate(age.int = replace(age.int, age.int == "85 - 89", "85 +")) %>%
  group_by(age.int, year) %>%
  summarize(total = sum(total), male = sum(male), female = sum(female)) %>%
  mutate(year = format(as.POSIXct(year, format="%d.%m.%Y"), format="%Y")) %>%
  filter(year < 2017)

stomach.cancer  <- read_excel("stomachCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  mutate(age = replace(age, age == "85", "85 +")) %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  mutate(t = as.integer(year)-1999) %>% 
  mutate(age.int = parse_number(age)) %>%
  mutate(x = parse_number(age)%/%5) %>% mutate(x.1 = x) %>%
  mutate(xt = (x*(2016-1998) +t)) %>% mutate(t.1 = t) %>%
  mutate(k = 5 * (max(x) - x) + t) %>%
  mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
  mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t")) %>%
  mutate("mortality rate" = total/total.t)

# read and format lung cancer data 
lung.cancer  <- read_excel("lungCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  mutate(age = replace(age, age == "85", "85 +")) %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  mutate(t = as.integer(year)-1999) %>% 
  mutate(age.int = parse_number(age)) %>%
  mutate(x = parse_number(age)%/%5) %>% mutate(x.1 = x) %>%
  mutate(xt = (x*(2016-1998) +t)) %>% mutate(t.1 = t) %>%
  mutate(k = 5 * (max(x) - x) + t) %>%
  mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
  mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t")) %>%
  mutate("mortality rate" = total/total.t)

pc.prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))

A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$x))) 
e.vec = 1

comp.lc = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) + 
  kappa(t.1, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma(k, model = "rw1", values =  unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xt, model = "iid", hyper = pc.prior)

# Lee-Carter
form.lc = total ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon
likelihood.lc.l = like(formula = form.lc, family = "poisson", data = lung.cancer, E = lung.cancer$total.t )
likelihood.lc.s = like(formula = form.lc, family = "poisson", data = stomach.cancer, E = stomach.cancer$total.t )

# cohort-extended Lee-Carter
form.lcc = total ~ -1 + Int + alpha + beta*phi + beta*kappa +  gamma + epsilon
likelihood.lcc.l = like(formula = form.lcc, family = "poisson", data = lung.cancer, E = lung.cancer$total.t)
likelihood.lcc.s = like(formula = form.lcc, family = "poisson", data = stomach.cancer, E = stomach.cancer$total.t)

# control compute
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

# Lee-Carter lung cancer
res.lc.l = bru(components = comp.lc,
               likelihood.lc.l, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 50
               )) 

# Lee-Carter stomach cancer
res.lc.s = bru(components = comp.lc,
               likelihood.lc.s, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 50
               )) 

# Cohort-extended Lee-Carter lung cancer
res.lcc.l = bru(components = comp.lc,
               likelihood.lcc.l, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 50
               )) 
bru_rerun(res.lcc.l)

# Cohort-extended Lee-Carter stomach cancer
res.lcc.s = bru(components = comp.lc,
                likelihood.lcc.s, 
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1),
                               bru_max_iter = 50
                )) 

# color palette:
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1',
                   '#5d8060', '#D7B36A', '#826133', '#A85150')

#   ----   plot effects from the different runs:   ----   

# Lee-Carter, lung
p.mu.lc.l <- ggplot(data.frame(res.lc.l$marginals.fixed$Int)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lc.l$summary.fixed, aes(xintercept = mean[1]), color = palette.basis[1]) + 
  labs(x = "Value of intercept", y = " ", title = "Intercept")

p.mu.lc.l

p.alpha.lc.l <- ggplot(data = res.lc.l$summary.random$alpha) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "alpha", title = "Alpha")

p.alpha.lc.l

p.beta.lc.l <- ggplot(res.lc.l$summary.random$beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "beta", title = "Beta")

p.beta.lc.l

p.kappa.lc.l <- ggplot(res.lc.l$summary.random$kappa) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa.lc.l

p.phi.lc.l <- ggplot(data.frame(res.lc.l$marginals.fixed$phi)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lc.l$summary.fixed, aes(xintercept = mean[2]), color = palette.basis[1]) + 
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi.lc.l

p.fitted.lc.l <- ggplot(res.lc.l$summary.fitted.values %>%
                          slice(1:324) %>% bind_cols(lung.cancer)) +
  geom_point(aes(x = `mortality rate`, y = mean), color = palette.basis[1]) +
  labs(x = "Observed", y = "Estimated", title = "Mortality rate") + 
  scale_x_continuous(breaks = c(0.000, 0.001, 0.002))
p.fitted.lc.l

p.prec.kappa.lc.l <- ggplot(data.frame(res.lc.l$marginals.hyperpar) %>%
                         filter(Precision.for.kappa.x < 200)) + 
  geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), fill = palette.basis[1], alpha = 0.4) + 
  geom_vline(data = res.lc.l$summary.hyperpar, aes(xintercept = mean[3]), color = palette.basis[1]) + 
  labs(x = "Value of precision", y = " ", title = "Precision, kappa")
p.prec.kappa.lc.l
# 11 values over 200


p.lc.l <- (p.mu.lc.l | p.alpha.lc.l | p.beta.lc.l)/(p.phi.lc.l | p.kappa.lc.l | p.prec.kappa.lc.l | p.fitted.lc.l ) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects for LC-model",
                  subtitle = "Lung cancer")
p.lc.l

ggsave('uv-full-data-lc-l.png',
       plot = p.lc.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# Lee-Carter, stomach

p.mu.lc.s <- ggplot(data.frame(res.lc.s$marginals.fixed$Int)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lc.s$summary.fixed, aes(xintercept = mean[1]), color = palette.basis[1]) + 
  labs(x = "Value of intercept", y = " ", title = "Intercept")

p.mu.lc.s

p.alpha.lc.s <- ggplot(data = res.lc.s$summary.random$alpha) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "alpha", title = "Alpha")

p.alpha.lc.s

p.beta.lc.s <- ggplot(res.lc.s$summary.random$beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "beta", title = "Beta")

p.beta.lc.s

p.kappa.lc.s <- ggplot(res.lc.s$summary.random$kappa) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa.lc.s

p.phi.lc.s <- ggplot(data.frame(res.lc.s$marginals.fixed$phi)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lc.s$summary.fixed, aes(xintercept = mean[2]), color = palette.basis[1]) + 
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi.lc.s

p.fitted.lc.s <- ggplot(res.lc.s$summary.fitted.values %>%
                          slice(1:324) %>% bind_cols(stomach.cancer)) +
  geom_point(aes(x = `mortality rate`, y = mean), color = palette.basis[1]) +
  labs(x = "Observed", y = "Estimated", title = "Mortality rate") + 
  scale_x_continuous(breaks = c(0.000, 0.0015, 0.003))
p.fitted.lc.s

p.prec.kappa.lc.s <- ggplot(data.frame(res.lc.s$marginals.hyperpar) %>%
                              filter(Precision.for.kappa.x < 400)) + 
  geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), fill = palette.basis[1], alpha = 0.4) + 
  geom_vline(data = res.lc.s$summary.hyperpar, aes(xintercept = mean[3]), color = palette.basis[1]) + 
  labs(x = "Value of precision", y = " ", title = "Precision, kappa")
p.prec.kappa.lc.s
# 11 values over 400


p.lc.s <- (p.mu.lc.s | p.alpha.lc.s | p.beta.lc.s)/(p.phi.lc.s | p.kappa.lc.s | p.prec.kappa.lc.s | p.fitted.lc.s ) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects for LC-model",
                  subtitle = "Stomach cancer")
p.lc.s

ggsave('uv-full-data-lc-s.png',
       plot = p.lc.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# LCC-model, lung cancer:

p.mu.lcc.l <- ggplot(data.frame(res.lcc.l$marginals.fixed$Int)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lcc.l$summary.fixed, aes(xintercept = mean[1]), color = palette.basis[1]) + 
  labs(x = "Value of intercept", y = " ", title = "Intercept") +
  scale_x_continuous(breaks = c(-9.25, -8.75))

p.mu.lcc.l

p.alpha.lcc.l <- ggplot(data = res.lcc.l$summary.random$alpha) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "alpha", title = "Alpha")

p.alpha.lcc.l

p.beta.lcc.l <- ggplot(res.lcc.l$summary.random$beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "beta", title = "Beta")

p.beta.lcc.l

p.kappa.lcc.l <- ggplot(res.lcc.l$summary.random$kappa) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa.lcc.l

p.phi.lcc.l <- ggplot(data.frame(res.lcc.l$marginals.fixed$phi)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lcc.l$summary.fixed, aes(xintercept = mean[2]), color = palette.basis[1]) + 
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi.lcc.l

p.gamma.lcc.l <- ggplot(res.lcc.l$summary.random$gamma) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "k", y = "gamma", title = "Gamma")
p.gamma.lcc.l

p.prec.kappa.lcc.l <- ggplot(data.frame(res.lcc.l$marginals.hyperpar) %>%
                              filter(Precision.for.kappa.x < 1000)) + 
  geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), fill = palette.basis[1], alpha = 0.4) + 
  geom_vline(data = res.lcc.l$summary.hyperpar, aes(xintercept = mean[3]), color = palette.basis[1]) + 
  labs(x = "Value of precision", y = " ", title = "Precision, kappa")
p.prec.kappa.lcc.l
# 14 values above 1000

p.fitted.lcc.l <- ggplot(res.lcc.l$summary.fitted.values %>%
                          slice(1:324) %>% bind_cols(lung.cancer)) +
  geom_point(aes(x = `mortality rate`, y = mean), color = palette.basis[1]) +
  labs(x = "Observed", y = "Estimated", title = "Mortality rate") + 
  scale_x_continuous(breaks = c(0.000, 0.001, 0.002))
p.fitted.lcc.l


p.lcc.l <- (p.mu.lcc.l | p.alpha.lcc.l | p.beta.lcc.l | p.phi.lcc.l)/(p.kappa.lcc.l | p.gamma.lcc.l | p.prec.kappa.lcc.l | p.fitted.lcc.l ) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects for LLC-model",
                  subtitle = "Stomach cancer")
p.lcc.l

ggsave('uv-full-data-lcc-l.png',
       plot = p.lcc.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# LCC-model, stomach cancer:

p.mu.lcc.s <- ggplot(data.frame(res.lcc.s$marginals.fixed$Int)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lcc.s$summary.fixed, aes(xintercept = mean[1]), color = palette.basis[1]) + 
  labs(x = "Value of intercept", y = " ", title = "Intercept") +
  scale_x_continuous(breaks = c(-8.85, -8.6, -8.35))

p.mu.lcc.s

p.alpha.lcc.s <- ggplot(data = res.lcc.s$summary.random$alpha) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "alpha", title = "Alpha")

p.alpha.lcc.s

p.beta.lcc.s <- ggplot(res.lcc.s$summary.random$beta) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "x", y = "beta", title = "Beta")

p.beta.lcc.s

p.kappa.lcc.s <- ggplot(res.lcc.s$summary.random$kappa) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`),fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa.lcc.s

p.phi.lcc.s <- ggplot(data.frame(res.lcc.s$marginals.fixed$phi)) + 
  geom_area(aes(x = x, y = y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.lcc.s$summary.fixed, aes(xintercept = mean[2]), color = palette.basis[1]) + 
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi.lcc.s

p.fitted.lcc.s <- ggplot(res.lcc.s$summary.fitted.values %>%
                           slice(1:324) %>% bind_cols(stomach.cancer)) +
  geom_point(aes(x = `mortality rate`, y = mean), color = palette.basis[1]) +
  labs(x = "Observed", y = "Estimated", title = "Mortality rate") + 
  scale_x_continuous(breaks = c(0.000, 0.0015, 0.003))
p.fitted.lcc.s

p.gamma.lcc.s <- ggplot(res.lcc.s$summary.random$gamma) + 
  geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = palette.basis[1], alpha = 0.4) +
  geom_point(aes(x = ID, y = mean), color = palette.basis[1]) + 
  labs(x = "k", y = "gamma", title = "Gamma")
p.gamma.lcc.s

p.prec.kappa.lcc.s <- ggplot(data.frame(res.lcc.s$marginals.hyperpar) %>%
                               filter(Precision.for.kappa.x < 200)) + 
  geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), fill = palette.basis[1], alpha = 0.4) + 
  geom_vline(data = res.lcc.s$summary.hyperpar, aes(xintercept = mean[3]), color = palette.basis[1]) + 
  labs(x = "Value of precision", y = " ", title = "Precision, kappa")
p.prec.kappa.lcc.s
# 11 values above 200


p.lcc.s <- (p.mu.lcc.s | p.alpha.lcc.s | p.beta.lcc.s  | p.phi.lcc.s)/( p.kappa.lcc.s | p.gamma.lcc.s | p.prec.kappa.lcc.s | p.fitted.lcc.s) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects for LLC-model",
                  subtitle = "Stomach cancer")
p.lcc.s

ggsave('uv-full-data-lcc-s.png',
       plot = p.lcc.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)



save.image("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Workspaces/ws_uv-fit-full-data.RData")





