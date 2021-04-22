# multivariate LCC - first attempt with the LCC model

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
  filter(year < 2017) %>%
  pivot_longer(c(total, male, female), names_to = "sex", values_to="population")

lung.cancer  <- read_excel("lungCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  mutate(age = replace(age, age == "85", "85 +")) %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  pivot_longer(c(total, male, female), names_to = "sex", values_to = "cases") %>%
  mutate(t = as.integer(year)-1999) %>% 
  mutate(age.int = parse_number(age)) %>%
  mutate(x = parse_number(age)%/%5) %>% 
  mutate(k = 5 * (max(x) - x) + t) %>%
  mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
  mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int", "sex" = "sex"), suffix = c("", ".t")) %>%
  filter(sex != "total") %>%
  #filter(sex == "total") %>%
  mutate(s = recode(sex, male = 0, female = 1))  %>%
  mutate(xts = (x*(2016-1998) + t + s*(max(x)*(2016-1998) + max(t)))) %>% 
  mutate("mortality rate" = cases/population)

#add helper-factors for multivariate analysis:
lung.cancer <- lung.cancer %>%
  mutate(x.c = x) %>% mutate(x.c2 = x) %>% mutate(x.c3 = x) %>%
  mutate(t.c = t) %>% mutate(t.c2 = t) %>% mutate(t.c3 = t) %>%
  mutate(k.c = k) %>% mutate(k.c2 = k) #%>%

# set data for years 2008-2016 ro NA
# these will be used as the observation data when running inlabru. 
lung.cancer.until2007 <- lung.cancer %>% 
  mutate(cases = replace(cases, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) 

lung.cancer.until2007.0 <- lung.cancer.until2007 %>% filter(s == 0)
lung.cancer.until2007.1 <- lung.cancer.until2007 %>% filter(s == 1)

pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.658, 0.5)))
pc.prior.beta <- list(prec = list(prior = "pc.prec", param = c(0.0837, 0.5)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))  # perhaps update this after run
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.00813, 0.5)))
pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.0316, 0.5)))

#A.mat.0 = matrix(1, nrow = 1, ncol = length(unique(lung.cancer.until2007.0$x))) 
#A.mat.1 = matrix(1, nrow = 1, ncol = length(unique(lung.cancer.until2007.1$x))) 
A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$x))) 
e.vec = 1

# common age effect:
comp.lc = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  beta(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t.c1, model = "linear", prec.linear = 1) +
  kappa0(t.c2, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  kappa1(t.c3, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  gamma1(k.c, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon)

# define two different likelihoods and formulas, one for male and one for female:
form.lc.0 = cases ~ -1 + mu + alpha + beta*phi0 + beta*kappa0 + gamma0 + epsilon0
likelihood.lc.0 = like(formula = form.lc.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.lc.1 = cases ~ -1 + mu + alpha + beta*phi1 + beta*kappa1 + gamma1 + epsilon1
likelihood.lc.1 = like(formula = form.lc.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res.lc.a = bru(components = comp.lc,
                  likelihood.lc.0,
                  likelihood.lc.1,
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
res.lc.a = bru_rerun(res.lc.a)

# comparable process with APC model - shared age effect:
pc.prior.rho <- list(prec = list(prior = "pc.prec"))
pc.prior.phi <- list(prec = list(prior = "pc.prec"))
pc.prior.epsilon.apc <- list(prec = list(prior = "pc.prec"))
pc.prior.psi <- list(prec = list(prior = "pc.prec"))

# common age effect:
comp.apc = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi0(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  phi1(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon.apc) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon.apc)

# define two different likelihoods and formulas, one for male and one for female:
form.apc.0 = cases ~ -1 + mu + rho + phi0 + psi0 + epsilon0
likelihood.apc.0 = like(formula = form.apc.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.apc.1 = cases ~ -1 + mu + rho + phi1 + psi1 + epsilon1
likelihood.apc.1 = like(formula = form.apc.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.apc.a = bru(components = comp.apc,
               likelihood.apc.0,
               likelihood.apc.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.apc.a = bru_rerun(res.apc.a)

summary(res.apc.a)

# constructed to ensure that observation data is in the same order as prediction data
lung.cancer.0 <- lung.cancer %>% filter(s == 0)
lung.cancer.1 <- lung.cancer %>% filter(s == 1)

# combine predictions and observations into one dataframe. 
data.pred.lc <- res.lc.a$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.apc <- res.apc.a$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred <- rbind(data.pred.lc, data.pred.apc) %>%
  mutate("method" = rep(c("LCC", "APC rw2"), each = 648))

pred.statistics.cutoff <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

# calculate and print cpo for the two:
cat("CPO, LCC model: "); -1*mean(log(res.lc.a$cpo$cpo[!is.na(res.lc.a$cpo$cpo)]))
cat("CPO, APC model: "); -1*mean(log(res.apc.a$cpo$cpo[!is.na(res.apc.a$cpo$cpo)]))

cat("\n Age <= 5 omitted");cat("\n Lung cancer data: ");pred.statistics.cutoff

# plot:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')

gg.pred <- ggplot(data.pred %>%
                      filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                    aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Shared age effect", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred

#   ----   Keep period effects common   ----:
# common age effect:
comp.lc.p = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  alpha1(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  beta0(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
  beta1(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
  phi(t, model = "linear", prec.linear = 1) +
  kappa(t.c2, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  gamma1(k.c, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon)

# define two different likelihoods and formulas, one for male and one for female:
form.lc.p.0 = cases ~ -1 + mu + alpha0 + beta0*phi + beta0*kappa + gamma0 + epsilon0
likelihood.lc.p.0 = like(formula = form.lc.p.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.lc.p.1 = cases ~ -1 + mu + alpha1 + beta1*phi + beta1*kappa + gamma1 + epsilon1
likelihood.lc.p.1 = like(formula = form.lc.p.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.lc.p = bru(components = comp.lc.p,
               likelihood.lc.p.0,
               likelihood.lc.p.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.lc.p = bru_rerun(res.lc.p)

# comparable process with APC model - shared period effect:

# common period effect:
comp.apc.p = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  rho1(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon.apc) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon.apc)

# define two different likelihoods and formulas, one for male and one for female:
form.apc.p.0 = cases ~ -1 + mu + rho0 + phi + psi0 + epsilon0
likelihood.apc.p.0 = like(formula = form.apc.p.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.apc.p.1 = cases ~ -1 + mu + rho1 + phi + psi1 + epsilon1
likelihood.apc.p.1 = like(formula = form.apc.p.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.apc.p = bru(components = comp.apc.p,
                likelihood.apc.p.0,
                likelihood.apc.p.1,
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1)
                )) 
res.apc.p = bru_rerun(res.apc.p)

# combine predictions and observations into one dataframe. 
data.pred.lc.p <- res.lc.p$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.apc.p <- res.apc.p$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.p <- rbind(data.pred.lc.p, data.pred.apc.p) %>%
  mutate("method" = rep(c("LCC", "APC rw2"), each = 648))

pred.statistics.cutoff.p <- data.pred.p %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

# calculate and print cpo for the two:
cat("CPO, LCC model: "); -1*mean(log(res.lc.p$cpo$cpo[!is.na(res.lc.p$cpo$cpo)]))
cat("CPO, APC model: "); -1*mean(log(res.apc.p$cpo$cpo[!is.na(res.apc.p$cpo$cpo)]))

cat("\n Age <= 5 omitted");cat("\n Lung cancer data: ");pred.statistics.cutoff.p

# plot:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')

gg.pred.p <- ggplot(data.pred.p %>%
                    filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                  aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Shared period effect", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred.p

# plot cohortwise

ggplot(data.pred.p %>%
         filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
       aes(x = k)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Shared age effect", x = "Cohort", y = "Mortality rate") + 
  facet_wrap(~year)

#   ----   Keep cohort effects common   ----:
# common age effect:
comp.lc.c = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  alpha1(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  beta0(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
  beta1(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t, model = "linear", prec.linear = 1) +
  kappa0(t.c2, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  kappa1(t.c2, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  gamma(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon)

# define two different likelihoods and formulas, one for male and one for female:
form.lc.c.0 = cases ~ -1 + mu + alpha0 + beta0*phi0 + beta0*kappa0 + gamma + epsilon0
likelihood.lc.c.0 = like(formula = form.lc.c.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.lc.c.1 = cases ~ -1 + mu + alpha1 + beta1*phi1 + beta1*kappa1 + gamma + epsilon1
likelihood.lc.c.1 = like(formula = form.lc.c.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.lc.c = bru(components = comp.lc.c,
               likelihood.lc.c.0,
               likelihood.lc.c.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.lc.c = bru_rerun(res.lc.c)

# comparable process with APC model - shared period effect:

# common period effect:
comp.apc.c = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  rho1(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi0(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  phi1(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon.apc) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon.apc)

# define two different likelihoods and formulas, one for male and one for female:
form.apc.c.0 = cases ~ -1 + mu + rho0 + phi0 + psi + epsilon0
likelihood.apc.c.0 = like(formula = form.apc.c.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.apc.c.1 = cases ~ -1 + mu + rho1 + phi1 + psi + epsilon1
likelihood.apc.c.1 = like(formula = form.apc.c.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.apc.c = bru(components = comp.apc.c,
                likelihood.apc.c.0,
                likelihood.apc.c.1,
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1)
                )) 
res.apc.c = bru_rerun(res.apc.c)

# combine predictions and observations into one dataframe. 
data.pred.lc.c <- res.lc.c$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.apc.c <- res.apc.c$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.c <- rbind(data.pred.lc.c, data.pred.apc.c) %>%
  mutate("method" = rep(c("LCC", "APC rw2"), each = 648))

pred.statistics.cutoff.c <- data.pred.c %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

# calculate and print cpo for the two:
cat("CPO, LCC model: "); -1*mean(log(res.lc.c$cpo$cpo[!is.na(res.lc.c$cpo$cpo)]))
cat("CPO, APC model: "); -1*mean(log(res.apc.c$cpo$cpo[!is.na(res.apc.c$cpo$cpo)]))

cat("\n Age <= 5 omitted");cat("\n Lung cancer data: ");pred.statistics.cutoff.c

# plot:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')

gg.pred.c <- ggplot(data.pred.c %>%
                      filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                    aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Shared cohort effect", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred.c

# plot cohortwise

ggplot(data.pred.p %>%
         filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
       aes(x = k)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Shared age effect", x = "Cohort", y = "Mortality rate") + 
  facet_wrap(~year)
