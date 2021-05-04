# compare different configurations of multivariate LCC models

# load workspace
load("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Lung cancer/Workspaces/ws_multivariate-LCC-lung.RData")


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

pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))

# pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.658, 0.5)))
# pc.prior.beta <- list(prec = list(prior = "pc.prec", param = c(0.0837, 0.5)))
# pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))  # perhaps update this after run
# pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.00813, 0.5)))
# pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.0316, 0.5)))

A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$x))) 
e.vec = 1

# common age effect: ABkg
comp.ABkg = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t, model = "linear", prec.linear = 1) +
  kappa0(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  kappa1(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma1(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.ABkg.0 = cases ~ -1 + mu + alpha + beta*phi0 + beta*kappa0 + gamma0 + epsilon
likelihood.ABkg.0 = like(formula = form.ABkg.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.ABkg.1 = cases ~ -1 + mu + alpha + beta*phi1 + beta*kappa1 + gamma1 + epsilon
likelihood.ABkg.1 = like(formula = form.ABkg.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res.ABkg = bru(components = comp.ABkg,
                  likelihood.ABkg.0,
                  likelihood.ABkg.1,
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1),
                                 bru_max_iter = 100
                  )) 

# constructed to ensure that observation data is in the same order as prediction data
lung.cancer.0 <- lung.cancer %>% filter(s == 0)
lung.cancer.1 <- lung.cancer %>% filter(s == 1)

# combine predictions and observations into one dataframe. 
data.pred.ABkg <- res.ABkg$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Keep period effects common : abKg   ----:

comp.abKg = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  alpha1(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta0(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  beta1(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi(t, model = "linear", prec.linear = 1) +
  kappa(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma1(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.abKg.0 = cases ~ -1 + mu + alpha0 + beta0*phi + beta0*kappa + gamma0 + epsilon
likelihood.abKg.0 = like(formula = form.abKg.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.abKg.1 = cases ~ -1 + mu + alpha1 + beta1*phi + beta1*kappa + gamma1 + epsilon
likelihood.abKg.1 = like(formula = form.abKg.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.abKg = bru(components = comp.abKg,
               likelihood.abKg.0,
               likelihood.abKg.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 100
               )) 

# combine predictions and observations into one dataframe. 
data.pred.abKg <- res.abKg$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Keep cohort effects common : abkG   ----:
# common age effect:
comp.abkG = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  alpha1(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta0(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  beta1(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t, model = "linear", prec.linear = 1) +
  kappa0(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  kappa1(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior) 

# define two different likelihoods and formulas, one for male and one for female:
form.abkG.0 = cases ~ -1 + mu + alpha0 + beta0*phi0 + beta0*kappa0 + gamma + epsilon
likelihood.abkG.0 = like(formula = form.abkG.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.abkG.1 = cases ~ -1 + mu + alpha1 + beta1*phi1 + beta1*kappa1 + gamma + epsilon
likelihood.abkG.1 = like(formula = form.abkG.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.abkG = bru(components = comp.abkG,
               likelihood.abkG.0,
               likelihood.abkG.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 100
               )) 

# combine predictions and observations into one dataframe. 
data.pred.abkG <- res.abkG$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   No common effects : abkg   ----:
comp.abkg = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  alpha1(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta0(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  beta1(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t, model = "linear", prec.linear = 1) +
  kappa0(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  kappa1(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma1(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.abkg.0 = cases ~ -1 + mu + alpha0 + beta0*phi0 + beta0*kappa0 + gamma0 + epsilon
likelihood.abkg.0 = like(formula = form.abkg.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.abkg.1 = cases ~ -1 + mu + alpha1 + beta1*phi1 + beta1*kappa1 + gamma1 + epsilon
likelihood.abkg.1 = like(formula = form.abkg.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.abkg = bru(components = comp.abkg,
               likelihood.abkg.0,
               likelihood.abkg.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 100
               )) 

# combine predictions and observations into one dataframe. 
data.pred.abkg <- res.abkg$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   All effects common : ABKG   ----

comp.ABKG = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi(t, model = "linear", prec.linear = 1) +
  kappa(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.ABKG = cases ~ -1 + mu + alpha + beta*phi + beta*kappa + gamma + epsilon
likelihood.ABKG = like(formula = form.ABKG, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$population)

#   NOTE: This does not converge!!! reformulate? 
res.ABKG = bru(components = comp.ABKG,
               likelihood.ABKG,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 50
               )) 

# combine predictions and observations into one dataframe. 
data.pred.ABKG <- res.ABKG$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(lung.cancer) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

#   ----   Shared period and cohort effect : abKG   ----

comp.abKG = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  alpha1(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta0(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  beta1(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi(t, model = "linear", prec.linear = 1) +
  kappa(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.abKG.0 = cases ~ -1 + mu + alpha0 + beta0*phi + beta0*kappa + gamma + epsilon
likelihood.abKG.0 = like(formula = form.abKG.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.abKG.1 = cases ~ -1 + mu + alpha1 + beta1*phi + beta1*kappa + gamma + epsilon
likelihood.abKG.1 = like(formula = form.abKG.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

#  Note: This also really struggles to converge!!! 

res.abKG = bru(components = comp.abKG,
               likelihood.abKG.0,
               likelihood.abKG.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 50
               )) 

# combine predictions and observations into one dataframe. 
data.pred.abKG <- res.abKG$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Common age and period effects : ABKg   ----


comp.ABKg = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi(t, model = "linear", prec.linear = 1) +
  kappa(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma1(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.ABKg.0 = cases ~ -1 + mu + alpha + beta*phi + beta*kappa + gamma0 + epsilon
likelihood.ABKg.0 = like(formula = form.ABKg.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.ABKg.1 = cases ~ -1 + mu + alpha + beta*phi + beta*kappa + gamma1 + epsilon
likelihood.ABKg.1 = like(formula = form.ABKg.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

res.ABKg = bru(components = comp.ABKg,
               likelihood.ABKg.0,
               likelihood.ABKg.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 100
               )) 

# combine predictions and observations into one dataframe. 
data.pred.ABKg <- res.ABKg$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Common age and cohort effects : ABkG   ----


comp.ABkG = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t, model = "linear", prec.linear = 1) +
  kappa0(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  kappa1(t.c, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.ABkG.0 = cases ~ -1 + mu + alpha + beta*phi0 + beta*kappa0 + gamma + epsilon
likelihood.ABkG.0 = like(formula = form.ABkG.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.ABkG.1 = cases ~ -1 + mu + alpha + beta*phi1 + beta*kappa1 + gamma + epsilon
likelihood.ABkG.1 = like(formula = form.ABkG.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

# note: inlabru does not crash, but very far from convergence...
res.ABkG = bru(components = comp.ABkG,
               likelihood.ABkG.0,
               likelihood.ABkG.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),
                              bru_max_iter = 100
               )) 

# combine predictions and observations into one dataframe. 
data.pred.ABkG <- res.ABkG$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

cat("CPO, LCC model: "); -1*mean(log(res.lc.n$cpo$cpo[!is.na(res.lc.n$cpo$cpo)]))


#   ----   TODO: Compare all methods   ----

data.pred <- rbind(data.pred.abkg, data.pred.abKg, data.pred.abkG, data.pred.ABkg, data.pred.ABKg, data.pred.ABkG, data.pred.ABKG) %>%
  mutate("method" = rep(c("No common", "Common period", "Common cohort", "Common age", "Common age & period", "Common age & cohort", "All common"), each = 648))

# display four significant digits 
options(pillar.sigfig = 4)

pred.statistics.included <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 included");cat("\n Lung cancer data - LCC: ");pred.statistics.included

pred.statistics.cutoff <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 omitted");cat("\n Lung cancer data - LCC: ");pred.statistics.cutoff

# plot:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#1C84BB', '#A85150', '#DA871F',
                   '#4C7246', '#D7B36A', '#FB5E4E', '#696B8D', '#76A7A6', '#826133')

gg.pred <- ggplot(data.pred %>%
                    filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                    filter(method %in% c("No common", "Common period")),
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
  labs(title = "Multivariate LCC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred

ggsave('multivariate-LCC-by-age-lung.png',
       plot = gg.pred,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)


# plot cohortwise

gg.pred.cohort <- ggplot(data.pred %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                           filter(method %in% c("No common", "Common period")),
                         aes(x = k)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  labs(title = "Multivariate LCC models - lung cancer", x = "Cohort", y = "Mortality rate") + 
  facet_wrap(~year)

gg.pred.cohort

ggsave('multivariate-LCC-by-cohort-lung.png',
       plot = gg.pred.cohort,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot along years - for different age groups:
gg.pred.period <- ggplot(data.pred %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                           filter(x > 5) %>%
                           filter(method %in% c("No common", "Common period")),
                         aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = interaction(method, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1, group = interaction(method, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  labs(title = "Multivariate LCC models - lung cancer", x = "Year", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)

gg.pred.period

ggsave('multivariate-LCC-by-period-lung.png',
       plot = gg.pred.period,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot effects of abkg and abKg models:

# abkg - no common effects:

p.mu.abkg <- ggplot(data.frame(res.abkg$marginals.random$mu)) + 
  geom_area(aes(x = index.1.x, y = index.1.y, fill = "Male"), alpha = 0.4) + 
  geom_area(aes(x = index.2.x, y = index.2.y, fill = "Female"), alpha = 0.4) + 
  geom_vline(data = res.abkg$summary.random$mu, aes(xintercept = mean[1], color = "Male")) + 
  geom_vline(data = res.abkg$summary.random$mu, aes(xintercept = mean[2], color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of intercept", y = " ", title = "Mu")

p.mu.abkg

p.alpha.abkg <- ggplot() + 
  geom_ribbon(data = res.abkg$summary.random$alpha0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abkg$summary.random$alpha1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abkg$summary.random$alpha0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abkg$summary.random$alpha1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "x", y = "alpha", title = "Alpha")

p.alpha.abkg

p.beta.abkg <- ggplot() + 
  geom_ribbon(data = res.abkg$summary.random$beta0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abkg$summary.random$beta1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abkg$summary.random$beta0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abkg$summary.random$beta1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "x", y = "beta", title = "Beta")

p.beta.abkg

p.kappa.abkg <- ggplot() + 
  geom_ribbon(data = res.abkg$summary.random$kappa0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abkg$summary.random$kappa1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abkg$summary.random$kappa0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abkg$summary.random$kappa1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  geom_vline(data = res.abkg$summary.random$kappa1, aes(xintercept = 9), color = palette.basis[5]) +
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa.abkg

p.phi.abkg <- ggplot(data.frame(res.abkg$marginals.fixed)) + 
  geom_area(aes(x = phi0.x, y = phi0.y, fill = "Male"), alpha = 0.4) + 
  geom_area(aes(x = phi1.x, y = phi1.y, fill = "Female"), alpha = 0.4) + 
  geom_vline(data = res.abkg$summary.fixed, aes(xintercept = mean[1], color = "Male")) + 
  geom_vline(data = res.abkg$summary.fixed, aes(xintercept = mean[2], color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of Phi", y = " ", title = "Phi")

p.phi.abkg

p.gamma.abkg <- ggplot() +
  geom_ribbon(data = res.abkg$summary.random$gamma0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abkg$summary.random$gamma1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abkg$summary.random$gamma0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abkg$summary.random$gamma1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "t", y = "gamma", title = "Gamma")

p.gamma.abkg

p.abkg <- (p.mu.abkg | p.alpha.abkg | p.beta.abkg)/(p.phi.abkg | p.kappa.abkg | p.gamma.abkg) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects for multivariate LCC model with no common effects",
                  subtitle = "Lung cancer")

p.abkg

ggsave('effects-LCC-no-common-lung.png',
       plot = p.abkg,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# abKg - common period effects:

p.mu.abKg <- ggplot(data.frame(res.abKg$marginals.random$mu)) + 
  geom_area(aes(x = index.1.x, y = index.1.y, fill = "Male"), alpha = 0.4) + 
  geom_area(aes(x = index.2.x, y = index.2.y, fill = "Female"), alpha = 0.4) + 
  geom_vline(data = res.abKg$summary.random$mu, aes(xintercept = mean[1], color = "Male")) + 
  geom_vline(data = res.abKg$summary.random$mu, aes(xintercept = mean[2], color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of intercept", y = " ", title = "Mu")

p.mu.abKg

p.alpha.abKg <- ggplot() + 
  geom_ribbon(data = res.abKg$summary.random$alpha0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abKg$summary.random$alpha1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abKg$summary.random$alpha0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abKg$summary.random$alpha1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "x", y = "alpha", title = "Alpha")

p.alpha.abKg

p.beta.abKg <- ggplot() + 
  geom_ribbon(data = res.abKg$summary.random$beta0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abKg$summary.random$beta1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abKg$summary.random$beta0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abKg$summary.random$beta1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "x", y = "beta", title = "Beta")

p.beta.abKg

p.kappa.abKg <- ggplot() + 
  geom_ribbon(data = res.abKg$summary.random$kappa, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Common'), alpha = 0.4) +
  #geom_ribbon(data = res.abkg$summary.random$kappa1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abKg$summary.random$kappa, aes(x = ID, y = mean, color = "Common")) + 
  #geom_point(data = res.abkg$summary.random$kappa1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis[3:length(palette.basis)]) + 
  scale_fill_manual(name = " ", values = palette.basis[3:length(palette.basis)]) +
  geom_vline(data = res.abKg$summary.random$kappa, aes(xintercept = 9), color = palette.basis[5]) +
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa.abKg

p.phi.abKg <- ggplot(data.frame(res.abKg$marginals.fixed)) + 
  geom_area(aes(x = phi.x, y = phi.y, fill = "Common"), alpha = 0.4) + 
  #geom_area(aes(x = phi1.x, y = phi1.y, fill = "Female"), alpha = 0.4) + 
  geom_vline(data = res.abKg$summary.fixed, aes(xintercept = mean, color = "Male")) + 
  #geom_vline(data = res.abKg$summary.fixed, aes(xintercept = mean[2], color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis[3:length(palette.basis)]) + 
  scale_fill_manual(name = " ", values = palette.basis[3:length(palette.basis)]) +
  labs(x = "Value of Phi", y = " ", title = "Phi")

p.phi.abKg

p.gamma.abKg <- ggplot() +
  geom_ribbon(data = res.abKg$summary.random$gamma0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
  geom_ribbon(data = res.abKg$summary.random$gamma1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
  geom_point(data = res.abKg$summary.random$gamma0, aes(x = ID, y = mean, color = "Male")) + 
  geom_point(data = res.abKg$summary.random$gamma1, aes(x = ID, y = mean, color = "Female")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "t", y = "gamma", title = "Gamma")

p.gamma.abKg

p.abKg <- (p.mu.abKg | p.alpha.abKg | p.beta.abKg)/(p.phi.abKg | p.kappa.abKg | p.gamma.abKg) +
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Estimated random effects for multivariate LCC model with common period effects",
                  subtitle = "Lung cancer")

p.abkg

ggsave('effects-LCC-common-period-lung.png',
       plot = p.abKg,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)


# save workspace
save.image("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Lung cancer/Workspaces/ws_multivariate-LCC-lung.RData")
