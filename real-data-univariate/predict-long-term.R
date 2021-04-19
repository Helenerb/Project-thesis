# script for comparison of long-term predictions using different APC and LC models
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
  mutate(cohort = 5 * (max(x) - x) + t) %>%
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
  mutate(cohort = 5 * (max(x) - x) + t) %>%
  mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
  mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t")) %>%
  mutate("mortality rate" = total/total.t)

# set data for years 2008-2016 ro NA
# these will be used as the observation data when running inlabru. 
lung.cancer.until2007 <- lung.cancer %>% 
  mutate(total = replace(total, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(male = replace(male, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA))

stomach.cancer.until2007 <- stomach.cancer %>% 
  mutate(total = replace(total, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(male = replace(male, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA))

# for APC model:
# based on hyperparameters from previous runs
#pc.prior.rho <- list(prec = list(prior = "pc.prec", param = c(0.172, 0.5)))
#pc.prior.phi <- list(prec = list(prior = "pc.prec", param = c(0.00430, 0.5)))
# pc.prior.epsilon.apc <- list(prec = list(prior = "pc.prec", param = c(0.00929, 0.5)))
# pc.prior.psi <- list(prec = list(prior = "pc.prec", param = c(0.0154, 0.5)))
pc.prior.rho <- list(prec = list(prior = "pc.prec"))
pc.prior.phi <- list(prec = list(prior = "pc.prec"))
pc.prior.epsilon.apc <- list(prec = list(prior = "pc.prec"))
pc.prior.psi <- list(prec = list(prior = "pc.prec"))

# for LC model:
#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$x))) 
e.vec = 1

# priors based on hyperparameters from previous runs
pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.658, 0.5)))
pc.prior.beta <- list(prec = list(prior = "pc.prec", param = c(0.0837, 0.5)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))  # perhaps update this after run
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.00813, 0.5)))
pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.0316, 0.5)))
# pc.prior.alpha <- list(prec = list(prior = "pc.prec"))
# pc.prior.beta <- list(prec = list(prior = "pc.prec"))
# pc.prior.kappa <- list(prec = list(prior = "pc.prec"))
# pc.prior.gamma <- list(prec = list(prior = "pc.prec"))
# pc.prior.epsilon <- list(prec = list(prior = "pc.prec"))

# note: this way of defining the x, t and k - values are only possible knowing that 
# we do not have any missing data points in stomach.cancer (we can, however, have zero-data.)
# also, we know that the lung cancer data set and the stomach cancer data set contain
# the same covariate values for x, t, cohort. 
comp.apc.2 = ~ -1 + 
  Int(1) + 
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi(cohort, model = "rw2", values = unique(stomach.cancer$cohort), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon.apc)

comp.apc.1 = ~ -1 + 
  Int(1) + 
  rho(x, model = "rw1", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi(t, model = "rw1", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi(cohort, model = "rw1", values = unique(stomach.cancer$cohort), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon.apc)

# define the formula for the predictor for the APC model:
form.apc = total ~ -1 + Int + rho + phi + psi + epsilon

# add data as lung.cancer and offset as population$total, which is the number of people at-risk
likelihood.apc.l = like(formula = form.apc, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t)
likelihood.apc.s = like(formula = form.apc, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t)

# define the formula for the predictor for the AC model - APC (and LC) model without period effect
form.ac = total ~ -1 + Int + rho + psi + epsilon
likelihood.ac.l = like(formula = form.ac, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t)
likelihood.ac.s = like(formula = form.ac, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t)

# test also a model with only the age-effect - to ensure that period does influence the data
form.a = total ~ -1 + Int + rho + epsilon
likelihood.a.l = like(formula = form.a, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t)
likelihood.a.s = like(formula = form.a, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t)

# LC:
comp.lc = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) + 
  kappa(t.1, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  gamma(cohort, model = "rw1", values =  unique(lung.cancer$cohort), constr = TRUE, hyper = pc.prior.gamma) +
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

form.lc = total ~ -1 + Int + alpha + beta*phi + beta*kappa +  gamma + epsilon
likelihood.lc.l = like(formula = form.lc, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t)
likelihood.lc.s = like(formula = form.lc, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t)

form.lc.lin = total ~ -1 + Int + alpha + beta*phi + gamma + epsilon
likelihood.lc.lin.l = like(formula = form.lc.lin, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t )
likelihood.lc.lin.s = like(formula = form.lc.lin, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t )

form.lc.basic = total ~ -1 + Int + alpha + beta*phi + epsilon
likelihood.lc.basic.l = like(formula = form.lc.basic, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t )
likelihood.lc.basic.s = like(formula = form.lc.basic, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t )


# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

# run inlabru for all models:

# apc model, rw1, lung cancer
res.apc.1.l = bru(components = comp.apc.1,
                likelihood.apc.l, 
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1)
                )) 
#res.apc.1.l = bru_rerun(res.apc.1.l)

# apc model, rw1, stomach cancer
res.apc.1.s = bru(components = comp.apc.1,
                  likelihood.apc.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
res.apc.1.s = bru_rerun(res.apc.1.s)

# apc model, rw2, lung cancer
res.apc.2.l = bru(components = comp.apc.2,
                  likelihood.apc.l, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
#res.apc.2.l = bru_rerun(res.apc.2.l)

# apc model, rw2, stomach cancer
res.apc.2.s = bru(components = comp.apc.2,
                  likelihood.apc.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
#res.apc.s.s = bru_rerun(res.apc.s.s)

# ac model, rw1, lung cancer
res.ac.l = bru(components = comp.apc.1,
                  likelihood.ac.l, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
#res.apc.2.l = bru_rerun(res.apc.2.l)

# ac model, rw1, stomach cancer
res.ac.s = bru(components = comp.apc.1,
                  likelihood.ac.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 

# a model, rw1, lung cancer
res.a.l = bru(components = comp.apc.1,
               likelihood.a.l, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
#res.apc.2.l = bru_rerun(res.apc.2.l)

# a model, rw1, stomach cancer
res.a.s = bru(components = comp.apc.1,
               likelihood.a.s, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 

# lc model, both period effects and cohort, lung cancer
res.lc.l = bru(components = comp.lc,
              likelihood.lc.l, 
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1)
              )) 
res.lc.l = bru_rerun(res.lc.l)

# lc model, both period effects and cohort, stomach cancer
res.lc.s = bru(components = comp.lc,
              likelihood.lc.s, 
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1)
              )) 
res.lc.s = bru_rerun(res.lc.s)

# lc model, only linear period effect and cohort, lung cancer
res.lc.lin.l = bru(components = comp.lc,
               likelihood.lc.lin.l, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.lc.lin.l = bru_rerun(res.lc.lin.l)

# lc model, only linear period effect and cohort, stomach cancer
res.lc.lin.s = bru(components = comp.lc,
               likelihood.lc.lin.s, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.lc.lin.s = bru_rerun(res.lc.lin.s)

# lc model, only linear period effect and cohort, lung cancer
res.lc.basic.l = bru(components = comp.lc,
                   likelihood.lc.basic.l, 
                   options = list(verbose = F,
                                  bru_verbose = 1, 
                                  num.threads = "1:1",
                                  control.compute = c.c,
                                  control.predictor = list(link = 1)
                   )) 
res.lc.basic.l = bru_rerun(res.lc.basic.l)

# lc model, only linear period effect and cohort, stomach cancer
res.lc.basic.s = bru(components = comp.lc,
                   likelihood.lc.basic.s, 
                   options = list(verbose = F,
                                  bru_verbose = 1, 
                                  num.threads = "1:1",
                                  control.compute = c.c,
                                  control.predictor = list(link = 1)
                   )) 
res.lc.basic.s = bru_rerun(res.lc.basic.s)

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')
palette.light <- c('#ABC9E6', '#F3DC90', '#86A46F', '#40BBFD', '#C38281')


data.apc.1.l <- res.apc.1.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.apc.1.s <- res.apc.1.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.apc.2.l <- res.apc.2.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.apc.2.s <- res.apc.2.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.ac.l <- res.ac.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.ac.s <- res.ac.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.a.l <- res.a.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.a.s <- res.a.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.l <- res.lc.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.s <- res.lc.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.lin.l <- res.lc.lin.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.lin.s <- res.lc.lin.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.basic.l <- res.lc.basic.l$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.basic.s <- res.lc.basic.s$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(stomach.cancer %>% select(- c("x.1", "t.1", "xt")))

data.pred.l <- rbind(data.apc.1.l, data.apc.2.l, data.ac.l, data.a.l, data.lc.l, data.lc.lin.l, data.lc.basic.l) %>%
  mutate("method" = rep(c("APC rw1", "APC rw2", "AC", "A", "LCC", "LCC linear", "LC"), each = 324)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


data.pred.s <- rbind(data.apc.1.s, data.apc.2.s, data.ac.s, data.a.s, data.lc.s, data.lc.lin.s, data.lc.basic.s) %>%
  mutate("method" = rep(c("APC rw1", "APC rw2", "AC", "A", "LCC", "LCC linear", "LC"), each = 324)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

#   ----   Objective measures of prediction performance   ----
pred.statistics.l <- data.pred.l %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

pred.statistics.s <- data.pred.s %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))
cat("\n Lung cancer data: ");pred.statistics.l;cat("\n Stomach cancer data: ");pred.statistics.s

# compute statistics with cutoff at age group 5
pred.statistics.cutoff.l <- data.pred.l %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

pred.statistics.cutoff.s <- data.pred.s %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))
cat("\n Age <= 5 omitted");cat("\n Lung cancer data: ");pred.statistics.cutoff.l;cat("\n Stomach cancer data: ");pred.statistics.cutoff.s

#   ----   plot APC and LCC models for comparison   ----

gg.pred.s <- ggplot(data.pred.s %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                           filter(method %in% c("APC rw2", "LCC")),
                         aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Stomach cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred.s

gg.pred.l <- ggplot(data.pred.l %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                           filter(method %in% c("APC rw2", "LCC")),
                         aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred.l

