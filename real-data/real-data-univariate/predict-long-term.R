# script for comparison of long-term predictions using different APC and LC models

# load work space:
load("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Workspaces/ws_predict-long-term.RData")

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
# pc.prior.rho <- list(prec = list(prior = "pc.prec"))
# pc.prior.phi <- list(prec = list(prior = "pc.prec"))
# pc.prior.epsilon.apc <- list(prec = list(prior = "pc.prec"))
# pc.prior.psi <- list(prec = list(prior = "pc.prec"))

pc.prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))

# for LC model:
#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$x))) 
e.vec = 1

# priors based on hyperparameters from previous runs
# pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.658, 0.5)))
# pc.prior.beta <- list(prec = list(prior = "pc.prec", param = c(0.0837, 0.5)))
# pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))  # perhaps update this after run
# pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.00813, 0.5)))
# pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.0316, 0.5)))
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
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xt, model = "iid", hyper = pc.prior)

comp.apc.1 = ~ -1 + 
  Int(1) + 
  rho(x, model = "rw1", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "rw1", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi(k, model = "rw1", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xt, model = "iid", hyper = pc.prior)

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
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) + 
  kappa(t.1, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  gamma(k, model = "rw1", values =  unique(lung.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xt, model = "iid", hyper = pc.prior)

form.lc = total ~ -1 + Int + alpha + beta*phi + beta*kappa +  gamma + epsilon
likelihood.lc.l = like(formula = form.lc, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t)
likelihood.lc.s = like(formula = form.lc, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t)

form.lc.lin = total ~ -1 + Int + alpha + beta*phi + gamma + epsilon
likelihood.lc.lin.l = like(formula = form.lc.lin, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$total.t )
likelihood.lc.lin.s = like(formula = form.lc.lin, family = "poisson", data = stomach.cancer.until2007, E = stomach.cancer.until2007$total.t )

form.lc.basic = total ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon
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
                               control.predictor = list(link = 1),
                               bru_max_iter = 1
                )) 

# apc model, rw1, stomach cancer
res.apc.1.s = bru(components = comp.apc.1,
                  likelihood.apc.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1),
                                 bru_max_iter = 1
                  )) 

# apc model, rw2, lung cancer
res.apc.2.l = bru(components = comp.apc.2,
                  likelihood.apc.l, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1),
                                 bru_max_iter = 1
                  )) 

# apc model, rw2, stomach cancer
res.apc.2.s = bru(components = comp.apc.2,
                  likelihood.apc.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1),
                                 bru_max_iter = 1
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

# color palette:
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1',
                   '#5d8060', '#D7B36A', '#826133', '#A85150')


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

# combine the different predictions:

data.pred.l <- rbind(data.apc.1.l, data.apc.2.l, data.lc.l, data.lc.lin.l, data.lc.basic.l) %>%
  mutate("method" = rep(c("APC rw1", "APC rw2", "LCC", "LCC linear", "LC"), each = 324)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.s <- rbind(data.apc.1.s, data.apc.2.s, data.lc.s, data.lc.lin.s, data.lc.basic.s) %>%
  mutate("method" = rep(c("APC rw1", "APC rw2", "LCC", "LCC linear", "LC"), each = 324)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

#   ----   Objective measures of prediction performance   ----

# display four significant digits 
options(pillar.sigfig = 4)

# compute statistics for all age groups
pred.statistics.include.l <- data.pred.l %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

pred.statistics.include.s <- data.pred.s %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 included");cat("\n Lung cancer data: ");pred.statistics.include.l;cat("\n Stomach cancer data: ");pred.statistics.include.s

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

# plot LCC-results for stomach cancer

# by-age:
gg.LCC.a.s <- ggplot(data.pred.s %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                       filter(method %in% c("LC", "LCC", "LCC linear")),
                         aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "LCC models - stomach cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.LCC.a.s

ggsave('univariate-LCC-by-age-stomach.png',
       plot = gg.LCC.a.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# by-period
gg.LCC.p.s <- ggplot(data.pred.s %>%
                       filter(x > 5) %>%
                       filter(method %in% c("LC", "LCC", "LCC linear")),
                     aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  geom_vline(aes(xintercept = "2007"), color = palette.basis[length(palette.basis)]) +
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "LCC models - stomach cancer", x = "Age groups", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)
gg.LCC.p.s

ggsave('univariate-LCC-by-period-stomach.png',
       plot = gg.LCC.p.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot results for APC models for stomach cancer:
# by-age:
gg.APC.a.s <- ggplot(data.pred.s %>%
                       filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                       filter(method %in% c("APC rw1", "APC rw2")),
                     aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "APC models - stomach cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.APC.a.s

ggsave('univariate-APC-by-age-stomach.png',
       plot = gg.APC.a.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# by-period
gg.APC.p.s <- ggplot(data.pred.s %>%
                       filter(x > 5) %>%
                       filter(method %in% c("APC rw1", "APC rw2")),
                     aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  geom_vline(aes(xintercept = "2007"), color = palette.basis[length(palette.basis)]) +
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "APC models - stomach cancer", x = "Age groups", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)
gg.APC.p.s

ggsave('univariate-APC-by-period-stomach.png',
       plot = gg.APC.p.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot results for LCC-models for lung cancer

# by-age:
gg.LCC.a.l <- ggplot(data.pred.l %>%
                       filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                       filter(method %in% c("LC", "LCC", "LCC linear")),
                     aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "LCC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.LCC.a.l

ggsave('univariate-LCC-by-age-lung.png',
       plot = gg.LCC.a.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# by-period
gg.LCC.p.l <- ggplot(data.pred.l %>%
                       filter(x > 5) %>%
                       filter(method %in% c("LC", "LCC", "LCC linear")),
                     aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  geom_vline(aes(xintercept = "2007"), color = palette.basis[length(palette.basis)]) +
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "LCC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)
gg.LCC.p.l

ggsave('univariate-LCC-by-period-lung.png',
       plot = gg.LCC.p.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot results for APC-models for lung cancer

# by-age:
gg.APC.a.l <- ggplot(data.pred.l %>%
                       filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                       filter(method %in% c("APC rw1", "APC rw2")),
                     aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "APC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.APC.a.l

ggsave('univariate-APC-by-age-lung.png',
       plot = gg.APC.a.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# by-period
gg.APC.p.l <- ggplot(data.pred.l %>%
                       filter(x > 5) %>%
                       filter(method %in% c("APC rw1", "APC rw2")),
                     aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  geom_vline(aes(xintercept = "2007"), color = palette.basis[length(palette.basis)]) +
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "APC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)
gg.APC.p.l

ggsave('univariate-APC-by-period-lung.png',
       plot = gg.APC.p.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)


# plot results for best APC model and best LCC for stomach cancer

# by-age:
gg.compare.a.s <- ggplot(data.pred.s %>%
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
  labs(title = "APC rw2 and LCC models - stomach cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.compare.a.s

ggsave('univariate-comparison-by-age-stomach.png',
       plot = gg.compare.a.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# by-period
gg.compare.p.s <- ggplot(data.pred.s %>%
                       filter(x > 5) %>%
                       filter(method %in% c("APC rw2", "LCC")),
                     aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  geom_vline(aes(xintercept = "2007"), color = palette.basis[length(palette.basis)]) +
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "APC rw2 and LCC models - stomach cancer", x = "Age groups", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)
gg.compare.p.s

ggsave('univariate-comparison-by-period-stomach.png',
       plot = gg.compare.p.s,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot results for best APC model and best LCC for lung cancer

# by-age:
gg.compare.a.l <- ggplot(data.pred.l %>%
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
  labs(title = "APC rw2 and LCC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.compare.a.l

ggsave('univariate-comparison-by-age-lung.png',
       plot = gg.compare.a.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# by-period
gg.compare.p.l <- ggplot(data.pred.l %>%
                       filter(x > 5) %>%
                       filter(method %in% c("LCC", "APC rw2")),
                     aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`, group = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  geom_vline(aes(xintercept = "2007"), color = palette.basis[length(palette.basis)]) +
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "APC rw2 and LCC models - lung cancer", x = "Age groups", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)
gg.compare.p.l

ggsave('univariate-comparison-by-period-lung.png',
       plot = gg.compare.p.l,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# save work space image:
save.image("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Workspaces/ws_predict-long-term.RData")

