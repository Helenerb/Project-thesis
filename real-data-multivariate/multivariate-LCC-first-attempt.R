# first attempt at brute-force model for prediction using a multivariate LCC model

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(lubridate)
library(readxl)

# load and format data

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

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

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

# start fitting multivariate model:

lung.cancer.mv <- lung.cancer.until2007 %>%
  left_join(res.lc.l$summary.random$alpha %>% select(c("ID", "mean")), by = c("x" = "ID"), suffix = c("", ".alpha")) %>%
  mutate("male offset" = male.t * exp(mean)) %>%
  mutate("female offset" = female.t * exp(mean)) %>%
  rename(mean.alpha = mean)

# formula for fitting male effects:
form.lc.m = male ~ -1 + Int + beta*phi + beta*kappa + gamma + epsilon
likelihood.lc.m = like(formula = form.lc.m, family = "poisson", data = lung.cancer.mv, E = lung.cancer.mv$`male offset`)

res.lc.m = bru(components = comp.lc,
               likelihood.lc.m, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.lc.m = bru_rerun(res.lc.m)

# formula for fitting female effects:
form.lc.f = female ~ -1 + Int + beta*phi + beta*kappa + gamma + epsilon
likelihood.lc.f = like(formula = form.lc.f, family = "poisson", data = lung.cancer.mv, E = lung.cancer.mv$`female offset`)

res.lc.f = bru(components = comp.lc,
               likelihood.lc.f, 
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 
res.lc.f = bru_rerun(res.lc.f)

# present results:
data.lc.m <- res.lc.m$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.lc.f <- res.lc.f$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt")))

data.pred <- rbind(data.lc.m, data.lc.f) %>%
  mutate("sex" = rep(c("male", "female"), each = 324)) %>%
  mutate("female mortality rate" = female/female.t) %>%
  mutate("male mortality rate" = male/male.t) %>%
  mutate(SE.m = (mean - `male mortality rate`)^2) %>%
  mutate(SE.f = (mean - `female mortality rate`)^2) %>%
  mutate(DSS.m = ((`male mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(DSS.f = ((`female mortality rate` - mean)/sd)^2 + 2*log(sd)) 

data.pred.f <- data.pred %>% filter(sex == "female") %>%
  select(-c("SE.m", "male", "male.t", "male mortality rate"))

data.pred.m <- data.pred %>% filter(sex == "male") %>%
  select(-c("SE.f", "female", "female.t", "female mortality rate"))

#   ----   Objective measures of prediction performance   ----
pred.statistics.f <- data.pred.f %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  summarise(MSE = mean(SE.f), MDSS = mean(DSS.f))

pred.statistics.m <- data.pred.m %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  summarise(MSE = mean(SE.m), MDSS = mean(DSS.m))

cat("\n Female: \n");pred.statistics.f;cat("\n Male: \n");pred.statistics.f

# compute statistics with cutoff at age group 5
pred.statistics.cutoff.f <- data.pred.f %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  summarise(MSE = mean(SE.f), MDSS = mean(DSS.f))

pred.statistics.cutoff.m <- data.pred.m %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  summarise(MSE = mean(SE.m), MDSS = mean(DSS.m))

cat("\n Age <= 5 omitted");cat("\n Female: \n");pred.statistics.cutoff.m;cat("\n Male: \n");pred.statistics.cutoff.m

# plot the two sexes:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')
palette.light <- c('#ABC9E6', '#F3DC90', '#86A46F', '#40BBFD', '#C38281')


gg.pred.mv <- ggplot(data.pred %>%
                       filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                     aes(x = x)) +
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `sex`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `sex`, group = 1), shape = 19) + 
  geom_point(aes(y = `female mortality rate`, color = "Observed female", fill = "Observed female"), shape = 4, size = 3) + 
  geom_point(aes(y = `male mortality rate`, color = "Observed female", fill = "Observed female"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  labs(title = "Multivariate prediction, lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)

gg.pred.mv





