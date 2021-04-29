# preliminary sensitivity analysis when updating priors for the LC model to some that 
# are less informative. 

# use the multivariate LC with shared age effect to get an impression of the sensitivity to priors

# Simple implementation of multivariate APC model with shared age effect:
library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(lubridate)
library(readxl)

#   ----   Syntax should be common for inla and inlabru   ----

# read and format population data
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

# read and format lung cancer data
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
# note: not actually needed for inlabru
lung.cancer <- lung.cancer %>%
  mutate(x.c = x) %>% mutate(x.c2 = x) %>% mutate(x.c3 = x) %>%
  mutate(t.c = t) %>% mutate(t.c2 = t) %>% mutate(t.c3 = t) %>%
  mutate(k.c = k) %>% mutate(k.c2 = k)

# set data for years 2008-2016 ro NA
# these will be used as the observation data when running inlabru. 
lung.cancer.until2007 <- lung.cancer %>% 
  mutate(cases = replace(cases, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) 

# split lung cancer into two data sets - male (0) and female (1)
lung.cancer.until2007.0 <- lung.cancer.until2007 %>% filter(s == 0)
lung.cancer.until2007.1 <- lung.cancer.until2007 %>% filter(s == 1)

pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.658, 0.5)))
pc.prior.beta <- list(prec = list(prior = "pc.prec", param = c(0.0837, 0.5)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))  # perhaps update this after run
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.00813, 0.5)))
pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.0316, 0.5)))

pc.prior.alpha.ui <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.beta.ui <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.kappa.ui <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))  # perhaps update this after run
pc.prior.epsilon.ui <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.gamma.ui <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))

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

comp.lc.ui = ~ -1 +
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha.ui, scale.model = TRUE) +
  beta(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta.ui) +
  phi0(t, model = "linear", prec.linear = 1) +
  phi1(t.c1, model = "linear", prec.linear = 1) +
  kappa0(t.c2, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa.ui, scale.model = TRUE) +
  kappa1(t.c3, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa.ui, scale.model = TRUE) +
  gamma0(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma.ui, scale.model = TRUE) +
  gamma1(k.c, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma.ui, scale.model = TRUE) +
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon.ui) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon.ui)

# define two different likelihoods and formulas, one for male and one for female:
form.lc.0 = cases ~ -1 + mu + alpha + beta*phi0 + beta*kappa0 + gamma0 + epsilon0
likelihood.lc.0 = like(formula = form.lc.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.lc.1 = cases ~ -1 + mu + alpha + beta*phi1 + beta*kappa1 + gamma1 + epsilon1
likelihood.lc.1 = like(formula = form.lc.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res.lc = bru(components = comp.lc,
               likelihood.lc.0,
               likelihood.lc.1,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1)
               )) 

res.lc = bru_rerun(res.lc)

res.lc.ui = bru(components = comp.lc.ui,
             likelihood.lc.0,
             likelihood.lc.1,
             options = list(verbose = F,
                            bru_verbose = 1, 
                            num.threads = "1:1",
                            control.compute = c.c,
                            control.predictor = list(link = 1)
             )) 

res.lc.ui = bru_rerun(res.lc.ui)

# constructed to ensure that observation data is in the same order as prediction data
lung.cancer.0 <- lung.cancer %>% filter(s == 0)
lung.cancer.1 <- lung.cancer %>% filter(s == 1)

# print CPO of non-predicted values:
cat("CPO, APC - common age: \n"); -1*mean(log(res.apc.a$cpo$cpo[!is.na(res.apc.a$cpo$cpo)]))
cat("\n CPO, APC - common age - default priors: \n"); -1*mean(log(res.apc.a.def$cpo$cpo[!is.na(res.apc.a.def$cpo$cpo)]))

# calculate score statistics
data.pred.lc <- res.lc$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.lc.ui <- res.lc.ui$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred = rbind(data.pred.lc, data.pred.lc.ui) %>%
  mutate("method" = rep(c("informative", "less informative"), each = 648))

# extract predicted years, and ages above x = 5 (mortality rate very close to zero)
pred.statistics.cutoff <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age: x <= 5 omitted");cat("\n Lung cancer data: \n");pred.statistics.cutoff

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

