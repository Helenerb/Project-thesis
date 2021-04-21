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
  #mutate(x0 = replace(x, s == 1, NA)) %>% mutate(x0.c = x0) %>%
  #mutate(x1 = replace(x, s == 0, NA)) %>% mutate(x1.c = x1) %>%
  mutate(t.c = t) %>% mutate(t.c2 = t) %>% mutate(t.c3 = t) %>%
  #mutate(t0 = replace(t, s == 1, NA)) %>% mutate(t0.c = t0) %>%
  #mutate(t1 = replace(t, s == 0, NA)) %>% mutate(t1.c = t1) %>%
  mutate(k.c = k) %>% mutate(k.c2 = k) #%>%
  #mutate(k0 = replace(k, s == 1, NA)) %>%
  #mutate(k1 = replace(k, s == 0, NA))

# lung.cancer <- lung.cancer %>%
#   mutate(x.c = x) %>% mutate(x.c2 = x) %>% mutate(x.c3 = x) %>%
#   mutate(x0 = replace(x, s == 1, -1)) %>% mutate(x0.c = x0) %>%
#   mutate(x1 = replace(x, s == 0, -1)) %>% mutate(x1.c = x1) %>%
#   mutate(t.c = t) %>% mutate(t.c2 = t) %>% mutate(t.c3 = t) %>%
#   mutate(t0 = replace(t, s == 1, -1)) %>% mutate(t0.c = t0) %>%
#   mutate(t1 = replace(t, s == 0, -1)) %>% mutate(t1.c = t1) %>%
#   mutate(k.c = k) %>% mutate(k.c2 = k) %>%
#   mutate(k0 = replace(k, s == 1, -1)) %>%
#   mutate(k1 = replace(k, s == 0, -1))

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

# comp.fix.alpha = ~ -1 +
#   mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
#   alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
#   beta0(x0, model = "rw1", values = unique(lung.cancer$x0), extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
#   beta1(x1, model = "rw1", values = unique(lung.cancer$x1), extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
#   phi0(t0, model = "linear", prec.linear = 1) +
#   phi1(t1, model = "linear", prec.linear = 1) +
#   kappa0(t0.c, model = "rw1", values = unique(lung.cancer$t0), constr = TRUE, hyper = pc.prior.kappa) +
#   kappa1(t1.c, model = "rw1", values = unique(lung.cancer$t1), constr = TRUE, hyper = pc.prior.kappa) +
#   gamma0(k0, model = "rw1", values = unique(lung.cancer$k0), constr = TRUE, hyper = pc.prior.gamma) +
#   gamma1(k1, model = "rw1", values = unique(lung.cancer$k1), constr = TRUE, hyper = pc.prior.gamma) +
#   epsilon(xts, model = "iid", hyper = pc.prior.epsilon)

comp.fix.alpha = ~ -1 +
  Int(1) + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  alpha0(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  alpha1(x.c, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) +
  beta(x.c2, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior.beta) +
  #beta1(x.c2, model = "rw1", values = unique(lung.cancer$x), extraconstr = list(A = A.mat.1, e = e.vec), hyper = pc.prior.beta) +
  phi(t, model = "linear", prec.linear = 1) +
  #phi1(t.c1, model = "linear", prec.linear = 1) +
  kappa(t.c2, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  #kappa1(t.c3, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  gamma(k, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  #gamma1(k.c, model = "rw1", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.gamma) +
  epsilon(xts, model = "iid", hyper = pc.prior.epsilon)

# define two different likelihoods and formulas, one for male and one for female:
form.0 = cases ~ -1 + mu + alpha0 + beta*phi + beta*kappa + gamma + epsilon
likelihood.0 = like(formula = form.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.1 = cases ~ -1 + mu + alpha1 + beta*phi + beta*kappa + gamma + epsilon
likelihood.1 = like(formula = form.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

# test with "normal" version - single likelihood, both male and female observations:
#form.common <- cases ~ -1 + Int + alpha + beta*phi + beta*kappa + gamma + epsilon
#likelihood.common = like(formula = form.common, family = "poisson", data = lung.cancer.until2007, E = lung.cancer.until2007$population)


c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res.fix.alpha = bru(components = comp.fix.alpha,
                  likelihood.0,
                  likelihood.1,
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
res.fix.alpha = bru_rerun(res.fix.alpha)

summary(res.fix.alpha)

# constructed to ensure that observation data is in the same order as prediction data
lung.cancer.0 <- lung.cancer %>% filter(s == 0)
lung.cancer.1 <- lung.cancer %>% filter(s == 1)

# combine predictions and observations into one dataframe. 
data.pred <- res.fix.alpha$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

pred.statistics.cutoff <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(sex) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 omitted");cat("\n Lung cancer data: ");pred.statistics.cutoff

# plot:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')

gg.pred <- ggplot(data.pred %>%
                      filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                    aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `sex`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `sex`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = `sex`, fill = `sex`), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Sex-specific alpha", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred

