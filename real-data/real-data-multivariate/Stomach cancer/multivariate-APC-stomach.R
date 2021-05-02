# compoarison of predictions for multivariate APC models, for stomach cancer
# compare different versions of multivariate APC-models.

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

stomach.cancer  <- read_excel("stomachCancer-germany.xls") %>%
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

# set data for years 2008-2016 ro NA
# these will be used as the observation data when running inlabru. 
stomach.cancer.until2007 <- stomach.cancer %>% 
  mutate(cases = replace(cases, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) 

stomach.cancer.until2007.0 <- stomach.cancer.until2007 %>% filter(s == 0)
stomach.cancer.until2007.1 <- stomach.cancer.until2007 %>% filter(s == 1)

pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))

# common age effect - Apc:
comp.Apc = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi0(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi1(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior) 

# define two different likelihoods and formulas, one for male and one for female:
form.Apc.0 = cases ~ -1 + mu + rho + phi0 + psi0 + epsilon
likelihood.Apc.0 = like(formula = form.Apc.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.Apc.1 = cases ~ -1 + mu + rho + phi1 + psi1 + epsilon
likelihood.Apc.1 = like(formula = form.Apc.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res.Apc = bru(components = comp.Apc,
              likelihood.Apc.0,
              likelihood.Apc.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

# constructed to ensure that observation data is in the same order as prediction data
stomach.cancer.0 <- stomach.cancer %>% filter(s == 0)
stomach.cancer.1 <- stomach.cancer %>% filter(s == 1)

data.pred.Apc <- res.Apc$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

#   ----   common period effect: aPc   ----:
comp.aPc = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.aPc.0 = cases ~ -1 + mu + rho0 + phi + psi0 + epsilon
likelihood.aPc.0 = like(formula = form.aPc.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.aPc.1 = cases ~ -1 + mu + rho1 + phi + psi1 + epsilon
likelihood.aPc.1 = like(formula = form.aPc.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.aPc = bru(components = comp.aPc,
              likelihood.aPc.0,
              likelihood.aPc.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.aPc <- res.aPc$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

#   ----   common cohort effect : apC   ----
comp.apC = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi0(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi1(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior) 

# define two different likelihoods and formulas, one for male and one for female:
form.apC.0 = cases ~ -1 + mu + rho0 + phi0 + psi + epsilon
likelihood.apC.0 = like(formula = form.apC.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.apC.1 = cases ~ -1 + mu + rho1 + phi1 + psi + epsilon
likelihood.apC.1 = like(formula = form.apC.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.apC = bru(components = comp.apC,
              likelihood.apC.0,
              likelihood.apC.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.apC <- res.apC$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

#   ----   no shared effects : apc   ----
comp.apc = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi0(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi1(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.apc.0 = cases ~ -1 + mu + rho0 + phi0 + psi0 + epsilon
likelihood.apc.0 = like(formula = form.apc.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.apc.1 = cases ~ -1 + mu + rho1 + phi1 + psi1 + epsilon
likelihood.apc.1 = like(formula = form.apc.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.apc = bru(components = comp.apc,
              likelihood.apc.0,
              likelihood.apc.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.apc <- res.apc$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   All effects shared : APC   ----

comp.APC = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.APC.0 = cases ~ -1 + mu + rho + phi + psi + epsilon
likelihood.APC.0 = like(formula = form.APC.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.APC.1 = cases ~ -1 + mu + rho + phi + psi + epsilon
likelihood.APC.1 = like(formula = form.APC.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.APC = bru(components = comp.APC,
              likelihood.APC.0,
              likelihood.APC.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.APC <- res.APC$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Shared age and period effects : APc   ----

comp.APc = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.APc.0 = cases ~ -1 + mu + rho + phi + psi0 + epsilon
likelihood.APc.0 = like(formula = form.APc.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.APc.1 = cases ~ -1 + mu + rho + phi + psi1 + epsilon
likelihood.APc.1 = like(formula = form.APc.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.APc = bru(components = comp.APc,
              likelihood.APc.0,
              likelihood.APc.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.APc <- res.APc$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Shared age and cohort effects : ApC   ----

comp.ApC = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi0(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi1(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.ApC.0 = cases ~ -1 + mu + rho + phi0 + psi + epsilon
likelihood.ApC.0 = like(formula = form.ApC.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.ApC.1 = cases ~ -1 + mu + rho + phi1 + psi + epsilon
likelihood.ApC.1 = like(formula = form.ApC.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.ApC = bru(components = comp.ApC,
              likelihood.ApC.0,
              likelihood.ApC.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.ApC <- res.ApC$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   shared period and cohort effects : aPC   ----

comp.aPC = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  psi(k, model = "rw2", values = unique(stomach.cancer$k), constr = TRUE, hyper = pc.prior, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior)

# define two different likelihoods and formulas, one for male and one for female:
form.aPC.0 = cases ~ -1 + mu + rho0 + phi + psi + epsilon
likelihood.aPC.0 = like(formula = form.aPC.0, family = "poisson", data = stomach.cancer.until2007.0, E = stomach.cancer.until2007.0$population)

form.aPC.1 = cases ~ -1 + mu + rho1 + phi + psi + epsilon
likelihood.aPC.1 = like(formula = form.aPC.1, family = "poisson", data = stomach.cancer.until2007.1, E = stomach.cancer.until2007.1$population)

res.aPC = bru(components = comp.aPC,
              likelihood.aPC.0,
              likelihood.aPC.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.aPC <- res.aPC$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(stomach.cancer.0, stomach.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))


#   ----   Plot and compare results   ----

data.pred <- rbind(data.pred.APC, data.pred.APc, data.pred.ApC, data.pred.Apc, data.pred.aPC, data.pred.aPc, data.pred.apC, data.pred.apc) %>%
  mutate("method" = rep(c("APC", "APc", "ApC", "Apc", "aPC", "aPc", "apC", "apc"), each = 648))

pred.statistics.cutoff <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(method) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 omitted");cat("\n Stomach cancer data: ");pred.statistics.cutoff

# plot:

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#1C84BB', '#A85150', '#DA871F',
                   '#4C7246', '#D7B36A', '#FB5E4E', '#696B8D', '#76A7A6', '#826133')

gg.pred <- ggplot(data.pred %>%
                    filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"))%>%
                    filter(method %in% c("apc", "aPc")),
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
  labs(title = "Multivariate APC models - Stomach cancer data", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred

ggsave('multivariate-APC-by-age-stomach.png',
       plot = gg.pred,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures'
)

# plot cohortwise

ggplot(data.pred %>%
         filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
         filter(method %in% c("apc", "aPc", "apC", "Apc")),
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
  labs(title = "Multivariate APC models", x = "Cohort", y = "Mortality rate") + 
  facet_wrap(~year)







