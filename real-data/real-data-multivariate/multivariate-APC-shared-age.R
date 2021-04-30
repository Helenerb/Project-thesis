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


#   ----   Inlabru syntax:   ----



# uninformative priors
pc.prior.rho <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.phi <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.epsilon.apc <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.psi <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))


# common age effect:
comp.apc = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho, scale.model = TRUE) + 
  phi0(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi, scale.model = TRUE) + 
  phi1(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi, scale.model = TRUE) + 
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon.apc) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon.apc)

# define two different likelihoods and formulas, one for male and one for female:
form.apc.0 = cases ~ -1 + mu + rho + phi0 + psi0 + epsilon0
likelihood.apc.0 = like(formula = form.apc.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.apc.1 = cases ~ -1 + mu + rho + phi1 + psi1 + epsilon1
likelihood.apc.1 = like(formula = form.apc.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

res.apc.a = bru(components = comp.apc,
                likelihood.apc.0,
                likelihood.apc.1,
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1),
                               bru_max_iter=1
                )) 

#default priors:
# comparable process with APC model - shared age effect:
pc.prior.rho.def <- list(prec = list(prior = "pc.prec"))
pc.prior.phi.def <- list(prec = list(prior = "pc.prec"))
pc.prior.epsilon.apc.def <- list(prec = list(prior = "pc.prec"))
pc.prior.psi.def <- list(prec = list(prior = "pc.prec"))

# common age effect:
comp.apc.def = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.rho.def, scale.model = TRUE) + 
  phi0(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi.def, scale.model = TRUE) + 
  phi1(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.phi.def, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi.def, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.psi.def, scale.model = TRUE) + 
  epsilon0(xts, model = "iid", hyper = pc.prior.epsilon.apc.def) + 
  epsilon1(xts, model = "iid", hyper = pc.prior.epsilon.apc.def)

res.apc.a.def = bru(components = comp.apc.def,
                likelihood.apc.0,
                likelihood.apc.1,
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1),
                               bru_max_iter=1
                )) 


#   ----   Display results   ----

# constructed to ensure that observation data is in the same order as prediction data
lung.cancer.0 <- lung.cancer %>% filter(s == 0)
lung.cancer.1 <- lung.cancer %>% filter(s == 1)

# print CPO of non-predicted values:
cat("CPO, APC - common age: \n"); -1*mean(log(res.apc.a$cpo$cpo[!is.na(res.apc.a$cpo$cpo)]))
cat("\n CPO, APC - common age - default priors: \n"); -1*mean(log(res.apc.a.def$cpo$cpo[!is.na(res.apc.a.def$cpo$cpo)]))

# calculate score statistics
data.pred.apc <- res.apc.a$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred.apc.def <- res.apc.a.def$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred = rbind(data.pred.apc, data.pred.apc.def) %>%
  mutate("method" = rep(c("uninformative priors", "default priors"), each = 648))

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
