# First attempt at prediction with inlabru. 

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(lubridate)
library(readxl)

#   ----   read data from excel files   ----

# read population data
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
  mutate(age.int = replace(age.int, age.int == "85 - 89", "85")) %>%
  group_by(age.int, year) %>%
  summarize(total = sum(total), male = sum(male), female = sum(female)) %>%
  mutate(year = format(as.POSIXct(year, format="%d.%m.%Y"), format="%Y")) %>%
  filter(year < 2017)

lung.cancer  <- read_excel("lungCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  mutate(t = as.integer(year)-1999) %>% mutate(t.1 = t) %>%
  mutate(x = parse_number(age)) %>% mutate(x.1 = x) %>%
  mutate(xt = ((x%/%5)*(2016-1998) +t)) %>%
  mutate(cohort = t - x) %>%
  mutate(birth.year = 1999 + cohort) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t"))

# extract 2016-observations:
obs.2016 <- lung.cancer %>% filter(year == "2016")

# set 2016-observations to NA in the lung.cancer data
# use this when running inlabru - should produce predictive results for the 2016-observations
lung.cancer.wo2016 <- lung.cancer %>% 
  mutate(total = replace(total, year == "2016", NA)) %>%
  mutate(male = replace(male, year == "2016", NA)) %>%
  mutate(female = replace(female, year == "2016", NA))

# attempt first with basic LC-model:
#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$x)))  #  not sure if you did this correctly
e.vec = 1

# The following priors are not adjusted to the data:
pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.4)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))
#pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.5)))

# this is just how we define our model
comp = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t.1, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  #gamma(cohort, model = "rw1", values =  unique(lung.cancer$cohort), constr = TRUE, hyper = pc.prior.gamma) +
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)


# Simplest LC model:
form.1 = total ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon

# add data as lung.cancer and offset as population$total, which is the number of people at-risk
likelihood.1 = like(formula = form.1, family = "poisson", data = lung.cancer.wo2016, E = lung.cancer.wo2016$total.t)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha.c = alpha, beta.c = beta, kappa.c = kappa, phi.t = phi*(1:nt))

res = bru(components = comp,
          likelihood.1, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c,
                         control.predictor = list(link = 1)
          )) 
res = bru_rerun(res)

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#A85150', '#607A4D', '#026AA1')
palette.light <- c('#ABC9E6', '#F3DC90', '#C38281', '#86A46F', '#40BBFD')

summary(res)
plot(res$summary.fitted.values$mean[307:324], obs.2016$total/obs.2016$total.t)

data.pred <- res$summary.fitted.values %>%
  slice(1:324) %>%
  bind_cols(lung.cancer %>% select(- c("x.1", "t.1", "xt"))) %>%
  mutate("mortality rate" = total/total.t) %>%
  mutate(age = replace(age, age == "85", "85 + "))

gg.pred.age.facet <- ggplot(data = data.pred, aes(x = year)) + 
  geom_vline(xintercept = "2015", color = palette.light[3]) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) +
  geom_errorbar(aes(min = `0.025quant`, ymax = `0.975quant`), color = palette.light[1], position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean, color = "LC")) + 
  scale_color_manual(name = "Prediction method",
                     breaks = c("Observed", "LC"),
                     values = c("LC" = palette.basis[1], "Observed" = palette.basis[2]) ) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  #theme(axis.text.x = element_text(angle = -45, vjust = -1, hjust=1)) + 
  theme(panel.spacing.x = unit(1, "lines")) + 
  facet_wrap(~age)

gg.pred.age.facet

gg.pred.2016 <- ggplot(data.pred %>% filter(year == "2016") %>%
                         mutate(id.order = factor(age, levels=age))) + 
  geom_point(aes(x = id.order, y = `mortality rate`, color = "Observed")) + 
  geom_errorbar(aes(x = id.order, min = `0.025quant`, ymax = `0.975quant`),
                color = palette.light[1], position=position_dodge(width=0.5)) +
  geom_point(aes(x = id.order, y = mean, color = "LC")) + 
  scale_color_manual(name = "Prediction method",
                     breaks = c("Observed", "LC"),
                     values = c("LC" = palette.basis[1], "Observed" = palette.basis[2])) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Prediction for 2016", x = "Age groups", y = "Mortality rate")
  
gg.pred.2016

