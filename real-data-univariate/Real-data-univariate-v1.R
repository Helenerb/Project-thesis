# This is the first attempt at running inlabru inference on real data. 
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
  filter(year != age) %>%
  filter(age != "Insgesamt") %>%
  mutate(age.number =  parse_number(age)) %>%
  mutate(age.number = replace(age.number, age == "unter 1 Jahr", 0)) %>%
  mutate(age.int = 5*(age.number%/%5)) %>%
  mutate(age.int = paste(age.int, "-", age.int + 4)) %>%
  mutate(age.int = replace(age.int, age.int == "85 - 89", "85")) %>%
  group_by(age.int, year) %>%
  summarize(total = sum(total), male = sum(male), female = sum(female)) %>%
  mutate(year = format(as.POSIXct(year, format="%d.%m.%Y"), format="%Y"))

# read lung cancer data
lung.cancer  <- read_excel("lungCancer-germany.xls") %>%
  rename(sex = "...1") %>%
  rename(age = "...2") %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  #mutate(age.1 = age) %>%
  #mutate(year.1 = year) %>%
  #mutate(xt = paste(age,year))
  mutate(t = as.integer(year)-1998) %>% mutate(t.1 = t) %>%
  mutate(x = parse_number(age)) %>% mutate(x.1 = x) %>%
  mutate(xt = )


#   ----  Start defining the inlabru model components  ----   
#  this first attempt is based on the L-C-cohort-v2

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = length(unique(lung.cancer$age)))  #  not sure if you did this correctly
e.vec = 1

# keeping the "less informative priors" that seemed to work well
pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.4)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.5)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))
pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.5)))

# this is just how we define our model
comp = ~ -1 + 
  Int(1) + 
  alpha(age, model = "rw1", values = unique(lung.cancer$age), constr = TRUE, hyper = pc.prior.alpha) + 
  phi(year, model = "linear", mean.linear = -0.5, prec.linear = 0.25) +
  beta(age.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(year.1, model = "rw1", values = unique(lung.cancer$year), constr = TRUE, hyper = pc.prior.kappa) +
  #gamma(cohort, model = "rw1", values = cohort.min:cohort.max, constr = TRUE, hyper = pc.prior.gamma) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

# here total will refer to the total deaths for each age-year
#form.1 = total ~ -1 + Int + alpha + beta*phi + beta*kappa + gamma + epsilon

# first: without cohort
form.1 = total ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon 

# add data as lung.cancer and offset as population$total, which is the number of people at-risk
likelihood.1 = like(formula = form.1, family = "poisson", data = lung.cancer, E = population$total)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha.c = alpha, beta.c = beta, kappa.c = kappa, phi.t = phi*(1:nt))

res = bru(components = comp,
          likelihood.1, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c
          )) 
