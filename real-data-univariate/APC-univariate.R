# first run with APC model, check that it works. Using real lung and stomach 
# cancer data. 

# Try to introduce cohort effect and see how that works:) 
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


# read stomach cancer data
stomach.cancer  <- read_excel("stomachCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  mutate(t = as.integer(year)-1999) %>% 
  mutate(age.int = parse_number(age)) %>%
  mutate(x = parse_number(age)%/%5) %>% 
  mutate(xt = (x*(2016-1998) +t)) %>%
  mutate(cohort = 5 * (max(stomach.cancer$x) - x) + t) %>%
  mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
  mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t"))


#   ----  Start defining the inlabru model components  ----   
# implement basic APC model

# use same priors as for LC in first step
pc.prior.rho <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.4)))
pc.prior.phi <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))
pc.prior.psi <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.5)))

# note: this way of defining the x, t and k - values are only possible knowing that 
# we do not have any missing data points in stomach.cancer (we can, however, have zero-data.)
comp.apc = ~ -1 + 
  Int(1) + 
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi(cohort, model = "rw2", values = unique(stomach.cancer$cohort), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

# define the formula for the predictor for the APC model:
form.apc = total ~ -1 + Int + rho + phi + psi + epsilon

# add data as lung.cancer and offset as population$total, which is the number of people at-risk
likelihood.apc = like(formula = form.apc, family = "poisson", data = stomach.cancer, E = stomach.cancer$total.t)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha.c = alpha, beta.c = beta, kappa.c = kappa, phi.t = phi*(1:nt))

res.apc = bru(components = comp.apc,
          likelihood.apc, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c
          )) 


res.apc$summary.fixed
res.apc$summary.hyperpar

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#A85150', '#607A4D', '#026AA1')
palette.light <- c('#ABC9E6', '#F3DC90', '#C38281', '#86A46F', '#40BBFD')

data.rho = res.apc$summary.random$rho
gg.rho <- ggplot(data.rho, aes(x = ID)) + 
  geom_errorbar(aes(ID, min = `0.025quant`, ymax =`0.975quant`), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean), fill = palette.basis[1]) + 
  labs(title = "Age effect", x = "x", y = "rho(x)")

data.phi = res.apc$summary.random$phi
gg.phi <- ggplot(data.phi, aes(x = ID)) + 
  geom_errorbar(aes(ID, min = `0.025quant`, ymax =`0.975quant`), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean), fill = palette.basis[1]) + 
  labs(title = "Period effect", x = "t", y = "phi(t)")

data.psi = res.apc$summary.random$psi
gg.psi <- ggplot(data.psi, aes(x = ID)) + 
  geom_errorbar(aes(ID, min = `0.025quant`, ymax =`0.975quant`), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean), fill = palette.basis[1]) + 
  labs(title = "Cohort effect", x = "k", y = "psi(k)")

data.eta <- data.frame(eta = res.apc$summary.linear.predictor$mean[1:324]) %>%
  mutate(t = stomach.cancer$t, x = stomach.cancer$x, cohort = stomach.cancer$cohort)
gg.eta <- ggplot(data = data.eta, aes(x=t, y=x, fill = eta)) + geom_tile() + 
  scale_fill_gradient(high = palette.basis[1], low = palette.basis[2]) + 
  labs(x = "Time: 1999 - 2016", y = "Age: 0 - 85+")
gg.eta

data.eta.t <- data.eta %>%
  group_by(t) %>% summarize(eta.t = mean(eta))
gg.eta.t <- ggplot(data = data.eta.t, aes(x = t)) + 
  geom_point(aes(y = eta.t), fill = palette.basis[1]) + 
  labs(title = "Predictor eta, grouped by period", x = "t")

data.eta.x <- data.eta %>% 
  group_by(x) %>% summarize(eta.x = mean(eta))
gg.eta.x <- ggplot(data = data.eta.x, aes(x = x)) + 
  geom_point(aes(y = eta.x), fill = palette.basis[1]) + 
  labs(title = "Predictor eta, grouped by age group", x = "x")

data.eta.k <- data.eta %>% 
  group_by(cohort) %>% summarize(eta.k = mean(eta))
gg.eta.k <- ggplot(data = data.eta.k, aes(x = cohort)) + 
  geom_point(aes(y = eta.k), fill = palette.basis[1]) + 
  labs(title = "Predictor eta, grouped by cohort", x = "k")

(gg.rho | gg.phi / gg.psi)/(gg.eta.x | gg.eta.t | gg.eta.k)/(gg.eta) + 
  plot_annotation(title = "Results from runninge inlabru on APC model")

