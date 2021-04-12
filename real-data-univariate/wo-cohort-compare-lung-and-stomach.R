# This script will combine the univariate no-cohort LC model for the stomach cancer data
# and the lung cancer data. 
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

# read lung cancer data
lung.cancer  <- read_excel("lungCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  mutate(t = as.integer(year)-1999) %>% mutate(t.1 = t) %>%
  mutate(x = parse_number(age)) %>% mutate(x.1 = x) %>%
  mutate(xt = ((x%/%5)*(2016-1998) +t))

# read stomach cancer data
stomach.cancer  <- read_excel("stomachCancer-germany.xls") %>%
  rename(sex = "...1") %>% rename(age = "...2") %>%
  pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
  mutate(sex = replace(sex, sex == "männlich", "male")) %>%
  mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(total = male + female) %>%
  mutate(t = as.integer(year)-1999) %>% mutate(t.1 = t) %>%
  mutate(x = parse_number(age)) %>% mutate(x.1 = x) %>%
  mutate(xt = ((x%/%5)*(2016-1998) +t))

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


#   ----  Start defining the inlabru model components  ----   
#  this first attempt is based on the L-C-cohort-v2

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = length(unique(stomach.cancer$age)))  #  not sure if you did this correctly
e.vec = 1

# keeping the "less informative priors" that seemed to work well
pc.prior.alpha <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.4)))
pc.prior.kappa <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))
pc.prior.gamma <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.5)))

# this is just how we define our model
comp.sto = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior.alpha) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t.1, model = "rw1", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

comp.lun = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.alpha) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.1, model = "iid", extraconstr = list(A = A.mat, e = e.vec)) + 
  kappa(t.1, model = "rw1", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.kappa) +
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

# here total will refer to the total deaths for each age-year
#form.1 = total ~ -1 + Int + alpha + beta*phi + beta*kappa + gamma + epsilon

# first: without cohort
form = total ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon 

# add data as stomach.cancer and offset as population$total, which is the number of people at-risk
likelihood.sto = like(formula = form, family = "poisson", data = stomach.cancer, E = population$total)
likelihood.lun = like(formula = form, family = "poisson", data = lung.cancer, E = population$total)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha.c = alpha, beta.c = beta, kappa.c = kappa, phi.t = phi*(1:nt))

res.sto = bru(components = comp.sto,
          likelihood.sto, 
          options = list(verbose = F,
                         bru_verbose = 1, 
                         num.threads = "1:1",
                         control.compute = c.c
          )) 

res.sto = bru_rerun(res.sto)

res.lun = bru(components = comp.lun,
              likelihood.lun, 
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c
              )) 

res.lun = bru_rerun(res.lun)

cat("Fixed effects and hyper parameters for stomach cancer data:")
res.sto$summary.fixed
res.sto$summary.hyperpar

rbind(res.sto$summary.random$alpha, res.lun$summary.random$alpha) %>%
  mutate(type = rep(c("Stomach cancer", "Lung cancer"), each=length(res.lun$summary.random$alpha$ID))) %>%
  ggplot(aes(x = ID)) + 
  geom_errorbar(aes(ID, min = `0.025quant`, ymax =`0.975quant`), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean, color = type)) + 
  ggtitle("Alpha")

rbind(res.sto$summary.random$beta, res.lun$summary.random$beta) %>%
  mutate(type = rep(c("Stomach cancer", "Lung cancer"), each=length(res.lun$summary.random$beta$ID))) %>%
  ggplot(aes(x = ID)) + 
  geom_errorbar(aes(ID, min = `0.025quant`, ymax =`0.975quant`), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean, color = type)) + 
  ggtitle("Beta")

rbind(res.sto$summary.random$kappa, res.lun$summary.random$kappa) %>%
  mutate(type = rep(c("Stomach cancer", "Lung cancer"), each=length(res.lun$summary.random$kappa$ID))) %>%
  ggplot(aes(x = ID)) + 
  geom_errorbar(aes(ID, min = `0.025quant`, ymax =`0.975quant`), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean, color = type)) + 
  ggtitle("Kappa")

data.phi.sto = data.frame(cbind(ID = unique(stomach.cancer$t), 
                            mean = res.sto$summary.fixed$mean[2]*unique(stomach.cancer$t),
                            X0.025quant = res.sto$summary.fixed$`0.025quant`[2]*unique(stomach.cancer$t),
                            X0.975quant = res.sto$summary.fixed$`0.975quant`[2]*unique(stomach.cancer$t)))

data.phi.lun = data.frame(cbind(ID = unique(lung.cancer$t), 
                                mean = res.lun$summary.fixed$mean[2]*unique(lung.cancer$t),
                                X0.025quant = res.lun$summary.fixed$`0.025quant`[2]*unique(lung.cancer$t),
                                X0.975quant = res.lun$summary.fixed$`0.975quant`[2]*unique(lung.cancer$t)))

rbind(data.phi.sto, data.phi.lun) %>%
  mutate(type = rep(c("Stomach cancer", "Lung cancer"), each=length(data.phi.lun$ID))) %>%
  ggplot(aes(x = ID)) + 
  geom_errorbar(aes(ID, min = X0.025quant, ymax =X0.975quant), position=position_dodge(width=0.5)) +
  geom_point(aes(y = mean, color = type)) + 
  ggtitle("Phi")

data.eta.sto <- data.frame(eta = res.sto$summary.linear.predictor$mean[1:324]) %>%
  mutate(t = stomach.cancer$t, x = stomach.cancer$x)
data.eta.lun <- data.frame(eta = res.lun$summary.linear.predictor$mean[1:324]) %>%
  mutate(t = lung.cancer$t, x = lung.cancer$x)

ggplot(data = data.eta.sto, aes(x=t, y=x, fill = eta)) + geom_tile() + 
  xlab("Time: 1999 - 2016") + 
  ylab("Age: 0 - 85+") + 
  ggtitle("Stomach cancer")

ggplot(data = data.eta.lun, aes(x=t, y=x, fill = eta)) + geom_tile() + 
  xlab("Time: 1999 - 2016") + 
  ylab("Age: 0 - 85+") + 
  ggtitle("Lung cancer")

data.eta.t <- rbind(data.eta.sto, data.eta.lun) %>%
  mutate(type = rep(c("Stomach cancer", "Lung cancer"), each=length(data.eta.lun$eta))) %>%
  group_by(t, type) %>% summarize(eta.t = mean(eta))
ggplot(data = data.eta.t, aes(x = t)) + 
  geom_point(aes(y = eta.t, color = type)) + 
  ggtitle("Eta for t - real data") + xlab("t") + ylab("Predictor")

data.eta.x <- rbind(data.eta.sto, data.eta.lun) %>%
  mutate(type = rep(c("Stomach cancer", "Lung cancer"), each=length(data.eta.lun$eta))) %>%
  group_by(x, type) %>% summarize(eta.x = mean(eta))
ggplot(data = data.eta.x, aes(x = x)) + 
  geom_point(aes(y = eta.x, color = type)) + 
  ggtitle("Eta for x - real data") + xlab("x") + ylab("Predictor")

