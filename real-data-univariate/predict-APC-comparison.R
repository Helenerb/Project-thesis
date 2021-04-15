# Compare predictions for stomach and lung cancer data in year 2015 and 1016
# using different versions of the APC

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
  mutate(age.int = replace(age.int, age.int == "85 - 89", "85 +")) %>%
  group_by(age.int, year) %>%
  summarize(total = sum(total), male = sum(male), female = sum(female)) %>%
  mutate(year = format(as.POSIXct(year, format="%d.%m.%Y"), format="%Y")) %>%
  filter(year < 2017)

# read and format stomach cancer data
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
  mutate(cohort = 5 * (max(x) - x) + t) %>%
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
  mutate(cohort = 5 * (max(x) - x) + t) %>%
  mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
  mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
  left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t")) %>%
  mutate("mortality rate" = total/total.t)

# extract 2015 and 2016-observations:
lung.20162015 <- lung.cancer %>% filter(year == "2016" | year == "2015")
stomach.20162015 <- stomach.cancer %>% filter(year == "2016" | year == "2015")

# set 2015- and 2016 data as missing:
# these will be used as the observation data when running inlabru. 
lung.cancer.wo1516 <- lung.cancer %>% 
  mutate(total = replace(total, year == "2016" | year == "2015", NA)) %>%
  mutate(male = replace(male, year == "2016" | year == "2015", NA)) %>%
  mutate(female = replace(female, year == "2016" | year == "2015", NA))

stomach.cancer.wo1516 <- stomach.cancer %>% 
  mutate(total = replace(total, year == "2016" | year == "2015", NA)) %>%
  mutate(male = replace(male, year == "2016" | year == "2015", NA)) %>%
  mutate(female = replace(female, year == "2016" | year == "2015", NA))

#   ----   define inlabru model components:

pc.prior.rho <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.4)))
pc.prior.phi <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.6)))
pc.prior.epsilon <- list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))
pc.prior.psi <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.5)))

# note: this way of defining the x, t and k - values are only possible knowing that 
# we do not have any missing data points in stomach.cancer (we can, however, have zero-data.)
comp.apc.1 = ~ -1 + 
  Int(1) + 
  rho(x, model = "rw2", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi(t, model = "rw2", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi(cohort, model = "rw2", values = unique(stomach.cancer$cohort), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

comp.apc.2 = ~ -1 + 
  Int(1) + 
  rho(x, model = "rw1", values = unique(stomach.cancer$x), constr = TRUE, hyper = pc.prior.rho) + 
  phi(t, model = "rw1", values = unique(stomach.cancer$t), constr = TRUE, hyper = pc.prior.phi) + 
  psi(cohort, model = "rw1", values = unique(stomach.cancer$cohort), constr = TRUE, hyper = pc.prior.psi) + 
  epsilon(xt, model = "iid", hyper = pc.prior.epsilon)

# define the formula for the predictor for the APC model:
form.apc = total ~ -1 + Int + rho + phi + psi + epsilon

# add data as lung.cancer and offset as population$total, which is the number of people at-risk
likelihood.apc.l = like(formula = form.apc, family = "poisson", data = lung.cancer.wo1516, E = lung.cancer.wo1516$total.t)
likelihood.apc.s = like(formula = form.apc, family = "poisson", data = stomach.cancer.wo1516, E = stomach.cancer.wo1516$total.t)

# the same control compute as in Sara's first example 
c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

#initial.values = list(alpha.c = alpha, beta.c = beta, kappa.c = kappa, phi.t = phi*(1:nt))

res.apc.1.l = bru(components = comp.apc.1,
              likelihood.apc.l, 
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1)
              )) 
res.apc.1.l = bru_rerun(res.apc.1.l)

res.apc.1.s = bru(components = comp.apc.1,
                  likelihood.apc.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
res.apc.1.s = bru_rerun(res.apc.1.s)

res.apc.2.l = bru(components = comp.apc.2,
                  likelihood.apc.l, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 control.predictor = list(link = 1)
                  )) 
res.apc.2.l = bru_rerun(res.apc.2.l)

res.apc.2.s = bru(components = comp.apc.2,
                  likelihood.apc.s, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c, 
                                 control.predictor = list(link = 1)
                  )) 
res.apc.2.s = bru_rerun(res.apc.2.s)

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#607A4D', '#026AA1', '#A85150')
palette.light <- c('#ABC9E6', '#F3DC90', '#86A46F', '#40BBFD', '#C38281')


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

data.pred.apc.l <- rbind(data.apc.1.l, data.apc.2.l) %>%
  mutate("method" = rep(c("rw2", "rw1"), each = 324))

data.pred.apc.s <- rbind(data.apc.1.s, data.apc.2.s) %>%
  mutate("method" = rep(c("rw2", "rw1"), each = 324))

gg.pred.apc.age.facet.l <- ggplot(data = data.pred.apc.l, aes(x = year)) + 
  geom_vline(xintercept = "2014", color = palette.light[5]) + 
  geom_errorbar(aes(min = `0.025quant`, ymax = `0.975quant`, color = `method`), position=position_dodge(width=0.5)) +
  #geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`)) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) +
  geom_point(aes(y = mean, color = `method`)) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  #theme(axis.text.x = element_text(angle = -45, vjust = -1, hjust=1)) + 
  theme(panel.spacing.x = unit(1, "lines")) + 
  facet_wrap(~age) + 
  labs(title = "Mortality rate for lung cancer by calendar year", y = "Mortality rate", x = "Year")

gg.pred.apc.age.facet.l

gg.pred.apc.age.facet.s <- ggplot(data = data.pred.apc.s, aes(x = year)) + 
  geom_vline(xintercept = "2014", color = palette.light[5]) + 
  geom_errorbar(aes(min = `0.025quant`, ymax = `0.975quant`, color = `method`), position=position_dodge(width=0.5)) +
  #geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`)) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) +
  geom_point(aes(y = mean, color = `method`)) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  #theme(axis.text.x = element_text(angle = -45, vjust = -1, hjust=1)) + 
  theme(panel.spacing.x = unit(1, "lines")) + 
  facet_wrap(~age) + 
  labs(title = "Mortality rate for stomach cancer by calendar year", y = "Mortality rate", x = "Year")

gg.pred.apc.age.facet.s

gg.pred.2016.apc.l <- ggplot(data.pred.apc.l %>% filter(year == "2016"), aes(x = fct_inorder(age))) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_point(aes(y = mean, color = `method`, group = 1)) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2016", x = "Age groups", y = "Mortality rate")

gg.pred.2015.apc.l <- ggplot(data.pred.apc.l %>% filter(year == "2015"), aes(x = fct_inorder(age))) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_point(aes(y = mean, color = `method`, group = 1)) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2015", x = "Age groups", y = "Mortality rate")

(gg.pred.2015.apc.l |gg.pred.2016.apc.l)  + 
  plot_annotation(title = "Predicted values for lung cancer by age groups") + 
  plot_layout(guides = "collect")

# stomach cancer plots:
gg.pred.2016.apc.s <- ggplot(data.pred.apc.s %>% filter(year == "2016"), aes(x = fct_inorder(age))) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_point(aes(y = mean, color = `method`, group = 1)) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2016", x = "Age groups", y = "Mortality rate")

gg.pred.2015.apc.s <- ggplot(data.pred.apc.s %>% filter(year == "2015"), aes(x = fct_inorder(age))) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_point(aes(y = mean, color = `method`, group = 1)) + 
  geom_point(aes(y = `mortality rate`, color = "Observed")) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2015", x = "Age groups", y = "Mortality rate")

(gg.pred.2015.apc.s |gg.pred.2016.apc.s)  + 
  plot_annotation(title = "Predicted values for stomach cancer by age groups") + 
  plot_layout(guides = "collect")

# plots using just x on the x-axis 
gg.pred.2016.apc.s <- ggplot(data.pred.apc.s %>% filter(year == "2016"), aes(x = x)) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  #geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2016", x = "Age groups", y = "Mortality rate")

gg.pred.2015.apc.s <- ggplot(data.pred.apc.s %>% filter(year == "2015"), aes(x = x)) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  #geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2016", x = "Age groups", y = "Mortality rate")

(gg.pred.2015.apc.s |gg.pred.2016.apc.s)  + 
  plot_annotation(title = "Predicted values for stomach cancer by age groups") + 
  plot_layout(guides = "collect")

gg.pred.2016.apc.l <- ggplot(data.pred.apc.l %>% filter(year == "2016"), aes(x = x)) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  #geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2016", x = "Age groups", y = "Mortality rate")

gg.pred.2015.apc.l <- ggplot(data.pred.apc.l %>% filter(year == "2015"), aes(x = x)) + 
  #geom_errorbar(aes(x = x, min = `0.025quant`, ymax = `0.975quant`, color = `method`)) +
  #geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = `method`), alpha = 0.5) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `method`), alpha = 0.5) +
  geom_point(aes(y = mean, color = `method`, group = 1), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed"), shape = 4, size = 3) + 
  scale_color_manual(name = "Prediction method",
                     values = palette.basis) +
  scale_fill_manual(name = "Prediction method",
                    values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "2016", x = "Age groups", y = "Mortality rate")

(gg.pred.2015.apc.l |gg.pred.2016.apc.l)  + 
  plot_annotation(title = "Predicted values for lung cancer by age groups") + 
  plot_layout(guides = "collect")

