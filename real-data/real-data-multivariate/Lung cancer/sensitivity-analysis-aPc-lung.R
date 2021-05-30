# sensitivity analysis of APC model with shared period effect:

load("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Lung cancer/Workspaces/ws_sensitivity-analysis-aPc-lung.RData")


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
  mutate(s = recode(sex, male = 0, female = 1))  %>%
  mutate(xts = (x*(2016-1998) + t + s*(max(x)*(2016-1998) + max(t)))) %>% 
  mutate("mortality rate" = cases/population)

# set data for years 2008-2016 ro NA
# these will be used as the observation data when running inlabru. 
lung.cancer.until2007 <- lung.cancer %>% 
  mutate(cases = replace(cases, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) 

lung.cancer.until2007.0 <- lung.cancer.until2007 %>% filter(s == 0)
lung.cancer.until2007.1 <- lung.cancer.until2007 %>% filter(s == 1)

# constructed to ensure that observation data is in the same order as prediction data
lung.cancer.0 <- lung.cancer %>% filter(s == 0)
lung.cancer.1 <- lung.cancer %>% filter(s == 1)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

pc.prior.1 <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc.prior.2 <- list(prec = list(prior = "pc.prec", param = c(3, 0.05)))
pc.prior.3 <- list(prec = list(prior = "pc.prec", param = c(1, 0.95)))

# using prior 1
comp.aPc.1 = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.1, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.1, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.1, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.1, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.1, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior.1)

# using prior 2
comp.aPc.2 = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.2, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.2, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.2, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.2, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.2, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior.2)

# using prior 3
comp.aPc.3 = ~ -1 + 
  mu(s, model = "iid", hyper = list(prec = list(fixed = TRUE, initial = log(0.001)))) +
  rho0(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.3, scale.model = TRUE) + 
  rho1(x, model = "rw2", values = unique(lung.cancer$x), constr = TRUE, hyper = pc.prior.3, scale.model = TRUE) + 
  phi(t, model = "rw2", values = unique(lung.cancer$t), constr = TRUE, hyper = pc.prior.3, scale.model = TRUE) + 
  psi0(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.3, scale.model = TRUE) + 
  psi1(k, model = "rw2", values = unique(lung.cancer$k), constr = TRUE, hyper = pc.prior.3, scale.model = TRUE) + 
  epsilon(xts, model = "iid", hyper = pc.prior.3)

# define two different likelihoods and formulas, one for male and one for female:
form.aPc.0 = cases ~ -1 + mu + rho0 + phi + psi0 + epsilon
likelihood.aPc.0 = like(formula = form.aPc.0, family = "poisson", data = lung.cancer.until2007.0, E = lung.cancer.until2007.0$population)

form.aPc.1 = cases ~ -1 + mu + rho1 + phi + psi1 + epsilon
likelihood.aPc.1 = like(formula = form.aPc.1, family = "poisson", data = lung.cancer.until2007.1, E = lung.cancer.until2007.1$population)

# results with prior 1
res.aPc.1 = bru(components = comp.aPc.1,
              likelihood.aPc.0,
              likelihood.aPc.1,
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             control.predictor = list(link = 1),
                             bru_max_iter = 1
              )) 

data.pred.aPc.1 <- res.aPc.1$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

# results with prior 2
res.aPc.2 = bru(components = comp.aPc.2,
                likelihood.aPc.0,
                likelihood.aPc.1,
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1),
                               bru_max_iter = 1
                )) 

data.pred.aPc.2 <- res.aPc.2$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

# results with prior 3
res.aPc.3 = bru(components = comp.aPc.3,
                likelihood.aPc.0,
                likelihood.aPc.1,
                options = list(verbose = F,
                               bru_verbose = 1, 
                               num.threads = "1:1",
                               control.compute = c.c,
                               control.predictor = list(link = 1),
                               bru_max_iter = 1
                )) 

data.pred.aPc.3 <- res.aPc.3$summary.fitted.values %>%
  slice(1:648) %>%
  bind_cols(rbind(lung.cancer.0, lung.cancer.1)) %>%
  mutate(SE = (mean - `mortality rate`)^2) %>%
  mutate(DSS = ((`mortality rate` - mean)/sd)^2 + 2*log(sd)) %>%
  mutate(contained = as.integer((`mortality rate` >= `0.025quant` & `mortality rate` <= `0.975quant`)))

data.pred <- rbind(data.pred.aPc.1, data.pred.aPc.2, data.pred.aPc.3) %>%
  mutate("prior" = rep(c("P(sd > 1) = 0.05", "P(sd > 3) = 0.05", "P(sd > 1) = 0.95"), each = 648))

# display four significant digits 
options(pillar.sigfig = 4)

pred.statistics.cutoff <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  filter(x > 5) %>%
  group_by(prior) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 omitted");cat("\n Lung cancer data - senstivity analysis: ");pred.statistics.cutoff

pred.statistics.included <- data.pred %>% 
  filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>% 
  group_by(prior) %>%
  summarise(MSE = mean(SE), MDSS = mean(DSS), contained = mean(contained))

cat("\n Age <= 5 included");cat("\n Lung cancer data - sensitivity analysis: ");pred.statistics.included

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#1C84BB', '#A85150', '#DA871F',
                   '#4C7246', '#D7B36A', '#FB5E4E', '#696B8D', '#76A7A6', '#826133')

gg.pred <- ggplot(data.pred %>%
                    filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                  aes(x = x)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `prior`, group = interaction(prior, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `prior`, group = 1, group = interaction(prior, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prior",
                     values = palette.basis) +
  scale_fill_manual(name = "Prior",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  labs(title = "Sensitivity analysis on precision priors - aPc model, lung cancer", x = "Age groups", y = "Mortality rate") + 
  facet_wrap(~year)
gg.pred

ggsave('sensitivity-analysis-aPcby-age-lung.png',
       plot = gg.pred,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot cohortwise

gg.pred.cohort <- ggplot(data.pred %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")),
                         aes(x = k)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `prior`, group = interaction(prior, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `prior`, group = 1, group = interaction(prior, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prior",
                     values = palette.basis) +
  scale_fill_manual(name = "Prior",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  labs(title = "Sensitivity analysis on precision priors - aPc model, lung cancer", x = "Cohort", y = "Mortality rate") + 
  facet_wrap(~year)

gg.pred.cohort

ggsave('sensitivity-analysis-aPc-by-cohort-lung.png',
       plot = gg.pred.cohort,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

# plot along years - for different age groups:
gg.pred.period <- ggplot(data.pred %>%
                           filter(year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) %>%
                           filter(x > 5),
                         aes(x = year)) + 
  geom_ribbon(aes(min = `0.025quant`, ymax = `0.975quant`, fill = `prior`, group = interaction(prior, sex)), alpha = 0.5) +
  geom_point(aes(y = mean, color = `prior`, group = 1, group = interaction(prior, sex)), shape = 19) + 
  geom_point(aes(y = `mortality rate`, color = "Observed", fill = "Observed", shape = `sex`), size = 2) + 
  scale_color_manual(name = "Prior",
                     values = palette.basis) +
  scale_fill_manual(name = "Prior",
                    values = palette.basis) +
  scale_shape_manual(values = c(3,2)) + 
  labs(title = "Sensitivity analysis on precision priors - aPc model, lung cancer", x = "Year", y = "Mortality rate") + 
  theme(axis.text.x = element_text(angle = -30, hjust=0)) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) + 
  facet_wrap(~age)

gg.pred.period

ggsave('sensitivity-analysis-aPc-by-period-lung.png',
       plot = gg.pred.period,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

save.image("/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-multivariate/Lung cancer/Workspaces/ws_sensitivity-analysis-aPc-lung.RData")


