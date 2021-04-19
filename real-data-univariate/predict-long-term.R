# script for comparison of long-term predictions using different APC and LC models
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
  filter(year < 2017)

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

# set 2015- and 2016 data as missing:
# these will be used as the observation data when running inlabru. 
lung.cancer.until2007 <- lung.cancer %>% 
  mutate(total = replace(total, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(male = replace(male, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA))

stomach.cancer.until2007 <- stomach.cancer %>% 
  mutate(total = replace(total, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(male = replace(male, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA))


