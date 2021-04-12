# script used to display the data sets that we will use, to get an overview of 
# the tendencies we might find. 

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(lubridate)
library(readxl)

# color palette.
palette.basis <- c('#70A4D4', '#ECC64B', '#A85150', '#607A4D', '#026AA1')
palette.light <- c('#ABC9E6', '#F3DC90', '#C38281', '#86A46F', '#40BBFD')

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

# read and format stomach cancer data
stomach.cancer  <- read_excel("stomachCancer-germany.xls") %>%
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

# read and format lung cancer data:
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


#  plot percentwise occurrance:



#  aggregated by age:
stomach.cancer.age = stomach.cancer %>%
  group_by(x) %>% 
  summarize(male = sum(male), female = sum(female), male.t = sum(male.t), female.t = sum(female.t)) %>%
  mutate(female.p = female/female.t, male.p = male/male.t) %>%
  pivot_longer(!x, names_to = "sex", values_to = "occurrances")

lung.cancer.age = lung.cancer %>%
  group_by(x) %>% 
  summarize(male = sum(male), female = sum(female), male.t = sum(male.t), female.t = sum(female.t)) %>%
  mutate(female.p = female/female.t, male.p = male/male.t) %>%
  pivot_longer(!x, names_to = "sex", values_to = "occurrances")

# plot total occurrances:

p.stomach.age <- ggplot(stomach.cancer.age %>% filter(sex == "male" | sex == "female")) + 
  geom_point(aes(x = x, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Stomach cancer", x = "Age group", y = "Total occurrances")

p.lung.age <- ggplot(lung.cancer.age %>% filter(sex == "male" | sex == "female")) + 
  geom_point(aes(x = x, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Lung cancer", x = "Age group", y = "Total occurrances")

(p.stomach.age | p.lung.age) + 
  plot_annotation(title = "Total occurrances of cases by age group") +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# plot percentwise occurances:
p.stomach.age.p <- ggplot(stomach.cancer.age %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = x, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Stomach cancer", x = "Age group", y = "Percent-wise occurrances")

p.lung.age.p <- ggplot(lung.cancer.age %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = x, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Lung cancer", x = "Age group", y = "Percent-wise occurrances")

(p.stomach.age.p | p.lung.age.p) + 
  plot_annotation(title = "Percent-wise occurrances of cases by age group") +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# aggregated by year:
stomach.cancer.year = stomach.cancer %>%
  group_by(year) %>% 
  summarize(total = sum(total), male = sum(male), female = sum(female), male.t = sum(male.t), female.t = sum(female.t)) %>%
  mutate(female.p = female/female.t, male.p = male/male.t) %>%
  pivot_longer(!year, names_to = "sex", values_to = "occurrances")

lung.cancer.year = lung.cancer %>%
  group_by(year) %>% summarize(total = sum(total), male = sum(male), female = sum(female), male.t = sum(male.t), female.t = sum(female.t)) %>%
  mutate(female.p = female/female.t, male.p = male/male.t) %>%
  pivot_longer(!year, names_to = "sex", values_to = "occurrances")

# plot total cases: 
p.stomach.year <- ggplot(stomach.cancer.year %>% filter(sex == "male" | sex == "female")) + 
  geom_point(aes(x = year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  labs(title="Stomach cancer", x = "Calendar year", y = "Total occurrances")

p.lung.year <- ggplot(lung.cancer.year %>% filter(sex == "male" | sex == "female")) + 
  geom_point(aes(x = year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  labs(title="Lung cancer", x = "Calendar year", y = "Total occurrances")

(p.stomach.year | p.lung.year) + 
  plot_annotation(title = "Total occurrances of cases by calendar year") + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# plot percent-wise cases:
p.stomach.year.p <- ggplot(stomach.cancer.year %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  labs(title="Stomach cancer", x = "Calendar year", y = "Total occurrances")

p.lung.year.p <- ggplot(lung.cancer.year %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  labs(title="Lung cancer", x = "Calendar year", y = "Total occurrances")

(p.stomach.year.p | p.lung.year.p) + 
  plot_annotation(title = "Percent-wise occurrances of cases by calendar year") + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')


# aggregated by cohort:
stomach.cancer.cohort = stomach.cancer %>%
  group_by(birth.year) %>% 
  summarize(male = sum(male), female = sum(female), male.t = sum(male.t), female.t = sum(female.t)) %>%
  mutate(male.p = male/male.t, female.p = female/female.t) %>%
  pivot_longer(!birth.year, names_to = "sex", values_to = "occurrances")
  
lung.cancer.cohort = lung.cancer %>%
  group_by(birth.year) %>% 
  summarize(male = sum(male), female = sum(female), male.t = sum(male.t), female.t = sum(female.t)) %>%
  mutate(male.p = male/male.t, female.p = female/female.t) %>%
  pivot_longer(!birth.year, names_to = "sex", values_to = "occurrances")

# plot total cases:

p.stomach.cohort <- ggplot(stomach.cancer.cohort %>% filter(sex == "male" | sex == "female")) + 
  geom_point(aes(x = birth.year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Stomach cancer", x = "Birth year", y = "Total occurrances")

p.lung.cohort <- ggplot(lung.cancer.cohort %>% filter(sex == "male" | sex == "female")) + 
  geom_point(aes(x = birth.year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Lung cancer", x = "Birth year", y = "Total occurrances")

(p.stomach.cohort | p.lung.cohort) + 
  plot_annotation(title = "Total occurrances of cases by birth year") +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# plot percent-wise cases:

p.stomach.cohort.p <- ggplot(stomach.cancer.cohort %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = birth.year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Stomach cancer", x = "Birth year", y = "Total occurrances")

p.lung.cohort.p <- ggplot(lung.cancer.cohort %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = birth.year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Lung cancer", x = "Birth year", y = "Total occurrances")

(p.stomach.cohort.p | p.lung.cohort.p) + 
  plot_annotation(title = "Percent-wise occurrances of cases by birth year") + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')


