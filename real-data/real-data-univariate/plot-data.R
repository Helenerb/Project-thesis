# script used to display the data sets that we will use, to get an overview of 
# the tendencies we might find. 

library(INLA)
library(inlabru)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(lubridate)
library(readxl)

# color palette:
palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

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

pop.by.year <- population %>%
  group_by(year) %>%
  summarize(total = sum(total), male = sum(male), female = sum(female))

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
  mutate("mortality rate" = total/total.t)%>%
  mutate(male.p = male/male.t) %>%
  mutate(female.p = female/female.t)

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
  mutate("mortality rate" = total/total.t) %>%
  mutate(male.p = male/male.t) %>%
  mutate(female.p = female/female.t)

# most recent plotting scheme of mortality rate:

p.lung.age.rate = ggplot(data = lung.cancer %>%
                           filter(year %in% c("1999", "2005", "2011" ,"2016"))) +
  geom_line(aes(x = age.int, y = male.p, color = year)) + 
  geom_line(aes(x = age.int, y = female.p, color = year)) + 
  geom_point(aes(x = age.int, y = male.p, shape = "Male", color = year)) + 
  geom_point(aes(x = age.int, y = female.p, shape = "Female", color = year)) + 
  scale_color_manual(name = "Year", values = palette.basis) +
  scale_shape_manual(name = "Sex", values  = c(3,2)) + 
  labs(x = "Age group", y = "Mortality rate", title = "Lung cancer")

p.stomach.age.rate = ggplot(data = stomach.cancer %>%
                           filter(year %in% c("1999", "2005", "2011" ,"2016"))) +
  geom_line(aes(x = age.int, y = male.p, color = year)) + 
  geom_line(aes(x = age.int, y = female.p, color = year)) + 
  geom_point(aes(x = age.int, y = male.p, shape = "Male", color = year)) + 
  geom_point(aes(x = age.int, y = female.p, shape = "Female", color = year)) + 
  scale_color_manual(name = "Year", values = palette.basis) +
  scale_shape_manual(name = "Sex", values  = c(3,2)) + 
  labs(x = "Age group", y = "Mortality rate", title = "Stomach cancer")

p.data.age.rate <- (p.lung.age.rate | p.stomach.age.rate) +
  plot_annotation(title = "Observed mortality rates") +
  plot_layout(guides = "collect") & theme(legend.position = 'right')
p.data.age.rate
  

ggsave('data-age-rate',
       plot = p.data.age.rate,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

p.lung.cohort.rate = ggplot(data = lung.cancer %>%
                           filter(year %in% c("1999", "2005", "2011" ,"2016"))) +
  geom_line(aes(x = cohort, y = male.p, color = year)) + 
  geom_line(aes(x = cohort, y = female.p, color = year)) + 
  geom_point(aes(x = cohort, y = male.p, shape = "Male", color = year)) + 
  geom_point(aes(x = cohort, y = female.p, shape = "Female", color = year)) + 
  scale_color_manual(name = "Year", values = palette.basis) +
  scale_shape_manual(name = "Sex", values  = c(3,2)) + 
  labs(x = "Cohort", y = "Mortality rate", title = "Lung cancer")
p.lung.cohort.rate

p.stomach.cohort.rate = ggplot(data = stomach.cancer %>%
                              filter(year %in% c("1999", "2005", "2011" ,"2016"))) +
  geom_line(aes(x = cohort, y = male.p, color = year)) + 
  geom_line(aes(x = cohort, y = female.p, color = year)) + 
  geom_point(aes(x = cohort, y = male.p, shape = "Male", color = year)) + 
  geom_point(aes(x = cohort, y = female.p, shape = "Female", color = year)) + 
  scale_color_manual(name = "Year", values = palette.basis) +
  scale_shape_manual(name = "Sex", values  = c(3,2)) + 
  labs(x = "Cohort", y = "Mortality rate", title = "Stomach cancer")

p.data.cohort.rate <- (p.lung.cohort.rate | p.stomach.cohort.rate) +
  plot_annotation(title = "Observed mortality rates") +
  plot_layout(guides = "collect") & theme(legend.position = 'right')
p.data.cohort.rate


ggsave('data-cohort-rate',
       plot = p.data.cohort.rate,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)


# old plotting scheme: 

#  aggregated by age:
stomach.cancer.age = stomach.cancer %>%
  group_by(x) %>% 
  summarize(male = sum(male), female = sum(female),
            male.t = sum(male.t), female.t = sum(female.t),
            male.p = mean(male.p), female.p = mean(female.p)) %>%
  pivot_longer(!x, names_to = "sex", values_to = "occurrances")

lung.cancer.age = lung.cancer %>%
  group_by(x) %>% 
  summarize(male = sum(male), female = sum(female),
            male.t = sum(male.t), female.t = sum(female.t),
            male.p = mean(male.p), female.p = mean(female.p)) %>%
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

# use the values of the total population from the lung cancer data set, could have been from the stomach cancer data set as well
p.total.age <- ggplot(lung.cancer.age %>% filter(sex == "male.t" | sex == "female.t") %>%
                        mutate(sex = replace(sex, sex == "male.t", "male")) %>%
                        mutate(sex = replace(sex, sex == "female.t", "female"))) + 
  geom_point(aes(x = x, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="German population", x = "Age group", y = "Total occurrances")

p.age.total <- (p.total.age | p.stomach.age | p.lung.age) + 
  plot_annotation(title = "Total occurrances of cases by age group") +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
p.age.total

ggsave('data-age-total',
       plot = p.age.total,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

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
  summarize(male = sum(male), female = sum(female),
            male.t = sum(male.t), female.t = sum(female.t),
            male.p = mean(male.p), female.p = mean(female.p)) %>%
  mutate(female.p = female/female.t, male.p = male/male.t) %>%
  pivot_longer(!year, names_to = "sex", values_to = "occurrances")

lung.cancer.year = lung.cancer %>%
  group_by(year) %>% 
  summarize(male = sum(male), female = sum(female),
            male.t = sum(male.t), female.t = sum(female.t),
            male.p = mean(male.p), female.p = mean(female.p)) %>%
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

p.total.year <- ggplot(stomach.cancer.year %>% filter(sex == "male.t" | sex == "female.t") %>%
                        mutate(sex = replace(sex, sex == "male.t", "male")) %>%
                        mutate(sex = replace(sex, sex == "female.t", "female"))) + 
  geom_point(aes(x = year, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  labs(title="German population", x = "Calendar year", y = "Total occurrances")

p.year.total <- (p.total.year | p.stomach.year | p.lung.year) + 
  plot_annotation(title = "Total occurrances of cases by calendar year") + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
p.year.total

ggsave('data-year-total',
       plot = p.year.total,
       device = "png",
       path = '/Users/helen/OneDrive - NTNU/Vår 2021/Project-thesis/real-data/real-data-univariate/Figures',
       height = 5, width = 8, 
       dpi = "retina"
)

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
  group_by(cohort) %>% 
  summarize(male = sum(male), female = sum(female), male.t = sum(male.t),
            female.t = sum(female.t)) %>%
  mutate(male.p = male/male.t, female.p = female/female.t) %>%
  pivot_longer(!cohort, names_to = "sex", values_to = "occurrances")
  
lung.cancer.cohort = lung.cancer %>%
  group_by(cohort) %>% 
  summarize(male = sum(male), female = sum(female), male.t = sum(male.t),
            female.t = sum(female.t)) %>%
  mutate(male.p = male/male.t, female.p = female/female.t) %>%
  pivot_longer(!cohort, names_to = "sex", values_to = "occurrances")

# plot total cases:

p.stomach.cohort <- ggplot(stomach.cancer.cohort %>% filter(sex == "male" | sex == "female")) + 
  geom_col(aes(x = cohort, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  labs(title="Stomach cancer", x = "Cohort", y = "Total occurrances")

p.lung.cohort <- ggplot(lung.cancer.cohort %>% filter(sex == "male" | sex == "female")) + 
  geom_col(aes(x = cohort, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  labs(title="Lung cancer", x = "Cohort", y = "Total occurrances")

p.total.cohort <- ggplot(stomach.cancer.cohort %>% filter(sex == "male.t" | sex == "female.t") %>%
                         mutate(sex = replace(sex, sex == "male.t", "male")) %>%
                         mutate(sex = replace(sex, sex == "female.t", "female"))) + 
  geom_col(aes(x = cohort, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="German population", x = "Cohort", y = "Total occurrances")

(p.total.cohort | p.stomach.cohort | p.lung.cohort) + 
  plot_annotation(title = "Total occurrances of cases by birth year") +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# plot percent-wise cases:

p.stomach.cohort.p <- ggplot(stomach.cancer.cohort %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = cohort, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Stomach cancer", x = "Cohort", y = " ")

p.lung.cohort.p <- ggplot(lung.cancer.cohort %>% filter(sex == "male.p" | sex == "female.p")) + 
  geom_point(aes(x = cohort, y = occurrances, color = sex)) + 
  scale_color_manual(name = "Sex", values = palette.basis) +
  labs(title="Lung cancer", x = "Cohort", y = " ")

(p.stomach.cohort.p | p.lung.cohort.p) + 
  plot_annotation(title = "Percent-wise occurrances of cases by birth year") + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')


