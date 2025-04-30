#################################################
## CAUSAL INFERENCE AND POLICY EVALUATION      ##
## -- Assignment 3                             ##
##################################################

#############
## Setup   ##
#############

# Empty working space
rm(list=ls())

#Load Packages
# Define packages that you need,
packages_vector <- c( "haven", "dplyr",  "sandwich",  "jtools", "fBasics",
                      "xtable",  "stargazer", "data.table", "tidyverse", "ggplot2",
                      "AER", # AER package for ivreg command
                      "causalweight") # for semiparametric LATE
lapply(packages_vector, require, character.only = TRUE)

##################
## Question 1   ##
##################

### (a)
# Only include rows which describe a married couple
data_m <- data %>% filter(msample == 1)

print("Average Age and Income of Mothers and Fathers When First Child was Born")
age_income_df <- data.frame(Parent = c("Mother", "Father"), Age = c(data_m$agem %>% mean(), data_m$aged %>% mean()),
                            Income = c(data_m$incomem %>% mean(), data_m$incomed %>% mean())) %>%
  mutate(across(c(Age, Income), ~ round(.x, 2)))

### (b)
d1b <- data_m$morekids # binary: had more than two kids
y1b <- data_m$hourswd # mens hours worked per week
z1b <- data_m$samesex # binary: first two children same sex

os.1b <- lm(d1b ~ z1b)
summ(os.1b, robust = "HC1")

ss.1b <- ivreg(y1b ~ d1b  | z1b)
summary(ss.1b, vcov = sandwich, diagnostics = TRUE)

### (c)
d1c <- data_m$morekids # binary: had more than two kids
y1c <- data_m$hourswd # mens hours worked per week
z1c_boys <- data_m$boys2 # binary: first two births boys
z1c_girls <- data_m$girls2 # binary: first two births girls

os.1c <- lm(d1c ~ z1c_boys + z1c_girls)
summ(os.1c, robust = "HC1")

ss.1c <- ivreg(y1c ~ d1c  | z1c_boys + z1c_girls)
summary(ss.1c, vcov = sandwich, diagnostics = TRUE)

### (d) See PDF

### (e)
d1e <- data_m$morekids # binary: had more than two kids
y1e <- data_m$hourswd # mens hours worked per week

ols.1e <- lm(y1e ~ d1e)
summ(ols.1e, robust = "HC1")

##################
## Question 2   ##
##################
### (a)
median_age_m <- median(data$agem)
data_om <- data %>% filter(agem >= median_age_m) # households where woman is over median age
data_um <- data %>% filter(agem < median_age_m) # households where woman is under median age

data_om$morekids %>% mean() %>% round(3)
data_um$morekids %>% mean() %>% round(3)

### (b)
y2b <- data_um$weeksm
d2b <- data_um$morekids
x2b <- as.matrix(data_um[,c("agem", "agefstm", "blackm", "hispm", "othracem")])

ols.2b <- lm(y2b ~ d2b + x2b)
summ(ols.2b, robust = "HC1")


### (c)
ss_2c <- ivreg(weeksm ~ morekids + agem + agefstm + blackm + hispm + othracem | samesex + agem + agefstm + blackm + hispm + othracem, data = data_um)
summary(ss_2c, vcov = sandwich, diagnostics = TRUE)

### (d) See PDF

##################
## Question 3   ##
##################
### (a)
# Are some households missing mothers or fathers?
data$hourswm %>% is.na() %>% sum()
data$hourswd %>% is.na() %>% sum()

# Check employment rates for men and women
women_er <- mean(na.omit(data$hourswm) > 0) %>% round(3) # women
men_er <- mean(na.omit(data$hourswd) > 0) %>% round(3) # men

# Check how common part time work is for men and women
women_pt <- mean((na.omit(data$hourswm) > 0) & (na.omit(data$hourswm) < 40)) %>% round(3) # women
men_pt <- mean((na.omit(data$hourswd) > 0) & (na.omit(data$hourswd) < 40)) %>% round(3) # men

women_men_df <- data.frame(Measure = c("Emp. Rate", "P.T. rate"),
                           Women = c(women_er, women_pt),
                           Men = c(men_er, men_pt))

print(women_men_df)

### (b) and (c) See PDF


##################
## Question 4   ##
##################
# See PDF

##################
## Question 5   ##
##################
# See PDF
