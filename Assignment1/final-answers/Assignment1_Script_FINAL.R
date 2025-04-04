###############################################################
## -- Causal Inference: Policy Evaluation — Assignment 1' -- ##
###############################################################

########################################################
# Authors:                                            ##
## - Amit Aryal, Nicholas Putney, Dominik Schawrzkopf ##
########################################################


# Define packages that you need
packages_vector <- c("tidyverse",    # package collection for data management
                      "dplyr",       # set of functions to work with data
                      "foreign",     # loading data (e.g. data from old Stata versions)
                      "haven",       # loading data (e.g. data from newer Stata versions)
                      "knitr",       # integrates computing and reporting
                      "fastDummies", # dummy creation
                      "sandwich",    # robust standard errors
                      "lmtest",      # robust standard errors
                      "jtools",      # summary statistics (robust standard errors)
                      "fBasics",     # summary statistics
                      "arsenal",     # summary statistics
                      "data.table",  # tables
                      "stargazer",   # tables
                      "xtable",      # tables
                      "expss",       # tables, labels
                      "sjlabelled",  # labels
                      "lubridate",   # dates
                      "mfx",         # marginal effects for nonlinear models, e.g. probit
                      "ggplot2",
                     "causalweight",
                     "gt")     # plots and graphs

lapply(packages_vector, require, character.only = TRUE)

# Loading given dataset
load("G:/My Drive/SwissTPH/Spring2025ClassesBasel/CausalInferencePolicyEvaluation/04_Small_Assignments/Assignment_1/Assignment_1.RData")

glimpse(raw.data)

## Question 1

# Creating a dummy variable indicating if someone is under 40 years old

# Getting an idea of current variables
raw.data[c("agegr_2", "agegr_3", "agegr_4")] %>% head()

# Seeing if there are rows with 0s for each dummy
raw.data[c("agegr_2", "agegr_3", "agegr_4")] %>% rowSums()

# Finding rows which correspond to those under 40 (note all those under 18 and above 60 were dropped)
raw.data$under_40yo <- ifelse(rowSums(raw.data[c("agegr_3", "agegr_4")]) == 0, 1, 0)

# Getting an idea of current variables
raw.data[c("agegr_2", "agegr_3", "agegr_4", "under_40yo")] %>% head()

# Examining the balance of this dummy variable across treatment groups
cro(raw.data$treat, raw.data$under_40yo)
cro_cpct_responses(raw.data$under_40yo, raw.data$treat)
cro_cpct_responses(raw.data$treat, raw.data$under_40yo)

attach(raw.data)

# Calculating the standardize bias and difference in means

# a selection of controls
x_desc <- cbind(under_40yo)
x_desc_names <- colnames(x_desc)

# Define a function estimating the differences in variables across D
balance_check.model <- function(x){

  # Conditional means
  mean_d0 <- mean(x[treat==0])
  mean_d1 <- mean(x[treat==1])

  # Variances in subsamples
  var_d0 <- var(x[treat==0])
  var_d1 <- var(x[treat==1])

  # Difference in means
  diff_d <- lm(x ~ treat)
  cov <- vcovHC(diff_d, type = "HC")
  robust.se <- sqrt(diag(cov))

  # Absolute standardized bias (difference in means over average std dev - assuming similar size of the groups)
  sb <- abs((mean_d1-mean_d0)/sqrt((var_d0+var_d1)/2))*100

  # Store output as a list
  list(mean_d0 = mean_d0,
       mean_d1 = mean_d1,
       diff_d = diff_d$coefficients[2],
       robust.se = robust.se[2],
       pval = 2*pnorm(-abs(diff_d$coefficients[2]/robust.se[2])),
       SB = sb)
}

diff_output <- apply(X = x_desc, MARGIN = 2, FUN = balance_check.model)

# Convert output in list format into a data frame
diff_output <-as.data.frame(rbindlist(diff_output))

# Add number of observations (don't forget the comma!)
raw.data$treat <- as.vector(raw.data$treat)
obs <- c(nrow(raw.data[raw.data$treat==0,]),
         nrow(raw.data[raw.data$treat==1,]),
         NA, NA, NA, NA)
diff_output <- rbind(diff_output, obs)

# Display in desired format
rownames(diff_output)<- c(x_desc_names, "Observations")
colnames(diff_output)<- c("E(X|D=0)", "E(X|D=1)", "Difference", "s.e.",
                          "p-value", "Abs. SB")
print("Difference in means by treatment status and standardized bias")
print(round(diff_output,digits = 3))

# Convert rownames into a column for better display
diff_output$Variable <- rownames(diff_output)
diff_output <- diff_output %>%
  dplyr::select(Variable, everything())

# Format as a gt table
diff_output_gt <- diff_output %>%
  gt() %>%
  tab_header(
    title = "Balance Check by Treatment Status",
    subtitle = "Difference in Means and Absolute Standardized Bias"
  ) %>%
  fmt_number(
    columns = c(`E(X|D=0)`, `E(X|D=1)`, Difference, `s.e.`, `p-value`, `Abs. SB`),
    decimals = 3
  ) %>%
  cols_label(
    Variable = "Variable",
    `E(X|D=0)` = "E(X | D = 0)",
    `E(X|D=1)` = "E(X | D = 1)",
    Difference = "Difference",
    `s.e.` = "Std. Error",
    `p-value` = "p-value",
    `Abs. SB` = "Abs. Std. Bias"
  ) %>%
  tab_options(
    table.font.size = "small",
    heading.title.font.size = "medium",
    heading.subtitle.font.size = "small"
  )

# Display the table
diff_output_gt

# Dummy variable under_40yo seems well balanced between treatment and control groups

# Question 2a
# Set the horizon
maxdur <- 24

# Create a matrix with n rows and 24 columns
emp <- matrix(0, nrow = nrow(raw.data), ncol = maxdur)

# Create column names like emp_1 to emp_24
emp_list <- paste("emp", 1:maxdur, sep = "_")
colnames(emp) <- emp_list

# Fill in the matrix
for (i in 1:maxdur) {
  emp[, i] <- ifelse(
    raw.data$date_end < raw.data$date_start + 30 * i,
    1,  # employed in month i
    0   # still unemployed in month i
  )
}

# Convert to data frame and combine with original data
emp <- as.data.frame(emp)
raw.data <- cbind(raw.data, emp)

# Create the outcome variable - "monthly employment probability"

raw.data[c("id", "idobs", "date_start", "date_end", "treat", paste0("emp_", 1:24))] %>% head()
raw.data$month_employ_prob <- raw.data %>% dplyr::select(paste0("emp_", 1:24)) %>% rowMeans()

# Inspect distribution of monthly employment probability (for each spell)
hist(raw.data$month_employ_prob, main = "Monthly Employment Probability", xlab = "Probability")
hist(subset(raw.data$month_employ_prob, treat == 0), main = "Monthly Employment Probability - Untreated", xlab = "Probability")
hist(subset(raw.data$month_employ_prob, treat == 1), main = "Monthly Employment Probability - Treated", xlab = "Probability")

# Sanity check for employment duration length
raw.data$duration <- raw.data$date_end - raw.data$date_start
hist(as.numeric(raw.data$duration), xlab = "days", main = "Length of Unemployment Spells")
hist(subset(as.numeric(raw.data$duration), treat == 0), xlab = "days", main = "Length of Unemployment Spells - Untreated")
hist(subset(as.numeric(raw.data$duration), treat == 1), xlab = "days", main = "Length of Unemployment Spells - Treated")

# Question 2b

# Splitting data into under 40 years old and over 40 years old
data_u40 <- raw.data %>% filter(under_40yo == 1)
data_o40 <- raw.data %>% filter(under_40yo == 0)

treat_u40 <- data_u40$treat
treat_o40 <- data_o40$treat

x_u40 <- as.matrix(dplyr::select(data_u40,
                              covs))
x_o40 <- as.matrix(dplyr::select(data_o40,
                              covs))
x_names <- colnames(x_u40)

# Estimate the p-score model for under 40
pscore.model_u40 <- glm(treat_u40 ~ x_u40, family = binomial(link = "probit"))
summ(pscore.model_u40, robust = "HC1")

# Estimate the p-score model for over 40
pscore.model_o40 <- glm(treat_o40 ~ x_o40, family = binomial(link = "probit"))
summ(pscore.model_o40, robust = "HC1")

data_u40$pscore <- pscore.model_u40$fitted.values
summary(data_u40$pscore)

data_o40$pscore <- pscore.model_o40$fitted.values
summary(data_o40$pscore)

# Create a factor variable for treatment status for the plot
data_u40$treat_f <- factor(treat_u40,
                           levels = c(0,1),
                           label = c("D=0", "D=1"))

# Density plot for the propensity score by treatment status
ggplot(data_u40, aes(x = pscore, fill = treat_f)) +
    geom_density(alpha=0.4) +
    scale_fill_grey()+
    theme_bw(base_size = 20) +
    xlim(0, 1) + ggtitle("P-scores for under 40 years of age")


# Create a factor variable for treatment status for the plot
data_o40$treat_f <- factor(treat_o40,
                           levels = c(0,1),
                           label = c("D=0", "D=1"))

# Density plot for the propensity score by treatment status
ggplot(data_o40, aes(x = pscore, fill = treat_f)) +
    geom_density(alpha=0.4) +
    scale_fill_grey()+
    theme_bw(base_size = 20) +
    xlim(0, 1) + ggtitle("P-scores for over 40 years of age")

# As it is the ATE, the min and max pscores for both treated and not treated should be similar - i.e. "common support"
min(data_u40$pscore[data_u40$treat == 1])
max(data_u40$pscore[data_u40$treat == 1])

min(data_u40$pscore[data_u40$treat == 0])
max(data_u40$pscore[data_u40$treat == 0])

trim_pscores <- function(data) {
  # Compute min and max pscore for treated and control groups
  treated_min <- min(data$pscore[data$treat == 1], na.rm = TRUE)
  treated_max <- max(data$pscore[data$treat == 1], na.rm = TRUE)
  control_min <- min(data$pscore[data$treat == 0], na.rm = TRUE)
  control_max <- max(data$pscore[data$treat == 0], na.rm = TRUE)

  # Define common support
  min_overlap <- max(treated_min, control_min)
  max_overlap <- min(treated_max, control_max)

  # Trim the dataset
  trimmed_data <- data[data$pscore >= min_overlap & data$pscore <= max_overlap, ]

  return(trimmed_data)
}

data_u40_trimmed <- trim_pscores(data_u40)
data_o40_trimmed <- trim_pscores(data_o40)

# Only 1 observation removed from each
nrow(data_u40) - nrow(data_u40_trimmed)
nrow(data_o40) - nrow(data_o40_trimmed)

# Checking for aligment
min(data_u40_trimmed$pscore[data_u40_trimmed$treat == 1])
max(data_u40_trimmed$pscore[data_u40_trimmed$treat == 1])

min(data_u40_trimmed$pscore[data_u40_trimmed$treat == 0])
max(data_u40_trimmed$pscore[data_u40_trimmed$treat == 0])

# Density plot for the propensity score by treatment status
ggplot(data_u40_trimmed, aes(x = pscore, fill = treat_f)) +
    geom_density(alpha=0.4) +
    scale_fill_grey()+
    theme_bw(base_size = 20) +
    xlim(0, 1) + ggtitle("P-scores for under 40 years of age - Trimmed")

# Density plot for the propensity score by treatment status
ggplot(data_o40_trimmed, aes(x = pscore, fill = treat_f)) +
    geom_density(alpha=0.4) +
    scale_fill_grey()+
    theme_bw(base_size = 20) +
    xlim(0, 1) + ggtitle("P-scores for over 40 years of age - Trimmed")

# Question 2c
# Define the outcome - "monthly employment probability"
y1_u40 <- data_u40_trimmed$month_employ_prob
y1_o40 <- data_o40_trimmed$month_employ_prob
x_u40 <- as.matrix(data_u40_trimmed[covs])
x_o40 <- as.matrix(data_o40_trimmed[covs])
treat_u40 <- data_u40_trimmed$treat
treat_o40 <- data_o40_trimmed$treat

s_boot <- 100
length(y1_u40) == length(treat_u40)
length(y1_u40) == dim(x_u40)[1]

length(y1_o40) == length(treat_o40)
length(y1_o40) == dim(x_o40)[[1]]

# Estimating the ATE

# Create function to do IPW for each month
reg_monthly_ipw <- function(y, d, x) {
  ipw <- treatweight(y = y, d = d, x = x, ATET = FALSE, trim = 0, boot = s_boot)
  list(effect = ipw$effect, se = ipw$se)
}


# For over 40 years of age

# Create a matrix for employment status for each month
emp_o40 <- data_o40_trimmed %>% dplyr::select(colnames(emp)) %>% as.matrix()

# Apply IPW each month
ipw_monthly_o40 <- apply(emp_o40, 2, function(y) reg_monthly_ipw(y, treat_o40, x_o40))

# Look at results
ipw_monthly_o40_df <- rbindlist(ipw_monthly_o40)
ipw_monthly_o40_df$month <- 1:maxdur

# Add confidence intervals
ipw_monthly_o40_df <- ipw_monthly_o40_df %>%
  mutate(
    cil = effect - 1.96 * se,
    cih = effect + 1.96 * se,
    sig = ifelse(abs(effect / se) > 1.64, effect, NA)
  )

# For under 40 years of age

# Create a matrix for employment status for each month
emp_u40 <- data_u40_trimmed %>% dplyr::select(colnames(emp)) %>% as.matrix()

# Apply IPW each month
ipw_monthly_u40 <- apply(emp_u40, 2, function(y) reg_monthly_ipw(y, treat_u40, x_u40))

# Look at results
ipw_monthly_u40_df <- rbindlist(ipw_monthly_u40)
ipw_monthly_u40_df$month <- 1:maxdur

# Add confidence intervals
ipw_monthly_u40_df <- ipw_monthly_u40_df %>%
  mutate(
    cil = effect - 1.96 * se,
    cih = effect + 1.96 * se,
    sig = ifelse(abs(effect / se) > 1.64, effect, NA)
  )

ipw_monthly_u40_df$Group <- "Under 40"
ipw_monthly_o40_df$Group <- "Over 40"
ipw_monthly_df <- rbind(ipw_monthly_u40_df, ipw_monthly_o40_df)
ipw_monthly_df <- ipw_monthly_df %>% dplyr::select(Group, month, everything())

# Create tables by group
ipw_monthly_u40_gt <- ipw_monthly_df %>%
  filter(Group == "Under 40") %>%
  gt() %>%
  tab_header(
    title = "Monthly IPW Estimates: Under 40",
    subtitle = "Inverse Probability Weighting (ATE Estimates per Month)"
  ) %>%
  fmt_number(columns = c(effect, se, cil, cih, sig), decimals = 3) %>%
  cols_label(
    month = "Month",
    effect = "ATE",
    se = "Std. Error",
    cil = "CI Lower",
    cih = "CI Upper",
    sig = "Significant Effects"
  ) %>%
  tab_options(
    table.font.size = px(10),
    data_row.padding = px(1))

ipw_monthly_o40_gt <- ipw_monthly_df %>%
  filter(Group == "Over 40") %>%
  gt() %>%
  tab_header(
    title = "Monthly IPW Estimates: Over 40",
    subtitle = "Inverse Probability Weighting (ATE Estimates per Month)"
  ) %>%
  fmt_number(columns = c(effect, se, cil, cih, sig), decimals = 3) %>%
  cols_label(
    month = "Month",
    effect = "ATE",
    se = "Std. Error",
    cil = "CI Lower",
    cih = "CI Upper",
    sig = "Significant Effects"
  ) %>%
  tab_options(
    table.font.size = px(10),
    data_row.padding = px(1))

ipw_monthly_u40_gt
ipw_monthly_o40_gt

gtsave(ipw_monthly_u40_gt, "ipw_under40_table.pdf")
gtsave(ipw_monthly_o40_gt, "ipw_over40_table.pdf")


## Question 3

# Plotting results by month by age group
ipw_monthly_u40_df$Group <- "Under 40"
ipw_monthly_o40_df$Group <- "Over 40"
ipw_monthly_all <- bind_rows(ipw_monthly_u40_df, ipw_monthly_o40_df)

ggplot(ipw_monthly_all, aes(x = month, y = effect, color = Group, fill = Group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = cil, ymax = cih), alpha = 0.2, color = NA) +
  geom_point(aes(y = sig), shape = 18, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 18) +
  labs(
    title = "Monthly ATEs on Employment Probability by Age Group",
    x = "Months after program start",
    y = "Average Treatment Effect"
  )

