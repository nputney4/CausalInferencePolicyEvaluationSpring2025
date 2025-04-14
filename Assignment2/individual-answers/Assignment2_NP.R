####################################################################
#   PHASE 0: LOAD PACKAGES AND READ IN DATA
####################################################################

rm(list=ls())


packages_vector <- c("haven", "dplyr", "tidyr", "sandwich", "expss",
                     "fBasics", "xtable", "data.table", "stargazer", "mfx", "jtools", "ggplot2")
# install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE)


# DiD-specific packages
packaged_vector_did <- c("fixest")
# install.packages(packaged_vector_did)
lapply(packaged_vector_did, require, character.only = TRUE)

# List loaded packages
# (.packages())

# Working Directory
work_dir <- "C:/Users/putne/OneDrive/Documents/CausalInferencePolicyEvaluationSpring2025/Assignment2/instructions-and-setup"
setwd(work_dir)

# Read in the data
load("data_assignment.RData")

# Brief examination of the data
glimpse(data)
Hmisc::describe(data)

treat_df <- data %>% dplyr::select(time, treatment, distance, connected, submarines, employed) %>% filter(treatment == 1)
View(treat_df)
treat_df$distance %>% max()
treat_df$employed %>% mean()

untreat_df <- data %>% dplyr::select(time, treatment, distance, connected, submarines, employed) %>% filter(treatment == 0)
View(untreat_df)
untreat_df$distance %>% max()
untreat_df$employed %>% mean()

####################
## Question 1     ##
####################

## (a)
# Point estimate
ATET_np_fn <- function(data){
  treated_before <- data %>% filter(connected == 1, submarines == 0) %>% dplyr::select(employed) %>% unlist() %>% mean()
  treated_after <- data %>% filter(connected == 1, submarines == 1) %>% dplyr::select(employed) %>% unlist() %>% mean()
  untreated_before <- data %>% filter(connected == 0, submarines == 0) %>% dplyr::select(employed) %>% unlist() %>% mean()
  untreated_after <- data %>% filter(connected == 0, submarines == 1) %>% dplyr::select(employed) %>% unlist() %>% mean()

  ATET_np <- (treated_after - treated_before) - (untreated_after - untreated_before)

  return(ATET_np)
}
ATET_np_est <- ATET_np_fn(data)

# Standard errors
B = 100
ATET_boot_sample <- replicate(B, {
  boot_data <- data %>% slice_sample(n = nrow(data), replace = TRUE)
  ATET_np_fn(boot_data)
})

ATET_np_se <- sd(ATET_boot_sample)

cat("ATET Point Estimate:", round(ATET_np_est,5))
cat("ATET Standard Error:", round(ATET_np_se, 5))
cat("ATET 2.5% and 97.5% Quantiles: ", quantile(ATET_boot_sample, c(0.025, 0.975)))


## (b)
did1 <- feols(employed ~ connected + submarines + treatment, data)
summary(did1)

### (i)
# They are similar, it is expected because they are doing the same calculations

### (ii)
# The standard error is also similar

## (c)
did1_clust_se <- summary(did1, cluster = "location")
cat("Not Clustered SE:", round(did1$se[["treatment"]], 4))
cat("Clustered SE:", round(did1_clust_se$se[["treatment"]], 4))

# The standard error doubles from about 0.0051 to 0.0101, implies employment (really residuals) within locations are correlated

## (d)
did2 <- feols(employed ~ treatment | time + location, data)
summary(did2, cluster = "location")

