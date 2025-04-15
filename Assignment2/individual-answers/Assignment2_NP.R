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

treat_df <- data %>% dplyr::select(time, treatment, distance, connected, submarines, employed, skilled) %>% filter(treatment == 1)
View(treat_df)
treat_df$distance %>% max()
treat_df$employed %>% mean()

untreat_df <- data %>% dplyr::select(time, treatment, distance, connected, submarines, employed, skilled) %>% filter(treatment == 0)
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
did1_clust_summary <- summary(did1, cluster = "location")
cat("Not Clustered SE:", round(did1$se[["treatment"]], 4))
cat("Clustered SE:", round(did1_clust_summary$se[["treatment"]], 4))

# The standard error doubles from about 0.0051 to 0.0101, implies employment (really residuals) within locations are correlated

## (d)
did2 <- feols(employed ~ treatment | time + location, data)
did2_summary <- summary(did2, cluster = "location")
round(did2_summary$se[[1]], 4)
cat("Location and Time Fixed Effects ATET:", round(did2_summary$coefficients, 4))
cat("Location and Time Fixed Effects ATET SE:", round(did2_summary$se[[1]], 4))

# The point estimate corresponds to a change in employment of about 2.2 percentage points, compared to
# about 0.83 percentage points when location and time effects were not accounted for.
# The difference could be due to the fact that the baseline level of employment for each location is different
# and/or the the difference between the treated and untreated are different at different time points (due to e.g. "shocks"
# that affect all locations the same).

####################
## Question 2     ##
####################
## (a)
# If it is believed that fast internet might improve skill acquisition and skill acquisition improve
# employment outcomes, it should not be included. This is because in such a case, skill acquisition is acting as a mediator,
# meaning part of the effect of fast internet on employment is through this variable. As such, by controlling for it,
# we are removing part of the impact of fast internet

## (b)
did3 <- feols(employed ~ treatment + skilled | time + location, data)
did3_summary <- summary(did3, cluster = "location")
did3_summary
cat("Location and Time Fixed Effects with Skilled ATET:", round(did3_summary$coefficients[1], 4))
cat("Location and Time Fixed Effects with Skilled ATET SE:", round(did3_summary$se[[1]], 4))

# The effect did decrease, supporting the hypothesis that part of the effect of fast internet on employment
# is through skill acquisition, suggesting it may be best not to adjust our estimates for it

####################
## Question 3     ##
####################
## (a)
common_trends <- data %>%
  group_by(time, connected) %>%
  summarise(edu_high_mean = mean(educ_high, na.rm = TRUE))

common_trends$connected <- factor(common_trends$connected,
                                  levels = c(0,1),
                                  labels = c("Connected = 0", "Connected = 1"))

View(common_trends)

ggplot(data = common_trends,
       aes(x = time, y = edu_high_mean,
           group = connected,
           color = connected)) +
  geom_line(size = 1.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed") + # why at 0.5?
  scale_x_continuous(breaks = seq(-5, 4, by = 1)) +
  theme_bw(base_size = 14) +
  labs(y = "Average high education rate", x = "Time", colour = "Group")

# Prior to the introduction of submarine cables, the difference in education between the connected and unconnected
# groups were much smaller than after the introduction of submarine cables. However, the difference is not
# constant over time. There could be many reasons for the different trajectories prior to the arrival of high-speed internet -
# anything that affects education and changed over time differently between treated and untreated could explain
# this difference. After the onset of high-speed internet, we see that this difference widens, possibly indicating that
# high-speed internet leads to higher education attainment. Another possibility, however, is that
# highly educated workers moved from the unconnected to the connected locations, where jobs requiring
# higher education became more available due to onset of high-speed internet.


## (b)
did4 <- feols(educ_high ~ i(time, connected, ref = 0) | location + time,
                     cluster = ~location,
                     data = data)
summary(did4)

iplot(did4,
      ref.line = 0,
      xlab = "Time to treatment",
      ylab = "Effect on high education rate",
      main = "Event Study: Fast Internet and High Education")

# We see that other than at time point -5, the ATET is not significant prior to the introduction of submarine cables
# but becomes significant afterwards. Additionally, we see that there is a greater effect at time 2 than at time 1 before
# stabilizing, which is also plausible given that it may take time for firms to adapt to high-speed internet, people to find
# employment, etc

## (c)
# Given the common trends plot, a DiD strategy would likely overestimate the effect. This is because we are assuming that the
# trend observed in the untreated is the trend that the treated group would have had in the absence of treatment. However,
# based on the plot, it appears that the difference in outcomes between the treated and untreated would have continued to
# increase over time, even in the absence of the intervention. This would lead to an overestimation of the impact of high-speed
# internet, as only part of the observed difference after the arrival of submarine cables
# can likely be attributed to high-speed internet.

####################
## Question 4     ##
####################
## (a)

# We have to sample at the cluster level, not individual level, and then take all individuals within a cluster
# It is how you would calculate variance if you really used a clustered sampling strategy
ATET_np_cluster_bootstrap <- function(data, B = 100, seed = 23) {
  set.seed(seed)

  # Find all locations
  cluster_ids <- unique(data[["location"]])
  G <- length(cluster_ids)

  # Initialize vector for storing estimates
  atet_boot <- numeric(B)

  # Run bootstraps B times
  for (b in 1:B) {

    # Resample clusters with replacement
    sampled_clusters <- sample(cluster_ids, size = G, replace = TRUE)

    # Initialize data frame to store data from each cluster
    boot_data <- data.frame()

    # Add data from each cluster to the dataframe
    for (clust_id in sampled_clusters) {
      cluster_rows <- data[data[["location"]] == clust_id, ]
      boot_data <- rbind(boot_data, cluster_rows)
    }

    # Compute ATET for the boostrap sample
    atet_boot[b] <- ATET_np_fn(boot_data)
  }

  return(atet_boot)
}

cluster_bootstrap_sample <- ATET_np_cluster_bootstrap(data, B = 100, seed = 23)

####################
## Question 5     ##
####################
