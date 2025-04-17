####################################
## LOAD PACKAGES AND READ IN DATA ##
####################################
rm(list=ls())

packages_vector <- c("haven", "dplyr", "tidyr", "sandwich", "expss",
                     "fBasics", "xtable", "data.table", "stargazer", "mfx", "jtools", "ggplot2")
# install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE)

# DiD-specific packages
packaged_vector_did <- c("fixest")
# install.packages(packaged_vector_did)
lapply(packaged_vector_did, require, character.only = TRUE)

# Working directory
work_dir <- "C:/Users/putne/OneDrive/Documents/CausalInferencePolicyEvaluationSpring2025/Assignment2/instructions-and-setup"
setwd(work_dir)

# Read in the data
load("data_assignment.RData")

############################
## BRIEF DATA EXPLORATION ##
############################
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

################################################################################
## QUESTION 1 - ESTIMATING THE EFFECT OF FAST INTERNET ON EMPLOYMENT          ##
################################################################################

###############################################################
## (a) - Nonparametric estimate with boostrap standard error ##
###############################################################
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

###############################
## (b) - Parametric estimate ##
###############################
did1 <- feols(employed ~ connected + submarines + treatment, data)
summary(did1) # estimates are identical/se are similar - same means are being compared

##############################################################
## (c) - Parametric estimate with clustered standard errors ##
##############################################################
did1_clust_summary <- summary(did1, cluster = "location")
cat("Not Clustered SE:", round(did1$se[["treatment"]], 4))
cat("Clustered SE:", round(did1_clust_summary$se[["treatment"]], 4)) # se increases

#########################################################################
## (d) - Parametric estimate including location and time fixed effects ##
#########################################################################
did2 <- feols(employed ~ treatment | time + location, data)
did2_summary <- summary(did2, cluster = "location")
round(did2_summary$se[[1]], 4)
cat("Location and Time Fixed Effects ATET:", round(did2_summary$coefficients, 4)) # much larger point estimate
cat("Location and Time Fixed Effects ATET SE:", round(did2_summary$se[[1]], 4))


################################################################################
## QUESTION 2 - CONTROLLING FOR BEING SKILLED                                 ##
################################################################################

#########################################################
## (a) - Assessing if skilled should be controlled for ##
#########################################################
# => See written answers in PDF

########################################################
## (b) - Parametric estimate, controlling for skilled ##
########################################################
# Still including time and location fixed effects
did3 <- feols(employed ~ treatment + skilled | time + location, data)
did3_summary <- summary(did3, cluster = "location")
did3_summary
cat("Location and Time Fixed Effects with Skilled ATET:", round(did3_summary$coefficients[1], 4))
cat("Location and Time Fixed Effects with Skilled ATET SE:", round(did3_summary$se[[1]], 4))


################################################################################
## QUESTION 3 - EFFECT OF FAST INTERNET ON EDUCATION                          ##
################################################################################

###############################################
## (a) - Examining common trends assumption  ##
###############################################
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

##################################
## (b) - Conducting event study ##
##################################
# Again still location and time fixed effects but now an ATET is estimated
# for each time point
did4 <- feols(educ_high ~ i(time, connected, ref = 0) | location + time,
              cluster = ~location,
              data = data)
summary(did4)

iplot(did4,
      ref.line = 0,
      xlab = "Time to treatment",
      ylab = "Effect on high education rate",
      main = "Event Study: Fast Internet and High Education")

##########################################
## (c) - Potential bias of DiD estimate ##
##########################################
# => See written answers in PDF

################################################################################
## QUESTION 4 - NON-PARAMETRIC CLUSTERED BOOSTRAP SE                          ##
################################################################################

#####################################################################
## (a) - Adjusting boostrap se function to account for clustering  ##
#####################################################################

# We have to sample at the cluster level, not individual level, and then take all individuals within a cluster
# It is how you would calculate variance if you really used a clustered sampling strategy
ATET_np_cluster_bootstrap <- function(data, B = 100, seed = 23) {
  set.seed(seed)

  cluster_list <- split(data, data$location)
  cluster_ids <- names(cluster_list)
  G <- length(cluster_ids)

  atet_boot <- numeric(B)

  for (b in 1:B) {
    # Sampling clusters with replacement
    sampled_clusters <- sample(cluster_ids, size = G, replace = TRUE)

    boot_data <- do.call(rbind, cluster_list[sampled_clusters])

    atet_boot[b] <- ATET_np_fn(boot_data)
  }

  return(atet_boot)
}

cluster_bootstrap_sample <- ATET_np_cluster_bootstrap(data, B = 100, seed = 23)

cluster_boostrap_se <- sd(cluster_bootstrap_sample)
cat("Clustered SE Bootstrap:", round(cluster_boostrap_se, 6))

################################################
## Question 5 - Next lab session preparation ##
###############################################
# => See written answers in PDF

################################################################################
## SUMMARY TABLE OF ATETS AND STANDARD ERRORS                                 ##
################################################################################
ATET_summary <- data.frame(
  Question = c("1a", "1b", "1c", "1d", "2b", "4a"),
  Description = c(
    "Nonparametric, not clustered",
    "Parametric, not clustered",
    "Parametric, clustered SE",
    "Parametric, location & time FE",
    "Parametric, location & time FE + skilled",
    "Nonparametric, clustered SE"
  ),
  ATET_Estimate = c(
    round(ATET_np_est, 6),
    round(coef(did1)["treatment"], 6),
    round(coef(did1)["treatment"], 6),
    round(coef(did2)["treatment"], 6),
    round(coef(did3)["treatment"], 6),
    round(ATET_np_est, 6)
  ),
  SE_Estimate = c(
    round(ATET_np_se, 6),
    round(did1$se["treatment"], 6),
    round(did1_clust_summary$se["treatment"], 6),
    round(did2_summary$se["treatment"], 6),
    round(did3_summary$se["treatment"], 6),
    round(cluster_boostrap_se, 6)
  )
)
print(ATET_summary)

# Saving summary table to latex to then be put in final PDF
stargazer(ATET_summary, type="latex", summary = FALSE, digits = 4, rownames = FALSE)
