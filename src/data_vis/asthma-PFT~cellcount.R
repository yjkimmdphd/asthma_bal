########
# data visualization for asthma phenotype and bal cell count
########
library(tidyverse)
library(gridExtra)
library(purrr)
bphen_path<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")

bphen<-read.csv(bphen_path,row.names = NULL)
variables_to_test <- c("asthma_phen_ACT.score", "asthma_phen_FEV1_perc", "asthma_phen_FVC_perc", 
                       "asthma_phen_FEV1.FVC", "asthma_phen_Post.alb_FEV1_perc", 
                       "asthma_phen_Post.alb_FVC_perc", "asthma_phen_Post.alb_FEV1.FVC")

# find how many people have the phenotype data
print(apply(!sapply(bphen[,variables_to_test],is.na),2,sum))

### # ACT,PFT ~ BAL Eos% 1% cutoff
eos_p_1<-bphen%>%mutate(BAL_eos_p_m_1 = BAL_eos_p>1)%>%group_by(BAL_eos_p_m_1)
      eos_p_1<-eos_p_1%>%arrange(by=BAL_eos_p_m_1)%>%select(SampleID,a_phen_cat)
      
      # Ensure BAL_eos_p_m_1 is a factor with exactly two levels
      if(!is.factor(eos_p_1$BAL_eos_p_m_1)) {
        eos_p_1$BAL_eos_p_m_1 <- as.factor(eos_p_1$BAL_eos_p_m_1)
      }
      
      if(length(levels(eos_p_1$BAL_eos_p_m_1)) != 2) {
        stop("BAL_eos_p_m_1 must have exactly two levels for t-tests.")
      }
      
      # Function to perform t-test for a given variable
      perform_t_test <- function(data, variable) {
        t_test_result <- t.test(data[[variable]] ~ data$BAL_eos_p_m_1)
        return(t_test_result$p.value)
      }
      
      # List of variables to test
      variables_to_test <- c("asthma_phen_ACT.score", "asthma_phen_FEV1_perc", "asthma_phen_FVC_perc", 
                             "asthma_phen_FEV1.FVC", "asthma_phen_Post.alb_FEV1_perc", 
                             "asthma_phen_Post.alb_FVC_perc", "asthma_phen_Post.alb_FEV1.FVC")
      
      # Apply the function to each variable and combine results
      t_test_results_all <- map_dfr(variables_to_test, function(var) {
        data.frame(
          Variable = var,
          P_Value = perform_t_test(eos_p_1, var)
        )
      })
      
      
      # Function to create boxplot for a given variable
      create_boxplot <- function(data, variable) {
        ggplot(data, aes(x = BAL_eos_p_m_1, y = .data[[variable]])) +
          geom_boxplot() +
          labs(title = paste("Boxplot of", variable, "by BAL_eos_p_m_1"),
               x = "BAL_eos_p_m_1",
               y = variable)
      }
      
      # Create a list of plots
      boxplots <- map(variables_to_test, ~ create_boxplot(eos_p_1, .x))
      
      # Display the plots together
      do.call("grid.arrange", c(boxplots, ncol = 2))
      
      t_test_results_all


# ACT,PFT ~ BAL Eos% 3% cutoff
bphen_more_than_3<-bphen$BAL_eos_p>3


      # Create the new eos_p_3 dataframe
      eos_p_3 <- bphen %>%
        mutate(BAL_eos_p_m_3 = BAL_eos_p > 3) %>%
        arrange(BAL_eos_p_m_3) %>%
        select(SampleID, BAL_eos_p_m_3, asthma_phen_ACT.score, 
               asthma_phen_FEV1_perc, asthma_phen_FVC_perc, asthma_phen_FEV1.FVC, 
               asthma_phen_Post.alb_FEV1_perc, asthma_phen_Post.alb_FVC_perc, asthma_phen_Post.alb_FEV1.FVC)
      
      # Ensure BAL_eos_p_m_3 is a factor
      eos_p_3$BAL_eos_p_m_3 <- as.factor(eos_p_3$BAL_eos_p_m_3)
      
      # List of variables to test
      variables_to_test <- c("asthma_phen_ACT.score", "asthma_phen_FEV1_perc", "asthma_phen_FVC_perc", 
                             "asthma_phen_FEV1.FVC", "asthma_phen_Post.alb_FEV1_perc", 
                             "asthma_phen_Post.alb_FVC_perc", "asthma_phen_Post.alb_FEV1.FVC")
      
      # Function to perform t-test for a given variable
      perform_t_test <- function(data, variable) {
        t_test_result <- t.test(data[[variable]] ~ data$BAL_eos_p_m_3)
        return(t_test_result$p.value)
      }
      
      # Apply the function to each variable and combine results
      t_test_results_all <- map_dfr(variables_to_test, function(var) {
        data.frame(
          Variable = var,
          P_Value = perform_t_test(eos_p_3, var)
        )
      })
      
      # View t-test results
      print(t_test_results_all)
      
      # Function to create boxplot for a given variable
      create_boxplot <- function(data, variable) {
        ggplot(data, aes(x = BAL_eos_p_m_3, y = .data[[variable]])) +
          geom_boxplot() +
          labs(title = paste("Boxplot of", variable, "by BAL_eos_p_m_3"),
               x = "BAL_eos_p_m_3",
               y = variable)
      }
      
      # Create a list of plots
      boxplots <- map(variables_to_test, ~ create_boxplot(eos_p_3, .x))
      
      # Display the plots together
      do.call("grid.arrange", c(boxplots, ncol = 2))

# ACT,PFT~BAL AEC 1000 above or below
      # Create the new eos_aec_1 dataframe
      eos_aec_1 <- bphen %>%
        mutate(BAL_eos_aec_m_1 = BAL_eos_ct > 1) %>%
        arrange(BAL_eos_aec_m_1) %>%
        select(SampleID, BAL_eos_aec_m_1, asthma_phen_ACT.score, 
               asthma_phen_FEV1_perc, asthma_phen_FVC_perc, asthma_phen_FEV1.FVC, 
               asthma_phen_Post.alb_FEV1_perc, asthma_phen_Post.alb_FVC_perc, asthma_phen_Post.alb_FEV1.FVC)
      
      # Ensure BAL_eos_aec_m_1 is a factor
      eos_aec_1$BAL_eos_aec_m_1 <- as.factor(eos_aec_1$BAL_eos_aec_m_1)
      
      # List of variables to test
      variables_to_test <- c("asthma_phen_ACT.score", "asthma_phen_FEV1_perc", "asthma_phen_FVC_perc", 
                             "asthma_phen_FEV1.FVC", "asthma_phen_Post.alb_FEV1_perc", 
                             "asthma_phen_Post.alb_FVC_perc", "asthma_phen_Post.alb_FEV1.FVC")
      
      # Function to perform t-test for a given variable
      perform_t_test <- function(data, variable) {
        t_test_result <- t.test(data[[variable]] ~ data$BAL_eos_aec_m_1)
        return(t_test_result$p.value)
      }
      
      # Apply the function to each variable and combine results
      t_test_results_all <- map_dfr(variables_to_test, function(var) {
        data.frame(
          Variable = var,
          P_Value = perform_t_test(eos_aec_1, var)
        )
      })
      
      # View t-test results
      print(t_test_results_all)
      
      # Function to create boxplot for a given variable
      create_boxplot <- function(data, variable) {
        ggplot(data, aes(x = BAL_eos_aec_m_1, y = .data[[variable]])) +
          geom_boxplot() +
          labs(title = paste("Boxplot of", variable, "by BAL_eos_aec_m_1"),
               x = "BAL_eos_aec_m_1",
               y = variable)
      }
      
      # Create a list of plots
      boxplots <- map(variables_to_test, ~ create_boxplot(eos_aec_1, .x))
      
      # Display the plots together
      do.call("grid.arrange", c(boxplots, ncol = 2))

      

# ACT PFT ~ blood AEC 300
      # Create the new blood_eos_aec_300 dataframe
      # Create the new blood_eos_aec_300 dataframe
      blood_eos_aec_300 <- bphen %>%
        mutate(blood_eos_aec_m_300 = blood_eos > 300) %>%
        arrange(blood_eos_aec_m_300) %>%
        select(SampleID, blood_eos_aec_m_300, asthma_phen_ACT.score, 
               asthma_phen_FEV1_perc, asthma_phen_FVC_perc, asthma_phen_FEV1.FVC, 
               asthma_phen_Post.alb_FEV1_perc, asthma_phen_Post.alb_FVC_perc, asthma_phen_Post.alb_FEV1.FVC) %>%
        filter(!is.na(blood_eos_aec_m_300))  # Filter out NA values
      
      # Ensure blood_eos_aec_m_300 is a factor
      blood_eos_aec_300$blood_eos_aec_m_300 <- as.factor(blood_eos_aec_300$blood_eos_aec_m_300)
      
      # List of variables to test
      variables_to_test <- c("asthma_phen_ACT.score", "asthma_phen_FEV1_perc", "asthma_phen_FVC_perc", 
                             "asthma_phen_FEV1.FVC", "asthma_phen_Post.alb_FEV1_perc", 
                             "asthma_phen_Post.alb_FVC_perc", "asthma_phen_Post.alb_FEV1.FVC")
      
      # Function to perform t-test for a given variable
      perform_t_test <- function(data, variable) {
        t_test_result <- t.test(data[[variable]] ~ data$blood_eos_aec_m_300, na.action = na.omit)
        return(t_test_result$p.value)
      }
      
      # Apply the function to each variable and combine results
      t_test_results_all <- map_dfr(variables_to_test, function(var) {
        data.frame(
          Variable = var,
          P_Value = perform_t_test(blood_eos_aec_300, var)
        )
      })
      
      # View t-test results
      print(t_test_results_all)
      
      # Function to create boxplot for a given variable
      create_boxplot <- function(data, variable) {
        ggplot(data, aes(x = blood_eos_aec_m_300, y = .data[[variable]])) +
          geom_boxplot(na.rm = TRUE) +
          labs(title = paste("Boxplot of", variable, "by blood_eos_aec_m_300"),
               x = "blood_eos_aec_m_300",
               y = variable)
      }
      
      # Create a list of plots
      boxplots <- map(variables_to_test, ~ create_boxplot(blood_eos_aec_300, .x))
      
      # Display the plots together
      do.call("grid.arrange", c(boxplots, ncol = 2))
      
# ACT~BAL ANC 1000 above or below
neut_anc_1<-bphen%>%mutate(BAL_neut_anc_m_1 = BAL_neut_ct>1)%>%group_by(BAL_neut_anc_m_1)
neut_anc_1<-neut_anc_1%>%arrange(by=BAL_neut_anc_m_1)%>%select(SampleID,BAL_neut_anc_m_1,asthma_phen_ACT.score)

t_result<-t.test(neut_anc_1$BAL_neut_anc_m_1,neut_anc_1$asthma_phen_ACT.score)
boxplot(asthma_phen_ACT.score~BAL_neut_anc_m_1,data=neut_anc_1,
        xlab="BAL ANC >1000",
        ylab="ACT score",
        main=paste("ACT ~ BAL ANC above 1000; t-test p-value",signif(t_result$p.value, digits=2),sep=":"))
# Add the p-value to the plot
text(x = 1.5, 
     y = max(neut_anc_1$asthma_phen_ACT.score) * 0.9, 
     labels = paste("p-value =", format(t_result$p.value, digits = 2)), 
     pos = 4, 
     col = "red")
neut_anc_1%>%summarize(median_act=median(asthma_phen_ACT.score,na.rm=TRUE))


##########correlations 

# correlation ACT ~ AEC (all counts including 0)
par(mfrow=c(2,2))

library(dplyr)

  # Extract the relevant vectors
  x <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(BAL_eos_ct_log)
  y <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal AEC,all)",ylab="ACT score")
  
  result <- cor.test(bphen$BAL_eos_ct_log, bphen$asthma_phen_ACT.score, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_eos_ct_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2

  # Display the correlation coefficient on the plot
  text(x = min(x), 
       y = max(y), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
# correlation ACT ~ AEC (>0)
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(BAL_eos_ct>0,asthma_phen_ACT.score > 0) %>% pull(BAL_eos_ct_log)
  y <- bphen %>% filter(BAL_eos_ct>0,asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal AEC, >0)",ylab="ACT score")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_eos_ct_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  
  # Display the correlation coefficient on the plot
  text(x = min(x), 
       y = max(y), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
# correlation ACT ~ Eos P
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(BAL_eos_p_log)
  y <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal Eos %,all)",ylab="ACT score")
  
  result <- cor.test(bphen$BAL_eos_p_log, bphen$asthma_phen_ACT.score, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_eos_p_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x), 
       y = max(y), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
# correlation ACT ~ Eos %>0
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(BAL_eos_p>0,asthma_phen_ACT.score > 0) %>% pull(BAL_eos_p_log)
  y <- bphen %>% filter(BAL_eos_p>0,asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal Eos %, >0)",ylab="ACT score")
  
  result <- cor.test(bphen$BAL_eos_p_log, bphen$asthma_phen_ACT.score, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_eos_p_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x), 
       y = max(y), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")

  
  
# ACT~BAL ANC

########## correlation ACT ~ ANC (all counts including 0)
  par(mfrow=c(2,2))
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(BAL_neut_ct_log)
  y <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal ANC,all)",ylab="ACT score")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_neut_ct_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  
  # Display the correlation coefficient on the plot
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  # correlation ACT ~ ANC (>0)
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(BAL_neut_ct>0,asthma_phen_ACT.score > 0) %>% pull(BAL_neut_ct_log)
  y <- bphen %>% filter(BAL_neut_ct>0,asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal ANC, >0)",ylab="ACT score")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_neut_ct_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  
  # Display the correlation coefficient on the plot
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  # correlation ACT ~ neut P
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(BAL_neut_p_log)
  y <- bphen %>% filter(asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal neut %,all)",ylab="ACT score")
  
  result <- cor.test(bphen$BAL_neut_p_log, bphen$asthma_phen_ACT.score, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_neut_p_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  # correlation ACT ~ neut %>0
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(BAL_neut_p>0,asthma_phen_ACT.score > 0) %>% pull(BAL_neut_p_log)
  y <- bphen %>% filter(BAL_neut_p>0,asthma_phen_ACT.score > 0) %>% pull(asthma_phen_ACT.score)
  
  # Plot the data
  plot(x, y,xlab="log(bal neut %, >0)",ylab="ACT score")
  
  result <- cor.test(bphen$BAL_neut_p_log, bphen$asthma_phen_ACT.score, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_ACT.score ~ BAL_neut_p_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  

# correlation PFT ~ ANC 
#########here
  # correlation PFT ~ ANC (all counts including 0)
  par(mfrow=c(2,2))
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(asthma_phen_FEV1_perc> 0) %>% pull(BAL_neut_ct_log)
  y <- bphen %>% filter(asthma_phen_FEV1_perc> 0) %>% pull(asthma_phen_FEV1_perc)
  
  # Plot the data
  plot(x, y,xlab="log(bal ANC,all)",ylab="FEV1_%")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_FEV1_perc ~ BAL_neut_ct_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  
  # Display the correlation coefficient on the plot
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  # correlation ACT ~ ANC (>0)
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(BAL_neut_ct>0,asthma_phen_FEV1_perc > 0) %>% pull(BAL_neut_ct_log)
  y <- bphen %>% filter(BAL_neut_ct>0,asthma_phen_FEV1_perc > 0) %>% pull(asthma_phen_FEV1_perc)
  
  # Plot the data
  plot(x, y,xlab="log(bal ANC, >0)",ylab="FEV1_perc")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_FEV1_perc ~ BAL_neut_ct_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  
  # Display the correlation coefficient on the plot
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  # correlation ACT ~ neut P
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(asthma_phen_FEV1_perc > 0) %>% pull(BAL_neut_p_log)
  y <- bphen %>% filter(asthma_phen_FEV1_perc > 0) %>% pull(asthma_phen_FEV1_perc)
  
  # Plot the data
  plot(x, y,xlab="log(bal neut %,all)",ylab="ACT score")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_FEV1_perc ~ BAL_neut_p_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  # correlation ACT ~ neut %>0
  
  
  library(dplyr)
  
  # Extract the relevant vectors
  x <- bphen %>% filter(BAL_neut_p>0,asthma_phen_FEV1_perc > 0) %>% pull(BAL_neut_p_log)
  y <- bphen %>% filter(BAL_neut_p>0,asthma_phen_FEV1_perc > 0) %>% pull(asthma_phen_FEV1_perc)
  
  # Plot the data
  plot(x, y,xlab="log(bal neut %, >0)",ylab="ACT score")
  
  result <- cor.test(x, y, use = "complete.obs")
  
  # Extract the correlation coefficient
  correlation_coefficient <- result$estimate
  
  # Calculate r-squared
  r_squared <- correlation_coefficient^2
  
  # Extract the p-value
  p_value <- result$p.value
  
  # Print the results
  cat("Correlation coefficient:", correlation_coefficient, "\n")
  cat("R-squared:", r_squared, "\n")
  cat("P-value:", p_value, "\n")
  
  model <- lm(asthma_phen_FEV1_perc ~ BAL_neut_p_log, data = bphen)
  
  # Add the regression line
  abline(model, col = "red", lwd = 2)  # Red line with width of 2
  text(x = min(x,na.rm=TRUE), 
       y = max(y,na.rm=TRUE), 
       labels = paste("r =", round(correlation_coefficient, 2),"R-squared:", round(r_squared,2),"P-value:", round(p_value,2)), 
       pos = 4, 
       col = "black")
  
  
  
  
  
  
  
  
  
  