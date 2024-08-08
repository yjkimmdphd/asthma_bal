########
# data visualization for asthma phenotype and bal cell count
########
bphen_path<-file.path("./resources/processed_data/scaled_phenotype_studyID_asthmaPhenotype_batch_cellCount_20240731.csv")

bphen<-read.csv(bphen_path,row.names = NULL)


### BAL Eos% 1% cutoff
bphen_more_than_1<-bphen$BAL_eos_p>1

eos_p_1<-bphen%>%mutate(BAL_eos_p_m_1 = BAL_eos_p>1)%>%group_by(BAL_eos_p_m_1)
eos_p_1<-eos_p_1%>%arrange(by=BAL_eos_p_m_1)%>%select(SampleID,BAL_eos_p_m_1,asthma_phen_ACT.score)

t.test(eos_p_1$BAL_eos_p_m_1,eos_p_1$asthma_phen_ACT.score)
boxplot(asthma_phen_ACT.score~BAL_eos_p_m_1,data=eos_p_1)

eos_p_1%>%summarize(median_act=median(asthma_phen_ACT.score,na.rm=TRUE))

### BAL Eos% 3% cutoff
bphen_more_than_3<-bphen$BAL_eos_p>3

eos_p_3<-bphen%>%mutate(BAL_eos_p_m_3 = BAL_eos_p>3)%>%group_by(BAL_eos_p_m_3)
eos_p_3<-eos_p_3%>%arrange(by=BAL_eos_p_m_3)%>%select(SampleID,BAL_eos_p_m_3,asthma_phen_ACT.score)

t.test(eos_p_3$BAL_eos_p_m_3,eos_p_3$asthma_phen_ACT.score)
boxplot(asthma_phen_ACT.score~BAL_eos_p_m_3,data=eos_p_3)

eos_p_3%>%summarize(median_act=median(asthma_phen_ACT.score,na.rm=TRUE))

# linear regression 
lm_bal_eos_p<-lm(asthma_phen_ACT.score~BAL_eos_p_log, data=bphen)
summary(lm_bal_eos_p)


# ACT~BAL AEC 1000 above or below
eos_aec_1<-bphen%>%mutate(BAL_eos_aec_m_1 = BAL_eos_ct>1)%>%group_by(BAL_eos_aec_m_1)
eos_aec_1<-eos_aec_1%>%arrange(by=BAL_eos_aec_m_1)%>%select(SampleID,BAL_eos_aec_m_1,asthma_phen_ACT.score)

t.test(eos_aec_1$BAL_eos_aec_m_1,eos_aec_1$asthma_phen_ACT.score)
boxplot(asthma_phen_ACT.score~BAL_eos_aec_m_1,data=eos_aec_1)

eos_aec_1%>%summarize(median_act=median(asthma_phen_ACT.score,na.rm=TRUE))


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

#########here
  # correlation ACT ~ ANC (all counts including 0)
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  