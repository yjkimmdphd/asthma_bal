# previous attepmpt 

  # sample_list<-vector(mode="list",length=2)
  # 
  # n.samples<-grepl("^N",colnames(counts))
  # n.samples<-which(n.samples==TRUE)
  # n.counts<-counts[,c(1,n.samples)]
  # 
  # colnames(n.counts)[2:length(n.samples)]<-substr(colnames(n.counts)[2:length(n.samples)],1,4)
  # head(n.counts)
  # 
  # sample_list[[1]]<-colnames(n.counts)[2:length(n.counts)]
  # 
  # b.samples<-grepl("^B",colnames(counts))
  # b.samples<-which(b.samples==TRUE)
  # b.counts<-counts[,c(1,b.samples)]
  # 
  # colnames(b.counts)[2:length(b.samples)]<-substr(colnames(b.counts)[2:length(b.samples)],1,4)
  # head(b.counts)
  # 
  # sample_list[[2]]<-colnames(b.counts)[2:length(b.counts)]
  # 
  # names(sample_list)<-c("nasal_samples","bronch_samples")
  # 
  # 
  # # asthma biomarker phenotype file saved in  'phenotype'
  # phenotype<-file.path("./resources/processed_data/nasal_biomarker_phenotype_batch12346_merged.csv")
  # phenotype<-if(file.exists(phenotype)){read.csv(phenotype, row.names = NULL)}
  # 
  # phen_n_b<-vector(mode="list",length=2)
  # 
  # subsetPhen<-function(d,index,sample.list,old_prefix,new_prefix){
  #   phen<-d[gsub(old_prefix, new_prefix, d$SampleID)%in%sample.list[[index]],]
  #   phen$SampleID<-gsub(old_prefix, new_prefix, phen$SampleID)
  #   phen<-mutate_at(phen,vars(all_of(source.cell.log)),scale)
  #   return(phen)
  # }
  # 
  # phen_n_b[[1]]<-subsetPhen(phenotype,1,sample_list,"B","N")
  # phen_n_b[[2]]<-subsetPhen(phenotype,2,sample_list,"B","B")
  # 
  # phen_n_b<-lapply(phen_n_b,function(data){
  #   mutate(data,swab_date=Date.of.nasal.swabs.%>%as.Date(format="%m/%d/%Y"), 
  #          Blood_draw_date=Date.of.blood.draw%>%as.Date(format="%m/%d/%Y"), 
  #          BAL_date=BAL_date%>%as.Date(format="%m/%d/%Y"), 
  #          swab_blood_delay=as.Date(swab_date,format="%m/%d/%Y")-as.Date(Blood_draw_date,format="%m/%d/%Y"),
  #          swab_bal_delay=as.Date(swab_date,format="%m/%d/%Y")-as.Date(BAL_date,format="%m/%d/%Y"))
  # })
  # ################################################
  # BAL_asthma_sampling_dates<-vector(mode="list",length=2)
  # BAL_asthma_sampling_dates<-lapply(phen_n_b,function(data){
  #   data[,c("Study.ID","SampleID","swab_date","Blood_draw_date","BAL_date","swab_blood_delay","swab_bal_delay")]
  # })
  # 
  # write.csv(rbind(BAL_asthma_sampling_dates[[1]],BAL_asthma_sampling_dates[[2]]),"./resources/processed_data/nasal-swab_blood-draw_BAL-date_delay.csv")

# new attempt
  # Load necessary libraries
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  
  # Load your data
  data <- read.csv("./resources/processed_data/sampling_dates/MS_nasalswb-bal-blood-draw_dates_filtered.csv")
  
  
  
  # Convert date columns to Date type with error handling
  data <- data %>%
    mutate(
      Date_of_nasal_swabs = mdy(Date_of_nasal_swabs._, quiet = TRUE),
      BAL_Date = mdy(BAL_Date, quiet = TRUE),
      Blood_Draw_Date_CBC = mdy(Blood_Draw_Date_.CBC., quiet = TRUE),
      Date_of_blood_draw_IgE = mdy(Date_of_blood_draw_.IgE., quiet = TRUE)
    ) %>%
    select(-Date_of_nasal_swabs._, -Blood_Draw_Date_.CBC., -Date_of_blood_draw_.IgE.) # Drop original character columns
  
  # Print column names to ensure correctness
  print(colnames(data))
  # Calculate specific pairwise differences
  data <- data %>%
    mutate(
      bal_nasal = as.numeric(difftime(BAL_Date, Date_of_nasal_swabs, units = "days")),
      blood_nasal = as.numeric(difftime( Blood_Draw_Date_CBC,Date_of_nasal_swabs, units = "days")),
      blood_bal = as.numeric(difftime(Blood_Draw_Date_CBC, BAL_Date, units = "days"))
    )
  
  # Reshape the data for plotting
  plot_data <- data %>%
    select(Study_ID, bal_nasal, blood_nasal, blood_bal) %>%
    pivot_longer(cols = c(bal_nasal, blood_nasal, blood_bal), names_to = "Comparison", values_to = "Difference")
  
  # Sort the data by differences
  plot_data <- plot_data %>%
    arrange(Difference)%>%
    mutate(Study_ID=factor(Study_ID,levels=unique(Study_ID)))
  
  # Plot the differences as a dot plot, sorted by differences
  ggplot(plot_data, aes(x = reorder(Study_ID, Difference), y = Difference/365, color = Comparison)) +
    geom_point(size=3) +
    geom_hline(yintercept = c(1,-1), color='red')+
    geom_label_repel(data=plot_data%>%filter(Difference/365>1),aes(label=Study_ID),nudge_y = 1, show.legend = FALSE)+
    geom_label_repel(data=plot_data%>%filter(Difference/365<(-1)),aes(label=Study_ID),nudge_y = -1,show.legend = FALSE)+
    scale_y_continuous(breaks = seq(-3,3, by=0.5)) +
    labs(title = "Pairwise Differences Between Sampling Dates (nasal swab, BAL, CBC) \nfor Each Study_ID",
         x = "Study_ID (sorted by difference)",
         y = "Difference in Years",
         color = "Comparison") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  write.table(plot_data,"./resources/processed_data/sampling_dates/swab-bal-cbc_differences_in_days.txt", sep="\t",row.names = FALSE)
