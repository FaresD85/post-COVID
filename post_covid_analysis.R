### 31/01/2024
#### post covid biomarkers study   
###Fares Darawshy


### set wd 
setwd("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/post_covid")

library(patchwork)
library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(ade4)
library("matrixStats")
library(RColorBrewer)
library(data.table)
library(knitr)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(decontam)
library(DESeq2)
library(microbiomeMarker)
library(stringr)
library(ggprism)
library(forcats)
library(rstatix)
library(patchwork)
library(rmarkdown)
library(ggtext)
library(glue)
library(tidyverse)
library(edgeR)
library(gplots)
library(randomForest)
library(caret)
library(table1)
library(mdatools)
library(glmnet)

################################################################################

#read count data 
mycounts <- read.delim2(file = "biomarkers_res.txt", row.names = 1)

#read metadata
metadata <- read.csv(file = "post_covid_metadata_upd.csv")

#update dsiease severity to create three groups only 
metadata <- metadata %>% 
  mutate(disease_severity_upd=case_when(
    disease_severity == "mild" ~ "mild_mod",
    disease_severity == "moderate" ~ "mild_mod",
    disease_severity=="severe" ~"severe", 
    disease_severity=="critical" ~"critical"
  ))

#update dsiease severity to create two groups only 
metadata <- metadata %>% 
  mutate(disease_severity_two=case_when(
    disease_severity == "mild" ~ "mild_mod",
    disease_severity == "moderate" ~ "mild_mod",
    disease_severity=="severe" ~"severe", 
    disease_severity=="critical" ~"severe"
  ))

rownames(metadata) <- metadata$sample

#transpose mycounts to check levels across samples types 
mycounts_t <- t(mycounts)
mycounts_t <- data.frame(mycounts_t)

#add sample type variable 
mycounts_t$sample <- rownames(mycounts_t)

mycounts_t[,2:10] <- lapply(mycounts_t[,2:10], as.numeric)
#plot comparison between samples types for each molecule 

setdiff(rownames(mycounts_t), rownames(metadata))



###### start with building table 1 for the cohort 

table_1_dat <- metadata %>% filter(sample_type=="baseline")


#add statistical analysis column 
#define p value function 

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times = sapply(x, length)))
  
  if (is.numeric(y)) {
    # For numeric variables, perform an analysis of variance (ANOVA)
    p <- "NA"
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # Format the p-value
  formatted_p <- ifelse(is.numeric(p), sub("<", "&lt;", format.pval(p, digits = 3, eps = 0.001)), NA)
  
  # Return the formatted p-value
  c("", formatted_p)
}

#work on baseline table for all patients 

colnames(table_1_dat)
#create table 1 using table 1 package 
table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         as.numeric(Age)+as.factor(Gender)+as.factor(Origin)+as.factor(smoking)+as.numeric(PY)+
         as.numeric(time_to_fu)+as.numeric(duration_of_admisson)+ as.factor(Death)+as.factor(HTN)+as.factor(DM)+as.factor(HLP)+as.factor(IHD)+as.factor(CHF)+
         as.factor(CKD)+as.factor(transplant)+as.factor(cancer_y_n)+as.factor(Liver.Disease)+as.factor(COPD)+as.factor(Asthma)+as.factor(other_lung_disease)+
         as.factor(Vaccination.status)+as.factor(disease_severity)+as.factor(Intubated)+
         as.numeric(sao2_ra)+as.numeric(Crea_baseline)+as.numeric(Crea_peak)+as.numeric(Crea_discharge)+as.numeric(LDH_peak)+
         as.numeric(WBC_peak)+as.numeric(DDM_peak)+as.numeric(ferritin_peak)+as.numeric(crp_peak)+as.factor(ICU_y_n)
         
       #set data 
       ,data=table_1_dat, 
       #set stat display options for continous variables (options from stat.default)
       render.continuous = c((.="Median [Q1, Q3]")))

### get follow up table 
table_1_dat_fu <- metadata %>% filter(sample_type=="clinic")

colnames(table_1_dat_fu)
#create table 1 using table 1 package 
table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         as.numeric(Age)+as.factor(Gender)+as.factor(Origin)+as.factor(smoking)+as.numeric(PY)+
         as.numeric(time_to_fu)+as.numeric(duration_of_admisson)+ as.factor(Death)+as.factor(HTN)+as.factor(DM)+as.factor(HLP)+as.factor(IHD)+as.factor(CHF)+
         as.factor(CKD)+as.factor(transplant)+as.factor(cancer_y_n)+as.factor(Liver.Disease)+as.factor(COPD)+as.factor(Asthma)+as.factor(other_lung_disease)+
         as.factor(Vaccination.status)+as.factor(disease_severity)+as.factor(Intubated)+
         as.numeric(sao2_ra)+as.numeric(Crea_baseline)+as.numeric(Crea_peak)+as.numeric(Crea_discharge)+as.numeric(LDH_peak)+
         as.numeric(WBC_peak)+as.numeric(DDM_peak)+as.numeric(ferritin_peak)+as.numeric(crp_peak)+
         as.factor(Fatigue_fu)+as.factor(Fever_fu)+as.factor(Headache_fu)+as.factor(Dyspnea_fu)+as.factor(Cough_fu)+
         as.factor(Chest_pain_fu)+as.factor(Nausea_fu)+as.factor(Diarrhea_fu)+as.factor(Sore_throat_fu)+as.factor(Nasal_congestion_fu)+as.factor(Red_eyes_fu)+
         as.factor(Vision_loss_fu)+as.factor(Concentraion_difficulties_fu)+as.factor(Memory_decline_fu)+as.factor(Weight_loss_fu)+as.factor(Arthalgia_fu)+
         as.factor(Parasthesia_Persistent_pain_fu)+ as.factor(Hair_Loss_fu)+as.factor(Aguesia_Dysguesia_fu)+
         as.numeric(Total_Symptoms_score_fu)++as.numeric(sao2_ra_fu)+
         as.numeric(FEV1_lit)+as.numeric(FEV1_p)+as.numeric(FVC_lit)+as.numeric(FVC_p)+as.numeric(FEV1_FVC)+as.numeric(TLC_lit)+
         as.numeric(TLC_p)+as.numeric(RV_lit)+as.numeric(RV_p)+as.numeric(RV_TLC)+as.numeric(KCO)+as.numeric(TLC_lit)+as.factor(abnormal_CT_fu_y_n)
       
       
       #set data 
       ,data=table_1_dat_fu, 
       #set stat display options for continous variables (options from stat.default)
       render.continuous = c((.="Mean (SD)")))


#plot prevelance of symptoms in follow up 
table_1_dat_fu$


# Assuming 'symptoms_df' is your data frame with symptoms columns

# Convert 'yes' and 'no' to 1 and 0
symptoms_df_numeric <- table_1_dat_fu %>% select(ends_with("fu")) %>% select(-c(time_to_fu, blood_fu, Total_Symptoms_score_fu, sao2_ra_fu, 
                                                                                CXR_fu, CT_fu))
colnames(symptoms_df_numeric) <- gsub("_fu", "", colnames(symptoms_df_numeric))
colnames(symptoms_df_numeric) <- gsub("_", " ", colnames(symptoms_df_numeric))


# Calculate the prevalence of each symptom
symptoms_prevalence <- colMeans(symptoms_df_numeric, na.rm = TRUE)

# Create a data frame for plotting
plot_data <- data.frame(symptom = names(symptoms_prevalence), prevalence = symptoms_prevalence)
plot_data <- plot_data[order(plot_data$prevalence, decreasing = TRUE), ]

#export table 
write_csv(plot_data, file = "Results/symptoms_prevalence_at_follow_up.csv")

# Plot the barplot with symptoms ordered by prevalence
pdf(file = "Figures/symptoms_at_fu.pdf", height = 12, width = 16)
ggplot(plot_data, aes(x = reorder(symptom, -prevalence), y = prevalence*100)) +
  geom_bar(stat = "identity", fill = "lightgrey", color = "black") +
  ylim(c(0,100))+
  labs( x="", y = "Prevalence (%)") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, face = "bold"), 
        axis.title.y = element_text(size = 22, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"))
dev.off()


#compare symptoms score 
plot_data <- table_1_dat_fu

plot_data$disease_severity_upd <- factor(plot_data$disease_severity_upd, levels = c("mild_mod", "severe", "critical"))


pdf(file = "Figures/symptoms_at_fu_score.pdf", height = 10, width = 8)
ggplot(plot_data, aes(x = disease_severity_upd, y = Total_Symptoms_score_fu)) +
  geom_boxplot(width=0.6)+
  geom_jitter(alpha=0.5, width = 0.2)+
  stat_compare_means(comparisons = list(c("mild_mod", "severe"), 
                                        c("severe", "critical"), 
                                        c("mild_mod", "critical")))+
  ylab("Total Symptoms Score")+
  xlab("")+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe", "Critical"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
dev.off()



#compare symptoms score between two groups 

plot_data <- table_1_dat_fu

plot_data$disease_severity_two <- factor(plot_data$disease_severity_two, levels = c("mild_mod", "severe"))


pdf(file = "Figures/symptoms_at_fu_score_two_groups.pdf", height = 10, width = 8)
ggplot(plot_data, aes(x = disease_severity_two, y = Total_Symptoms_score_fu)) +
  geom_boxplot(width=0.6)+
  geom_jitter(alpha=0.5, width = 0.2)+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
  ylab("Total Symptoms Score")+
  xlab("")+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
dev.off()



#compare symptoms score according to inubation and ICU vs. not 

plot_data <- table_1_dat_fu

plot_data$ICU_y_n <- factor(plot_data$ICU_y_n)


pdf(file = "Figures/symptoms_at_fu_ICU_vs_not.pdf", height = 10, width = 8)
ggplot(plot_data, aes(x = ICU_y_n, y = Total_Symptoms_score_fu)) +
  geom_boxplot(width=0.6)+
  geom_jitter(alpha=0.5, width = 0.2)+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  ylab("Total Symptoms Score")+
  xlab("")+
  scale_x_discrete(labels=c("non-ICU", "ICU"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
dev.off()

plot_data$Intubated

pdf(file = "Figures/symptoms_at_fu_Intubated.pdf", height = 10, width = 8)
ggplot(plot_data, aes(x =Intubated , y = Total_Symptoms_score_fu)) +
  geom_boxplot(width=0.6)+
  geom_jitter(alpha=0.5, width = 0.2)+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  ylab("Total Symptoms Score")+
  xlab("")+
  scale_x_discrete(labels=c("non-ICU", "ICU"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"))
dev.off()




##### PFT at folloow up comparison 
plot_data <- table_1_dat_fu

plot_data$disease_severity_upd <- factor(plot_data$disease_severity_upd, levels = c("mild_mod", "severe", "critical"))

#compare PFTs according to disease severity 

PFT_var <- colnames(plot_data)
PFT_var <- PFT_var[120:130]

plot_data <- plot_data %>% mutate_at(., .vars = PFT_var, .funs = as.numeric)

summary_data <- plot_data %>% dplyr::select(PFT_var)

# Function to calculate median, Q1, and Q3 for a vector (handling NA values)
calculate_summary <- function(x) {
  median_val <- median(x, na.rm = TRUE)
  q1_val <- quantile(x, 0.25, na.rm = TRUE)
  q3_val <- quantile(x, 0.75, na.rm = TRUE)
  return(c(Median = median_val, Q1 = q1_val, Q3 = q3_val))
}

# Apply the function to each column of the dataframe
summary_stats <- sapply(summary_data, calculate_summary)

# Convert to dataframe and transpose to have variables as rows
summary_df <- data.frame(t(summary_stats))

# Export to CSV
write.csv(summary_df, "Results/summary_statistics_PFT_at_FU.csv", row.names = TRUE)



#go for a loop for PFT according to symptoms score and according to disease severity 
for (i in PFT_var) {
  p<- ggplot(plot_data, aes_string(x=plot_data$disease_severity_upd, y=i))+
    geom_boxplot(color="black", width=0.6)+
    geom_jitter(alpha=0.5, width = 0.2)+
    stat_compare_means(comparisons = list(c("mild_mod", "severe"), 
                                          c("severe", "critical"), 
                                          c("mild_mod", "critical")))+
    ylab("")+
    xlab("")+
    scale_x_discrete(labels=c("Mild-Moderate", "Severe", "Critical"))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"))
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_disease_severity_at_fu.pdf")))
  pdf(file=pdf_output,height = 10, width = 8)
  show(p)
  dev.off()
}

#repeat for total symtpoms score 
plot_data$symptom_score_categ <- ifelse(plot_data$Total_Symptoms_score_fu<3, "low", ifelse(plot_data$Total_Symptoms_score_fu>6, "high", "mod"))

plot_data$symptom_score_categ <- factor(plot_data$symptom_score_categ, levels = c("low", "mod", "high"))


summary_data <- plot_data %>% dplyr::select(PFT_var, symptom_score_categ)


# Calculate summary statistics for each variable by symptom_score_categ
summary_stats <- summary_data %>%
  group_by(symptom_score_categ) %>%
  summarize_all(funs(median = median(., na.rm = TRUE),
                     q1 = quantile(., 0.25, na.rm = TRUE),
                     q3 = quantile(., 0.75, na.rm = TRUE),
                     IQR = IQR(., na.rm = TRUE))) %>%
  ungroup()

# Save the results to a CSV file
write.csv(combined_results, "summary_and_comparison.csv", row.names = FALSE)


# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_PFT_by_category.csv", row.names = TRUE)


####### PFT comparison according to intubation and ICU

#intubation 
plot_data$Intubated <- factor(plot_data$Intubated)

for (i in PFT_var) {
  p<- ggplot(plot_data, aes_string(x=plot_data$Intubated, y=i))+
    geom_boxplot(color="black", width=0.6)+
    geom_jitter(alpha=0.5, width = 0.2)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    ylab("")+
    xlab("")+
    scale_x_discrete(labels=c("Non-Intubated", "Intubated"))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"))
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_intubation.pdf")))
  pdf(file=pdf_output,height = 10, width = 8)
  show(p)
  dev.off()
}


summary_data <- plot_data %>% dplyr::select(PFT_var, Intubated)


# Calculate summary statistics for each variable by symptom_score_categ
summary_stats <- summary_data %>%
  group_by(Intubated) %>%
  summarize_all(funs(median = median(., na.rm = TRUE),
                     q1 = quantile(., 0.25, na.rm = TRUE),
                     q3 = quantile(., 0.75, na.rm = TRUE),
                     IQR = IQR(., na.rm = TRUE))) %>%
  ungroup()

# Export to CSV
write.csv(summary_stats, "Results/summary_stats_PFT_by_intubated.csv", row.names = TRUE)


#intubation 
plot_data$ICU_y_n <- factor(plot_data$ICU_y_n)

for (i in PFT_var) {
  p<- ggplot(plot_data, aes_string(x=plot_data$ICU_y_n, y=i))+
    geom_boxplot(color="black", width=0.6)+
    geom_jitter(alpha=0.5, width = 0.2)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    ylab("")+
    xlab("")+
    scale_x_discrete(labels=c("Non-ICU", "ICU"))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"))
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_ICU_y_n.pdf")))
  pdf(file=pdf_output,height = 10, width = 8)
  show(p)
  dev.off()
}


summary_data <- plot_data %>% dplyr::select(PFT_var, ICU_y_n)


# Calculate summary statistics for each variable by symptom_score_categ
summary_stats <- summary_data %>%
  group_by(ICU_y_n) %>%
  summarize_all(funs(median = median(., na.rm = TRUE),
                     q1 = quantile(., 0.25, na.rm = TRUE),
                     q3 = quantile(., 0.75, na.rm = TRUE),
                     IQR = IQR(., na.rm = TRUE))) %>%
  ungroup()

# Export to CSV
write.csv(summary_stats, "Results/summary_stats_PFT_by_ICU_y_n.csv", row.names = TRUE)




########## general analysis looking at molecules in all patients while comparing baseline to follow up to control 

col_names <- colnames(mycounts_t)[2:10]

mycounts_t$sample_type <- factor(mycounts_t$sample_type)

#sstats 
summary_stats <- mycounts_t %>%
  dplyr::select(-c(sample)) %>% 
  group_by(sample_type) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_by_SAMPLE_TYPE.csv", row.names = TRUE)


for (i in col_names) {
  p <- ggplot(data = mycounts_t, aes_string(x=mycounts_t$sample_type, y=i, color=mycounts_t$sample_type))+
    geom_boxplot(outlier.shape = NA)+
    scale_color_manual(values=c("red", "orange", "blue"))+
    scale_x_discrete(labels=c("Baseline", "Follow Up", "Control"))+
    geom_jitter(alpha=0.5, width = 0.2)+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("baseline", "clinic"), 
                                          c("clinic", "control"), 
                                          c("baseline", "control")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_sample_type.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}


#specefically repeat IL-6 while limitng to 500 
temp_dat <- mycounts_t 
temp_dat <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>500, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat, aes(x=sample_type, y=IL6, color=sample_type))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values=c("red", "orange", "blue"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up", "Control"))+
  geom_jitter(alpha=0.5, width = 0.2)+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic"), 
                                        c("clinic", "control"), 
                                        c("baseline", "control")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_sample_type.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat <- mycounts_t %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat, aes(x=sample_type, y=CA15_3, color=sample_type))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values=c("red", "orange", "blue"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up", "Control"))+
  geom_jitter(alpha=0.5, width = 0.2)+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic"), 
                                        c("clinic", "control"), 
                                        c("baseline", "control")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_no_outliers_", paste0("_by_sample_type.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()







#patired analysis between baseline and foloow up 

col_names <- colnames(mycounts_t)[2:10]

#put data in TEMP dataframe 
temp_dat <- mycounts_t
#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
temp_dat$sample_type <- factor(temp_dat$sample_type)

#leave only those with paired samples 
temp_dat<- temp_dat %>% group_by(subject_ID) %>% filter(n()>1) %>% ungroup()
temp_dat$subject_ID <- factor(temp_dat$subject_ID)

ggplot(data = temp_dat, aes(x=sample_type, y=MMP9, color=sample_type))+
  geom_boxplot()+
  geom_point()+
  scale_color_manual(values=c("red", "orange"))+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  geom_line(aes(group = subject_ID), color = "grey")
  



#sstats 
summary_stats <- temp_dat %>%
  dplyr::select(-c(sample, cluster, disease_severity, subject_ID)) %>% 
  group_by(sample_type) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_by_SAMPLE_TYPE_baseline_clinic.csv", row.names = TRUE)



for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$sample_type, y=i, color=temp_dat$sample_type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(alpha=0.5, width = 0.2)+
    geom_line(aes(group = subject_ID), color = "grey")+
    scale_color_manual(values=c("red", "orange"))+
    scale_x_discrete(labels=c("Baseline", "Follow Up"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("baseline", "clinic")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_sample_type_paired_baseline_follow_up.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}


#repeat for IL-6 and KL-6 while removing outliers 



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=IL6, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_sample_type_paired.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=CA15_3, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_no_outliers_", paste0("_by_sample_type_paired.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()



###########compare baseline samples according to disease severity ###### 
temp_dat <- mycounts_t

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$disease_severity_two <- metadata$disease_severity_two

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$disease_severity_two, y=IL6, color=temp_dat$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("KL-6")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

temp_dat$disease_severity_two <- factor(temp_dat$disease_severity_two)


col_names <- colnames(mycounts_t)[2:10]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$disease_severity_two, y=i, color=as.factor(temp_dat$disease_severity_two)))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(alpha=0.5, width = 0.2)+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
      theme_prism()+
      theme(text = element_text(size = 16), 
            axis.title.x = element_text(size = 24, face = "bold"), 
            axis.text.x = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 16, face = "bold"), 
            axis.title.y = element_text(size = 24, face = "bold"))+
      guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_disease_severity.pdf")))
  pdf(file=pdf_output,height = 6,  width = 5)
  show(p)
  dev.off()
}



#repeat for IL-6 and KL-6 while removing outliers 



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$disease_severity_two, y=IL6, color=temp_dat_2$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_disease_severity_baseline.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$disease_severity_two, y=CA15_3, color=temp_dat_2$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_by_disease_severity_baseline.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()




#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID)) %>% 
  group_by(disease_severity_two) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_baseline_by_disease_severity.csv", row.names = TRUE)



#repeat comparison only for follow up time point
temp_dat <- mycounts_t

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$disease_severity_two <- metadata$disease_severity_two

temp_dat <- temp_dat %>% filter(sample_type!="baseline")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$disease_severity_two, y=IL6, color=temp_dat$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("KL-6")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

temp_dat$disease_severity_two <- factor(temp_dat$disease_severity_two)


col_names <- colnames(mycounts_t)[2:10]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$disease_severity_two, y=i, color=as.factor(temp_dat$disease_severity_two)))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(alpha=0.5, width = 0.2)+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_follow_up_levels_by_disease_severity.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}


#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$disease_severity_two, y=IL6, color=temp_dat_2$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_disease_severity_follow_up.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$disease_severity_two, y=CA15_3, color=temp_dat_2$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_by_disease_severity__follow_up.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID)) %>% 
  group_by(disease_severity_two) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_follow_up_by_disease_severity.csv", row.names = TRUE)






###### compare at baseline between ICU vs. non ICU #######
temp_dat <- mycounts_t

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$ICU_y_n <- metadata$ICU_y_n

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$ICU_y_n, y=IL6, color=temp_dat$ICU_y_n))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("non-ICU", "ICU"))+
  xlab("")+ylab("IL-6")+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

temp_dat$ICU_y_n <- factor(temp_dat$ICU_y_n)


col_names <- colnames(mycounts_t)[2:10]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$ICU_y_n, y=i, color=as.factor(temp_dat$ICU_y_n)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Non ICU", "ICU"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_ICU.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$ICU_y_n, y=IL6, color=temp_dat_2$ICU_y_n))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Non ICU", "ICU"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_ICU_baseline.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$ICU_y_n, y=CA15_3, color=temp_dat_2$ICU_y_n))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Non ICU", "ICU"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_by_ICU_baseline.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()






#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID, disease_severity)) %>% 
  group_by(ICU_y_n) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_baseline_by_ICU_Y_N.csv", row.names = TRUE)



### repeat for follow up 
temp_dat <- mycounts_t

#filter for follow_up and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$ICU_y_n <- metadata$ICU_y_n

temp_dat <- temp_dat %>% filter(sample_type!="baseline")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$ICU_y_n, y=IL6, color=temp_dat$ICU_y_n))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("non-ICU", "ICU"))+
  xlab("")+ylab("IL-6")+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

temp_dat$ICU_y_n <- factor(temp_dat$ICU_y_n)


col_names <- colnames(mycounts_t)[2:10]


#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$ICU_y_n, y=i, color=as.factor(temp_dat$ICU_y_n)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Non ICU", "ICU"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_follow_up_levels_by_ICU.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$ICU_y_n, y=IL6, color=temp_dat_2$ICU_y_n))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Non ICU", "ICU"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_ICU_follow_up.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$ICU_y_n, y=CA15_3, color=temp_dat_2$ICU_y_n))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Non ICU", "ICU"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_by_ICU_follow_up.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()



#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID, disease_severity)) %>% 
  group_by(ICU_y_n) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_follow_up_by_ICU_Y_N.csv", row.names = TRUE)



######## baseline and follow up levels comparison according to intubation ######

temp_dat <- mycounts_t

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$Intubated <- metadata$Intubated

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$Intubated, y=IL6, color=temp_dat$Intubated))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("No Intubation", "Intubation"))+
  xlab("")+ylab("IL-6")+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

temp_dat$Intubated <- factor(temp_dat$Intubated)


col_names <- colnames(mycounts_t)[2:10]

for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$Intubated, y=i, color=as.factor(temp_dat$Intubated)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("No Intubation", "Intubation"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_intubation.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$Intubated, y=IL6, color=temp_dat_2$Intubated))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("No Intubation", "Intubation"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_baseline_levels_by_intubation.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$Intubated, y=CA15_3, color=temp_dat_2$Intubated))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("No Intubation", "Intubation"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_baseline_levels_by_intubation.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID, disease_severity)) %>% 
  group_by(Intubated) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_baseline_by_Intubated.csv", row.names = TRUE)



### repeat for follow up 
temp_dat <- mycounts_t

#filter for follow_up and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$Intubated <- metadata$Intubated

temp_dat <- temp_dat %>% filter(sample_type!="baseline")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$Intubated, y=IL6, color=temp_dat$Intubated))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("non-ICU", "ICU"))+
  xlab("")+ylab("IL-6")+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

temp_dat$Intubated <- factor(temp_dat$Intubated)


col_names <- colnames(mycounts_t)[2:10]


for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$Intubated, y=i, color=as.factor(temp_dat$Intubated)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("No Intubation", "Intubation"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_follow_up_levels_by_intubation.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$Intubated, y=IL6, color=temp_dat_2$Intubated))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("No Intubation", "Intubation"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_follow_up_levels_by_intubation.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$Intubated, y=CA15_3, color=temp_dat_2$Intubated))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("No Intubation", "Intubation"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_follow_up_levels_by_intubation.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID, disease_severity)) %>% 
  group_by(Intubated) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_follow_up_by_Intubated.csv", row.names = TRUE)




#### levels of baseline only comparing dead to alive 


temp_dat <- mycounts_t

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
#temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$Death <- metadata$Death

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$Death, y=IL6, color=temp_dat$Death))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Alive", "Dead"))+
  xlab("")+ylab("IL-6")+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

temp_dat$Death <- factor(temp_dat$Death)


col_names <- colnames(mycounts_t)[2:10]



for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$Death, y=i, color=as.factor(temp_dat$Death)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Alive", "Dead"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_follow_up_levels_by_Death.pdf")))
  pdf(file=pdf_output,height = 6, width = 5)
  show(p)
  dev.off()
}



#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>300, 500,temp_dat$IL6))


p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$Death, y=IL6, color=temp_dat_2$Death))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Alive", "Dead"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_follow_up_levels_by_Death.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% filter(CA15_3 < 1000)

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$Death, y=CA15_3, color=temp_dat_2$Death))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Alive", "Dead"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("N", "Y")), size=6)+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("KL_6_no_outliers_", paste0("_follow_up_levels_by_Death.pdf")))
pdf(file=pdf_output,height = 6, width = 5)
p
dev.off()


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(-c(sample_type, sample, subject_ID, disease_severity)) %>% 
  group_by(Death) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_at_baseline_by_death.csv", row.names = TRUE)











##### compare baseline and follow up samples between 2 categories of disease severity ##### 

temp_dat <- mycounts_t
#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
temp_dat$sample_type <- factor(temp_dat$sample_type)

#get disease severity of two categories 
temp_dat$disease_severity_two <- metadata$disease_severity_two

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$sample_type, y=IL6, color=temp_dat$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("KL-6")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~disease_severity_two)+
  theme_prism()+
  guides(color="none")


#sstats 
summary_stats <- temp_dat %>%
  dplyr::select(-c(sample, cluster, disease_severity, sample_type, subject_ID)) %>% 
  group_by(disease_severity_two) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_by_disease_severity_baseline_clinic_paired.csv", row.names = TRUE)


col_names <- colnames(mycounts_t)[2:10]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$sample_type, y=i, color=temp_dat$sample_type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    geom_line(aes(group = subject_ID), color = "grey")+
    scale_color_manual(values=c("red", "orange"))+
    scale_x_discrete(labels=c("Baseline", "Follow Up"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("baseline", "clinic")))+
    facet_wrap(~disease_severity_two, labeller = labeller(severe = "Severe", mild_mod = "Mild"), scales = "free")+
    theme_prism()+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_sample_type_paired_baseline_follow_up_by_disease_severity.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}

#again repeat for IL6 and KL6 while removing outliers 


#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>500, 500, temp_dat$IL6))

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=IL6, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("IL6")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~disease_severity_two, labeller = labeller(severe = "Severe", mild_mod = "Mild"), scales = "free")+
  theme_prism()+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_sample_type_paired_by_disease_Severity.pdf")))
pdf(file=pdf_output,height = 8, width = 6)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% mutate(CA15_3=ifelse(temp_dat$CA15_3>1000, 1000, temp_dat$CA15_3))

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=CA15_3, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("KL-6")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~disease_severity_two, labeller = labeller(severe = "Severe", mild_mod = "Mild"), scales = "free")+
  theme_prism()+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_no_outliers_", paste0("_by_sample_type_paired_by_disease_seveity.pdf")))
pdf(file=pdf_output,height = 8, width = 6)
p
dev.off()






####### comparison between intubated vs non-intubated patients ###### 

temp_dat$Intubated <- metadata$Intubated

#repeat the same loop again 



col_names <- colnames(mycounts_t)[2:10]


#sstats 
summary_stats <- temp_dat %>%
  dplyr::select(-c(sample, cluster, disease_severity, subject_ID, sample_type,disease_severity_two)) %>% 
  group_by(Intubated) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_by_intubation.csv", row.names = TRUE)


#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$sample_type, y=i, color=temp_dat$sample_type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    geom_line(aes(group = subject_ID), color = "grey")+
    scale_color_manual(values=c("red", "orange"))+
    scale_x_discrete(labels=c("Baseline", "Follow Up"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("baseline", "clinic")))+
    facet_wrap(~Intubated, scales = "free")+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_sample_type_paired_baseline_follow_up_by_intubation.pdf")))
  pdf(file=pdf_output,height = 8, width = 7.5)
  show(p)
  dev.off()
}

#again repeat for IL6 and KL6 while removing outliers 


#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>500, 500, temp_dat$IL6))

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=IL6, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~Intubated, scales = "free")+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_sample_type_paired_by_intubation.pdf")))
pdf(file=pdf_output,height = 8, width = 7.5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% mutate(CA15_3=ifelse(temp_dat$CA15_3>1000, 1000, temp_dat$CA15_3))

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=CA15_3, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("KL6 units/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~Intubated, scales = "free")+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_no_outliers_", paste0("_by_sample_type_paired_by_intubation.pdf")))
pdf(file=pdf_output,height = 8, width = 7.5)
p
dev.off()













#######comparison of baseline vs. follow up among ICU vs. non ICU ###############

temp_dat <- mycounts_t
#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type!="control")
#get subject IDs
temp_dat$subject_ID <- metadata$subject_ID
temp_dat$sample_type <- factor(temp_dat$sample_type)
temp_dat$ICU_y_n <- metadata$ICU_y_n

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$sample_type, y=IL6, color=temp_dat$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("IL-6")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~ICU_y_n)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")


#sstats 
summary_stats <- temp_dat %>%
  dplyr::select(-c(sample, sample_type, subject_ID)) %>% 
  group_by(ICU_y_n) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

# Export to CSV
write.csv(summary_stats, "Results/summary_statistics_cytokines_by_ICU_y_n_baseline_clinic_paired.csv", row.names = TRUE)


col_names <- colnames(mycounts_t)[2:10]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$sample_type, y=i, color=temp_dat$sample_type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    geom_line(aes(group = subject_ID), color = "grey")+
    scale_color_manual(values=c("red", "orange"))+
    scale_x_discrete(labels=c("Baseline", "Follow Up"))+
    xlab("")+ylab(paste0(i, paste0("  pg/mL")))+
    stat_compare_means(comparisons = list(c("baseline", "clinic")))+
    facet_wrap(~ICU_y_n, scales = "free")+
    theme_prism()+
    theme(text = element_text(size = 16), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_sample_type_paired_baseline_follow_up_by_ICU.pdf")))
  pdf(file=pdf_output,height = 8, width = 7.5)
  show(p)
  dev.off()
}




#specefically repeat IL-6 while limitng to 100 
temp_dat_2 <- temp_dat %>% mutate(IL6=ifelse(temp_dat$IL6>500, 500, temp_dat$IL6))

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=IL6, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("IL6  pg/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~ICU_y_n, scales = "free")+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("IL_6_no_outliers_", paste0("_by_sample_type_paired_by_ICU_y_n.pdf")))
pdf(file=pdf_output,height = 8, width = 7.5)
p
dev.off()

#same for CA15-3
temp_dat_2 <- temp_dat %>% mutate(CA15_3=ifelse(temp_dat$CA15_3>1000, 1000, temp_dat$CA15_3))

p <- ggplot(data = temp_dat_2, aes(x=temp_dat_2$sample_type, y=CA15_3, color=temp_dat_2$sample_type))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_line(aes(group = subject_ID), color = "grey")+
  scale_color_manual(values=c("red", "orange"))+
  scale_x_discrete(labels=c("Baseline", "Follow Up"))+
  xlab("")+ylab("KL6  units/mL")+
  stat_compare_means(comparisons = list(c("baseline", "clinic")))+
  facet_wrap(~ICU_y_n, scales = "free")+
  theme_prism()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_no_outliers_", paste0("_by_sample_type_paired_by_ICU_y_n.pdf")))
pdf(file=pdf_output,height = 8, width = 7.5)
p
dev.off()










####### K medoids clustering does not work #########
library(factoextra)
library(cluster)

temp_dat <- mycounts_t

temp_dat <- temp_dat %>% select(-c(sample_type, sample))

temp_dat <- na.omit(temp_dat)

temp_dat<- scale(temp_dat)

fviz_nbclust(temp_dat, pam, method = "wss")

gap_stat <- clusGap(temp_dat,
                    FUN = pam,
                    K.max = 10, #max clusters to consider
                    B = 50) #total bootstrapped iterations

#plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)



######## herarchial clustering - does not work in this data #####
library(factoextra)
library(cluster)

#load data
df <- temp_dat

#remove rows with missing values
df <- na.omit(df)

#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)

#perform hierarchical clustering using average
clust <- agnes(df, method = "ward")

#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 

#calculate gap statistic for each number of clusters (up to 10 clusters)
gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)

#produce plot of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

#compute distance matrix
d <- dist(df, method = "euclidean")

#perform hierarchical clustering using Ward's method
final_clust <- hclust(d, method = "ward.D2" )

#cut the dendrogram into 4 clusters
groups <- cutree(final_clust, k=4)

# Number of members in each cluster
table(groups)

#append cluster labels to original data
final_data <- cbind(temp_dat, cluster = groups)

#display first six rows of final data
head(final_data)

#find mean values for each cluster
aggregate(final_data, by=list(cluster=final_data$cluster), mean)


######### PCA of all samples ######

#first perform clustering or PCA 

#keep from mycounts_t only samples without control 

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  select(-sample_type, -sample) %>%
  scale(.) %>% 
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = FALSE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
pc_scores$sample <- mycounts_t$sample
pc_scores$sample_type <- mycounts_t$sample_type

# Plot PCA
ggplot(pc_scores, aes(x = PC1, y = PC2, color = sample_type, shape = sample_type)) +
  geom_point() +
  labs(title = "PCA Plot by sample type") +
  theme_minimal()


#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ sample_type,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="sample_type",suffixes=c("",".centroid"))


pdf(file = "Figures/all_samples_PCA.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= sample_type)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("red", "orange", "blue")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= sample_type), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= sample_type)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Baseline", "Follow Up", "Control")), size=10) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ sample_type, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ sample_type, data = pc_scores))

temp_dat <- pc_scores %>% filter(sample_type!= "control")
summary(manova(cbind(PC1, PC2) ~ sample_type, data = temp_dat))


#### PCA plot according to disease severity at baseline only####

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="baseline") %>% 
  select(-sample_type, -sample) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="baseline") %>% select(disease_severity_upd)
temp <- temp$disease_severity_upd
pc_scores$disease_severity_upd <- temp
pc_scores$disease_severity_upd <- factor(pc_scores$disease_severity_upd, levels = c("mild_mod", "severe", "critical"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ disease_severity_upd,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="disease_severity_upd",suffixes=c("",".centroid"))


pdf(file = "Figures/baseline_samples_PCA_by_disease_sev.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= disease_severity_upd)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("grey", "goldenrod", "orangered")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= disease_severity_upd), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= disease_severity_upd)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Mild Moderate", "Severe", "Critical")), size=10) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)) # p= 0.67




#####now repeat for mild-mod vs. severe and critical##### 

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="baseline") %>% 
  select(-sample_type, -sample) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="baseline") %>% select(disease_severity_upd)
temp$disease_severity_upd <- ifelse(temp$disease_severity_upd=="mild_mod", "mild_mod", "severe_critical")
temp <- temp$disease_severity_upd
pc_scores$disease_severity_upd <- temp
pc_scores$disease_severity_upd <- factor(pc_scores$disease_severity_upd, levels = c("mild_mod", "severe_critical"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ disease_severity_upd,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="disease_severity_upd",suffixes=c("",".centroid"))

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)
x <- summary(manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)) # p= 0.33


pdf(file = "Figures/baseline_samples_PCA_by_disease_sev_mild_mod_vs_sevecritic.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= disease_severity_upd)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("grey", "goldenrod")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= disease_severity_upd), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= disease_severity_upd)) + 
  #labels centroids should be same number of categories
  ggtitle(paste0("p=", paste0(x[4])))+
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Mild Moderate", "Severe Critical")), size=10) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()




###### cluster baseline data using k means silhoete method ##### 
baseline_mycounts <- mycounts_t %>% filter(sample_type=="baseline") %>% select(-c(sample, sample_type))
scaled_baseline_dat <- scale(baseline_mycounts)

library(factoextra)
library(fpc)

# Run k-means clustering for different values of k
k_values <- 2:10
sil_widths <- numeric(length(k_values))
#prepeare data, should be scaled before that 
temp_dat <- scaled_baseline_dat

for (k in k_values) {
  kmeans_results <- kmeans(temp_dat, centers = k, nstart = 10)
  sil_widths[k - 1] <- cluster.stats(dist(temp_dat), kmeans_results$cluster)$avg.silwidth
}

# Plot silhouette widths for different k values
plot(k_values, sil_widths, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters (K)", ylab = "Average Silhouette Width")

# Identify the optimal number of clusters with the maximum silhouette width
optimal_k <- k_values[which.max(sil_widths)]
cat("Optimal number of clusters (k) using silhouette analysis:", optimal_k, "\n")

##### here you get 2 clusters as optimal 

###now assign clusters 

set.seed(123)  # for reproducibility
kmeans_result <- kmeans(scaled_baseline_dat, centers = 2)

baseline_mycounts$cluster <- kmeans_result$cluster

# plot 
library(factoextra)
p <- fviz_cluster(kmeans_result, data = scaled_baseline_dat, stand = FALSE)
pdf(file = "Figures/cluster_plot_kmeans_baseline_samples.pdf", height = 8, width = 12)
p+
  scale_color_manual(values = c("goldenrod", "darkblue"))+
  theme_prism()
dev.off()

baseline_metadata <- metadata %>% filter(sample_type=="baseline")

baseline_metadata$cluster <- kmeans_result$cluster
baseline_metadata$cluster <- factor(baseline_metadata$cluster)
baseline_mycounts$cluster <- factor(baseline_mycounts$cluster)

#compare cytokiens levels across clustering gorups 

ggplot(baseline_mycounts, aes(x=as.factor(cluster), y=MMP7))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

for (i in colnames(baseline_mycounts)[1:9]){
  p <- ggplot(data = baseline_mycounts, aes_string(x=baseline_mycounts$cluster, y=i, color=baseline_mycounts$cluster))+
    geom_boxplot()+
    geom_jitter(color="black", alpha=0.6, width = 0.2)+
    scale_color_manual(values = c("goldenrod", "darkblue"))+
    stat_compare_means(comparisons = list(c("1", "2")))+
    scale_x_discrete(labels=c(",Cluster 1", "Cluster 2"))+
    xlab("")+ylab(i)+
    guides(color="none")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_cluster_kmeans_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#now write data with cluster assignment 
write.csv(baseline_metadata, "Results/cluster_assignment_baseline.csv")


##### compare clinical data ebtween clusters 
table_1_dat <- baseline_metadata

colnames(table_1_dat)
#create table 1 using table 1 package 
table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         as.numeric(Age)+as.factor(Gender)+as.factor(Origin)+as.factor(smoking)+as.numeric(PY)+
         as.numeric(time_to_fu)+as.numeric(duration_of_admisson)+ as.factor(Death)+as.factor(HTN)+as.factor(DM)+as.factor(HLP)+as.factor(IHD)+as.factor(CHF)+
         as.factor(CKD)+as.factor(transplant)+as.factor(cancer_y_n)+as.factor(Liver.Disease)+as.factor(COPD)+as.factor(Asthma)+as.factor(other_lung_disease)+
         as.factor(Vaccination.status)+as.factor(disease_severity)+as.factor(Intubated)+
         as.numeric(sao2_ra)+as.numeric(Crea_baseline)+as.numeric(Crea_peak)+as.numeric(Crea_discharge)+as.numeric(LDH_peak)+
         as.numeric(WBC_peak)+as.numeric(DDM_peak)+as.numeric(ferritin_peak)+as.numeric(crp_peak)+
         as.factor(Fatigue_fu)+as.factor(Fever_fu)+as.factor(Headache_fu)+as.factor(Dyspnea_fu)+as.factor(Cough_fu)+
         as.factor(Chest_pain_fu)+as.factor(Nausea_fu)+as.factor(Diarrhea_fu)+as.factor(Sore_throat_fu)+as.factor(Nasal_congestion_fu)+as.factor(Red_eyes_fu)+
         as.factor(Vision_loss_fu)+as.factor(Concentraion_difficulties_fu)+as.factor(Memory_decline_fu)+as.factor(Weight_loss_fu)+as.factor(Arthalgia_fu)+
         as.factor(Parasthesia_Persistent_pain_fu)+ as.factor(Hair_Loss_fu)+as.factor(Aguesia_Dysguesia_fu)+
         as.numeric(Total_Symptoms_score_fu)++as.numeric(sao2_ra_fu)+
         as.numeric(FEV1_lit)+as.numeric(FEV1_p)+as.numeric(FVC_lit)+as.numeric(FVC_p)+as.numeric(FEV1_FVC)+as.numeric(TLC_lit)+
         as.numeric(TLC_p)+as.numeric(RV_lit)+as.numeric(RV_p)+as.numeric(RV_TLC)+as.numeric(KCO)+as.numeric(TLC_lit)
       
       
       #set data 
       |cluster ,data=table_1_dat, 
       #set stat display options for continous variables (options from stat.default)
       overall = FALSE,
       render.continuous = c((.="Median [Q1, Q3]")),
       #add p value as extra column (defined in function above)
       extra.col = list("P-value" = pvalue))







#######clustering all smaples using k means #####


temp_dat <- mycounts_t %>% select(-c(sample_type, sample)) %>% scale(.)

library(factoextra)
library(fpc)

# Run k-means clustering for different values of k
k_values <- 2:10
sil_widths <- numeric(length(k_values))
#prepeare data, should be scaled before that 

for (k in k_values) {
  kmeans_results <- kmeans(temp_dat, centers = k, nstart = 10)
  sil_widths[k - 1] <- cluster.stats(dist(temp_dat), kmeans_results$cluster)$avg.silwidth
}

# Plot silhouette widths for different k values
plot(k_values, sil_widths, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters (K)", ylab = "Average Silhouette Width")

# Identify the optimal number of clusters with the maximum silhouette width
optimal_k <- k_values[which.max(sil_widths)]
cat("Optimal number of clusters (k) using silhouette analysis:", optimal_k, "\n")


#### here we also get 2 clusters 
###now assign clusters 

set.seed(123)  # for reproducibility
kmeans_result <- kmeans(temp_dat, centers = 2)

mycounts_t$cluster <- kmeans_result$cluster


# plot 
library(factoextra)
p <- fviz_cluster(kmeans_result, data = temp_dat, stand = FALSE)
pdf(file = "Figures/cluster_plot_kmeans_all_samples.pdf", height = 8, width = 12)
p+
  scale_color_manual(values = c("goldenrod", "darkblue"))+
  theme_prism()
dev.off()


metadata$cluster <- kmeans_result$cluster[1:91]
metadata$cluster <- factor(metadata$cluster)
mycounts_t$cluster <- factor(mycounts_t$cluster)

#compare cytokiens levels across clustering gorups 

ggplot(mycounts_t, aes(x=as.factor(cluster), y=MMP7))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

for (i in colnames(mycounts_t)[2:10]){
  p <- ggplot(data = mycounts_t, aes_string(x=mycounts_t$cluster, y=i, color=mycounts_t$cluster))+
    geom_boxplot()+
    geom_jitter(color="black", alpha=0.6, width = 0.2)+
    scale_color_manual(values = c("goldenrod", "darkblue"))+
    stat_compare_means(comparisons = list(c("1", "2")))+
    scale_x_discrete(labels=c(",Cluster 1", "Cluster 2"))+
    xlab("")+ylab(i)+
    guides(color="none")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_cluster_kmeans_all_Samples.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


for (i in colnames(mycounts_t)[2:10]){
  p <- ggplot(data = mycounts_t, aes_string(x=mycounts_t$cluster, y=i, color=mycounts_t$cluster))+
    geom_boxplot()+
    geom_jitter(color="black", alpha=0.6, width = 0.2)+
    scale_color_manual(values = c("goldenrod", "darkblue"))+
    stat_compare_means(comparisons = list(c("1", "2")))+
    scale_x_discrete(labels=c(",Cluster 1", "Cluster 2"))+
    xlab("")+ylab(i)+
    facet_wrap(~sample_type, scales = "free")+
    guides(color="none")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_cluster_kmeans_all_Samples_facet_by_sample_type.pdf")))
  pdf(file=pdf_output,height = 8, width = 10)
  show(p)
  dev.off()
}




#now write data with cluster assignment 
write.csv(metadata, "Results/cluster_assignment_all_samples.csv")


##### compare clinical data ebtween clusters 
table_1_dat <- metadata

colnames(table_1_dat)
#create table 1 using table 1 package 
table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         as.numeric(Age)+as.factor(Gender)+as.factor(Origin)+as.factor(smoking)+as.numeric(PY)+
         as.numeric(time_to_fu)+as.numeric(duration_of_admisson)+ as.factor(Death)+as.factor(HTN)+as.factor(DM)+as.factor(HLP)+as.factor(IHD)+as.factor(CHF)+
         as.factor(CKD)+as.factor(transplant)+as.factor(cancer_y_n)+as.factor(Liver.Disease)+as.factor(COPD)+as.factor(Asthma)+as.factor(other_lung_disease)+
         as.factor(Vaccination.status)+as.factor(disease_severity)+as.factor(Intubated)+
         as.numeric(sao2_ra)+as.numeric(Crea_baseline)+as.numeric(Crea_peak)+as.numeric(Crea_discharge)+as.numeric(LDH_peak)+
         as.numeric(WBC_peak)+as.numeric(DDM_peak)+as.numeric(ferritin_peak)+as.numeric(crp_peak)+
         as.factor(Fatigue_fu)+as.factor(Fever_fu)+as.factor(Headache_fu)+as.factor(Dyspnea_fu)+as.factor(Cough_fu)+
         as.factor(Chest_pain_fu)+as.factor(Nausea_fu)+as.factor(Diarrhea_fu)+as.factor(Sore_throat_fu)+as.factor(Nasal_congestion_fu)+as.factor(Red_eyes_fu)+
         as.factor(Vision_loss_fu)+as.factor(Concentraion_difficulties_fu)+as.factor(Memory_decline_fu)+as.factor(Weight_loss_fu)+as.factor(Arthalgia_fu)+
         as.factor(Parasthesia_Persistent_pain_fu)+ as.factor(Hair_Loss_fu)+as.factor(Aguesia_Dysguesia_fu)+
         as.numeric(Total_Symptoms_score_fu)++as.numeric(sao2_ra_fu)+
         as.numeric(FEV1_lit)+as.numeric(FEV1_p)+as.numeric(FVC_lit)+as.numeric(FVC_p)+as.numeric(FEV1_FVC)+as.numeric(TLC_lit)+
         as.numeric(TLC_p)+as.numeric(RV_lit)+as.numeric(RV_p)+as.numeric(RV_TLC)+as.numeric(KCO)+as.numeric(TLC_lit)
       
       
       #set data 
       |cluster ,data=table_1_dat, 
       #set stat display options for continous variables (options from stat.default)
       overall = FALSE,
       render.continuous = c((.="Median [Q1, Q3]")),
       #add p value as extra column (defined in function above)
       extra.col = list("P-value" = pvalue))


#stats of numerical variables 

chisq.test(table_1_dat$Age, table_1_dat$cluster)
chisq.test(table_1_dat$PY, table_1_dat$cluster)
chisq.test(table_1_dat$time_to_fu, table_1_dat$cluster)
chisq.test(table_1_dat$duration_of_admisson, table_1_dat$cluster)
chisq.test(table_1_dat$sao2_ra, table_1_dat$cluster)
chisq.test(table_1_dat$Crea_baseline, table_1_dat$cluster)
chisq.test(table_1_dat$Crea_peak, table_1_dat$cluster)
chisq.test(table_1_dat$Crea_discharge, table_1_dat$cluster)
chisq.test(table_1_dat$LDH_peak, table_1_dat$cluster)
chisq.test(table_1_dat$WBC_peak, table_1_dat$cluster)
chisq.test(table_1_dat$DDM_peak, table_1_dat$cluster)
chisq.test(table_1_dat$ferritin_peak, table_1_dat$cluster)
chisq.test(table_1_dat$crp_peak, table_1_dat$cluster)
chisq.test(table_1_dat$FEV1_lit, table_1_dat$cluster)
chisq.test(table_1_dat$FEV1_p, table_1_dat$cluster)
chisq.test(table_1_dat$FVC_lit, table_1_dat$cluster)
chisq.test(table_1_dat$FVC_p, table_1_dat$cluster)
chisq.test(table_1_dat$FEV1_FVC, table_1_dat$cluster)
chisq.test(table_1_dat$TLC_lit, table_1_dat$cluster)
chisq.test(table_1_dat$TLC_p, table_1_dat$cluster)
chisq.test(table_1_dat$RV_lit, table_1_dat$cluster)
chisq.test(table_1_dat$RV_p, table_1_dat$cluster)
chisq.test(table_1_dat$RV_TLC, table_1_dat$cluster)
chisq.test(table_1_dat$KCO, table_1_dat$cluster)








######## corelation plots between  counts and PFT at follow up#######
#merge data 
baseline_metadata
baseline_mycounts
merged_basseline <- merge(baseline_mycounts, baseline_metadata, by="row.names")
merged_basseline$FVC_lit <- as.numeric(merged_basseline$FVC_lit)
merged_basseline$FVC_p <- as.numeric(merged_basseline$FVC_p)
merged_basseline$TLC_lit <- as.numeric(merged_basseline$TLC_lit)
merged_basseline$TLC_p <- as.numeric(merged_basseline$TLC_p)
merged_basseline$KCO <- as.numeric(merged_basseline$KCO)

col_names <- colnames(baseline_mycounts)[1:9]

# FVC lit 
for (i in col_names) {
  p <- ggplot(data = merged_basseline, aes_string(x=i, y=merged_basseline$FVC_lit))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_FVC_lit_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#repeat for IL6 while removing outliers 
plot_data <- merged_basseline %>% filter(IL6<100)

p <- ggplot(data = plot_data, aes(x=IL6, y=FVC_lit))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("IL6_correl_FVC_lit_baseline_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()

#repeat for KL6 while removing outliers 
plot_data <- merged_basseline %>% filter(CA15_3<120)

p <- ggplot(data = plot_data, aes(x=CA15_3, y=FVC_lit))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_correl_FVC_lit_baseline_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()




# DLCO
for (i in col_names) {
  p <- ggplot(data = merged_basseline, aes_string(x=i, y=merged_basseline$KCO))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_KCO_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#same plots repetition for IL6 and KL6 
plot_data <- merged_basseline %>% filter(IL6<100)

p <- ggplot(data = plot_data, aes(x=IL6, y=KCO))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("IL6_correl_KCO_baseline_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()

#repeat for KL6 while removing outliers 
plot_data <- merged_basseline %>% filter(CA15_3<120)

p <- ggplot(data = plot_data, aes(x=CA15_3, y=KCO))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_correl_KCO_baseline_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()



#######correlaiton of PFT with follow up samples ######## 
fu_metadata <- metadata %>% filter(sample_type=="clinic")
fu_counts <- mycounts_t %>% filter(sample_type=="clinic")
merged_fu <- merge(fu_counts, fu_metadata, by="row.names")
merged_fu$FVC_lit <- as.numeric(merged_fu$FVC_lit)
merged_fu$FVC_p <- as.numeric(merged_fu$FVC_p)
merged_fu$TLC_lit <- as.numeric(merged_fu$TLC_lit)
merged_fu$TLC_p <- as.numeric(merged_fu$TLC_p)
merged_fu$KCO <- as.numeric(merged_fu$KCO)

col_names <- colnames(fu_counts)[2:10]

# FVC lit 
for (i in col_names) {
  p <- ggplot(data = merged_fu, aes_string(x=i, y=merged_fu$FVC_lit))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_FVC_lit_fu.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#repeat for IL6 and KL6 

plot_data <- merged_fu %>% filter(IL6<100)

p <- ggplot(data = plot_data, aes(x=IL6, y=FVC_lit))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("IL6_correl_FVC_lit_follow_up_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()

#repeat for KL6 while removing outliers 
plot_data <- merged_fu %>% filter(CA15_3<120)

p <- ggplot(data = plot_data, aes(x=CA15_3, y=FVC_lit))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_correl_FVC_lit_follow_up_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()










# DLCO
for (i in col_names) {
  p <- ggplot(data = merged_fu, aes_string(x=i, y=merged_fu$KCO))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_prism()
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_KCO_fu.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}



#repeat for IL6 and KL6 

plot_data <- merged_fu %>% filter(IL6<100)

p <- ggplot(data = plot_data, aes(x=IL6, y=KCO))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("IL6_correl_KCO_follow_up_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()

#repeat for KL6 while removing outliers 
plot_data <- merged_fu %>% filter(CA15_3<120)

p <- ggplot(data = plot_data, aes(x=CA15_3, y=KCO))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman")+
  theme_prism()

pdf_output <- paste0("Figures/post_covid_", paste0("CA15_3_correl_KCO_follow_up_no_outliers.pdf"))
pdf(file=pdf_output,height = 8, width = 6)
show(p)
dev.off()






####### correlation does not work, but does the change in cytokines levels correlated with FVC ?###### 
mycounts_t
temp <- mycounts_t
last_two_digits <- substr(temp$sample, nchar(temp$sample) - 1, nchar(temp$sample))

# Create the 'subj_id' column by concatenating 'subj' with the last two digits
temp$subj_id <- paste0("subj_", last_two_digits)

baseline_data <- temp %>% filter(sample_type == "baseline") %>% select(-sample_type)
clinic_data <- temp %>% filter(sample_type == "clinic") %>% select(-sample_type)

# Ensure the samples are matched properly by joining them on the 'sample' column
matched_data <- inner_join(baseline_data, clinic_data, by = "subj_id")

# Calculate absolute differences for each cytokine
matched_data$TNFa_diff <- abs(matched_data$TNFa.y - matched_data$TNFa.x)
matched_data$IL6_diff <- abs(matched_data$IL6.y - matched_data$IL6.x)
matched_data$MMP9_diff <- abs(matched_data$MMP9.y - matched_data$MMP9.x)
matched_data$MMP2_diff <- abs(matched_data$MMP2.y - matched_data$MMP2.x)
matched_data$MMP7_diff <- abs(matched_data$MMP7.y - matched_data$MMP7.x)
matched_data$Angiopoietin2_diff <- abs(matched_data$Angiopoietin2.y - matched_data$Angiopoietin2.x)
matched_data$Thrombomodulin_diff <- abs(matched_data$Thrombomodulin.y - matched_data$Thrombomodulin.x)
matched_data$CA15_3_diff <- abs(matched_data$CA15_3.y - matched_data$CA15_3.x)
matched_data$Hyaluronan_diff <- abs(matched_data$Hyaluronan.y - matched_data$Hyaluronan.x)

#get diffecnes data 
difference_data <- matched_data %>% dplyr::select(c(sample.x, sample.y, ends_with("diff")))
difference_data$sample <- difference_data$sample.y
# Print the updated difference_data to verify
print(difference_data)

#add FVC data 
temp <- metadata %>% dplyr::select(FVC_lit, FVC_p, FEV1_lit, FEV1_p, KCO, sample)

temp <- inner_join(difference_data, temp, by="sample")
temp$FVC_lit <- as.numeric(temp$FVC_lit)
temp$FVC_p <- as.numeric(temp$FVC_p)
temp$FEV1_lit <- as.numeric(temp$FEV1_lit)
temp$FEV1_p <- as.numeric(temp$FEV1_p)
temp$KCO <- as.numeric(temp$KCO)

col_names <- colnames(temp)[3:11]

# FVC lit 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$FVC_lit))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_FVC_lit_difference_fu_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


# FVC % 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$FVC_p))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_FVC_per_difference_fu_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}

#KCO %pred
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$KCO))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_correl_KCO_per_difference_fu_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}










#### follow up damples########## 
#### PCA plot according to disease severity at baseline only 

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  select(-sample_type, -sample) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% select(disease_severity_upd)
temp <- temp$disease_severity_upd
pc_scores$disease_severity_upd <- temp
pc_scores$disease_severity_upd <- factor(pc_scores$disease_severity_upd, levels = c("mild_mod", "severe", "critical"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ disease_severity_upd,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="disease_severity_upd",suffixes=c("",".centroid"))


pdf(file = "Figures/fu_samples_PCA_by_disease_sev.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= disease_severity_upd)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("grey", "goldenrod", "orangered")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= disease_severity_upd), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= disease_severity_upd)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Mild Moderate", "Severe", "Critical")), size=10) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)) # p= 0.64

#now repeat for mild-mod vs. severe and critical 

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  select(-sample_type, -sample) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% select(disease_severity_upd)
temp$disease_severity_upd <- ifelse(temp$disease_severity_upd=="mild_mod", "mild_mod", "severe_critical")
temp <- temp$disease_severity_upd
pc_scores$disease_severity_upd <- temp
pc_scores$disease_severity_upd <- factor(pc_scores$disease_severity_upd, levels = c("mild_mod", "severe_critical"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ disease_severity_upd,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="disease_severity_upd",suffixes=c("",".centroid"))

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)
x <- summary(manova(cbind(PC1, PC2) ~ disease_severity_upd, data = pc_scores)) # p= 0.33


pdf(file = "Figures/clinic_samples_PCA_by_disease_sev_mild_mod_vs_sevecritic.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= disease_severity_upd)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("grey", "goldenrod")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= disease_severity_upd), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= disease_severity_upd)) + 
  #labels centroids should be same number of categories
  ggtitle(paste0("p=", paste0(x[4])))+
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Mild Moderate", "Severe Critical")), size=10) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


















#### heatmap of all results 

#add disease seceity 
mycounts_t$disease_severity <- NA
mycounts_t$disease_severity[1:91] <- metadata$disease_severity_two


scaled_data <- mycounts_t %>% select(-c(sample, sample_type, disease_severity)) %>% scale(.)

max(scaled_data)
min(scaled_data)
#set colors 
library(circlize)
col_fun = colorRamp2(c(-3.252642, 0, 8.438242), c("blue", "white", "red"))
col_fun(seq(-3.252642, 0, 8.438242))


#plot heatmap using complex heatmaps 
set.seed(1234)
heatmap_mat <- as.matrix(scaled_data)

# Define colors for each levels of qualitative variables
col = list(Sample_Type = c("clinic"="red", "baseline"="orange", "control"="blue"), 
           Severity=c("mild_mod"="grey", "severe"="darkred", "NA"="white"))
#define annotations
ha <- HeatmapAnnotation(
  Sample_Type=mycounts_t$sample_type,
  Severity=mycounts_t$disease_severity,
  col = col)


# Create the heatmap annotation
library(ComplexHeatmap)


#plot
pdf(file = "Figures/heatmap_cytokines_all_Samples_clustering_rows.pdf", width = 16, height = 6)
ComplexHeatmap::Heatmap(t(heatmap_mat), 
                        column_dend_side = "top", 
                        row_dend_side = "left",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        column_title_side = "bottom", 
                        cluster_rows = TRUE, 
                        cluster_columns = FALSE, 
                        #clustering_method_columns = "aver", 
                        clustering_method_rows = "aver",
                        top_annotation = ha, 
                        row_names_gp = gpar(fontface="bold"), 
                        row_dend_width = unit(5, "cm"), 
                        column_dend_height = unit(5, "cm"))
dev.off()





####### comparison of cytokines levels according to disease severity for each timepoint 

#baseline
temp_dat <- mycounts_t %>%filter(sample_type=="baseline") %>% dplyr::select(-c(sample, sample_type))
baseline_metadata<- metadata %>% filter(sample_type=="baseline")

col_names <- colnames(temp_dat)

temp_dat$disease_severity_upd <- factor(baseline_metadata$disease_severity_upd, levels = c("mild_mod", "severe", "critical"))

for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$disease_severity_upd, y=i, color=temp_dat$disease_severity_upd))+
    geom_boxplot(outlier.shape = NA)+
    scale_color_manual(values=c("grey", "goldenrod", "orangered"))+
    scale_x_discrete(labels=c("Mild Moderate", "Severe", "Critical"))+
    geom_jitter(alpha=0.5, width = 0.2)+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("mild_mod", "severe"), 
                                          c("severe", "critical"), 
                                          c("mild_mod", "critical")))+
    theme_prism()+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_dise_sev_baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}

###fu
temp_dat <- mycounts_t %>%filter(sample_type=="clinic") %>% dplyr::select(-c(sample, sample_type))
fu_metadata<- metadata %>% filter(sample_type=="clinic")

col_names <- colnames(temp_dat)

temp_dat$disease_severity_upd <- factor(fu_metadata$disease_severity_upd, levels = c("mild_mod", "severe", "critical"))

for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$disease_severity_upd, y=i, color=temp_dat$disease_severity_upd))+
    geom_boxplot(outlier.shape = NA)+
    scale_color_manual(values=c("grey", "goldenrod", "orangered"))+
    scale_x_discrete(labels=c("Mild Moderate", "Severe", "Critical"))+
    geom_jitter(alpha=0.5, width = 0.2)+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("mild_mod", "severe"), 
                                          c("severe", "critical"), 
                                          c("mild_mod", "critical")))+
    theme_prism()+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_by_dise_sev_clinic.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}



############# random forest analysis of cytokines levels at two time points. correlating with FVC as outcome #########

####start with follow up samples 
fu_counts <- mycounts_t %>% filter(sample_type=="clinic")
temp_dat <- metadata %>% filter(sample_type=="clinic") %>% dplyr::select(FVC_lit)
fu_counts$FVC_lit <- temp_dat$FVC_lit

temp <- fu_counts %>% dplyr::select(-c(sample, sample_type, disease_severity))
temp <- temp %>% filter(!is.na(FVC_lit))
temp$FVC_lit <- as.numeric(temp$FVC_lit)
temp <- drop_na(temp)
temp_2 <- temp %>% dplyr::select(-FVC_lit)
rf_model <- randomForest(temp_2, as.numeric(temp$FVC_lit), ntree = 500)

### stratify 5-fold CV (80% discovery and 20 validation)
set.seed(1234)
train_index <- sample(nrow(temp), 0.8 * nrow(temp))
train_data <- temp[train_index, ]
test_data <- temp[-train_index, ]

#tune 
mtry <- sqrt(ncol(train_data))
tunegrid <- expand.grid(.mtry=mtry)
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, search = 'random', savePredictions = TRUE, sampling = NULL, classProbs = TRUE)

#run model 
fit <- randomForest(train_data[,-ncol(train_data)], as.numeric(train_data$FVC_lit), 
                    trControl=control, 
                    ntree=500, preProcess=c("BoxCox"), importance=TRUE, scale=TRUE)

#get variable importance
# Extract the mean decrease impurity (MDI) values for each feature
mdi <- randomForest::importance(fit) %>% as.data.frame()

mdi=mdi[order(mdi$IncNodePurity, decreasing = TRUE),]




















# Load required libraries
library(randomForest)
library(caret)

# Assuming 'mycounts_t' contains cytokine data, and 'fu_metadata' contains metadata including 'FVC_lit'

# Combine cytokine data and FVC_lit
fu_counts <- mycounts_t %>%
  filter(sample_type == "clinic") %>%
  dplyr::select(-c(sample, sample_type, disease_severity))

temp_dat <- metadata %>% filter(sample_type=="clinic") %>% dplyr::select(FVC_lit)
fu_counts$FVC_lit <- temp_dat$FVC_lit

# Convert FVC_lit to numeric
fu_counts$FVC_lit <- as.numeric(fu_counts$FVC_lit)

# Drop rows with missing FVC_lit values
fu_counts <- drop_na(fu_counts)
#sett the model and tune it 
#repeat for x times 

# Set the number of repetitions
num_repetitions <- 10  # Choose the number of repetitions you want

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  # Split data into training and test sets
  set.seed(1234 + i)
  train_index <- sample(nrow(fu_counts), 0.8 * nrow(fu_counts))
  train_data <- fu_counts[train_index, ]
  test_data <- fu_counts[-train_index, ]
  
  # Train the random forest model
  rf_model <- randomForest(FVC_lit ~ ., data = train_data, ntree = 500)
  

  # Get variable importance
  mdi <- data.frame(importance(rf_model))
  mdi <- mdi %>% arrange(desc(IncNodePurity))
  mdi$variable <- rownames(mdi)
  
  # Store results in the data frame
  results <- rbind(results, data.frame(Replication = i, mdi))
}

# Create a boxplot
pdf(file = "Figures/rf_improtance_fu_FVC.pdf", height = 10, width = 12)
ggplot(results, aes(x=IncNodePurity, y=reorder(variable, IncNodePurity))) +
  geom_boxplot(color="orange") +
  geom_jitter(alpha=0.5)+
  guides(fill="none", size="none")+
  ylab("")+xlab("Importance")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 24, face = "bold"))
dev.off()


#now repeat for KCO 

# Combine cytokine data and FVC_lit
fu_counts <- mycounts_t %>%
  filter(sample_type == "clinic") %>%
  dplyr::select(-c(sample, sample_type, disease_severity, cluster))

fu_counts$KCO <- fu_metadata$KCO

# Convert FVC_lit to numeric
fu_counts$KCO <- as.numeric(fu_counts$KCO)

# Drop rows with missing FVC_lit values
fu_counts <- drop_na(fu_counts)
#sett the model and tune it 
#repeat for x times 

# Set the number of repetitions
num_repetitions <- 10  # Choose the number of repetitions you want

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  # Split data into training and test sets
  set.seed(1234 + i)
  train_index <- sample(nrow(fu_counts), 0.8 * nrow(fu_counts))
  train_data <- fu_counts[train_index, ]
  test_data <- fu_counts[-train_index, ]
  
  # Train the random forest model
  rf_model <- randomForest(KCO ~ ., data = train_data, ntree = 500)
  
  
  # Get variable importance
  mdi <- data.frame(importance(rf_model))
  mdi <- mdi %>% arrange(desc(IncNodePurity))
  mdi$variable <- rownames(mdi)
  
  # Store results in the data frame
  results <- rbind(results, data.frame(Replication = i, mdi))
}

# Create a boxplot
pdf(file = "Figures/rf_improtance_fu_KCO.pdf", height = 10, width = 12)
ggplot(results, aes(x=IncNodePurity, y=reorder(variable, IncNodePurity))) +
  geom_boxplot(color="orange") +
  geom_jitter(alpha=0.5)+
  guides(fill="none", size="none")+
  ylab("")+xlab("Importance")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 24, face = "bold"))
dev.off()

  




#### random foerst for basleine smaples 

# Combine cytokine data and FVC_lit
baseline_counts <- mycounts_t %>%
  filter(sample_type == "baseline") %>%
  dplyr::select(-c(sample, sample_type, disease_severity, cluster))

baseline_counts$FVC_lit <- baseline_metadata$FVC_lit

# Convert FVC_lit to numeric
baseline_counts$FVC_lit <- as.numeric(baseline_counts$FVC_lit)

# Drop rows with missing FVC_lit values
baseline_counts <- drop_na(baseline_counts)
#sett the model and tune it 
#repeat for x times 

# Set the number of repetitions
num_repetitions <- 10  # Choose the number of repetitions you want

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  # Split data into training and test sets
  set.seed(1234 + i)
  train_index <- sample(nrow(baseline_counts), 0.8 * nrow(baseline_counts))
  train_data <- baseline_counts[train_index, ]
  test_data <- baseline_counts[-train_index, ]
  
  # Train the random forest model
  rf_model <- randomForest(FVC_lit ~ ., data = train_data, ntree = 500)
  
  
  # Get variable importance
  mdi <- data.frame(importance(rf_model))
  mdi <- mdi %>% arrange(desc(IncNodePurity))
  mdi$variable <- rownames(mdi)
  
  # Store results in the data frame
  results <- rbind(results, data.frame(Replication = i, mdi))
}

# Create a boxplot
pdf(file = "Figures/rf_improtance_baseline_FVC.pdf", height = 10, width = 12)
ggplot(results, aes(x=IncNodePurity, y=reorder(variable, IncNodePurity))) +
  geom_boxplot(color="red") +
  geom_jitter(alpha=0.5)+
  guides(fill="none", size="none")+
  ylab("")+xlab("Importance")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 24, face = "bold"))
dev.off()


#now repeat for KCO 


# Combine cytokine data and KCO
baseline_counts <- mycounts_t %>%
  filter(sample_type == "baseline") %>%
  dplyr::select(-c(sample, sample_type, disease_severity, cluster))

baseline_counts$KCO <- baseline_metadata$KCO

# Convert KCO to numeric
baseline_counts$KCO <- as.numeric(baseline_counts$KCO)

# Drop rows with missing KCO values
baseline_counts <- drop_na(baseline_counts)
#sett the model and tune it 
#repeat for x times 

# Set the number of repetitions
num_repetitions <- 10  # Choose the number of repetitions you want

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  # Split data into training and test sets
  set.seed(1234 + i)
  train_index <- sample(nrow(baseline_counts), 0.8 * nrow(baseline_counts))
  train_data <- baseline_counts[train_index, ]
  test_data <- baseline_counts[-train_index, ]
  
  # Train the random forest model
  rf_model <- randomForest(KCO ~ ., data = train_data, ntree = 500)
  
  
  # Get variable importance
  mdi <- data.frame(importance(rf_model))
  mdi <- mdi %>% arrange(desc(IncNodePurity))
  mdi$variable <- rownames(mdi)
  
  # Store results in the data frame
  results <- rbind(results, data.frame(Replication = i, mdi))
}

# Create a boxplot
pdf(file = "Figures/rf_improtance_baseline_KCO.pdf", height = 10, width = 12)
ggplot(results, aes(x=IncNodePurity, y=reorder(variable, IncNodePurity))) +
  geom_boxplot(color="red") +
  geom_jitter(alpha=0.5)+
  guides(fill="none", size="none")+
  ylab("")+xlab("Importance")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 24, face = "bold"))
dev.off()








############ Next: correlation of cytokines with CT at follow up ##########

temp <- mycounts_t
temp <- temp %>% filter(sample_type!="control")

temp$Age <- metadata$Age
temp$Sex <- metadata$Gender
temp$IHD <- metadata$IHD
temp$CKD <- metadata$CKD
metadata$lung_dis <- ifelse(metadata$Asthma=="Y", "Y", 
                            ifelse(metadata$COPD=="Y", "Y", 
                                   ifelse(metadata$other_lung_disease!="N", "Y", "N")))
temp$lung_dis <- metadata$lung_dis
temp$Death <- metadata$Death

temp_2 <- temp %>% filter(sample_type!="clinic")
temp_2 <- temp %>% dplyr::select(-c(sample_type, sample))
temp_2$Death <- ifelse(temp_2$Death=="Y", "1", "0")
temp_2$Death <- factor(temp_2$Death)

# Fit the model for baseline only and moretality as outcome 
library(broom)

# Fit the logistic regression model
model <- glm(Death ~ ., data = temp_2, family = binomial)

# Summarize the model with tidy and confint_tidy from the broom package
tidy_model <- tidy(model)
confint_model <- confint_tidy(model)

# Combine the results
summary_model <- cbind(tidy_model, confint_model)
summary_model <- summary_model[, c("term", "estimate", "conf.low", "conf.high", "p.value")]

# Rename columns for clarity
colnames(summary_model) <- c("Term", "Log(Odds Ratio)", "Conf.Low", "Conf.High", "p.value")

# Add Odds Ratios
summary_model$`Odds Ratio` <- exp(summary_model$`Log(Odds Ratio)`)

#add adj p value 
summary_model$adj.p <- p.adjust(summary_model$p.value, method = "BH")

# Print summary
print(summary_model)

#save 
write.csv(summary_model, file = "Results/OR_baseline_mortality_as_outcome.csv")






###### reapt model for follow up samples looking at abnormal CT or abnormal FVC as outcome 
#does not work due to small number of samples  
###ignore it 


temp <- mycounts_t
temp <- temp %>% filter(sample_type!="control")

temp$Age <- metadata$Age
temp$Sex <- metadata$Gender
temp$IHD <- metadata$IHD
temp$CKD <- metadata$CKD
metadata$lung_dis <- ifelse(metadata$Asthma=="Y", "Y", 
                            ifelse(metadata$COPD=="Y", "Y", 
                                   ifelse(metadata$other_lung_disease!="N", "Y", "N")))
temp$lung_dis <- metadata$lung_dis
temp$abnormal_CT_fu_y_n <- metadata$abnormal_CT_fu_y_n
temp$abnormal_fvc <- ifelse(metadata$FVC_p<80, "y", "n")

temp_2 <- temp %>% filter(sample_type!="baseline")
temp_2 <- temp %>% dplyr::select(-c(sample_type, sample, disease_severity))
temp_2$abnormal_CT_fu_y_n <- factor(temp_2$abnormal_CT_fu_y_n)
temp_2$abnormal_fvc <- factor(temp_2$abnormal_fvc)

# Fit the model for baseline only and moretality as outcome 
library(broom)

# Fit the logistic regression model
model <- glm(abnormal_CT_fu_y_n ~ ., data = temp_2, family = binomial)

# Summarize the model with tidy and confint_tidy from the broom package
tidy_model <- tidy(model)
confint_model <- confint_tidy(model)

# Combine the results
summary_model <- cbind(tidy_model, confint_model)
summary_model <- summary_model[, c("term", "estimate", "conf.low", "conf.high", "p.value")]

# Rename columns for clarity
colnames(summary_model) <- c("Term", "Log(Odds Ratio)", "Conf.Low", "Conf.High", "p.value")

# Add Odds Ratios
summary_model$`Odds Ratio` <- exp(summary_model$`Log(Odds Ratio)`)

#add adj p value 
summary_model$adj.p <- p.adjust(summary_model$p.value, method = "BH")

# Print summary
print(summary_model)

#save 
write.csv(summary_model, file = "Results/OR_baseline_mortality_as_outcome.csv")







##### build PCoA that differntiate baseline and follow up samples accoridng to : disease severity, ICU, Intubation, abnormal CT FU, abnormal FVC fu 


###### PCOAs of baseline samples only ####

#disease severity 

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="baseline") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="baseline") %>% dplyr::select(disease_severity)
temp <- temp$disease_severity
pc_scores$disease_severity <- temp
pc_scores$disease_severity <- ifelse(pc_scores$disease_severity=="mild", "mild_mod", 
                                     ifelse(pc_scores$disease_severity=="moderate", "mild_mod", 
                                            ifelse(pc_scores$disease_severity=="critical", "severe", "severe")))
pc_scores$disease_severity <- factor(pc_scores$disease_severity, levels = c("mild_mod", "severe"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ disease_severity,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="disease_severity",suffixes=c("",".centroid"))


pdf(file = "Figures/baseline_samples_PCA_by_disease_sev.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= disease_severity)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= disease_severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= disease_severity)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Mild Moderate", "Severe")), size=10) +
  ggtitle("p=0.33")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ disease_severity, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ disease_severity, data = pc_scores)) # p= 0.3373




#by intubation y_n 
# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="baseline") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="baseline") %>% dplyr::select(Intubated)
temp <- temp$Intubated
pc_scores$Intubated <- temp
pc_scores$Intubated <- factor(pc_scores$Intubated, levels = c("N", "Y"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Intubated,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="Intubated",suffixes=c("",".centroid"))


pdf(file = "Figures/baseline_samples_PCA_by_Intubated.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Intubated)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Intubated), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Intubated)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non-Intubated", "Intubated")), size=10) +
  ggtitle("p=0.819")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ Intubated, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ Intubated, data = pc_scores)) # p= 0.3373



#### by ICU _y_n 
# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="baseline") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="baseline") %>% dplyr::select(ICU_y_n)
temp <- temp$ICU_y_n
pc_scores$ICU_y_n <- temp
pc_scores$ICU_y_n <- factor(pc_scores$ICU_y_n, levels = c("N", "Y"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ ICU_y_n,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="ICU_y_n",suffixes=c("",".centroid"))


pdf(file = "Figures/baseline_samples_PCA_by_ICU_y_n.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= ICU_y_n)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= ICU_y_n), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= ICU_y_n)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non-ICU", "ICU")), size=10) +
  ggtitle("p=0.7765")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ ICU_y_n, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ ICU_y_n, data = pc_scores)) # p= 0.3373





###### PCOAs of follow up samples ########



#disease severity 

# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% dplyr::select(disease_severity)
temp <- temp$disease_severity
pc_scores$disease_severity <- temp
pc_scores$disease_severity <- ifelse(pc_scores$disease_severity=="mild", "mild_mod", 
                                     ifelse(pc_scores$disease_severity=="moderate", "mild_mod", 
                                            ifelse(pc_scores$disease_severity=="critical", "severe", "severe")))
pc_scores$disease_severity <- factor(pc_scores$disease_severity, levels = c("mild_mod", "severe"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ disease_severity,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="disease_severity",suffixes=c("",".centroid"))


pdf(file = "Figures/clinic_samples_PCA_by_disease_sev.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= disease_severity)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= disease_severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= disease_severity)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Mild Moderate", "Severe")), size=10) +
  ggtitle("p=0.8849")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ disease_severity, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ disease_severity, data = pc_scores)) # p= 0.3373




#by intubation y_n 
# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% dplyr::select(Intubated)
temp <- temp$Intubated
pc_scores$Intubated <- temp
pc_scores$Intubated <- factor(pc_scores$Intubated, levels = c("N", "Y"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Intubated,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="Intubated",suffixes=c("",".centroid"))


pdf(file = "Figures/clinic_samples_PCA_by_Intubated.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Intubated)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Intubated), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Intubated)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non-Intubated", "Intubated")), size=10) +
  ggtitle("p=0.7409")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ Intubated, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ Intubated, data = pc_scores)) # p= 0.3373



#### by ICU _y_n 
# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% dplyr::select(ICU_y_n)
temp <- temp$ICU_y_n
pc_scores$ICU_y_n <- temp
pc_scores$ICU_y_n <- factor(pc_scores$ICU_y_n, levels = c("N", "Y"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ ICU_y_n,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="ICU_y_n",suffixes=c("",".centroid"))


pdf(file = "Figures/clinic_samples_PCA_by_ICU_y_n.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= ICU_y_n)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= ICU_y_n), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= ICU_y_n)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non-ICU", "ICU")), size=10) +
  ggtitle("p=0.5413")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ ICU_y_n, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ ICU_y_n, data = pc_scores)) # p= 0.3373




#### adding abnormal CT y_n 
# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% dplyr::select(abnormal_CT_fu_y_n)
temp <- temp$abnormal_CT_fu_y_n
pc_scores$abnormal_CT_fu_y_n <- temp
pc_scores$abnormal_CT_fu_y_n <- factor(pc_scores$abnormal_CT_fu_y_n, levels = c("N", "Y"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ abnormal_CT_fu_y_n,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="abnormal_CT_fu_y_n",suffixes=c("",".centroid"))


pdf(file = "Figures/clinic_samples_PCA_by_abnormal_CT_fu_y_n.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= abnormal_CT_fu_y_n)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= abnormal_CT_fu_y_n), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= abnormal_CT_fu_y_n)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Normal", "Abnormal")), size=10) +
  ggtitle("p=0.1914")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ abnormal_CT_fu_y_n, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ abnormal_CT_fu_y_n, data = pc_scores)) # p= 0.3373



#### abnormal FVC y_n 
# Extract numerical columns for PCA
pca_data <- mycounts_t %>%
  filter(sample_type=="clinic") %>% 
  dplyr::select(-sample_type, -sample, -disease_severity) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PC scores
pc_scores <- as.data.frame(pca_result$x)

# Add sample information to PC scores
temp <- metadata %>%filter(sample_type=="clinic") %>% dplyr::select(FVC_p)
temp <- temp$FVC_p
pc_scores$FVC_p <- temp
pc_scores$FVC_p <- ifelse(pc_scores$FVC_p<80, "Y", "N")
pc_scores$FVC_p <- factor(pc_scores$FVC_p, levels = c("N", "Y"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ FVC_p,data= pc_scores, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(pc_scores,centroids,by="FVC_p",suffixes=c("",".centroid"))


pdf(file = "Figures/clinic_samples_PCA_by_abnormal_FVC_p.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= FVC_p)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green3", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= FVC_p), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= FVC_p)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Normal", "Abnormal")), size=10) +
  ggtitle("p=0.1682")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

# Assuming 'pc_scores' is your PCA scores data frame

# Perform MANOVA
manova_result <- manova(cbind(PC1, PC2) ~ FVC_p, data = pc_scores)

summary(manova(cbind(PC1, PC2) ~ FVC_p, data = pc_scores)) # p= 0.3373







































# Load necessary libraries
library(tidyverse)
library(glmnet)
library(caret)
library(boot)

temp <- mycounts_t
last_two_digits <- substr(temp$sample, nchar(temp$sample) - 1, nchar(temp$sample))

# Create the 'subj_id' column by concatenating 'subj' with the last two digits
temp$subj_id <- paste0("subj_", last_two_digits)

baseline_data <- temp %>% filter(sample_type == "baseline") %>% select(-sample_type)
clinic_data <- temp %>% filter(sample_type == "clinic") %>% select(-sample_type)

# Ensure the samples are matched properly by joining them on the 'sample' column
matched_data <- inner_join(baseline_data, clinic_data, by = "subj_id")

# Calculate absolute differences for each cytokine
matched_data$TNFa_diff <- matched_data$TNFa.y - matched_data$TNFa.x
matched_data$IL6_diff <- matched_data$IL6.y - matched_data$IL6.x
matched_data$MMP9_diff <- matched_data$MMP9.y - matched_data$MMP9.x
matched_data$MMP2_diff <- matched_data$MMP2.y - matched_data$MMP2.x
matched_data$MMP7_diff <- matched_data$MMP7.y - matched_data$MMP7.x
matched_data$Angiopoietin2_diff <- matched_data$Angiopoietin2.y - matched_data$Angiopoietin2.x
matched_data$Thrombomodulin_diff <- matched_data$Thrombomodulin.y - matched_data$Thrombomodulin.x
matched_data$CA15_3_diff <- matched_data$CA15_3.y - matched_data$CA15_3.x
matched_data$Hyaluronan_diff <- matched_data$Hyaluronan.y - matched_data$Hyaluronan.x

#get diffecnes data 
difference_data <- matched_data %>% dplyr::select(c(sample.x, sample.y, ends_with("diff")))
difference_data$sample <- difference_data$sample.y
# Print the updated difference_data to verify
print(difference_data)






######comaprre differnce data #######

temp <- inner_join(difference_data, metadata, by="sample")

#overall median and CI 
stats_data <- temp %>% 
  dplyr::select(ends_with("diff")) %>% 
  summarize_all(funs(Median = median(., na.rm = TRUE),
       Q1 = quantile(., 0.25, na.rm = TRUE),
       Q3 = quantile(., 0.75, na.rm = TRUE))) %>% 
  t(.)
stats_data
#save 
write.csv(stats_data, file="Results/difference_data_stats.csv")


###compare according to intubation 
temp$Intubated <- as.factor(temp$Intubated)

col_names <- temp %>% dplyr::select(ends_with("diff"))
col_names<- colnames(col_names)

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=temp$Intubated, y=i, color=as.factor(temp$Intubated)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Non-Intubated", "Intubated"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_difference_by_intubation.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}



####compare according to ICU y_n 
temp$ICU_y_n <- as.factor(temp$ICU_y_n)

col_names <- temp %>% dplyr::select(ends_with("diff"))
col_names<- colnames(col_names)

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=temp$ICU_y_n, y=i, color=as.factor(temp$ICU_y_n)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Non-Intubated", "Intubated"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_difference_by_ICU_y_n.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}












####### correlation btween cytokines data and FVC_lit########
library(MASS)
library(parameters)
library(GGally)
#correltion matrix of baseline wtih FVC_lit
temp <- inner_join(baseline_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(FVC_lit))
features <- temp %>% dplyr::select(c(TNFa, IL6, MMP2, MMP7, MMP9, Angiopoietin2, Thrombomodulin, CA15_3, Hyaluronan))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                     crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(FVC_lit)

skewness <- apply(features, 2, function(x) skewness(x, na.rm = TRUE))
print(skewness)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)

#extract the coef 
cor_matrix <- as.data.frame(cor_matrix)
#save
write.csv(cor_matrix, file = "Results/correlation_coef_baseline.csv")


library(corrplot)


# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_baseline_with_FVC_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients

dev.off()

### compute p values 
data$outcome <- data$FVC_lit
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
#save 
write.csv(correlation_results, file = "Results/correlation_cytokines_baseline_fvc_lit_at_fu.csv")


# get all p values and put into a table 
library(Hmisc)
# Combine the features and outcome into a single data frame
data <- cbind(features, outcome)

# Compute the correlation matrix along with p-values using rcorr
cor_results <- rcorr(as.matrix(data), type = "spearman")



# Extract the p-values matrix
p_matrix <- cor_results$P
print("P-values Matrix:")
print(p_matrix)
p_matrix <- as.data.frame(p_matrix)
write.csv(p_matrix, file = "Results/correlation_p_values_baseline.csv")


################ ridge moidel for baseline data 
# Load necessary library
library(glmnet)
library(boot)

# Prepare the data
x <- as.matrix(scale(features))
y <- outcome$FVC_lit

# Perform Ridge Regression with cross-validation to find the optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0)  # alpha = 0 for Ridge Regression
ridge_model <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)

# Function to fit ridge regression and extract coefficients
ridge_fit <- function(data, indices) {
  boot_x <- data[indices, -ncol(data)]
  boot_y <- data[indices, ncol(data)]
  model <- glmnet(boot_x, boot_y, alpha = 0, lambda = cv_ridge$lambda.min)
  return(as.vector(coef(model)))
}

# Combine features and outcome for bootstrapping
data <- data.frame(cbind(x, y))

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data, ridge_fit, R = 1000)

# Extract bootstrap estimates and calculate standard errors and confidence intervals
bootstrap_coefs <- t(sapply(1:ncol(data), function(i) boot::boot.ci(bootstrap_results, type = "perc", index = i)$percent))
coef_names <- rownames(coef(ridge_model))
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = as.vector(coef(ridge_model)),
  `2.5 %` = bootstrap_coefs[, 4],
  `97.5 %` = bootstrap_coefs[, 5]
)

print(bootstrap_summary)

#save 
write.csv(bootstrap_summary, file = "Results/Coef_summary_ridge_model_baseline_data_fvc_lit.csv")

# Plot the coefficients with error bars

p <- ggplot(bootstrap_summary[-1, ], aes(x = Estimate, y = Coefficient)) + 
  geom_point() +
  geom_errorbar(aes(xmin = X2.5.., xmax = X97.5..), width = 0.2) +
  theme_minimal() +
  labs(title = "Ridge Regression Coefficients with 95% Confidence Intervals",
       x = "Estimate", 
       y="")+
  theme(axis.text.y = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(size = 18, face = "bold"))

pdf(file = "Figures/Coef_summary_ridge_model_baseline_data_fvc_lit.pdf", height = 8, width = 10)
p
dev.off()


# LOOCV for Ridge Regression
loocv_ridge <- train(y ~ ., data = cbind(scale(features), y), method = "glmnet", trControl = loocv_results,
                     tuneGrid = expand.grid(alpha = 0, lambda = cv_ridge$lambda.min))
print(loocv_ridge)

#extract model performance 
model_perf <- loocv_ridge$results

#save
write.csv(model_perf, file = "Results/ridge_model_performance_baseline_data_fvc_lit.csv")






############ repeeat for follow up data 
temp <- inner_join(clinic_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(FVC_lit))
features <- temp %>% dplyr::select(c(TNFa, IL6, MMP2, MMP7, MMP9, Angiopoietin2, Thrombomodulin, CA15_3, Hyaluronan))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(FVC_lit)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


#extract the coef 
cor_matrix <- as.data.frame(cor_matrix)
#save
write.csv(cor_matrix, file = "Results/correlation_coef_follow_up.csv")



library(corrplot)

# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_clinic_with_FVC_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients
dev.off()

### compute p values 
data$outcome <- data$FVC_lit
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
write.csv(correlation_results, file = "Results/correlation_cytokines_clinic_fvc_lit_at_fu.csv")

# get all p values and put into a table 
library(Hmisc)
# Combine the features and outcome into a single data frame
data <- cbind(features, outcome)

# Compute the correlation matrix along with p-values using rcorr
cor_results <- rcorr(as.matrix(data), type = "spearman")

# Extract the p-values matrix
p_matrix <- cor_results$P
print("P-values Matrix:")
print(p_matrix)
p_matrix <- as.data.frame(p_matrix)
write.csv(p_matrix, file = "Results/correlation_p_values_clinic.csv")


################ ridge moidel for follow up  data 
# Load necessary library
# Load necessary library
library(glmnet)
library(boot)

# Prepare the data
x <- as.matrix(scale(features))
y <- outcome$FVC_lit

# Perform Ridge Regression with cross-validation to find the optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0)  # alpha = 0 for Ridge Regression
ridge_model <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)

# Function to fit ridge regression and extract coefficients
ridge_fit <- function(data, indices) {
  boot_x <- data[indices, -ncol(data)]
  boot_y <- data[indices, ncol(data)]
  model <- glmnet(boot_x, boot_y, alpha = 0, lambda = cv_ridge$lambda.min)
  return(as.vector(coef(model)))
}

# Combine features and outcome for bootstrapping
data <- data.frame(cbind(x, y))

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data, ridge_fit, R = 1000)

# Extract bootstrap estimates and calculate standard errors and confidence intervals
bootstrap_coefs <- t(sapply(1:ncol(data), function(i) boot::boot.ci(bootstrap_results, type = "perc", index = i)$percent))
coef_names <- rownames(coef(ridge_model))
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = as.vector(coef(ridge_model)),
  `2.5 %` = bootstrap_coefs[, 4],
  `97.5 %` = bootstrap_coefs[, 5]
)

print(bootstrap_summary)

#save 
write.csv(bootstrap_summary, file = "Results/Coef_summary_ridge_model_follow_up_data_fvc_lit.csv")

# Plot the coefficients with error bars
p <- ggplot(bootstrap_summary[-1, ], aes(x = Estimate, y = Coefficient)) + 
  geom_point() +
  geom_errorbar(aes(xmin = X2.5.., xmax = X97.5..), width = 0.2) +
  theme_minimal() +
  labs(title = "Ridge Regression Coefficients with 95% Confidence Intervals",
       x = "Estimate", 
       y="")+
  theme(axis.text.y = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(size = 18, face = "bold"))

pdf(file = "Figures/Coef_summary_ridge_model_follow_up_data_fvc_lit.pdf", height = 8, width = 10)
p
dev.off()


# LOOCV for Ridge Regression
loocv_ridge <- train(y ~ ., data = cbind(scale(features), y), method = "glmnet", trControl = loocv_results,
                     tuneGrid = expand.grid(alpha = 0, lambda = cv_ridge$lambda.min))
print(loocv_ridge)

#extract model performance 
model_perf <- loocv_ridge$results

#save
write.csv(model_perf, file = "Results/ridge_model_performance_follow_up_data_fvc_lit.csv")








#### repeat for differnce data 
temp <- inner_join(difference_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(FVC_lit))
features <- temp %>% dplyr::select(c(TNFa_diff, IL6_diff, MMP2_diff, MMP7_diff, MMP9_diff,
                                     Angiopoietin2_diff, Thrombomodulin_diff, CA15_3_diff, Hyaluronan_diff))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(FVC_lit)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


#extract the coef 
cor_matrix <- as.data.frame(cor_matrix)
#save
write.csv(cor_matrix, file = "Results/correlation_coef_difference.csv")



library(corrplot)

# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_difference_with_FVC_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients
dev.off()

### compute p values 
data$outcome <- data$FVC_lit
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
write.csv(correlation_results, file = "Results/correlation_cytokines_difference_fvc_lit_at_fu.csv")

# get all p values and put into a table 
library(Hmisc)
# Combine the features and outcome into a single data frame
data <- cbind(features, outcome)

# Compute the correlation matrix along with p-values using rcorr
cor_results <- rcorr(as.matrix(data), type = "spearman")

# Extract the p-values matrix
p_matrix <- cor_results$P
print("P-values Matrix:")
print(p_matrix)
p_matrix <- as.data.frame(p_matrix)
write.csv(p_matrix, file = "Results/correlation_p_values_difference.csv")

# Plotting the correlation matrix with corrplot
#corrplot::corrplot(cor_matrix, method = "circle", type = "upper", 
#                   tl.col = "black", p.mat = p_matrix, sig.level = 0.05)






















# Load necessary library
library(glmnet)
library(boot)

# Prepare the data
x <- as.matrix(scale(features))
y <- outcome$FVC_lit

# Perform Ridge Regression with cross-validation to find the optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0)  # alpha = 0 for Ridge Regression
ridge_model <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)

# Function to fit ridge regression and extract coefficients
ridge_fit <- function(data, indices) {
  boot_x <- data[indices, -ncol(data)]
  boot_y <- data[indices, ncol(data)]
  model <- glmnet(boot_x, boot_y, alpha = 0, lambda = cv_ridge$lambda.min)
  return(as.vector(coef(model)))
}

# Combine features and outcome for bootstrapping
data <- data.frame(cbind(x, y))

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data, ridge_fit, R = 1000)

# Extract bootstrap estimates and calculate standard errors and confidence intervals
bootstrap_coefs <- t(sapply(1:ncol(data), function(i) boot::boot.ci(bootstrap_results, type = "perc", index = i)$percent))
coef_names <- rownames(coef(ridge_model))
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = as.vector(coef(ridge_model)),
  `2.5 %` = bootstrap_coefs[, 4],
  `97.5 %` = bootstrap_coefs[, 5]
)

print(bootstrap_summary)

#save 
write.csv(bootstrap_summary, file = "Results/Coef_summary_ridge_model_differnce_data_fvc_lit.csv")

# Plot the coefficients with error bars
p <- ggplot(bootstrap_summary[-1, ], aes(x = Estimate, y = Coefficient)) + 
  geom_point() +
  geom_errorbar(aes(xmin = X2.5.., xmax = X97.5..), width = 0.2) +
  theme_minimal() +
  labs(title = "Ridge Regression Coefficients with 95% Confidence Intervals",
       x = "Estimate", 
       y="")+
  theme(axis.text.y = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(size = 18, face = "bold"))

pdf(file = "Figures/Coef_summary_ridge_model_differnce_data_fvc_lit.pdf", height = 8, width = 10)
p
dev.off()


# LOOCV for Ridge Regression
loocv_ridge <- train(y ~ ., data = cbind(scale(features), y), method = "glmnet", trControl = loocv_results,
                     tuneGrid = expand.grid(alpha = 0, lambda = cv_ridge$lambda.min))
print(loocv_ridge)

#extract model performance 
model_perf <- loocv_ridge$results

#save
write.csv(model_perf, file = "Results/ridge_model_performance_difference_data_fvc_lit.csv")
















#repeat for FVC % as outcome 
####### correlation btween cytokines data and FVC_p########

#correltion matrix of baseline wtih FVC_p
temp <- inner_join(baseline_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(FVC_p))
features <- temp %>% dplyr::select(c(TNFa, IL6, MMP2, MMP7, MMP9, Angiopoietin2, Thrombomodulin, CA15_3, Hyaluronan))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(FVC_p)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


library(corrplot)


# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_baseline_with_FVC_p_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients

dev.off()

### compute p values 
data$outcome <- data$FVC_p
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
#save 
write.csv(correlation_results, file = "Results/correlation_cytokines_baseline_FVC_p_at_fu.csv")




################ ridge moidel for baseline data 
# Load necessary library
library(glmnet)
library(boot)

# Prepare the data
x <- as.matrix(scale(features))
y <- outcome$FVC_p

# Perform Ridge Regression with cross-validation to find the optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0)  # alpha = 0 for Ridge Regression
ridge_model <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)

# Function to fit ridge regression and extract coefficients
ridge_fit <- function(data, indices) {
  boot_x <- data[indices, -ncol(data)]
  boot_y <- data[indices, ncol(data)]
  model <- glmnet(boot_x, boot_y, alpha = 0, lambda = cv_ridge$lambda.min)
  return(as.vector(coef(model)))
}

# Combine features and outcome for bootstrapping
data <- data.frame(cbind(x, y))

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data, ridge_fit, R = 1000)

# Extract bootstrap estimates and calculate standard errors and confidence intervals
bootstrap_coefs <- t(sapply(1:ncol(data), function(i) boot::boot.ci(bootstrap_results, type = "perc", index = i)$percent))
coef_names <- rownames(coef(ridge_model))
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = as.vector(coef(ridge_model)),
  `2.5 %` = bootstrap_coefs[, 4],
  `97.5 %` = bootstrap_coefs[, 5]
)

print(bootstrap_summary)

#save 
write.csv(bootstrap_summary, file = "Results/Coef_summary_ridge_model_baseline_data_FVC_p.csv")

# Plot the coefficients with error bars

p <- ggplot(bootstrap_summary[-1, ], aes(x = Estimate, y = Coefficient)) + 
  geom_point() +
  geom_errorbar(aes(xmin = X2.5.., xmax = X97.5..), width = 0.2) +
  theme_minimal() +
  labs(title = "Ridge Regression Coefficients with 95% Confidence Intervals",
       x = "Estimate", 
       y="")+
  theme(axis.text.y = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(size = 18, face = "bold"))

pdf(file = "Figures/Coef_summary_ridge_model_baseline_data_FVC_p.pdf", height = 8, width = 10)
p
dev.off()


# LOOCV for Ridge Regression
loocv_ridge <- train(y ~ ., data = cbind(scale(features), y), method = "glmnet", trControl = loocv_results,
                     tuneGrid = expand.grid(alpha = 0, lambda = cv_ridge$lambda.min))
print(loocv_ridge)

#extract model performance 
model_perf <- loocv_ridge$results

#save
write.csv(model_perf, file = "Results/ridge_model_performance_baseline_data_FVC_p.csv")






############ repeeat for follow up data 
temp <- inner_join(clinic_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(FVC_p))
features <- temp %>% dplyr::select(c(TNFa, IL6, MMP2, MMP7, MMP9, Angiopoietin2, Thrombomodulin, CA15_3, Hyaluronan))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(FVC_p)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


library(corrplot)

# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_clinic_with_FVC_p_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients
dev.off()

### compute p values 
data$outcome <- data$FVC_p
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
write.csv(correlation_results, file = "Results/correlation_cytokines_clinic_FVC_p_at_fu.csv")



################ ridge moidel for follow up  data 
# Load necessary library
# Load necessary library
library(glmnet)
library(boot)

# Prepare the data
x <- as.matrix(scale(features))
y <- outcome$FVC_p

# Perform Ridge Regression with cross-validation to find the optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0)  # alpha = 0 for Ridge Regression
ridge_model <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)

# Function to fit ridge regression and extract coefficients
ridge_fit <- function(data, indices) {
  boot_x <- data[indices, -ncol(data)]
  boot_y <- data[indices, ncol(data)]
  model <- glmnet(boot_x, boot_y, alpha = 0, lambda = cv_ridge$lambda.min)
  return(as.vector(coef(model)))
}

# Combine features and outcome for bootstrapping
data <- data.frame(cbind(x, y))

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data, ridge_fit, R = 1000)

# Extract bootstrap estimates and calculate standard errors and confidence intervals
bootstrap_coefs <- t(sapply(1:ncol(data), function(i) boot::boot.ci(bootstrap_results, type = "perc", index = i)$percent))
coef_names <- rownames(coef(ridge_model))
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = as.vector(coef(ridge_model)),
  `2.5 %` = bootstrap_coefs[, 4],
  `97.5 %` = bootstrap_coefs[, 5]
)

print(bootstrap_summary)

#save 
write.csv(bootstrap_summary, file = "Results/Coef_summary_ridge_model_follow_up_data_FVC_p.csv")

# Plot the coefficients with error bars

p <- ggplot(bootstrap_summary[-1, ], aes(x = Estimate, y = Coefficient)) + 
  geom_point() +
  geom_errorbar(aes(xmin = X2.5.., xmax = X97.5..), width = 0.2) +
  theme_minimal() +
  labs(title = "Ridge Regression Coefficients with 95% Confidence Intervals",
       x = "Estimate", 
       y="")+
  theme(axis.text.y = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(size = 18, face = "bold"))

pdf(file = "Figures/Coef_summary_ridge_model_follow_up_data_FVC_p.pdf", height = 8, width = 10)
p
dev.off()


# LOOCV for Ridge Regression
loocv_ridge <- train(y ~ ., data = cbind(scale(features), y), method = "glmnet", trControl = loocv_results,
                     tuneGrid = expand.grid(alpha = 0, lambda = cv_ridge$lambda.min))
print(loocv_ridge)

#extract model performance 
model_perf <- loocv_ridge$results

#save
write.csv(model_perf, file = "Results/ridge_model_performance_follow_up_data_FVC_p.csv")








#### repeat for differnce data 
temp <- inner_join(difference_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(FVC_p))
features <- temp %>% dplyr::select(c(TNFa_diff, IL6_diff, MMP2_diff, MMP7_diff, MMP9_diff,
                                     Angiopoietin2_diff, Thrombomodulin_diff, CA15_3_diff, Hyaluronan_diff))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(FVC_p)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


library(corrplot)

# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_difference_with_FVC_p_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients
dev.off()

### compute p values 
data$outcome <- data$FVC_p
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
write.csv(correlation_results, file = "Results/correlation_cytokines_difference_FVC_p_at_fu.csv")

# Load necessary library
library(glmnet)
library(boot)

# Prepare the data
x <- as.matrix(scale(features))
y <- outcome$FVC_p

# Perform Ridge Regression with cross-validation to find the optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0)  # alpha = 0 for Ridge Regression
ridge_model <- glmnet(x, y, alpha = 0, lambda = cv_ridge$lambda.min)

# Function to fit ridge regression and extract coefficients
ridge_fit <- function(data, indices) {
  boot_x <- data[indices, -ncol(data)]
  boot_y <- data[indices, ncol(data)]
  model <- glmnet(boot_x, boot_y, alpha = 0, lambda = cv_ridge$lambda.min)
  return(as.vector(coef(model)))
}

# Combine features and outcome for bootstrapping
data <- data.frame(cbind(x, y))

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data, ridge_fit, R = 1000)

# Extract bootstrap estimates and calculate standard errors and confidence intervals
bootstrap_coefs <- t(sapply(1:ncol(data), function(i) boot::boot.ci(bootstrap_results, type = "perc", index = i)$percent))
coef_names <- rownames(coef(ridge_model))
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = as.vector(coef(ridge_model)),
  `2.5 %` = bootstrap_coefs[, 4],
  `97.5 %` = bootstrap_coefs[, 5]
)

print(bootstrap_summary)

#save 
write.csv(bootstrap_summary, file = "Results/Coef_summary_ridge_model_differnce_data_FVC_p.csv")

# Plot the coefficients with error bars

p <- ggplot(bootstrap_summary[-1, ], aes(x = Estimate, y = Coefficient)) + 
  geom_point() +
  geom_errorbar(aes(xmin = X2.5.., xmax = X97.5..), width = 0.2) +
  theme_minimal() +
  labs(title = "Ridge Regression Coefficients with 95% Confidence Intervals",
       x = "Estimate", 
       y="")+
  theme(axis.text.y = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        axis.title.x = element_text(size = 18, face = "bold"))

pdf(file = "Figures/Coef_summary_ridge_model_differnce_data_FVC_p.pdf", height = 8, width = 10)
p
dev.off()


# LOOCV for Ridge Regression
loocv_ridge <- train(y ~ ., data = cbind(scale(features), y), method = "glmnet", trControl = loocv_results,
                     tuneGrid = expand.grid(alpha = 0, lambda = cv_ridge$lambda.min))
print(loocv_ridge)

#extract model performance 
model_perf <- loocv_ridge$results

#save
write.csv(model_perf, file = "Results/ridge_model_performance_difference_data_FVC_p.csv")









#random forest on didference data 
library(dplyr)
library(randomForest)
library(ggplot2)
library(caret)

#choose counts 
counts <- difference_data
counts <- counts %>% dplyr::select(-c(sample, sample.x, sample.y))

#choose metadata 
temp_dat <- metadata %>% filter(sample_type=="clinic") %>% dplyr::select(FVC_lit)
counts$FVC_lit <- temp_dat$FVC_lit

# Convert FVC_lit to numeric
counts$FVC_lit <- as.numeric(counts$FVC_lit)

# Drop rows with missing FVC_lit values
counts <- drop_na(counts)

# Set the number of repetitions
num_repetitions <- 10  # Choose the number of repetitions you want

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  # Split data into training and test sets
  set.seed(1234 + i)
  train_index <- sample(nrow(counts), 0.6 * nrow(counts))
  train_data <- counts[train_index, ]
  test_data <- counts[-train_index, ]
  
  control <- trainControl(method = 'repeatedcv', number = 10, repeats = 2, 
                          search = 'random', savePredictions = TRUE, 
                          sampling = NULL, classProbs = TRUE)
  
  # Train the random forest model
  rf_model <- randomForest(FVC_lit ~ ., data = train_data, ntree = 500, 
                           trControl=control, 
                           preProcess=c("BoxCox"), scale=TRUE)
  
  # Get variable importance
  mdi <- data.frame(importance(rf_model))
  mdi <- mdi %>% arrange(desc(IncNodePurity))
  mdi$variable <- rownames(mdi)
  
  # Store results in the data frame
  results <- rbind(results, data.frame(Replication = i, mdi))
}

# Summarize variable importance across repetitions
importance_summary <- results %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(IncNodePurity), sd_importance = sd(IncNodePurity))

# Define threshold for importance (e.g., mean importance > 1)
threshold <- median(importance_summary$mean_importance)
confirmed_biomarkers <- importance_summary %>% filter(mean_importance > threshold)
rejected_biomarkers <- importance_summary %>% filter(mean_importance <= threshold)

# Define a threshold based on mean + 1 standard deviation
threshold <- mean(importance_summary$mean_importance) + sd(importance_summary$mean_importance)
confirmed_biomarkers <- importance_summary %>% filter(mean_importance > threshold)
rejected_biomarkers <- importance_summary %>% filter(mean_importance <= threshold)

# Define a percentile threshold (e.g., top 25%)
percentile <- 0.75
threshold <- quantile(importance_summary$mean_importance, probs = percentile)
confirmed_biomarkers <- importance_summary %>% filter(mean_importance > threshold)
rejected_biomarkers <- importance_summary %>% filter(mean_importance <= threshold)

# Print confirmed and rejected biomarkers
print("Confirmed Biomarkers:")
print(confirmed_biomarkers)
print("Rejected Biomarkers:")
print(rejected_biomarkers)

# Create a boxplot
pdf(file = "Figures/rf_importance_biomarkers_differences_FVC_lit.pdf", height = 10, width = 12)
ggplot(results, aes(x=IncNodePurity, y=reorder(variable, IncNodePurity))) +
  geom_boxplot(color="red") +
  geom_jitter(alpha=0.5)+
  guides(fill="none", size="none")+
  ylab("")+xlab("Importance")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 24, face = "bold"))
dev.off()













######## correlaiton wit DLCO % pred as outcome ########


#repeat for DLCO % as outcome 
####### correlation btween cytokines data and DLCO_pred########

#correltion matrix of baseline wtih DLCO_pred
temp <- inner_join(baseline_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(KCO))
features <- temp %>% dplyr::select(c(TNFa, IL6, MMP2, MMP7, MMP9, Angiopoietin2, Thrombomodulin, CA15_3, Hyaluronan))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(KCO)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


library(corrplot)


# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_baseline_with_KCO_p_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients

dev.off()

### compute p values 
data$outcome <- data$KCO
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
#save 
write.csv(correlation_results, file = "Results/correlation_cytokines_baseline_KCO_p_at_fu.csv")




############ repeeat for follow up data 
temp <- inner_join(clinic_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(KCO))
features <- temp %>% dplyr::select(c(TNFa, IL6, MMP2, MMP7, MMP9, Angiopoietin2, Thrombomodulin, CA15_3, Hyaluronan))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(KCO)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


library(corrplot)

# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_clinic_with_KCO_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7)      # Size of the correlation coefficients
dev.off()

### compute p values 
data$outcome <- data$KCO
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
write.csv(correlation_results, file = "Results/correlation_cytokines_clinic_KCO_p_at_fu.csv")





#### repeat for differnce data 
temp <- inner_join(difference_data, metadata, by="sample")
temp <- temp %>% dplyr::filter(!is.na(KCO))
features <- temp %>% dplyr::select(c(TNFa_diff, IL6_diff, MMP2_diff, MMP7_diff, MMP9_diff,
                                     Angiopoietin2_diff, Thrombomodulin_diff, CA15_3_diff, Hyaluronan_diff))
clinical_variables <- temp %>% dplyr::select(c(Age, time_to_fu, Gender, IHD, immunsupression_y_n, Vaccination.status, 
                                               crp_peak,ICU_y_n, Intubated, sao2_ra, sao2_ra_fu))
outcome<- temp %>% dplyr::select(KCO)

#first analyze cytokines only 
data <- cbind(features, outcome)

# Correlation matrix
cor_matrix <- cor(data, method = "spearman")
print(cor_matrix)


library(corrplot)


# Plot the correlation matrix
pdf(file = "Figures/correlation_plot_difference_with_KCO_p_at_FU.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients on the plot
         number.cex = 0.7, 
         sig.level = 0.05, p.mat = p.mat)      # Size of the correlation coefficients
dev.off()



cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


p.mat <- cor.mtest(data)



### compute p values 
data$outcome <- data$KCO
# Initialize an empty data frame to store correlation coefficients and p-values
correlation_results <- data.frame(cytokine = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each cytokine
for (cytokine in colnames(data)[-ncol(data)]) {
  # Calculate Spearman correlation coefficient and p-value
  correlation_result <- cor.test(data[[cytokine]], data$outcome, method = "spearman")
  
  # Store cytokine name, correlation coefficient, and p-value
  correlation_results <- rbind(correlation_results, data.frame(cytokine = cytokine, correlation = correlation_result$estimate, p_value = correlation_result$p.value))
}

# Print the results
print(correlation_results)
write.csv(correlation_results, file = "Results/correlation_cytokines_difference_KCO_p_at_fu.csv")





############## random forest Boruta for FVC lit as outcome ###########

# Install and load necessary packages
library(Boruta)
library(dplyr)
library(ggplot2)

# Prepare your data as before
counts <- difference_data
counts <- counts %>% dplyr::select(-c(sample, sample.x, sample.y))

temp_dat <- metadata %>% filter(sample_type == "clinic") %>% dplyr::select(FVC_lit)
counts$FVC_lit <- temp_dat$FVC_lit
counts$FVC_lit <- as.numeric(counts$FVC_lit)
counts <- drop_na(counts)
atributes <- colnames(counts)
atributes <- atributes[-10]

# Assuming 'counts' and 'FVC_lit' are already prepared as per your previous examples

# Set the number of repetitions and folds for cross-validation
num_repetitions <- 10
num_folds <- 5  # Number of folds for cross-validation

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  
  # Set seed for reproducibility
  set.seed(1234 + i)
  
  # Run cross-validation
  folds <- createFolds(counts$FVC_lit, k = num_folds, list = TRUE, returnTrain = TRUE)
  
  for (fold in 1:num_folds) {
    # Split data into train and test for this fold
    train_index <- unlist(folds[fold])
    train_data <- counts[train_index, ]
    test_data <- counts[-train_index, ]
    
    # Run Boruta feature selection
    boruta_result <- Boruta(FVC_lit ~ ., data = train_data, doTrace = 0,
                            randomSeed = 1234 + i, ntree = 500)  # Adjust ntree as needed
    
    # Get the decision on each feature
    final_boruta <- TentativeRoughFix(boruta_result)
    feature_decisions <- attStats(final_boruta)
    
    # Store results in the data frame
    feature_decisions$Replication <- i
    results <- rbind(results, feature_decisions)
  }
}

results

# Prepare data for visualization (example assuming 'results' structure)
results$attr <- colnames(train_data)[-which(colnames(train_data) == "FVC_lit")]
results$color <- ifelse(results$decision == "Confirmed", "blue2", "red")

# Determine the decision for each feature based on the majority decision
feature_decision_summary <- results %>%
  group_by(attr) %>%
  summarise(decision = ifelse(sum(decision == "Confirmed") > sum(decision == "Rejected"), "Confirmed", "Rejected"))

# Merge the decision summary back into the results
results <- results %>%
  left_join(feature_decision_summary, by = "attr") %>%
  mutate(color = ifelse(decision.y == "Confirmed", "blue2", "red"))


# Plotting the results
p <- ggplot(results, aes(x = medianImp, y = reorder(attr, medianImp))) +
  geom_boxplot(aes(color = color)) +
  geom_jitter(alpha = 0.5) +
  scale_color_manual(values = c("blue", "red"))+
  guides(fill = "none", size = "none", color="none") +
  ylab("") +
  xlab("Importance") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 20, face = "bold"))

p

#save 
pdf(file = "Figures/RF_bouratta_diff_data_FVC_lit_outcome.pdf", height = 10, width = 14)
p
dev.off()




############# repeat random forest for  follow up data ######
# Install and load necessary packages
library(Boruta)
library(dplyr)
library(ggplot2)

# Prepare your data as before
counts <- followup_data
counts <- counts %>% dplyr::select(-c(sample_type, sample))

temp_dat <- metadata %>% filter(sample_type == "clinic") %>% dplyr::select(FVC_lit)
counts$FVC_lit <- temp_dat$FVC_lit
counts$FVC_lit <- as.numeric(counts$FVC_lit)
counts <- drop_na(counts)
atributes <- colnames(counts)
atributes <- atributes[-10]

# Assuming 'counts' and 'FVC_lit' are already prepared as per your previous examples

# Set the number of repetitions and folds for cross-validation
num_repetitions <- 10
num_folds <- 5  # Number of folds for cross-validation

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  
  # Set seed for reproducibility
  set.seed(1234 + i)
  
  # Run cross-validation
  folds <- createFolds(counts$FVC_lit, k = num_folds, list = TRUE, returnTrain = TRUE)
  
  for (fold in 1:num_folds) {
    # Split data into train and test for this fold
    train_index <- unlist(folds[fold])
    train_data <- counts[train_index, ]
    test_data <- counts[-train_index, ]
    
    # Run Boruta feature selection
    boruta_result <- Boruta(FVC_lit ~ ., data = train_data, doTrace = 0,
                            randomSeed = 1234 + i, ntree = 500)  # Adjust ntree as needed
    
    # Get the decision on each feature
    final_boruta <- TentativeRoughFix(boruta_result)
    feature_decisions <- attStats(final_boruta)
    
    # Store results in the data frame
    feature_decisions$Replication <- i
    results <- rbind(results, feature_decisions)
  }
}

# Prepare data for visualization (example assuming 'results' structure)
results$attr <- colnames(train_data)[-which(colnames(train_data) == "FVC_lit")]
results$color <- ifelse(results$decision == "Confirmed", "blue2", "red")

# Determine the decision for each feature based on the majority decision
feature_decision_summary <- results %>%
  group_by(attr) %>%
  summarise(decision = ifelse(sum(decision == "Confirmed") > sum(decision == "Rejected"), "Confirmed", "Rejected"))

# Merge the decision summary back into the results
results <- results %>%
  left_join(feature_decision_summary, by = "attr") %>%
  mutate(color = ifelse(decision.y == "Confirmed", "blue2", "red"))

results$decision.y

# Plotting the results
p <- ggplot(results, aes(x = medianImp, y = reorder(attr, medianImp))) +
  geom_boxplot(aes(color = color)) +
  geom_jitter(alpha = 0.5) +
  scale_color_manual(values = "red")+
  guides(fill = "none", size = "none", color="none") +
  ylab("") +
  xlab("Importance") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 20, face = "bold"))

p

#save 
pdf(file = "Figures/RF_bouratta_follow_up_data_FVC_lit_outcome.pdf", height = 10, width = 14)
p
dev.off()










############# repeat random forest for  baseline data ######
# Install and load necessary packages
library(Boruta)
library(dplyr)
library(ggplot2)

# Prepare your data as before
counts <- baseline_data
#counts <- counts %>% dplyr::select(-c(sample_type, sample))

temp_dat <- metadata %>% filter(sample_type == "clinic") %>% dplyr::select(FVC_lit)
row_names <- rownames(temp_dat)
last_two_digits <- substr(row_names, nchar(row_names) - 1, nchar(row_names))

temp_dat$subj_id <- paste0("subj_", paste0(last_two_digits))
counts <- inner_join(counts, temp_dat, by="subj_id")
counts <- counts %>% dplyr::select(-sample, -subj_id)
counts$FVC_lit <- as.numeric(counts$FVC_lit)
counts <- drop_na(counts)
atributes <- colnames(counts)
atributes <- atributes[-10]

# Assuming 'counts' and 'FVC_lit' are already prepared as per your previous examples

# Set the number of repetitions and folds for cross-validation
num_repetitions <- 10
num_folds <- 5  # Number of folds for cross-validation

# Create an empty data frame to store results
results <- data.frame()

# Repeat the process for the specified number of times
for (i in 1:num_repetitions) {
  
  # Set seed for reproducibility
  set.seed(1234 + i)
  
  # Run cross-validation
  folds <- createFolds(counts$FVC_lit, k = num_folds, list = TRUE, returnTrain = TRUE)
  
  for (fold in 1:num_folds) {
    # Split data into train and test for this fold
    train_index <- unlist(folds[fold])
    train_data <- counts[train_index, ]
    test_data <- counts[-train_index, ]
    
    # Run Boruta feature selection
    boruta_result <- Boruta(FVC_lit ~ ., data = train_data, doTrace = 0,
                            randomSeed = 1234 + i, ntree = 500)  # Adjust ntree as needed
    
    # Get the decision on each feature
    final_boruta <- TentativeRoughFix(boruta_result)
    feature_decisions <- attStats(final_boruta)
    
    # Store results in the data frame
    feature_decisions$Replication <- i
    results <- rbind(results, feature_decisions)
  }
}

# Prepare data for visualization (example assuming 'results' structure)
results$attr <- colnames(train_data)[-which(colnames(train_data) == "FVC_lit")]
results$color <- ifelse(results$decision == "Confirmed", "blue2", "red")

# Determine the decision for each feature based on the majority decision
feature_decision_summary <- results %>%
  group_by(attr) %>%
  summarise(decision = ifelse(sum(decision == "Confirmed") > sum(decision == "Rejected"), "Confirmed", "Rejected"))

# Merge the decision summary back into the results
results <- results %>%
  left_join(feature_decision_summary, by = "attr") %>%
  mutate(color = ifelse(decision.y == "Confirmed", "blue2", "red"))

results$decision.y
which(results$decision.y=="Confirmed")
# Plotting the results
p <- ggplot(results, aes(x = medianImp, y = reorder(attr, medianImp))) +
  geom_boxplot(aes(color = color)) +
  geom_jitter(alpha = 0.5) +
  scale_color_manual(values = "red")+
  guides(fill = "none", size = "none", color="none") +
  ylab("") +
  xlab("Importance") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 20, face = "bold"))

p

#save 
pdf(file = "Figures/RF_bouratta_baseline_data_FVC_lit_outcome.pdf", height = 10, width = 14)
p
dev.off()
























######## comparing baseline blood tests #########


############# according to disease severity 


temp_dat <- metadata

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type=="baseline")
#get subject IDs

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

temp_dat$disease_severity_two <- factor(temp_dat$disease_severity_two, levels = c("mild_mod", "severe"))

#example plot
ggplot(data = temp_dat, aes(x=temp_dat$disease_severity_two, y=DDM_peak, color=temp_dat$disease_severity_two))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
  xlab("")+ylab("KL-6")+
  stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")



col_names <- colnames(temp_dat)[55:62]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$disease_severity_two, y=i, color=as.factor(temp_dat$disease_severity_two)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Mild-Moderate", "Severe"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("mild_mod", "severe")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_disease_severity.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(c(col_names, disease_severity_two)) %>% 
  group_by(disease_severity_two) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

#save
write.csv(summary_stats, file = paste0("Results/blood_tests_comparison_baseline_by_", paste0("disease_severity", paste0(".csv"))))





############# according to intubation 


temp_dat <- metadata

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type=="baseline")
#get subject IDs

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

temp_dat$Intubated <- factor(temp_dat$Intubated, levels = c("N", "Y"))



col_names <- colnames(temp_dat)[55:62]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$Intubated, y=i, color=as.factor(temp_dat$Intubated)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Non-Intubated", "Intubated"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_Intubated.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(c(col_names, Intubated)) %>% 
  group_by(Intubated) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

#save
write.csv(summary_stats, file = paste0("Results/blood_tests_comparison_baseline_by_", paste0("Intubated", paste0(".csv"))))











############# according to ICU 


temp_dat <- metadata

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type=="baseline")
#get subject IDs

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

temp_dat$ICU_y_n <- factor(temp_dat$ICU_y_n, levels = c("N", "Y"))



col_names <- colnames(temp_dat)[55:62]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$ICU_y_n, y=i, color=as.factor(temp_dat$ICU_y_n)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Non-ICU", "ICU"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_ICU_y_n.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(c(col_names, ICU_y_n)) %>% 
  group_by(ICU_y_n) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

#save
write.csv(summary_stats, file = paste0("Results/blood_tests_comparison_baseline_by_", paste0("ICU_y_n", paste0(".csv"))))










############# according to death 


temp_dat <- metadata

#filter for baseline and follow up samples only 
temp_dat <- temp_dat %>% filter(sample_type=="baseline")
#get subject IDs

temp_dat <- temp_dat %>% filter(sample_type!="clinic")

temp_dat$Death <- factor(temp_dat$Death, levels = c("N", "Y"))



col_names <- colnames(temp_dat)[55:62]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp_dat, aes_string(x=temp_dat$Death, y=i, color=as.factor(temp_dat$Death)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Alive", "Dead"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_Death.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}


#export table if we are going to add it as supp table 
temp_dat

summary_stats <- temp_dat %>%
  dplyr::select(c(col_names, Death)) %>% 
  group_by(Death) %>%
  summarize_all(funs(Median = median(., na.rm = TRUE),
                     Q1 = quantile(., 0.25, na.rm = TRUE),
                     Q3 = quantile(., 0.75, na.rm = TRUE)))

# Print the summary statistics
print(summary_stats)

#save
write.csv(summary_stats, file = paste0("Results/blood_tests_comparison_baseline_by_", paste0("Death", paste0(".csv"))))





















####### repeat correaltion plots and save them into one folder for baseline and follow up #########

####### baseline 
baseline_data$sample
temp <- metadata %>% dplyr::select(FVC_lit, FVC_p, FEV1_lit, FEV1_p, KCO, sample)

temp <- inner_join(baseline_data, temp, by="sample")
temp$FVC_lit <- as.numeric(temp$FVC_lit)
temp$FVC_p <- as.numeric(temp$FVC_p)
temp$FEV1_lit <- as.numeric(temp$FEV1_lit)
temp$FEV1_p <- as.numeric(temp$FEV1_p)
temp$KCO <- as.numeric(temp$KCO)

temp$IL6 <- ifelse(temp$IL6 > 500, 500, temp$IL6)
temp$CA15_3 <- ifelse(temp$CA15_3 > 1000, 1000, temp$CA15_3)

col_names <- colnames(temp)[1:9]

# FVC lit 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$FVC_lit))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/correlation_plots/post_covid_", paste0(i, paste0("_correl_FVC_lit__baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}



#KCO %pred
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$KCO))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/correlation_plots/post_covid_", paste0(i, paste0("_correl_KCO_per__baseline.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}





########follow up 
followup_data$sample
temp <- metadata %>% dplyr::select(FVC_lit, FVC_p, FEV1_lit, FEV1_p, KCO, sample)

temp <- inner_join(followup_data, temp, by="sample")
temp$FVC_lit <- as.numeric(temp$FVC_lit)
temp$FVC_p <- as.numeric(temp$FVC_p)
temp$FEV1_lit <- as.numeric(temp$FEV1_lit)
temp$FEV1_p <- as.numeric(temp$FEV1_p)
temp$KCO <- as.numeric(temp$KCO)

temp$IL6 <- ifelse(temp$IL6 > 500, 500, temp$IL6)
temp$CA15_3 <- ifelse(temp$CA15_3 > 1000, 1000, temp$CA15_3)

col_names <- colnames(temp)[2:10]

# FVC lit 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$FVC_lit))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/correlation_plots/post_covid_", paste0(i, paste0("_correl_FVC_lit__follow_up.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}



#KCO %pred
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=i, y=temp$KCO))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman")+
    theme_classic()  
  
  pdf_output <- paste0("Figures/correlation_plots/post_covid_", paste0(i, paste0("_correl_KCO_per_follow_up.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}












######### check differnces according to normal abnormal CT #######

########baseline 

baseline_data$sample
temp <- metadata %>% dplyr::select(abnormal_CT_fu_y_n, sample)

temp <- inner_join(baseline_data, temp, by="sample")

temp <- temp %>% filter(!is.na(abnormal_CT_fu_y_n))
temp$abnormal_CT_fu_y_n <- factor(temp$abnormal_CT_fu_y_n, levels = c("N", "Y"))

temp$IL6 <- ifelse(temp$IL6 > 500, 500, temp$IL6)
temp$CA15_3 <- ifelse(temp$CA15_3 > 1000, 1000, temp$CA15_3)

ggplot(data = temp, aes(x=abnormal_CT_fu_y_n, y=IL6, color=as.factor(abnormal_CT_fu_y_n)))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=c("green3", "blue3"))+
  scale_x_discrete(labels=c("Normal", "Abnormal"))+
  xlab("")+ylab(i)+
  stat_compare_means(comparisons = list(c("N", "Y")))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"))+
  guides(color="none")

col_names <- colnames(temp)[1:9]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=temp$abnormal_CT_fu_y_n, y=i, color=as.factor(temp$abnormal_CT_fu_y_n)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Normal", "Abnormal"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_baseline_levels_by_abnormal_Ct.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}



###########follow up 

followup_data$sample
temp <- metadata %>% dplyr::select(abnormal_CT_fu_y_n, sample)

temp <- inner_join(followup_data, temp, by="sample")

temp <- temp %>% filter(!is.na(abnormal_CT_fu_y_n))
temp$abnormal_CT_fu_y_n <- factor(temp$abnormal_CT_fu_y_n, levels = c("N", "Y"))

temp$IL6 <- ifelse(temp$IL6 > 100, 100, temp$IL6)
temp$CA15_3 <- ifelse(temp$CA15_3 > 1000, 1000, temp$CA15_3)



col_names <- colnames(temp)[2:10]

#loop over al variables 
for (i in col_names) {
  p <- ggplot(data = temp, aes_string(x=temp$abnormal_CT_fu_y_n, y=i, color=as.factor(temp$abnormal_CT_fu_y_n)))+
    geom_boxplot(outlier.shape = NA)+
    geom_point()+
    scale_color_manual(values=c("green3", "blue3"))+
    scale_x_discrete(labels=c("Normal", "Abnormal"))+
    xlab("")+ylab(i)+
    stat_compare_means(comparisons = list(c("N", "Y")))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"), 
          axis.title.y = element_text(size = 20, face = "bold"))+
    guides(color="none")
  
  pdf_output <- paste0("Figures/post_covid_", paste0(i, paste0("_follow_up_levels_by_abnormal_Ct.pdf")))
  pdf(file=pdf_output,height = 8, width = 6)
  show(p)
  dev.off()
}











save.image("post_covid.RData")
