library(tidyverse)
library(ggrepel)
library(ggdark)
library(plotly)

recom.gene_summary <- read.delim("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/files_from_others/recombinant data/rra-High-Low-STMN2.gene_summary.txt")
recom.sgrna_summary <- read.delim("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/files_from_others/recombinant data/rra-High-Low-STMN2.sgrna_summary.txt")
recom.sgrna_summary$control_1 <- recom.sgrna_summary$control_count %>% str_split_i('/', 1) %>% as.numeric()
recom.sgrna_summary$control_2 <- recom.sgrna_summary$control_count %>% str_split_i('/', 2) %>% as.numeric() 
recom.sgrna_summary$control_3 <- recom.sgrna_summary$control_count %>% str_split_i('/', 3) %>% as.numeric() 
recom.sgrna_summary$treatment_1 <- recom.sgrna_summary$treatment_count %>% str_split_i('/', 1) %>% as.numeric()
recom.sgrna_summary$treatment_2 <- recom.sgrna_summary$treatment_count %>% str_split_i('/', 2) %>% as.numeric() 
recom.sgrna_summary$treatment_3 <- recom.sgrna_summary$treatment_count %>% str_split_i('/', 3) %>% as.numeric() 
recom.sgrna_summary <- recom.sgrna_summary %>% mutate(control_avg = (control_1 + control_2 + control_3)/3, treatment_avg = (treatment_1 + treatment_2 + treatment_3)/3)
recom_control_sd <- recom.sgrna_summary[c('control_1', 'control_2', 'control_3')] %>% apply(1, sd)
recom.sgrna_summary <- cbind(recom.sgrna_summary, control_sd = recom_control_sd)
recom_treatment_sd <- recom.sgrna_summary[c('treatment_1', 'treatment_2', 'treatment_3')] %>% apply(1, sd)
recom.sgrna_summary <- cbind(recom.sgrna_summary, treatment_sd = recom_treatment_sd)

hits.sgrna_summary <- read.delim("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_screenfiles/CRISPRi_alg_optimization/data_processed/rra.sgrna_summary.txt")
hits.sgrna_summary$control_1 <- hits.sgrna_summary$control_count %>% str_split_i('/', 1) %>% as.numeric()
hits.sgrna_summary$control_2 <- hits.sgrna_summary$control_count %>% str_split_i('/', 2) %>% as.numeric() 
hits.sgrna_summary$control_3 <- hits.sgrna_summary$control_count %>% str_split_i('/', 3) %>% as.numeric() 
hits.sgrna_summary$treatment_1 <- hits.sgrna_summary$treatment_count %>% str_split_i('/', 1) %>% as.numeric()
hits.sgrna_summary$treatment_2 <- hits.sgrna_summary$treatment_count %>% str_split_i('/', 2) %>% as.numeric() 
hits.sgrna_summary$treatment_3 <- hits.sgrna_summary$treatment_count %>% str_split_i('/', 3) %>% as.numeric() 
hits.sgrna_summary <- hits.sgrna_summary %>% mutate(control_avg = (control_1 + control_2 + control_3)/3, treatment_avg = (treatment_1 + treatment_2 + treatment_3)/3)
hits_control_sd <- hits.sgrna_summary[c('control_1', 'control_2', 'control_3')] %>% apply(1, sd)
hits.sgrna_summary <- cbind(hits.sgrna_summary, control_sd = hits_control_sd)
hits_treatment_sd <- hits.sgrna_summary[c('treatment_1', 'treatment_2', 'treatment_3')] %>% apply(1, sd)
hits.sgrna_summary <- cbind(hits.sgrna_summary, treatment_sd = hits_treatment_sd)

#sgrna compare lfc
sgrna_compare <- hits.sgrna_summary %>% 
  merge(recom.sgrna_summary, by="sgrna",#all=T
  )

#weight hits 2 times to recombinants
sgrna_compare <- sgrna_compare %>% mutate(combined_control_avg = (2*control_avg.x + control_avg.y)/3, combined_treatment_avg = (2*treatment_avg.x + treatment_avg.y)/3)
sgrna_compare <- sgrna_compare %>% mutate(avg_lfc.x = log2(treatment_avg.x / control_avg.x), avg_lfc.y = log2(treatment_avg.y / control_avg.y), avg_lfc_combined = log2(combined_treatment_avg / combined_control_avg))

#add columns for noise modeling
sgrna_compare <- sgrna_compare %>% mutate(recombinant_cell_count = control_avg.y+treatment_avg.y,hits_cell_count = control_avg.x + treatment_avg.x )
sgrna_compare <- sgrna_compare %>% mutate(cv.x = sqrt(control_sd.x^2 + treatment_sd.x^2)/(control_avg.x + treatment_avg.x))
sgrna_compare <- sgrna_compare %>% mutate(cv.y = sqrt(control_sd.y^2 + treatment_sd.y^2)/(control_avg.y + treatment_avg.y))

#hits: noise modeling
linear_model<-lm(formula = log2(cv.x) ~ hits_cell_count,
                 data = sgrna_compare)
sgrna_compare %>%
  ggplot(aes(hits_cell_count, cv.x))+
  geom_point(color = 'green', alpha = 0.3, size = 0.5)+
  labs(title = "Hits: Count vs CV, combined high and low groups", x = "Count", y = "Coefficient of Variation")+
  stat_function(fun=function(x) {2^((x*-0.000085)-2.143)}, color = 'orange')+
  geom_smooth( se = TRUE, color = 'blue', alpha = 0.4)+
  xlim(c(0,20000))+
  ylim(c(0,0.75))+
  theme_classic()
#hits: exponential regression
sgrna_compare %>%
  ggplot(aes(hits_cell_count, log2(cv.x)))+
  geom_point(color = 'green', alpha = 0.1, size = 0.5)+
  labs(title = "Hits: Count vs Log CV, combined high and low groups", x = "Count", y = "Log Coefficient of Variation")+
  geom_smooth(method = 'lm', se = TRUE, color = 'blue')
#hits: noise modeling with 1/x
sgrna_compare %>%
  ggplot(aes(hits_cell_count, cv.x))+
  geom_point(color = 'green', alpha = 0.3)+
  labs(title = "Hits: Count vs CV, combined high and low groups", x = "Count", y = "Coefficient of Variation")+
  stat_function(fun=function(x) {1/((0.0004722*x)+4.873)}, color = 'orange')+
  xlim(c(0,20000))+
  ylim(c(0,0.75))
sgrna_compare %>%
  ggplot(aes(hits_cell_count, 1/(cv.x)))+
  geom_point(color = 'green', alpha = 0.1, size = 0.2)+
  labs(title = "Hits: Count vs Inverse CV, combined high and low groups", x = "Count", y = "Inverse Coefficient of Variation")+
  geom_smooth(method = 'lm', se = TRUE, color = 'blue')

inverse_linear_model<-lm(formula = 1/cv.x ~ (hits_cell_count),
                         data = sgrna_compare)

#recombinants: noise modeling
recom_linear_model<-lm(formula = log2(cv.y) ~ recombinant_cell_count,
                       data = sgrna_compare)
sgrna_compare %>%
  ggplot(aes(recombinant_cell_count, cv.y))+
  geom_point(color = 'red', alpha = 0.1, size = 0.3)+
  labs(title = "Recombinants: Count vs CV, combined high and low groups", x = "Count", y = "Coefficient of Variation")+
  stat_function(fun=function(x) {2^((x*-0.0001488)-1.785)}, color = 'orange')+
  geom_smooth( se = TRUE, color = 'blue', alpha = 0.4)+
  xlim(c(0,20000))+
  ylim(c(0,0.75))+
  theme_classic()
sgrna_compare %>%
  ggplot(aes(recombinant_cell_count, log2(cv.y)))+
  geom_point(color = 'red', alpha = 0.1)+
  labs(title = "Recombinants: Count vs Log CV, combined high and low groups", x = "Count", y = "Log Coefficient of Variation")+
  geom_smooth(method = 'lm', se = TRUE, color = 'blue')