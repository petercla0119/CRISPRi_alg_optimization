library(tidyverse)
library(ggrepel)
library(ggdark)
library(plotly)

recom.gene_summary <- read.delim(
  "/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/files_from_others/recombinant data/rra-High-Low-STMN2.gene_summary.txt"
)
recom.sgrna_summary <- read.delim(
  "/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/files_from_others/recombinant data/rra-High-Low-STMN2.sgrna_summary.txt"
)
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
hits.sgrna_summary <- hits.sgrna_summary %>%
                        mutate(control_avg = (control_1 + control_2 + control_3)/3, 
                               treatment_avg = (treatment_1 + treatment_2 + treatment_3)/3
                               )
hits_control_sd <- hits.sgrna_summary[c('control_1', 'control_2', 'control_3')] %>% apply(1, sd)
hits.sgrna_summary <- cbind(hits.sgrna_summary, control_sd = hits_control_sd)
hits_treatment_sd <- hits.sgrna_summary[c('treatment_1', 'treatment_2', 'treatment_3')] %>% apply(1, sd)
hits.sgrna_summary <- cbind(hits.sgrna_summary, treatment_sd = hits_treatment_sd)

hits_missing <- hits.sgrna_summary %>% 
                  filter(control_count == '0/0/0' & treatment_count == '0/0/0')
recom_present <- recom.sgrna_summary %>% 
                    filter(control_count != '0/0/0' | treatment_count != '0/0/0')
hits_missing_genes <- hits_missing$Gene
recom_present_genes <- recom_present$Gene
unique_recom_genes <- intersect(hits_missing_genes, recom_present_genes)
unique_recom <- recom.sgrna_summary %>% 
                  filter(Gene %in% unique_recom_genes)

### Unique recombinants
unique_recom %>%
  ggplot(aes(control_mean, treat_mean, label=Gene)) +
  geom_point() +
  geom_point(data = unique_recom %>% 
               filter(Gene=="negative_control"), color="darkgrey", size=2, alpha=0.8) +
  #geom_text(hjust=0, vjust=0)+
  geom_point(data=unique_recom %>% filter(Gene=="OPTN"|Gene=="VCP"), color="red",size=2) +
  
  geom_label_repel(data=unique_recom %>% filter(control_mean > 100 | treat_mean > 100), box.padding = 0.05, max.overlaps = 15) +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = "Recombinant only STMN2 screen transcript count",
       x= "KD decreases STMN2 (median count)",
       y= "KD increases STMN2 (median count)",
       caption = "3 replicate whole genome CRISPRi KD screen in i3Neurons")+
  dark_theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
#sgrna compare lfc
sgrna_compare <- hits.sgrna_summary %>% 
                  merge(recom.sgrna_summary, by="sgrna",#all=TRUE
                        )

#weight hits 2 times to recombinants
sgrna_compare <- sgrna_compare %>%
                  mutate(combined_control_avg = (2 * control_avg.x + control_avg.y)/3, 
                         combined_treatment_avg = (2 * treatment_avg.x + treatment_avg.y)/3
                         )
sgrna_compare <- sgrna_compare %>% 
                  mutate(avg_lfc.x = log2(treatment_avg.x / control_avg.x), 
                         avg_lfc.y = log2(treatment_avg.y / control_avg.y), 
                         avg_lfc_combined = log2(combined_treatment_avg / combined_control_avg)
                         )
sgrna_compare <- sgrna_compare %>% 
                  mutate(recombinant_cell_count = control_avg.y + treatment_avg.y,
                         hits_cell_count = control_avg.x + treatment_avg.x
                         )

#adjusting minimum cell count
sgrna_compare_reduced <- sgrna_compare %>% 
                          filter((control_avg.x+treatment_avg.x) > 500 & (control_avg.y+treatment_avg.y) > 500)





# sgrna_compare_reduced %>%
#   ggplot(aes(log2(treatment_avg.x / control_avg.x), log2(treatment_avg.y / control_avg.y), label=Gene.x))+
#   labs(title = "Log 2 High / Low STMN2 Count Ratio in Hits vs Recombinants",
#        x= "Hits",
#        y= "Recombinants",
#        caption = "3 replicate whole genome CRISPRi KD screen in i3Neurons")+
#   theme(plot.title = element_text(hjust = 0.5))+
#   stat_bin2d(bins=50)+
#   scale_fill_gradient(low='lightblue', high = 'red', limits = c(0,100))+
#   geom_abline(intercept=0,slope=1,color='green')+
#   geom_smooth(method = 'lm', se = TRUE)

sgrna_compare_reduced %>%
  ggplot(aes(
    x = log2(treatment_avg.x / control_avg.x),
    y = log2(treatment_avg.y / control_avg.y)
  )) +
  geom_bin2d(bins = 50) +  # Binning the data
  scale_fill_gradient(low = 'lightblue', high = 'red', limits = c(0, 100)) +
  geom_abline(intercept = 0, slope = 1, color = 'green') +
  geom_smooth(method = 'lm', se = FALSE, color = 'black') +  # Linear model without confidence interval
  labs(
    title = "Log2 High / Low STMN2 Count Ratio in Hits vs Recombinants",
    x = "Hits (Log2 Ratio)",
    y = "Recombinants (Log2 Ratio)",
    caption = "3 replicate whole genome CRISPRi KD screen in i3Neurons"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.caption = element_text(size = 8, hjust = 0.5)
  ) +
  geom_text(
    aes(label = Gene.x),
    size = 3,
    check_overlap = TRUE,
    data = subset(sgrna_compare_reduced, !is.na(Gene.x))  # Ensure that labels are only plotted for valid data
  )


sgrna_compare_reduced %>%
  ggplot(aes(avg_lfc.x, avg_lfc_combined, label=Gene.x, color = recombinant_cell_count)) +
  geom_point(alpha=0.4) +
  #geom_point(data=sgrna_compare_reduced %>% filter(Gene.x=="STMN2"|Gene.x=="TARDBP"), color="blue",size=2, alpha =0.5)+
  #geom_label_repel(data = sgrna_compare_reduced %>% filter (Gene.x=="STMN2"|Gene.x=="TARDBP"))+
  geom_label_repel(data = sgrna_compare_reduced %>% filter (abs(avg_lfc.x - avg_lfc_combined) > 2 & recombinant_cell_count > 1500) )+
  scale_color_gradient(low='lightblue',high='red') +
  #scale_y_log10()+
  #scale_x_log10()+
  labs(title = "Log 2 High / Low STMN2 Count Ratio in Hits vs Combined Data",
       x= "Hits",
       y= "Hits + Recombinants",
       caption = "3 replicate whole genome CRISPRi KD screen in i3Neurons") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_abline(intercept=0,slope=1,color='green') +
  geom_smooth(method = 'lm', se = TRUE)
