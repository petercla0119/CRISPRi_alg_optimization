library(tidyverse)
library(ggrepel)
library(svglite)
library(plotly)
source("05_scripts/01_data_analysis/00_experiments/clean_preprocess_module.R")

# Import sgRNA summary file from MAGeCK output
rra.High.Low.STMN2.sgrna_summary <- read.delim("~/Desktop/Projects/dual_guide_optimization/01_sorting-based_screens/stmn2/analysis_method/mageck/initial_mageck_09122024/stmn2_unpaired/results/test/High_vs_Low.sgrna_summary.txt")

rra.High.Low.STMN2.sgrna_summary <- separate_sgrna_summary(rra.High.Low.STMN2.sgrna_summary)
# Convert columns to numeric
rra.High.Low.STMN2.sgrna_summary[, 3:ncol(rra.High.Low.STMN2.sgrna_summary)] <- lapply(rra.High.Low.STMN2.sgrna_summary[, 3:ncol(rra.High.Low.STMN2.sgrna_summary)], as.numeric)

##transcript plot - need to add noise modeling too
rra.High.Low.STMN2.sgrna_summary %>% 
  ggplot(aes(control_mean,treat_mean, label=Gene))+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
               filter(Gene!="negative_control",
                      (control_mean>=3000|
                        treat_mean>=3000)),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
               filter(Gene=="negative_control",
                      (control_mean>=3000|
                         treat_mean>=3000)),
             color="grey",
             alpha=0.7)+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
               filter(Gene=="STMN2"|
                        Gene=="TARDBP"),
             color="purple",
             alpha=0.7)+
  #geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$id),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$id),
  #           color="green",
  #           alpha=0.7)+
  geom_label_repel(data = rra.High.Low.STMN2.sgrna_summary %>%
                     dplyr::filter((control_mean>=3000 | treat_mean>=3000), abs(LFC)>4)
                              ,
                   box.padding = 0.5,
                   segment.color="black",
                   alpha=0.6,
                   max.overlaps = Inf
  ) +
  scale_x_log10()+
  scale_y_log10()+
  coord_fixed(1)+
  labs(
    title = "STMN2 Expression Screen",
    x="KD decreases STMN2 (median transcript count)",
    y="KD increases STMN2 (median transcript count)",
    caption = 
      "Grey: NT guides"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

