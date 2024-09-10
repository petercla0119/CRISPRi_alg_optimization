library(tidyverse)
library(ggrepel)
rra.High.Low.STMN2.sgrna_summary <- read.delim("~/Desktop/Projects/FACS_projects/STMN2/STMN2_mscarlet_221205/input/rra-High-Low-STMN2.sgrna_summary.txt")

rra.High.Low.STMN2.sgrna_summary<- rra.High.Low.STMN2.sgrna_summary%>%
  arrange(LFC)
rra.High.Low.STMN2.sgrna_summary<-rra.High.Low.STMN2.sgrna_summary %>% 
  mutate(
    zscore = scale(LFC),
    row = as.numeric(rownames(rra.High.Low.STMN2.sgrna_summary))
  )

rra.High.Low.STMN2.sgrna_summary %>%
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary %>% 
               filter(Gene=="TARDBP"|
                        Gene=="STMN2"|
                        Gene=="RANBP1"),
             color="purple",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
  geom_label_repel(data = rra.High.Low.STMN2.sgrna_summary %>% 
                     filter(Gene=="TARDBP"|
                              Gene=="STMN2"|
                              Gene=="RANBP1"),
                   fill="purple"
                   ,
                   box.padding = 0.5,
                   segment.color="black",
                   alpha=0.6,
                   max.overlaps = Inf)+
  labs(
    title = "STMN2-mscarlet rank plot",
    caption = 
      "Grey: NT guides
    Purple: genes of interest"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

rra.High.Low.STMN2.sgrna_summary_filtered<- rra.High.Low.STMN2.sgrna_summary%>%
  filter(control_mean>=3000|
           treat_mean>=3000) %>% 
  arrange(LFC)
rra.High.Low.STMN2.sgrna_summary_filtered<-rra.High.Low.STMN2.sgrna_summary_filtered %>% 
  mutate(
    zscore = scale(LFC),
    row = as.numeric(rownames(rra.High.Low.STMN2.sgrna_summary_filtered))
  )


rra.High.Low.STMN2.sgrna_summary_filtered %>%
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary_filtered %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary_filtered %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  geom_point(data = rra.High.Low.STMN2.sgrna_summary_filtered %>% 
               filter(Gene=="TARDBP"|
                        Gene=="STMN2"
                      #  Gene=="RANBP1"
               ),
             color="purple",
             alpha=0.7)+
  #geom_point(data = rra.High.Low.STMN2.sgrna_summary_filtered %>% 
  #             filter(Gene%in%mitochondrion$id
  #             ),
  #           color="red",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
#geom_point(data = both %>% 
#             filter((control_mean>=3000|
#                       treat_mean>=3000),
#                    Gene%in%translation$Gene),
#           color="green",
#           alpha=0.7)+
geom_label_repel(data = rra.High.Low.STMN2.sgrna_summary_filtered %>% 
                   filter(zscore < -5
                          # Gene=="RANBP1"
                   ),
                 box.padding = 0.5,
                 segment.color="black",
                 alpha=0.6,
                 max.overlaps = Inf)+
  labs(
    title = "STMN2-mScarlet rank plot filtered at 3000",
    caption = 
      "Grey: NT guides
    Purple: Genes of interest"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))




### Rank stats -#

rra.High.Low.STMN2.sgrna_summary_filtered %>% 
  filter(Gene!="negative_control") %>% 
  summary

rra.High.Low.STMN2.sgrna_summary_filtered %>% 
  filter(Gene=="negative_control") %>% 
  summarise(
    sd = sd(zscore),
    mean = mean(zscore)
  )

rra.High.Low.STMN2.sgrna_summary_filtered %>% 
  filter(Gene=="negative_control") %>% 
  summarise(
    sd = sd(zscore),
    mean = mean(zscore)
  )
up = -0.1179674+ 4*0.5991452
down = -0.1179674- 4*0.5991452

rra.High.Low.STMN2.sgrna_summary_filtered %>% 
  filter(Gene!="negative_control",
         zscore > up |
           zscore < down) %>% 
  summarise(
    n=n()
  )

