library(tidyverse)
library(ggplot2)

# Have not been able to run as of 9/18/24
STMN2 <- read.delim("~/Desktop/Projects/dual_guide_optimization/01_sorting-based_screens/stmn2/analysis_method/mageck/initial_mageck_09122024/stmn2_unpaired/results/test/High_vs_Low.gene_summary.txt")

negative<- STMN2 %>% filter(neg.lfc<0) %>% 
  select(id,lfc=neg.lfc, score = neg.score,pvalue =neg.p.value,fdr = neg.fdr)

positive<- STMN2 %>% filter(pos.lfc>0) %>% 
  select(id,lfc=pos.lfc, score = pos.score,pvalue =pos.p.value,fdr = pos.fdr)

both<-rbind(negative, positive) %>% 
  mutate(
    qvalue = -log10(pvalue)
  )
  
# test <-hgtRes1_up@result
# hgt_combined <- rbind(hgtRes1@result, hgtRes1_up@result)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(shiny)
library(ggiraph)
library(shinyWidgets)


###this is the custom UI
###Omic x is proteomics and y is transcriptomics
ui <- fluidPage(
  tags$h1("Volcano plot of all targets"),
  fluidRow(
    column(width = 3,
           numericInput("p", label="-log10(padj)", 
                        value = "8"),
           numericInput("lfc_z", label="LFC", 
                        value = "5"),
           multiInput(inputId = "Target",
                      label = "Select Target:",
                      choices = both$id %>% unique()),
           #multiInput(inputId = "GO",
            #          label = "Select GO Term:",
             #         choices =  %>% unique()),
           multiInput(inputId = "HGT",
                      label = "Select HGT Term:",
                      choices =  hgt_combined$Description %>% unique(),
                      selected = "SUMO E3 ligases SUMOylate target proteins"),
           #multiInput(inputId = "CORUM",
             #         label = "Select CORUM Term:",
            #          choices =  %>% unique()),
           downloadButton("downloadData", "Download")),
    column(width = 7,
           actionButton("reset", label = "Reset selection"),
           ggiraph::girafeOutput("plot")
    ),
    column(width = 10,
           h4("Selected targets"),
           # Button
           #downloadButton("downloadData", "Download"),
           tableOutput("datatab")
    )
  )
)

server<- function(input, output, session) {
  HGT<-reactive(hgt_combined %>% 
                    filter(Description %in% input$HGT)%>% 
                    pull(geneName))
  output$plot <- renderGirafe({
    x <- girafe(code = print(both %>% dplyr::mutate(
      diffex = case_when(
        both$lfc >= input$lfc_z & both$qvalue >= input$p ~ "UP",
        both$lfc <= -input$lfc_z & both$qvalue >= input$p ~ "DOWN",
        TRUE~"NO",
        both$lfc  < abs(input$lfc_z) | both$qvalue  < input$p ~ "NO"),
      delabel = case_when(
        both$lfc  >= input$lfc_z & both$qvalue  >= input$p ~ both$id,
        both$lfc  <= -input$lfc_z & both$qvalue  >= input$p ~ both$id)) %>% 
        ggplot(aes(x=lfc ,y=qvalue,col=diffex, label=delabel)) +
        geom_point_interactive(aes(color = diffex,tooltip = id, data_id = id)) +
        scale_color_manual(breaks = c("UP", "NO", "DOWN"),
                           values=c("red", "lightgrey", "blue"))+
        theme_minimal() +
        theme(legend.position = "none")+
        geom_label_repel(max.overlaps = Inf, stat = "unique")
        #geom_point_interactive(data = both %>% 
         #            filter(grepl(HGT(),id)))
        #geom_label_repel(both %>% filter(id %in% HGT()$geneName),label = id, max.overlaps = Inf)
      ),
      options = list(
        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),
        opts_selection(
          type = "multiple", css = "fill:#FF3333;stroke:black;"))
      #### add y lim buffer
    )
    x
  })
  selected <- reactive({
    input$plot_selected
  })
  targets <- reactive({
    input$Target
  })
  output$console <- renderPrint({
    input$plot_hovered
  })
  observeEvent(input$reset, {
    session$sendCustomMessage(type = 'plot_set', message = character(0))
  })
  
  data <- reactive({
    out <- both[both$id %in% selected() | 
                  both$id %in% targets() ]
    #out <- both %>% filter(both$id %in% selected() | both$id %in% Target())
    if( nrow(out) < 1 ) return(NULL)
    row.names(out) <- NULL
    out %>% arrange(lfc)  #%>% 
      #pivot_wider(names_from = Cells, values_from = c(lfc,qvalue)) %>% 
      #group_by(id) #%>%
      #summarise_each(funs(sum(., na.rm = TRUE)))
    #%>% 
    #select(Target=X, Cells, contains("log10")| contains("zscore"))
  })
  
  output$datatab <- renderTable(data())
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$Target,"_",input$p,"_",list(selected()),".csv", sep = "")
    },
    content = function(file) {
      write.csv(data(), file, row.names = FALSE)
    }
  )
  
}
shinyApp(ui, server)







borc <- c("SNAPIN","BLOC1S1","BLOC1S2","LOH12CR1","C10orf32","MEF2BNB","C17orf59")
tdp_stmn <-c("STMN2","TARDBP")
both %>%
  ggplot(aes(lfc, qvalue, label=id))+
  geom_point()+
  geom_point(data=both %>% filter(id %in% borc), color="blue",size=2)+
  geom_point(data=both %>% filter(id %in% tdp_stmn), color="red",size=2)+
  geom_label_repel(data=both %>% filter(id=="STMN2"|id=="TARDBP"))+
  geom_label_repel(data=both %>% filter(id %in% borc))+
  theme_classic()


test <-hgt_combined %>% 
  filter(Description %in% c("RNA Polymerase II Pre-transcription Events")) %>% 
  pull(geneName)
