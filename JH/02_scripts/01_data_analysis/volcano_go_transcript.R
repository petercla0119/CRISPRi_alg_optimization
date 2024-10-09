library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(shiny)
library(ggiraph)
library(shinyWidgets)
library(scales)
library(GO.db)
library(org.Hs.eg.db)
library(GOSim)
library(biomaRt)

##
#BiocManager::install("org.Hs.eg.db")

## load gene summary and guide summary
STMN2 <- read.delim("~/Desktop/Projects/FACS_projects/STMN2/STMN2_mscarlet_221205/input/rra-High-Low-STMN2.gene_summary.txt")
rra.High.Low.STMN2.sgrna_summary <- read.delim("~/Desktop/Projects/FACS_projects/STMN2/STMN2_mscarlet_221205/input/rra-High-Low-STMN2.sgrna_summary.txt")
## you just need to change these paths for the gene_summary and sgrna_summary files


## combine p values 
negative<- STMN2 %>% 
  filter(neg.lfc<0) %>% 
  dplyr::select(id,lfc=neg.lfc, score = neg.score,pvalue =neg.p.value,fdr = neg.fdr)

positive<- STMN2 %>% 
  filter(pos.lfc>0) %>% 
  dplyr::select(id,lfc=pos.lfc, score = pos.score,pvalue =pos.p.value,fdr = pos.fdr)

both<-rbind(negative, positive) %>% 
  mutate(
    qvalue = -log10(pvalue)
  )

# get all go terms via biomaRt
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
Both_list<- getBM(attributes=c("hgnc_symbol",'entrezgene_id',"go_id", "name_1006" ,"definition_1006"), filters="hgnc_symbol",values=both$id,
                  mart=mart, uniqueRows=T)
Both_list<- Both_list %>% 
  mutate(go_id = dplyr::na_if(go_id,y="")) %>% drop_na()
both_go <- both %>% merge(Both_list %>% 
                            mutate(id=hgnc_symbol), by = "id")

both_go_tran <- rra.High.Low.STMN2.sgrna_summary %>% 
  dplyr::select(Gene, control_mean, treat_mean) %>% 
  mutate(id=Gene) %>% merge(both_go, by = "id")

## nuclear import genes
transcription<-both_go %>% 
  filter(grepl("transcription",name_1006)) %>% 
  dplyr::select(id) %>% 
  unique()
translation<-both_go %>% 
  filter(grepl("translation",name_1006)) %>% 
  dplyr::select(id) %>% 
  unique()
nuclear_import<-both_go %>% 
  filter(grepl("nuclear import",name_1006)) %>% 
  dplyr::select(id) %>% 
  unique()
nuclear_export<-both_go %>% 
  filter(grepl("nuclear export",name_1006)) %>% 
  dplyr::select(id) %>% 
  unique()
nucleocytoplasmic_transport<-both_go %>% 
  filter(grepl("nucleocytoplasmic transport",name_1006)) %>% 
  dplyr::select(id) %>% 
  unique()
nuclear_envelope<-both_go %>% 
  filter(grepl("nuclear envelope",name_1006)) %>% 
  dplyr::select(id) %>% 
  unique()



## can use this for the kegg pathway
#library(KEGGREST)
#kegg<-left_join(keggLink("pathway", "hsa") %>% 
#                 tibble(pathway = ., eg = sub("hsa:", "", names(.))) %>%
#                mutate(
#                 symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
#                ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
#             ), keggList("pathway", "hsa") %>% 
#            tibble(pathway = names(.), description = .))

#both_go_tran_kegg <- kegg %>% mutate(
# id=symbol) %>% 
#select(id, description) %>% 
#merge(both_go_tran %>% 
#       select(id, control_mean,treat_mean, lfc, score,pvalue, go_id,
#             fdr,qvalue,name_1006,definition_1006), by="id")


###this is the custom UI
###Omic x is proteomics and y is transcriptomics
ui <- fluidPage(
  tags$h1("STMN2 screen hits GO"),
  fluidRow(
    column(width = 3,
           numericInput("low", label="Low mean count", 
                        value = "600"),
           numericInput("high", label="High mean count", 
                        value = "600"),
           numericInput("p", label="qvalue", 
                        value = "6"),
           numericInput("lfc_z", label="LFC", 
                        value = "5"),
           multiInput(inputId = "Target",
                      label = "Select Target:",
                      choices = both_go_tran$id %>% unique(),
           ),
           multiInput(inputId = "GO_Pathway",
                      label = "Select GO term:",
                      choices = both_go_tran$name_1006 %>% unique(),
           ),
           downloadButton("downloadData", "Download")),
    column(width = 7,
           actionButton("reset", label = "Reset selection"),
           ggiraph::girafeOutput("plot")
    ),
    column(width = 10,
           h4("Selected targets"),
           tableOutput("datatab")
    )
  )
)

server<- function(input, output, session) {
  output$plot <- renderGirafe({
    x <- girafe(code = print(
      both_go_tran %>%
        distinct(id, .keep_all = T) %>% 
        filter(control_mean>=input$low| 
                 treat_mean>=input$high) %>% 
        ggplot(aes(x=lfc ,y=qvalue, label=id
        )) +
        geom_point_interactive(aes(tooltip = id, data_id = id), color="grey", alpha=0.6) +
        geom_point_interactive(data=both_go_tran %>% 
                                 dplyr::filter(name_1006 %in% input$GO_Pathway),
                               aes(tooltip=id, data_id=id, color=name_1006), alpha = 0.6)+
        geom_point_interactive(data=both_go_tran %>% 
                                 dplyr::filter(id %in% input$Target),
                               aes(tooltip=id, data_id=id),color= "red", alpha = 0.6)+
        geom_label_repel_interactive(data=both_go_tran %>% 
                                       dplyr::filter(
                                         abs(lfc) >= input$lfc_z,
                                         qvalue >= input$p) %>% 
                                       dplyr::filter(
                                         name_1006 %in% input$GO_Pathway
                                       ) %>% 
                                       distinct(id, name_1006, .keep_all = T),
                                     aes(color=name_1006),
                                     box.padding = 0.7, max.overlaps = Inf)+
        ylim(NA,8)+
        theme_minimal()+
        labs(y="-log10pvalue",
             x="LFC of STMN2 expression")+
        theme(legend.position="bottom")+
        guides(colour = guide_legend(nrow = 6))
    ),
    options = list(
      opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),
      opts_selection(
        type = "multiple", css = "fill:#FF3333;stroke:black;"))
    )
    x
  })
  selected <- reactive({
    input$plot_selected
  })
  targets <- reactive({
    input$Target
  })
  GO <- reactive({
    input$GO_Pathway
  })
  output$console <- renderPrint({
    input$plot_hovered
  })
  observeEvent(input$reset, {
    session$sendCustomMessage(type = 'plot_set', message = character(0))
  })
  
  data <- reactive({
    out <- both_go_tran[both_go_tran$id %in% selected() | both_go_tran$id %in% targets()
                        | both_go_tran$id %in% (both_go_tran %>% 
                                                  dplyr::filter(
                                                    abs(lfc) >= input$lfc_z,
                                                    qvalue >= input$p,
                                                    name_1006 %in% input$GO_Pathway
                                                  ) %>% 
                                                  distinct(id) %>% 
                                                  pull), ]
    if( nrow(out) < 1 ) return(NULL)
    row.names(out) <- NULL
    out %>% 
      arrange(id) %>%
      dplyr::select(!definition_1006) %>% 
      group_by(id) %>% 
      mutate(GO_id = paste0(go_id, collapse = " "),
             GO_term = paste0(name_1006, collapse = " ")) %>% 
      distinct(GO_id, .keep_all = T) 
  })
  
  output$datatab <- renderTable(data())
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$GO_Pathway,"_",list(selected()),".csv", sep = "")
    },
    content = function(file) {
      write.csv(data(), file, row.names = FALSE)
    }
  )
  
}
shinyApp(ui, server)

