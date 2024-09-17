#BiocManager::install("tools")

if(!"MAGeCKFlute" %in% installed.packages()) BiocManager::install("MAGeCKFlute")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)
library(msigdbr)
library(patchwork)
library(ggnewscale)

#genelist
#file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
 #                 "testdata/rra.gene_summary.txt")
#gdata = ReadRRA(file1)


rra.High.Low.STMN2.sgrna_summary <- read.delim("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_screenfiles/rra-High-Low-STMN2.sgrna_summary.txt") 

rra.High.Low.STMN2.sgrna_summary_filtered<-rra.High.Low.STMN2.sgrna_summary%>% 
  dplyr::filter(control_mean>3000|treat_mean>3000) %>% 
  dplyr::select(Gene)

##test with GRN vs control microglia
gdata = ReadRRA("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_screenfiles/rra-High-Low-STMN2.gene_summary.txt")

#------------------------------------------


## filter gdata based on filtered sgrna file - TDP43 and counterscreens
##  ---- ignore next few lines (JH)
  # list_filtered<-Combined_STMN2_mean %>% 
  #   dplyr::filter(!Gene%in%Halo_high$Gene,
  #                 !Gene%in%Halo_low$Gene) %>% 
  #   dplyr::select(Gene)

gdata<-gdata %>% dplyr::filter(id%in%list_filtered$Gene)

genelist = gdata$Score
#genelist = sort(genelist, decreasing = TRUE)
#genelist_up = sort(genelist, decreasing = F)

names(genelist) = gdata$id
genelist[1:10]
names(genelist)


###Hypergeometric test
# Alternative functions EnrichAnalyzer and enrich.HGT.
hgtRes1 = EnrichAnalyzer(genelist[genelist< -1], method = "HGT")
head(hgtRes1@result)
# hgtRes2 = enrich.HGT(genelist[genelist< -1])


hgtRes1_up = EnrichAnalyzer(genelist[genelist> 1], method = "HGT")
head(hgtRes1_up@result)
# hgtRes2 = enrich.HGT(genelist[genelist< -1])

###over-representation test
# Alternative functions EnrichAnalyzer and enrich.ORT.
ortRes1 = EnrichAnalyzer(genelist[genelist< -1], method = "ORT")
head(ortRes1@result)
# ortRes2 = enrich.ORT(genelist[genelist< -1])


###gene set enrichment test
# Alternative functions EnrichAnalyzer and enrich.GSE.
gseRes1 = EnrichAnalyzer(genelist[genelist< -1], method = "GSEA")
head(gseRes1@result)
# gseRes2 = enrich.GSE(genelist)


#####visualize
#require(ggplot2)
#df = hgtRes1@result
#df$logFDR = -log10(df$p.adjust)
#p = BarView(df[1:5,], "Description", 'logFDR')
#p = p + labs(x = NULL) + coord_flip()
#p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 5, title = "STMN2 Meta filtered negative modulating genes enrichplot")
barplot(hgtRes1_up, showCategory = 5, title = "STMN2 Meta filtered positive modulating genes enrichplot")


###dot plot
## top: up-regulated pathways; 
## bottom: down-regulated pathways
EnrichedView(hgtRes1, top = 0, bottom = 5, mode = 1)
EnrichedView(hgtRes1_up, top = 10, bottom = 0, mode = 1)
dotplot(hgtRes1, showCategory = 5,title="HGT STMN2 Meta filtered negative modulating genes")
dotplot(hgtRes1_up, showCategory = 5, title="HGT STMN2 Meta filtered positive modulating genes")


###other plots
hgtRes1@result$geneID = hgtRes1@result$geneName
hgtRes1_up@result$geneID = hgtRes1_up@result$geneName

cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 5, foldChange=genelist)
#emapplot(hgtRes1)



#gseaplot
gseaplot(gseRes1, geneSetID = 1, title = gseRes1$Description[1])
gseaplot(gseRes1, geneSetID = 1, by = "runningScore", title = gseRes1$Description[1])
gseaplot(gseRes1, geneSetID = 1, by = "preranked", title = gseRes1$Description[1])
#or
gseaplot2(gseRes1, geneSetID = 1:5)


#####go terms and pathways
## KEGG and REACTOME pathways
enrich = EnrichAnalyzer(geneList = genelist[genelist< -1], type = "KEGG+REACTOME")
EnrichedView(enrich, bottom = 10)


enrich_up = EnrichAnalyzer(geneList = genelist[genelist>1], type = "KEGG+REACTOME")
EnrichedView(enrich_up, top = 10)

## Only KEGG pathways
enrich = EnrichAnalyzer(geneList = genelist[genelist< -1], type = "KEGG")
EnrichedView(enrich, bottom = 10)
## Gene ontology
enrichGo1_down = EnrichAnalyzer(genelist[genelist< -1], type = "GOBP+GOMF")
EnrichedView(enrichGo1_down, bottom = 10)


enrichGo1_up = EnrichAnalyzer(genelist[genelist>1], type = "GOBP+GOMF")
EnrichedView(enrichGo1_up, top = 10)



####protein complex analysis
enrichPro = EnrichAnalyzer(genelist[genelist< -1], type = "CORUM")
enrichPro2 = EnrichAnalyzer(genelist[genelist> 1], type = "CORUM")

EnrichedView(enrichPro, bottom = 5) ### issue as argument 1 is not a vector
EnrichedView(enrichPro2, top = 5)+
  xlim(0.8,4)+
  labs(title="STMN2 Meta filtered positive modulating genes CORUM analysis")+
  theme_classic()

###Enrichment analysis on the combination of the gene sets
enrichComb = EnrichAnalyzer(genelist[genelist< -1], type = "GOBP+KEGG")
EnrichedView(enrichComb, bottom = 5)

####Limit the size of gene sets for testing
#enrich = EnrichAnalyzer(genelist[genelist< -1], type = "GOMF+GOBP+GOCC", limit = c(50, 500))
#EnrichedView(enrich, bottom = 10)

###Remove redundant results using EnrichedFilter.
enrich1 = EnrichAnalyzer(genelist[genelist< -1], type = "GOBP")
EnrichedView(enrich1, bottom = 20)+
  labs(title="STMN2 Meta filtered negative modulating genes GOMF+GOBP+GOCC")+
  theme_classic()

enrich2 = EnrichAnalyzer(genelist[genelist> 1], type = "GOMF+GOBP+GOCC")
EnrichedView(enrich2, top = 20,
             mode = 1)+
  labs(title="STMN2 Meta filtered positive modulating genes GOMF+GOBP+GOCC")+
  theme_classic()
?EnrichedView


### patchwork
a1<-barplot(hgtRes1, showCategory = 10, title = "STMN2 Meta filtered negative modulating genes enrichplot")
a2<-barplot(hgtRes1_up, showCategory = 10, title = "STMN2 Meta filtered positive modulating genes enrichplot")
b1<-dotplot(hgtRes1, showCategory = 10)+
  labs(title = "STMN2 Meta filtered negative modulating HGT analyis")
b2<-dotplot(hgtRes1_up, showCategory = 10)+
  labs(title = "STMN2 Meta filtered positive modulating HGT analyis")

c1<-cnetplot(hgtRes1, 5, 
             foldChange= genelist, 
             circular = TRUE, colorEdge = TRUE,
             max.overlaps = Inf)+
  labs(title = "STMN2 Meta filtered negative modulating HGT analyis")
c2<-cnetplot(hgtRes1_up, 5, 
             foldChange= genelist, 
             circular = TRUE, colorEdge = TRUE,
             max.overlaps=Inf
             )+
  labs(title = "STMN2 Meta filtered positive modulating HGT analyis")
?cnetplot

d1<-EnrichedView(enrichPro, bottom = 10)+
  labs(title = "STMN2 Meta filtered negative modulating CORUM analyis")
d2<-EnrichedView(enrichPro2, top = 10)+
  labs(title = "STMN2 Meta filtered positive modulating genes CORUM analyis")

p<-a1+a2+b1+b2+
  plot_layout(ncol = 2)

?patchwork
p
ggsave("meta_filtered.tiff",p,width = 20,height = 15,units = "in")

a1
a2
b1
b2
c1
c2
d1
d2

ggsave("cnet_down.tiff",c1,width = 15,height = 10,units = "in",bg="white")
ggsave("cnet_up.tiff",c2,width = 15,height = 10,units = "in",bg="white")

ggsave("corum_down.tiff",d1+theme_classic(),width = 10,height = 5,units = "in",bg="white")
ggsave("corum_up.tiff",d2+theme_classic(),width = 10,height = 5,units = "in",bg="white")








###check the genes:
phos <- enrichGo_RNAseq_up@result %>% filter(Description=="nucleosome")
phos[["geneName"]]
telomerase <- hgtRes1 %>% filter(Description=="negative regulation of telomere maintenance via telomerase")
#hgtRes1 %>% filter(grepl('telomerase',Description))
telomerase@result[["geneName"]]



GOtable <-enrichGo_RNAseq_up@result

enrichPro@result %>% 
  filter(Description=="BORC complex")

####Session info
sessionInfo()



