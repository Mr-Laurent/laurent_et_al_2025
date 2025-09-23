#### app_Macs_hm_v3.R : added annotation selection, change cluster order etc 
## 3.3 : corrected gene log heatmap to normalize first by cpm the gene expression by metacell
## 3.5 : updated "Modscores by group" panel with genes and perso signature + colors + order
## 3.6 : updated "Modscores by group" panel with stats and personalization
## 3.7 : updated "Modscores by group" panel with more stats and personalization


setwd("~/laurent_et_al_2025/Viz_Apps/")

library(shiny)
library(ggplot2)
library(magick)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(raster)
library(seriation)
library(dendsort)
library(matrixStats)
library(datasets)
library(circlize)
library(heatmaply)
library(plotly)
library(scales)
library(Matrix)
library(shinyWidgets)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(rstatix)

cell_type="Macs"

load("~/laurent_et_al_2025/Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog.rd")
load("~/laurent_et_al_2025/Grouped_objects/genelist_16nov22_Macs.rd")
load("~/laurent_et_al_2025/Grouped_objects/Macs_dgCmc_16nov22_counts.rd")
hm_zmod<-as.matrix(metadf_z[,3:152])
load("~/laurent_et_al_2025/Grouped_objects/ht_Enrich_Macs.rd")


#---------------------------------------#

colnames(genelist)<-c("Mods","Gens")
hm_zmod<-as.matrix(metadf_z[,c(154:303)])
df_zmod<-as.data.frame(metadf_z[,c(2:3,154:303)])
colnames(df_zmod)<-sub("([a-zA-Z_]+)(\\d+)([a-zA-Z]+)", "Module_\\2", colnames(df_zmod))  
preloco_nc<-as.matrix(log(1+dgcmtxcounts))
allsum<-Matrix::colSums(preloco_nc) 
row_mod_ord<- as.vector(unlist(row_order(ht)))
data.cpm<-1e+06*sweep(dgcmtxcounts,2,colSums(dgcmtxcounts),`/`)  # if you don't want another package
preloco<-as.matrix(log2(1+(data.cpm/10)))  # Divide by 10 the CPM for more readability

###----------------------------------------------------------### 
#### Function to transform list of genes to vector readable #### 
#### by the program, with exclusion of unknown genes        ####
###----------------------------------------------------------### 
gene_input2vect<-function(genelist,dgcmtx){
  gnshow<-unlist(strsplit(genelist,split=","))
  gnshow<-gsub(" ","",gnshow)
  notfound<-setdiff(gnshow,rownames(dgcmtx))
  print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  return(setdiff(gnshow,notfound) )
}
###----------------------------------------------------------### 
####    Function to have multiple pvalue levels as stars     ### 
get_significance_label <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )}
###----------------------------------------------------------### 
#### Function to get all the genes as vector for each list  #### 
findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}  
###----------------------------------------------------------### 
gnlist<-"CXCL8,SERPINB2,MFSD2A,F3,IL1A,TNF,IL7R,IDO1,INHBA,CD274,IL3RA,P2RX7,SLAMF1,CRIM1,SLC39A8,SLC7A11,ADAM19,SPP1,DNAJB4,HSPA6,FCGR2B,CA2,SCD,BNIP3,RRAD,SLC2A1,IFI44L,CXCL10,CXCL9,IFIT2,OAS2,CXCL11,IL23A,CD40,SLAMF7,CCL20,CXCL3,CXCL2,
S100A8,S100A9,VCAN,FCN1,RIPOR2,S100A12,APOBEC3A,SELL,FCAR,ALOX5AP,RETN,MCEMP1,OSM,IL10,CDK1,RAB7B,TPRA1,LY6K,BCL11A,PHLDB1,MMP19,EREG,ATP5PB,MRPS21,RNF7,NDUFB5,ATP5ME,CSF1R,CD72,CD101,IL18,C1QA,VSIG4,DNASE1L3,DAB2,PDK4,SLCO2B1,CMKLR1,CD209,FUCA1,MERTK,SELENOP,FOLR2,IGF1,AXL,LYVE1"
#gnlist<-"CD7,CD2,CD3D,STMN1,PCNA,TNFRSF4,TNFRSF18,CTLA4,IL2RA,FOXP3,CCR7,SELL,LEF1,CST7,GZMA,GZMK,CD8A,CD8B,CD69,ID2,ANKRD28,GPR171,TNFAIP3,ANXA1,IL7R,TRDC,TYROBP,FCER1G,PTGDR,KRT81,KRT86,KIT,PCDH9,ALDOC,LINC00299,CMC1, XCL2, IL2RB, KLRD1,KLRF1, CLIC3, MCTP2,MZB1,TNFRSF17,SEC11C,DERL3,XBP1,IGHA2,IGHM,IGHG1,IGHG3,BANK1, VPREB3, CD24, ARHGAP24, FCRLA, RALGPS2, TNFRSF13C, SPIB,MS4A1, CD79B, CD19,IGHD, FCER2,TCL1A, CD72,CD27,CLECL1,TNFRSF13B,FCGR3A,CD14,CD68, SLC40A1, STAB1, CSF1R, MS4A4A, SLCO2B1, MAFB,MS4A7,C1QA,C1QB,C1QC,IL1RN,S100A8,S100A9,TREM1,KYNU,C1orf54,CLEC9A,CADM1,CPNE3,XCR1,IDO1,CLEC10A,FCER1A,CD1C,CD1D,GPR157,LAMP3,DAPP1,FSCN1,CCL19,CCL22,EBI3,GZMB,TCF4,LILRA4,IRF7,CLEC4C,IL1RL1,CPA3,MS4A2,TPSB2,TPSAB1,ADCYAP1,GATA2,HPGDS,HPGD,VWF, RAMP3, NPDC1, JAM2, PLVAP, NOTCH4, HSPG2, ESAM, CYYR1, CD93, ICAM2, S1PR1, RAMP2,CD36, CA4, TMEM88, FLT1,ACKR1, SELP,DUSP23,PROX1,MMRN1,CCL21,IGF2,LYVE1,RGS5,NDUFA4L2,C11orf96,ACTG2,MYH11,CXCL14,ADH1B,CTSK,MMP2,LUM,PTGDS,TCF21,ADAMDEC1,CCL13,CCL8,ADAM28,HAPLN1, ABCA8,CFD,THY1,CHI3L1,BGN,PDPN,TNFRSF12A,MPZ,NRXN1,SCN7A,SEMA3B,SOX2,MYOT,GFRA3,TUBB2B,GPM6B,PLP1,XKR4,ALDH1A1,CRYAB,NTM,ANK3,LGI4,S100B"
ann_df<-data.frame(metacell=rownames(metadf_z),annotation="undefined")

# ---- Define UI for app that draws a histogram ---- #
ui <- fluidPage(
  navbarPage("MCexplorer_HM",
             tabPanel("z-scored HM",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Heatmap of the Z-scaled module, displayed as mean in metacells"),
                                 verbatimTextOutput('prefoo')
                               )
                        ),
                        mainPanel(
                          plotOutput("hmzall",brush = "ht_brush", click = "ht_click", height = 500, width = 900),
                          tags$hr()
                        )
                      )
                      
             ),
             tabPanel("z-scored Heatmaply",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Heatmap of the Z-scaled modules enrichment by metacells"),
                                 textAreaInput("clsordhm", h4("reorder clusters (n?)"),height ='200%',width='100%',
                                               value = paste0(c(1:40),collapse=",")),
                                 selectInput("var", label = "Choose a module to display", choices=c(1:150),selected=1),
                                 textInput("geneask", h4("Looking for a gene ?"),placeholder = "no gene selected"),
                                 verbatimTextOutput("vgen"),
                                 # sliderInput(inputId = "HMwidz2", label = "Width", value=1000, min=0, max=1600,step=100),
                                 # sliderInput(inputId = "HMheign2", label = "Height", value=800, min=0, max=1600,step=100),
                                 sliderInput(inputId = "hmgnprop", label = "Gene heatmap size", value=0.4, min=0, max=1,step=0.05),
                                 sliderInput(inputId = "hmgnlbsz", label = "HM y-labels size", value=0.75, min=0, max=2,step=0.05),
                                 textAreaInput("geneshow", h4("Genes to plot"),height ='200%',width='100%',
                                               value = gnlist),
                                 actionButton("click", "Run Heatmap"),
                                 tags$hr(),
                                 verbatimTextOutput("notf_gen"),
                                 textInput("clusterschi", h4("Chi? test"),value = ""),
                                 actionButton("click_run_chi", "Run chi2 test"),
                                 verbatimTextOutput("chirep"),
                               )
                        ),
                        mainPanel(
                          textOutput("text"),
                          verbatimTextOutput("verb"),
                          tags$hr(),
                          plotlyOutput("hmzall2", height = 1500, width = 1000),
                          textInput("topmod", h4("Top modules (>0.8) for the cluster:"),value = ""),
                          verbatimTextOutput("topmodout"),
                        )
                      )   
             ),
             tabPanel("Modscores by clusters",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Heatmap of the Z-scaled modules enrichment by metacells"),
                                 selectInput("var2", label = "Choose a module to display", choices=c(1:150),selected=1),
                               )
                        ),
                        mainPanel(
                          textOutput("text2"),
                          verbatimTextOutput("verb2"),
                          plotOutput("mod_clust_bar", height = 500, width = 900),
                        )
                      )
             ),
             tabPanel("Mc annotation",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Annotation"),
                      fluidRow(column(4,
                                      wellPanel(
                                        h3("Cluster selection"),
                                        textAreaInput("ann_clu",label= NULL,height ='200%',width='100%',placeholder="Cluster numbers, comma-separated")
                                      )),
                               column(4,
                                      wellPanel(
                                        h3("Add metacells"),
                                        fluidRow(
                                          column(2, HTML('<h4><b>Mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_add_mc1",label= NULL,placeholder="XXX"),offset = 0, style='padding:0px;'),
                                          column(3, HTML('<h4><b>to mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_add_mc2",label= NULL,placeholder="XXX"), offset = 0,style='padding:0px;')
                                        )
                                      )),
                               column(4,
                                      wellPanel(
                                        h3("Remove metacells"),
                                        fluidRow(
                                          column(2, HTML('<h4><b>Mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_rem_mc1",label= NULL,placeholder="XXX"),offset = 0, style='padding:0px;'),
                                          column(3, HTML('<h4><b>to mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_rem_mc2",label= NULL,placeholder="XXX"), offset = 0,style='padding:0px;')
                                        )
                                      ))
                      ),
                      #br(),
                      #actionButton("click", "get metacell list"),
                      tags$hr(),
                      fluidRow(
                        column(8,
                               wellPanel(
                                 h3("Selected metacells to annotate"),
                                 verbatimTextOutput("out_mcs"),
                               )),
                        column(4,
                               wellPanel(
                                 h4("Name of the group"),
                                 textInput(inputId = "grp_name",label= NULL,value=""),
                                 actionButton("btn_save_name", "save name"),
                               ))
                      ),
                      fluidRow(
                        column(8,
                               wellPanel(
                                 h3("Table of proportions"),
                                 tableOutput("out_prop_cellannot")
                               )),
                        column(4,
                               wellPanel(
                                 textInput("file_name", "File Name (without extension):"),
                                 downloadButton("save_csv", "Save CSV")
                               ))),
                      fileInput("annot_csv", "Or choose an annotation file:"),
                      actionButton("btn_load_name", "load name")
             ),
             tabPanel("Modscores by group",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Boxplot of the Z-scaled modules enrichment or gene expression by cell category"),
                                 # selectInput("var3", label = "Choose a module to display", choices=c(1:150),selected=1),
                                 selectInput("name_fromlist", "Select data to plot",# \n(personal gene signatures are un-corrected for Tcells)"
                                             choices = c("Gene (log CPM) or Gene signature",colnames(df_zmod))),
                                 uiOutput("reactive_txt"),
                                 textInput("palette_1", "personnalized palette",value = "c('Mono-like 1'='#874037','Mono-like 2'='#C3A34B','Mono-like 3'='#f77774','Mono-like 4'='#ffa3a1','Mono-like 5'='#77c799','Macrophage-like 6'='#c967eb','Macrophage-like 7'='#5fc2ed','Macrophage-like 8'='#4F88B9','Macrophage-like 9'='#5C538B')" ),
                                 textInput("bx_order", "X axis order", value="c('Mono-like 1','Mono-like 2','Mono-like 3','Mono-like 4','Mono-like 5','Macrophage-like 6','Macrophage-like 7','Macrophage-like 8','Macrophage-like 9')"),
                                 # switchInput("toggle_legend", "Show Legend", value = TRUE),
                                 actionButton("toggle_legend", "Remove legend ?", icon("image"),  #icon can be found : https://fontawesome.com/start
                                              style="color: #000000; background-color: #f0dda8; border-color: #e39424"),
                                 # switchInput("toggle_stats", "Show Stats", value = FALSE),
                                 selectInput("choose_stats", "Show Stats", choices = c("No stats","Symbol stats","Numeric stats","Symbol stats - vs mean","Numeric stats - vs mean")),
                                 uiOutput("reactive_txt2"),
                                 uiOutput("reactive_slider"),
                                 sliderInput(inputId = "bx_hei", label = "Height", value=600, min=0, max=1200,step=50),
                                 sliderInput(inputId = "bx_wid", label = "Width", value=550, min=0, max=1200,step=50),
                                 sliderInput(inputId = "bx_siz", label = "Outline size", value=2, min=0, max=5,step=0.2),
                                 sliderInput(inputId = "vjust_stat", label = "vjust stat text", value=0.8, min=0, max=2,step=0.05),
                                 sliderInput(inputId = "y_txt_size", label = "y axis size", value=3, min=0, max=10,step=0.2),
                                 downloadButton('export_barplt'),
                                 # switchInput("toggle_code", "Show Code", value = TRUE)
                               )
                        ),
                        mainPanel(
                          textOutput("text3"),
                          verbatimTextOutput("verb3"),
                          plotOutput("mod_group_bar", height = 500, width = 900,  plotOutput(outputId = "plt1")),
                          
                        )
                      )
                      #  DEBUG :
                      # df<-read.csv2(file="G:/Mon Drive/UG_metacells/MyelEnr_v1/EnrMacs_Mono1merged_annot_june24.csv",header = TRUE,sep=",")
                      # input<-c()
                      # input$name_fromlist<-'Gene (log CPM) or Gene signature'
                      # input$Not_on_list  <-"CD14"
                      # input$palette_1="c('Mono-like 1'='#874037','Mono-like 2'='#C3A34B','Mono-like 3'='#f77774','Mono-like 4'='#ffa3a1','Mono-like 5'='#77c799','Macrophage-like 6'='#c967eb','Macrophage-like 7'='#5fc2ed','Macrophage-like 8'='#4F88B9','Macrophage-like 9'='#5C538B')"
                      # input$bx_order="c('Mono-like 2','Mono-like 3','Mono-like 4','Mono-like 5','Macrophage-like 6','Macrophage-like 7','Macrophage-like 8','Macrophage-like 9')"
                      # input$palette_1=""
                      # input$bx_order=""
                      # input$toggle_legend= TRUE
                      # input$choose_stats="No stats"
                      # input$bx_hei=600
                      # input$bx_wid=550
                      # input$bx_siz=1.3
                      # 
                      # input$compare_grp="all"
                      # input$bx_stat_siz=1.3
             )
             
  ))

hm_env = new.env()

server <- function(input, output, session) {
  ## Search and provide genes from modules:
  output$text <- renderText({paste0(genelist$Mod[as.numeric(input$var)]," genes: ")})
  output$verb <- renderText({genelist$Gens[as.numeric(input$var)] })
  output$text2 <- renderText({paste0(genelist$Mod[as.numeric(input$var2)]," genes: ")})
  output$verb2 <- renderText({genelist$Gens[as.numeric(input$var2)] })
  output$text3 <- renderText({paste0(genelist$Mod[as.numeric(input$var3)]," genes: ")})
  output$verb3 <- renderText({genelist$Gens[as.numeric(input$var3)] })
  output$vgen <- renderText({
    if(input$geneask%in%unlist(strsplit(genelist$Gens,split=","))){
      print(paste0("Module ",findgene(strsplit(genelist$Gens,split=","),input$geneask)))
    }else{print("Not found")}
  })
  bxhei<-function(){return(input$bx_hei)}
  bxwid<-function(){return(input$bx_wid)}
  bxsiz<-function(){return(input$bx_siz)}
  
  output$hmzall <- renderPlot({
    print("lezgo")
    ##Compute the heatmap figure (Macs/DC have different code as we already have a fixed order to match our previous analyses)
    set.seed(42)
    validate(need(exists("ht")==T,"HT OBJECT NOT FOUND / Please precompute heatmap \n\n(it allows a common ground of matrix order for all configs \nas setting seed doesn't always work)"))
    hm_env$ht = ht 
    ComplexHeatmap::draw(hm_env$ht)
    print("draw")
    
    
    metacells_kclust<- column_order(hm_env$ht) 
    hm_env$metacells_kclust<-metacells_kclust
    
    for (i in 1:length(metacells_kclust)){
      if (i == 1) {
        clu <- t(t(rownames(hm_zmod)[metacells_kclust[[i]]]))
        out <- cbind(clu, paste("cluster", i, sep=""))
        colnames(out) <- c("Metacell", "Cluster")
      } else {
        clu <- t(t(rownames(hm_zmod)[metacells_kclust[[i]]]))
        clu <- cbind(clu, paste("cluster", i, sep=""))
        out <- rbind(out, clu)
      }
    }
    rm(clu)
    out<-as.data.frame(out)
    hm_env$out<-out
    
    # Automatic Table to display on the shiny app (Later do a tool to select a k cluster (can interactive give which cluster?) and it shows the list the genes in Module positives) :
    mod_clu<-c()
    mcs_clu<-list()
    nmod_clu<-list()
    for(i in 1:40){
      clust=paste0("cluster",i)
      mcs_clu[[i]]<-out$Metacell[which(out$Cluster%in%clust)]
      mod_clu[[i]]<-if(length(mcs_clu[[i]])==1){colnames(hm_zmod)[which(hm_zmod[mcs_clu[[i]],]>0.8)]
      }else{colnames(hm_zmod)[which(colMeans2(as.matrix(hm_zmod[mcs_clu[[i]],]))>0.8)]}
      nmod_clu[[i]]<-length(mod_clu[[i]])
    }
    hm_env$mcs_clu<-mcs_clu
    hm_env$mod_clu<-mod_clu
    
    
  })
  
  # memus1<-mem_used() 
  # output$prefoo <- renderText({
  #   print(paste0("memusage: ",memus1))
  # })
  
  output$notf_gen = renderText({
    gnshow<-unlist(strsplit(input$geneshow,split=","))
    gnshow<-gsub(" ","",gnshow)
    notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
    print( paste0("[",paste0(notfound,collapse =","),"] not found") )
    ## ## gc()
  })
  
  
  output$hmzall2  = renderPlotly({
    input$click
    req(input$click) #to prevent print at first lauch
    isolate({
      
      clorhm<-as.numeric(unlist(strsplit(input$clsordhm,split=",")))
      htmxord<-t(hm_zmod)[ as.vector(unlist(row_order(hm_env$ht))) , as.vector(unlist(column_order(hm_env$ht)[clorhm])) ] 
      #htmxord<-t(hm_zmod)[ as.vector(unlist(row_order(hm_env$ht))) , as.vector(unlist(column_order(hm_env$ht))) ] 
      
      
      
      ## Prepare the same heatmap as in 1st window but not cut in clusters, and interactive (plotly) 
      hm_env$htpt1<-heatmaply(htmxord, dendrogram = "none",limits=range(-2,2), fontsize_row = 5,colors = "RdBu",Rowv=F,Colv=F,legendgroup="1st",showlegend = T,coloraxis = 'coloraxis',
                              scale_fill_gradient_fun =scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish) )
      
      
      ##  Second Heatmap displaying the selected genes (below the module heatmap)
      gnshow<-unlist(strsplit(input$geneshow,split=","))                                                                          
      gnshow<-gsub(" ","",gnshow)                                                                                                 
      notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
      loco<-preloco[setdiff(gnshow,notfound),colnames(htmxord)]

      pt2<-heatmaply(loco, fontsize_row = rel(input$hmgnlbsz),fontsize_col = 7,dendrogram = "none",Rowv=T,Colv=F,legendgroup="2nd",coloraxis = 'coloraxis2',scale_fill_gradient_fun =scale_fill_gradientn(
        colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis")) ) 
      
      
      rm(loco)
      rm(htmxord)
      gc()
      
      subplot(hm_env$htpt1, margin = 0.01) %>% 
        subplot(pt2,nrows=2, heights = c((1-input$hmgnprop),input$hmgnprop),shareX = TRUE,margin=0.01)
      
    })
    
  }) # ,height=hm_heign2, width=hm_widz2)   MARCHE PAS ??? Try to find a way to resize or at least change pt1/pt2 proportions ?
  
  gc()
  ## How much memory is in use ?
  # output$foo5 <- renderText({
  #   req(input$click2) #to prevent print at first lauch
  #   isolate({
  #     memus<-mem_used() 
  #     print(paste0("memusage: ",memus))
  #   })
  # })
  output$topmodout <- renderText({
    print(paste0(gsub("corr_z_","",hm_env$mod_clu[as.numeric(input$topmod)]),collapse = ","))
    
  })
  ## Chi test to get genes differential between clusters [Not working anymore? Adapt to the matrix format ?]
  output$chirep<- renderText({
    req(input$click_run_chi) #to prevent print at first lauch
    isolate({
      clusterschi<-unlist(strsplit(input$clusterschi,split=","))                                                                          
      clusterschi<-gsub(" ","",clusterschi) 
      cel_in_chi<-c()
      for(i in as.numeric(clusterschi)){
        cel_in_chi<-append(cel_in_chi,unlist(hm_env$mcs_clu[i])) }
      cel_in_chi
      counts=dgcmtxcounts[,cel_in_chi]
      
      gene_mask=apply(counts,1,max)>3
      counts=counts[gene_mask,,drop=F]
      metacell_tot=colSums(counts)
      testarray<-array(c(counts,matrix(metacell_tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts))
      
      testchi<-suppressWarnings({apply(counts,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}) })
      testchi<-t(testchi)
      testchi=testchi[!is.na(testchi[,1]),]
      
      res=cbind(testchi,adjp=p.adjust(testchi[,1],method="BH"))
      mask=rownames(res)%in%res[,3]<0.01&  rowSums(counts)>=10
      
      a=res[mask,]
      genes_to_show=head(rownames(a)[order(a[,2],decreasing=T)],30)
      
      genes_ts_good=paste(genes_to_show,collapse=",")
      rm(counts)
      rm(testarray)
      print(paste0(genes_ts_good))
    })
  })
  
  output$mod_clust_bar <- renderPlot({
    print("plotting")
    df_zmod_h<-as.data.frame(hm_zmod)
    df_zmod_h$cluster<-hm_env$out$Cluster[match(rownames(df_zmod_h),hm_env$out$Metacell )]
    df_zmod_h$cluster<-factor(df_zmod_h$cluster, levels = paste0("cluster", 1:length(hm_env$metacells_kclust)))
    # 
    print(ggplot(data = df_zmod_h, aes(x = cluster, y = !!sym(paste0("My_Mod_",input$var2,"scaled")) )) +
            geom_boxplot(width=0.8, fill="white",size=0.35,outlier.alpha = 0.5)+
            xlab("Cluster")+theme_bw()+
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.85,size = rel(1)),
                  axis.title.x =element_text(size=rel(1),face="bold"),
                  axis.text.y = element_text(size = rel(1.5)),
                  axis.title.y =element_text(size=rel(1),face="bold"),
                  plot.title = element_text(size = rel(1.5), face = "bold"),
                  plot.subtitle = element_text(size = rel(1)),
                  panel.border = element_blank(),
                  legend.position = "none",
                  axis.line = element_line(colour = "black",size = 0.6)) + geom_jitter(height = 0,size=0.3, width = 0.2,alpha=0.3))
    
  })
  
  # Show all metacells selected in tab Mc annotation
  output$out_mcs = renderText({
    ht_mc_ord<-rownames(hm_zmod)[as.vector(unlist(column_order(ht))) ]  
    
    if(!input$ann_clu=="" ){
      ann_clu_show<-unlist(strsplit(input$ann_clu,split=","))
      ann_clu_show<-gsub(" ","",ann_clu_show)
      #Get the clusters from the heatmap (as cluster 1 is not named "1", I use numeric to get the first one, whatever its name)
      # And affect it to rownames of  hm_zmod as metacell "42" in heatmap is named "mc41" in the end
      mc_part1<-rownames(hm_zmod)[as.vector(unlist(  column_order(ht)[as.numeric(ann_clu_show)]   )) ]
    }else{mc_part1=""}
    
    
    if(!input$ann_add_mc1=="" && !input$ann_add_mc2=="" ){
      mc_part2<-c( ht_mc_ord[c(which(ht_mc_ord==paste0("mc",input$ann_add_mc1)):which(ht_mc_ord==paste0("mc",input$ann_add_mc2)))])
    }else{mc_part2=""}
    if(!input$ann_rem_mc1=="" && !input$ann_rem_mc2=="" ){
      mc_part3<-c( ht_mc_ord[c(which(ht_mc_ord==paste0("mc",input$ann_rem_mc1)):which(ht_mc_ord==paste0("mc",input$ann_rem_mc2)))])
    }else{mc_part3=""}
    mc_selected<-setdiff( append(mc_part1,mc_part2), mc_part3 )
    hm_env$current_mc_selected<-mc_selected
    print(mc_selected)
  })
  
  
  env_ann_df <- reactiveValues(df = ann_df)
  observeEvent(input$btn_save_name, {
    
    annotation <- input$grp_name
    
    if (!is.null(hm_env$current_mc_selected) && !is.null(input$grp_name)) {
      # Update the annotation in the dataframe
      env_ann_df$df[env_ann_df$df$metacell %in% hm_env$current_mc_selected, "annotation"] <- input$grp_name
    }
  })
  
  # Display cell proportions from given names
  output$out_prop_cellannot<- renderTable({
    table(env_ann_df$df$annotation)
  })
  # Save the annotations in a csv, give file name and choose where to save
  output$save_csv <- downloadHandler(
    filename = function() {
      paste0(input$file_name, ".csv")
    },
    content = function(file) {
      write.csv(env_ann_df$df, file, quote = F,col.names=T,row.names=F,sep = ",")
    })
  loadednames <- reactiveVal(NULL)
  
  observeEvent(input$annot_csv, {
    if (!is.null(input$annot_csv)) {
      loadednames(read.csv(input$annot_csv$datapath))
    }
  })
  
  observeEvent(input$btn_load_name, {
    env_ann_df$df <- loadednames()
  })
  
  output$reactive_txt <- renderUI({
    if (!input$name_fromlist == 'Gene (log CPM) or Gene signature') return(NULL) else {
      textInput("Not_on_list", "Write the gene or create a gene signature:",value = "CD14")
    }
  })
  output$reactive_txt2 <- renderUI({
    if (input$choose_stats == "No stats") return(NULL) else {
      textInput("compare_grp", "Comparisons to do", value="all")
    }
  })
  output$reactive_slider <- renderUI({
    if (input$choose_stats == "No stats") return(NULL) else {
      sliderInput(inputId = "bx_stat_siz", label = "Stat size", value=1.3, min=0, max=3,step=0.1)
    }
  })
  
  # # Add or not the legend to the Bxplt
  # legend_state <- reactiveVal(TRUE)
  # observeEvent(input$toggle_legend, {
  #   legend_state(!legend_state())
  #   new_label <- ifelse(legend_state(), "Legend off", "Legend on")
  #   updateActionButton(session, "toggle_legend", label = new_label)
  # })
  
  legend_state <- reactiveVal("Remove legend?")
  observeEvent(input$toggle_legend, {
    # Cycle through the three states
    new_state <- switch(legend_state(),
                        "Remove legend?" = "Remove texts too?",
                        "Remove texts too?" = "Restore legends?",
                        "Restore legends?" = "Remove legend?") 
    legend_state(new_state) 
    # Update button label based on state
    new_label <- new_state
    updateActionButton(session, "toggle_legend", label = new_label)
  })
  
  
  
  output$mod_group_bar <- renderPlot({
    print("plotting")
    notfound = ""
    # df_zmod<-as.data.frame(hm_zmod)
    # df_zmod$annot<-env_ann_df$df$annotation[match(rownames(df_zmod),env_ann_df$df$metacell )]
    df_zmod2<-df_zmod
    df_zmod2$Annotation<-env_ann_df$df$annotation[match(rownames(df_zmod),env_ann_df$df$metacell )]
    #    df_zmod2$Annotation<-df$annotation[match(rownames(df_zmod),df$metacell )]
    # df_zmod2$Annotation<-factor(df_zmod2$Annotation, levels=eval(parse(text=input$compare_grp)))
    df_zmod2$currentscore<-NA
    brplt_title=""
    print("step1")
    # Personalized palette 
    if(input$palette_1 == ""){
      palette_en_cours<-rep("black",50)
    }else{eval(parse(text=paste0("palette_en_cours<-",input$palette_1))) }
    # Personalized order
    # if(input$bx_order == ""){
    if(input$bx_order==""){
      bx_ord<-unique(df_zmod2$Annotation)
      rn_slct_dfz2<-rownames(df_zmod2)  # allows to subset the object on which we do the total Kruskal-Wallis Test 
    }else{eval(parse(text=paste0("bx_ord<-",input$bx_order)))
      rn_slct_dfz2<-rownames(df_zmod2[which(df_zmod2$Annotation%in%bx_ord),])}   
    df_zmod2$Annotation<-factor(df_zmod2$Annotation,levels =bx_ord )
    # Show legend
    # if (input$toggle_legend==TRUE) {
    #   lgd_status="right"
    # }else{lgd_status="none"}
    print("step 2")
    if (input$name_fromlist == 'Gene (log CPM) or Gene signature'){
      var_select <- input$Not_on_list 
      ilist<-gene_input2vect(var_select,dgcmtxcounts)
      print(ilist)
      if(length(ilist)==1){
        print("ilist est de 1")
        df_zmod3<-df_zmod2[rn_slct_dfz2,] 
        eval(parse(text=paste0("df_zmod3$currentscore<-as.vector(preloco[ilist,rn_slct_dfz2,drop=F])")))
        brplt_title<-paste0("projection of ",ilist," expression in MoMacs metacells")
        stat_results <- compare_means(currentscore ~ Annotation, data=df_zmod3, method = "wilcox.test", p.adjust.method = "bonferroni")
        
        clr_ln_score=NA
        if (input$choose_stats!="No stats") {
          KW_rstatix<-df_zmod3 %>%rstatix::kruskal_test(currentscore ~ Annotation)
          pwc_label <- bquote(paste("pwc: ", bold("Wilcoxon test"),"; p.adjust: ", bold("Bonferroni")) )  #pwc= PairWise Comparison
          test_label <- get_test_label(KW_rstatix, detailed = TRUE)
          combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))
          
          results_wilcoxon<- df_zmod3 %>% group_by(Annotation) %>%
            summarise(p_value = wilcox.test(currentscore, mu = mean(df_zmod3$currentscore))$p.value, .groups = "drop") %>%
            mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
          
          if(input$compare_grp=="all"){   # If comparison=all, take all shown groups to compare them together
            # eval(parse(text=paste0("cmp_grp<-",input$bx_order)))
            # eval(parse(text=paste0("cmp_grp<-",bx_ord)))
            # cmp_grp<-factor(cmp_grp,levels=cmp_grp)
            # my_comparisons <- combn(cmp_grp, 2, simplify = FALSE)
            my_comparisons <- combn(bx_ord, 2, simplify = FALSE)
          }else if(unlist(strsplit(input$compare_grp,split=","))%in%bx_ord){
            # }else if(unlist(strsplit(input$compare_grp,split=","))%in%eval(parse(text=bx_ord))){
            # }else if(unlist(strsplit(input$compare_grp,split=","))%in%eval(parse(text=input$bx_order))){ # If comparison= written names of subtypes seen in shown groups, compare them together
            # eval(parse(text=paste0("cmp_grp<-",unlist(strsplit(input$compare_grp,split=",")) )))
            cmp_grp<-unlist(strsplit(input$compare_grp,split=",")) 
            my_comparisons <- combn(cmp_grp, 2, simplify = FALSE)
          }else{ # else, take it as a raw code text input, for instance:   list(c("Mono-like 1","Mono-like 2"),c("Mono-like 1","Mono-like 3"))
            # on.exit( eval(parse(text=paste0("my_comparisons<-",input$compare_grp)))  )
            eval(parse(text=paste0("my_comparisons<-",input$compare_grp)))
          }
          
          if (input$choose_stats=="Numeric stats") {
            grob_add <- grobTree(textGrob(label = paste("Kruskal-Wallis, p =", stat_results$p.format), x= 0.1,  y=1.95, hjust=0, #put too high to hide it but keep the code if needed
                                          gp=gpar(col="#222222", fontsize=rel(7)*input$bx_stat_siz, fontface="italic")))
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title),subtitle = combined_label)+annotation_custom(grob_add)+stat_compare_means(comparisons = my_comparisons,
                                                                                                                    size = 4*input$bx_stat_siz,bracket.size = rel(0.5)*input$bx_stat_siz,
                                                                                                                    vjust = input$vjust_stat,tip.length = 0)+scale_x_discrete(limits=bx_ord)
            print("p4 done")
          }else if (input$choose_stats=="Symbol stats") {
            grob_add <- grobTree(textGrob(label = paste("Kruskal-Wallis, p =", stat_results$p.format), x= 0.1,  y=1.95, hjust=0, #put too high to hide it but keep the code if needed
                                          gp=gpar(col="#222222", fontsize=rel(7)*input$bx_stat_siz, fontface="italic")))
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title),subtitle = combined_label)+annotation_custom(grob_add)+stat_compare_means(comparisons = my_comparisons,
                                                                                                                    size = 4*input$bx_stat_siz,bracket.size = rel(0.5)*input$bx_stat_siz,
                                                                                                                    vjust = input$vjust_stat,tip.length = 0, 
                                                                                                                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                                                                       symbols = c("****", "***", "**", "*", "ns")),
                                                                                                                    hide.ns=TRUE)+scale_x_discrete(limits=bx_ord)                                  # En réordonnant on corrige le bazar de stat_compare_means qui change les annotation de X
            print("p4 done")
          }else if (input$choose_stats=="Numeric stats - vs mean") {
            clr_ln_score="#444444"
            results_wilcoxon$p_adj<-format(results_wilcoxon$p_adj,scientific=TRUE, digits=3)
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title))+
              geom_text(data = results_wilcoxon, aes(x = Annotation, y = max(df_zmod3$currentscore) + 0.5*input$bx_stat_siz,
                                                     label = p_adj),
                        size = 2*input$bx_stat_siz, color="black", vjust = 0)+scale_x_discrete(limits=bx_ord)
            
            print("p4 done") 
          }else if (input$choose_stats=="Symbol stats - vs mean") {
            clr_ln_score="#444444"
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title))+
              geom_text(data = results_wilcoxon, aes(x = Annotation, y = max(df_zmod3$currentscore) + 0.5*input$bx_stat_siz,
                                                     label = get_significance_label(p_adj)),
                        size = 5*input$bx_stat_siz, color="black", vjust = 0)+scale_x_discrete(limits=bx_ord)
            
            print("p4 done")
          }
          
        }else{
          plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
            ggtitle(paste0(brplt_title))
          notfound = ""
          print("p4 done")}
        
      }else if(length(ilist)>=2){ 
        print("ilist is more than 1")
        values<- scale( Matrix::colSums(preloco_nc[ilist,rn_slct_dfz2,drop=F])/allsum[rn_slct_dfz2] )[,]  # Need [,] to drop the attributes of scale(x) so it works for after
        df_zmod3<-df_zmod2[rn_slct_dfz2,]  
        eval(parse(text=paste0("df_zmod3$currentscore[match(rownames(df_zmod3),names(values))]<-values")))
        brplt_title<-paste0("Current module score in MoMacs metacells")
        stat_results <- compare_means(currentscore ~ Annotation, data=df_zmod3, method = "wilcox.test", p.adjust.method = "bonferroni")
        
        clr_ln_score="black"
        if (input$choose_stats!="No stats") {
          KW_rstatix<-df_zmod3 %>%rstatix::kruskal_test(currentscore ~ Annotation)
          pwc_label <- bquote(paste("pwc: ", bold("Wilcoxon test"),"; p.adjust: ", bold("Bonferroni")) )  #pwc= PairWise Comparison
          test_label <- get_test_label(KW_rstatix, detailed = TRUE)
          combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))
          
          results_wilcoxon<- df_zmod3 %>% group_by(Annotation) %>%
            summarise(p_value = wilcox.test(currentscore, mu = mean(df_zmod3$currentscore))$p.value, .groups = "drop") %>%
            mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
          
          if(input$compare_grp=="all"){   # If comparison=all, take all shown groups to compare them together
            my_comparisons <- combn(bx_ord, 2, simplify = FALSE)
          }else if(unlist(strsplit(input$compare_grp,split=","))%in%bx_ord){
            cmp_grp<-unlist(strsplit(input$compare_grp,split=",")) 
            my_comparisons <- combn(cmp_grp, 2, simplify = FALSE)
          }else{ # else, take it as a raw code text input, for instance:   list(c("Mono-like 1","Mono-like 2"),c("Mono-like 1","Mono-like 3"))
            eval(parse(text=paste0("my_comparisons<-",input$compare_grp)))
          }
          
          if (input$choose_stats=="Numeric stats") {
            grob_add <- grobTree(textGrob(label = paste("K-W, p =", stat_results$p.format), x= 0.1,  y=1.95, hjust=0, #put too high to hide it but keep the code if needed
                                          gp=gpar(col="#222222", fontsize=rel(7)*input$bx_stat_siz, fontface="italic")))
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title),subtitle = combined_label)+annotation_custom(grob_add)+stat_compare_means(comparisons = my_comparisons,
                                                                                                                    size = 4*input$bx_stat_siz,bracket.size = rel(0.5)*input$bx_stat_siz,
                                                                                                                    vjust = input$vjust_stat,tip.length = 0)+scale_x_discrete(limits=bx_ord)
            print("p4 done")
          }else if (input$choose_stats=="Symbol stats") {
            grob_add <- grobTree(textGrob(label = paste("K-W, p =", stat_results$p.format), x= 0.1,  y=1.95, hjust=0, #put too high to hide it but keep the code if needed
                                          gp=gpar(col="#222222", fontsize=rel(7)*input$bx_stat_siz, fontface="italic")))
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title),subtitle = combined_label)+annotation_custom(grob_add)+stat_compare_means(comparisons = my_comparisons,
                                                                                                                    size = 4*input$bx_stat_siz,bracket.size = rel(0.5)*input$bx_stat_siz,
                                                                                                                    vjust = input$vjust_stat,tip.length = 0, 
                                                                                                                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                                                                       symbols = c("****", "***", "**", "*", "ns")),
                                                                                                                    hide.ns=TRUE)+scale_x_discrete(limits=bx_ord)                                  # En réordonnant on corrige la merde de stat_compare_means qui change les annotation de X
            print("p4 done")
          }else if (input$choose_stats=="Numeric stats - vs mean") {
            results_wilcoxon$p_adj<-format(results_wilcoxon$p_adj,scientific=TRUE, digits=3)
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title))+
              geom_text(data = results_wilcoxon, aes(x = Annotation, y = max(df_zmod3$currentscore) + 0.5*input$bx_stat_siz,
                                                     label = p_adj),
                        size = 2*input$bx_stat_siz, color="black", vjust = 0)+scale_x_discrete(limits=bx_ord)
            
            print("p4 done") 
          }else if (input$choose_stats=="Symbol stats - vs mean") {
            plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
              ggtitle(paste0(brplt_title))+
              geom_text(data = results_wilcoxon, aes(x = Annotation, y = max(df_zmod3$currentscore) + 0.5*input$bx_stat_siz,
                                                     label = get_significance_label(p_adj)),
                        size = 5*input$bx_stat_siz, color="black", vjust = 0)+scale_x_discrete(limits=bx_ord)
            
            print("p4 done")
            print(input$bx_stat_siz)
          }
        }else{
          plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
            ggtitle(paste0(brplt_title))
          notfound = ""
          print("p4 done")}
      }else{
        notfound = paste("Incorrect input")
        plot_w4<-ggplot() + 
          annotate("text", x = 4, y = 25, size=10, label = notfound) + 
          theme_void()
        print("p4 done")
      }
      
      
      
    }else{
      print("need var")
      var_select <- input$name_fromlist
      df_zmod3<-df_zmod2[rn_slct_dfz2,] 
      brplt_title<-paste0("projection of ",var_select," in MoMacs metacells")
      #NO STATS because it would be too complicated to show
      clr_ln_score=NA
      plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = !!sym(paste0(var_select)),fill=Annotation, color=Annotation ))+
        ggtitle(paste0(brplt_title)) 
      notfound = ""
      print("p4 done")
    }
    
    print("step end done, now print")
    
    
    
    
    if (notfound == "Incorrect input"){
      print(plot_w4)
      print(notfound)
      vals$plt1<-plot_w4
      
    }else{
      print(notfound)
      palette_en_cours2<-palette_en_cours[bx_ord]
      palette_en_cours2[which(is.na(palette_en_cours2))]<-"black"
      plot_w5<-plot_w4+
        geom_boxplot(width=0.8, size=0.95*input$bx_siz,outlier.alpha = 0.5)+ #fill="white",
        xlab("Cluster")+theme_bw()+
        scale_fill_manual(values = alpha(palette_en_cours2, .4)) +
        scale_color_manual(values = palette_en_cours2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.85,size = rel(1)),
              axis.title.x =element_text(size=rel(1),face="bold"),
              #axis.text.y = element_text(size = rel(1.5)),
              axis.text.y = element_text(size = rel(input$y_txt_size)),
              #axis.title.y =element_text(size=rel(1),face="bold"),
              axis.title.y=element_blank(),
              plot.title = element_text(size = rel(1.5), face = "bold"),
              plot.subtitle = element_text(size = rel(1)),
              panel.border = element_blank(),
              legend.key.size = unit(1*input$bx_siz,'cm'),
              legend.position = "right",
              legend.text = element_text(size=rel(1)*input$bx_siz),
              axis.line = element_line(colour = "black",size = 0.6)) +
        geom_jitter(height = 0,size=0.5*input$bx_siz, width = 0.2,alpha=0.5) +
        geom_hline(yintercept = mean(df_zmod3$currentscore), linetype = "dashed",color=clr_ln_score, size=rel(0.7)*input$bx_siz) 
      
      
      if (legend_state() == "Restore legends?") {
        plot_w5 <- plot_w5 + theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_blank(),
                                   plot.title = element_blank())
      } else if (legend_state() == "Remove texts too?") {
        plot_w5 <- plot_w5 + theme(legend.position = "none")
      }
      
      print(plot_w5) 
      vals$plt1<-plot_w5
      
    }
    
    
    
    
    
  },height=bxhei, width=bxwid)
  
  
  
  vals <- reactiveValues(plt1=NULL)
  
  output$export_barplt = downloadHandler(
    
    filename = function() {paste0("plots",Sys.Date(),".pdf")},
    content = function(file) {
      pdf(file, onefile = TRUE, width = 0.016*input$bx_wid, height = 0.016*input$bx_hei)
      print(vals$plt1)
      
      dev.off()
    })
  
  
  
}



shinyApp(ui = ui, server = server)





## Easy to copy cheatsheet:

# HM Cluster order:
# 6,40,2,33,39,38,37,34,35,36,3,9,11,25,4,13,7,8,5,10,24,26,27,28,12,14,16,18,1,19,20,23,29,30,32,15,17,21,22,31

# Boxplot color palette + order
# c('Mono-like 1A'='#874037','Mono-like 1B'='#A3672C','Mono-like 2'='#C3A34B','Mono-like 3'='#f77774','Mono-like 4'='#ffa3a1','Mono-like 5'='#B4DEC6','Macrophage-like 6'='#c967eb','Macrophage-like 7'='#5fc2ed','Macrophage-like 8'='#4F88B9','Macrophage-like 9'='#5C538B')
# c('Mono-like 1A','Mono-like 1B','Mono-like 2','Mono-like 3','Mono-like 4','Mono-like 5','Macrophage-like 6','Macrophage-like 7','Macrophage-like 8','Macrophage-like 9')

# Boxplot comparisons
# list(c('Mono-like 1','Mono-like 2'),c('Mono-like 1','Mono-like 3'),c('Mono-like 1','Mono-like 4'),c('Mono-like 1','Mono-like 5'), c('Mono-like 1','Macrophage-like 6'),c('Mono-like 1','Macrophage-like 7'),c('Mono-like 1','Macrophage-like 8'),c('Mono-like 1','Macrophage-like 9'))
# list(c('Mono-like 1','Mono-like 3'),c('Mono-like 1','Mono-like 4'))


# Regulons gene lists:
# regulons <- readRDS("./2.6_regulons_asGeneSet.Rds")
# paste0(regulons$HIF1A,collapse=",")

# SPI1 : AAGAB,AAK1,AAR2,ABCB7,ABCF1,ABHD17A,ABR,ABTB1,ACAA1,ACAP2,ACBD5,ACP5,ACTB,ACTG1,ACTR1A,ACVR1,ADGRG5,ADIPOR2,ADO,ADSL,AGPAT1,AGPAT3,AHCY,AIF1,AKIRIN2,AKT2,ALDOA,ANAPC15,ANGEL1,ANKRD13A,ANXA1,ANXA4,AP1M1,AP2A1,AP2B1,AP3D1,APEX2,APH1A,APLP2,APOM,APP,AQR,ARAP1,ARCN1,ARF1,ARF5,ARF6,ARFIP1,ARHGAP1,ARHGAP19,ARHGAP27,ARHGAP30,ARHGAP9,ARHGDIA,ARHGDIB,ARHGEF1,ARHGEF12,ARHGEF3,ARHGEF6,ARL6IP5,ARL6IP6,ARMCX5,ARPC1A,ARPC1B,ARPC2,ARPC3,ARPC4,ARRB2,ASCC2,ASGR1,ASGR2,ATF5,ATG12,ATG5,ATOH8,ATP11C,ATP13A1,ATP6V0D1,ATP6V0E1,ATP6V1A,ATP6V1E1,ATXN7L1,AURKAIP1,B3GAT3,B9D2,BAK1,BBIP1,BCAP31,BCL2L14,BCL3,BDH2,BECN1,BIN3,BOD1L1,BORCS8,BROX,BTBD10,BTF3,BTK,BUB3,BZW1,C14orf119,C16orf91,C17orf49,C18orf21,C19orf25,C19orf38,C1orf162,C1RL,C3AR1,C5AR1,C6orf89,CAB39,CACYBP,CALHM2,CALM2,CAMTA1,CAP1,CAPG,CAPN1,CAPN15,CAPZA2,CARD9,CARS,CASP4,CASP8,CASP8AP2,CAST,CBFA2T3,CBX3,CCAR2,CCDC12,CCDC126,CCDC61,CCNG1,CCNI,CCR1,CCRL2,CCT2,CCT7,CCZ1,CCZ1B,CD14,CD163,CD2BP2,CD300E,CD300LF,CD33,CD36,CD4,CD48,CD53,CD63,CD68,CDC123,CDC40,CDC42,CDH11,CDIPT,CDK10,CDK9,CELF1,CELF6,CENPQ,CERS6,CFL1,CFLAR,CFP,CGGBP1,CHCHD2,CHCHD5,CHID1,CHMP1A,CHMP1B,CHMP2A,CHMP2B,CHP1,CHPF,CITED2,CLEC4A,CLEC4E,CLIC1,CLINT1,CLIP4,CLN3,CLP1,CLPP,CLPTM1L,CLTA,CLTB,CLTC,CMAS,CMC1,CNDP2,CNOT6L,CNOT7,COA5,COA6,COMMD1,COMMD2,COMMD6,COPA,COPB2,COPE,COPS7A,COQ6,CORO1B,COX14,COX17,COX5B,COX6A1,COX7C,COX8A,CPSF3,CPT2,CREB3,CREM,CRIPT,CRTAM,CSK,CSNK1D,CSRP1,CSTA,CSTF2,CTBS,CTDSP1,CTNNA1,CTNNBL1,CTU1,CUL4B,CXCR5,CXXC5,CYB561A3,CYB5R4,CYBA,CYBB,CYSTM1,DAD1,DCAF11,DCAF12,DCAF13,DCAF7,DCTN2,DCTN3,DDAH2,DDT,DDX18,DDX21,DDX50,DEK,DENND1C,DENND4B,DERL1,DET1,DGUOK,DHDDS,DHRS4L2,DHRS7,DHX34,DHX8,DIAPH1,DIAPH2,DNAJC1,DNAJC10,DNAJC14,DNAJC18,DNASE2,DNM1,DOCK8,DOK3,DPH2,DPM1,DPM3,DPP8,DPYSL2,DR1,DTD2,DYNC1I2,DYNLL2,DYNLT1,E2F4,EBAG9,EDEM2,EFCAB14,EFHD2,EHBP1L1,EIF2AK1,EIF2B2,EIF2S2,EIF2S3,EIF3A,EIF3D,EIF3G,EIF3H,EIF4A1,EIF4E2,EIF4G2,EIF4H,EIF5AL1,EIF6,ELAC2,ELOVL1,ELP5,EMC4,EMC7,EMILIN2,ENG,ENO1,EPN1,EPS15,EPSTI1,ERBB2,ERMAP,ERP27,ERP44,ESRRA,ETF1,ETFB,ETV6,EVI2B,EVI5,EXOC6,EXOSC3,EXOSC4,F11R,F13A1,FAAP20,FAF2,FAM107A,FAM110A,FAM120AOS,FAM172A,FAM174A,FAM199X,FAM200A,FAM32A,FAM49B,FAM89B,FAM91A1,FARSA,FAU,FBP1,FBXL15,FBXL18,FBXO28,FBXO36,FBXW11,FBXW5,FCER1G,FCGR1A,FCGR1B,FCHO1,FCN1,FDPS,FERMT3,FES,FGD2,FGR,FHL3,FIBP,FKBP1A,FMN1,FMNL1,FOLR3,FOXJ3,FPGS,FRG1,FSIP1,FUBP3,FURIN,FUT4,FXYD5,FZR1,G6PD,GABARAP,GABARAPL2,GABPA,GALE,GAPDH,GAPT,GAS7,GATB,GBA,GDAP1,GDI2,GEMIN2,GEMIN8,GGA1,GIMAP4,GIN1,GIT2,GLMP,GLUL,GLYCTK,GMFG,GMPPB,GMPR2,GNA13,GNAI2,GNAI3,GNB2,GNB4,GNG11,GNGT2,GNL2,GNPAT,GOLGA7,GOLPH3,GORASP2,GOSR2,GPAA1,GPANK1,GPD2,GPR108,GPS2,GPSM3,GRB2,GRIPAP1,GSKIP,GSN,GSPT1,GTF2A2,GTF2F2,GTF3C5,GTPBP4,GTSF1,HARS2,HAT1,HCCS,HCK,HCLS1,HCST,HDAC2,HERPUD2,HEXB,HIST1H4E,HK1,HK3,HLA-B,HLA-F,HM13,HMGN4,HNRNPA3,HNRNPD,HNRNPF,HNRNPH2,HNRNPK,HOPX,HOXB-AS1,HPCAL1,HPF1,HRH2,HSBP1,HSP90B1,HSPA4,HVCN1,HYAL2,ICAM3,IDH3B,IDI1,IFI35,IFITM2,IFITM3,IGFLR1,IGSF6,IKBKG,IKZF2,IL16,IL17RA,IL1R2,IL4I1,ILK,IMP3,IMPDH1,ING4,INIP,INO80C,INTS1,INTS9,IPO4,IPO7,IQSEC2,IRF2,IRF5,ISCU,ISY1,ITGA4,ITGAL,ITGAM,ITGB1,ITGB1BP1,ITGB2,JADE1,JAK2,JAK3,JOSD2,KAT2B,KATNB1,KCMF1,KCNAB2,KCNK13,KCTD12,KCTD5,KIAA0040,KIF3B,KLF7,KLHL18,KPNB1,KRTCAP2,KTN1,L3MBTL2,LAIR1,LAMP1,LAMTOR2,LAMTOR3,LAPTM5,LARP7,LAS1L,LCP1,LEMD2,LEPROT,LEPROTL1,LGALS9,LILRA1,LILRA6,LILRB1,LILRB2,LILRB3,LIMK1,LIMS1,LIN7A,LINC00167,LLPH,LMNB1,LPXN,LRP10,LRRC25,LSM12,LSM3,LSM5,LST1,LTA4H,LTB4R,LTBR,LTV1,LUC7L2,LY96,LYL1,LYN,LYPLA2,LYRM1,LYZ,LZTFL1,M6PR,MAGED2,MAGEH1,MAGOH,MAN1A2,MANF,MAP3K11,MAP3K3,MAP3K8,MAP7,MAPKAP1,MAPKAPK3,MARCH5,MARK2,MBD1,MBNL1,MCRS1,MCTS1,MED10,MED23,MED27,MED31,MED7,MED8,MEPCE,METTL17,METTL2A,METTL3,METTL7A,MFAP1,MFF,MFSD11,MGAT1,MID1IP1,MIEF1,MIEN1,MIR497HG,MIR99AHG,MLF2,MLXIP,MMP2,MMP3,MMP7,MNDA,MOB1A,MOB2,MON1B,MORF4L1,MOSPD3,MPEG1,MPG,MPLKIP,MR1,MRPL18,MRPL27,MRPL28,MRPL3,MRPL33,MRPL43,MRPL51,MRPL52,MRPS10,MRPS21,MS4A7,MSANTD2,MSN,MTA2,MTIF3,MTMR11,MTMR14,MTPN,MVB12A,MVP,MYL12A,MYL12B,MYO1F,N4BP2L2,NAA20,NACA,NADK,NARS,NAT1,NCEH1,NCF2,NCF4,NCSTN,NDUFA13,NDUFA2,NDUFAF3,NDUFS3,NECAP2,NEK6,NELFB,NELFE,NFAM1,NFATC3,NFIX,NFS1,NFYC,NGRN,NIF3L1,NKAP,NKG7,NLK,NLRC4,NME1,NME2,NMI,NOMO1,NOP10,NOSIP,NPRL2,NPTN,NR1H2,NR2F6,NRBF2,NRDC,NRG1,NSMCE1,NTAN1,NUB1,NUBP2,NUDCD3,NUMA1,NUP214,OAS1,OCIAD1,OGFOD2,OGFR,OGG1,ORC2,ORC4,ORMDL2,OSBPL11,OSM,OST4,OSTM1,OTUD5,P2RX1,P4HB,PABPC1,PACSIN2,PARN,PARP10,PARVG,PCED1B-AS1,PCIF1,PCYT1A,PDAP1,PDCD10,PDCD6,PDE8B,PDIA6,PDLIM2,PDRG1,PECAM1,PEMT,PER1,PES1,PET100,PEX16,PEX26,PFDN5,PFN1,PGAM5,PGK1,PGRMC1,PHB2,PHF5A,PHKA2,PICALM,PIGC,PIGM,PIGT,PIH1D1,PIK3AP1,PIK3C3,PIK3CD,PIK3R2,PILRA,PIN4,PITPNM1,PLA2G7,PLRG1,PLSCR1,PLXNC1,PMP22,PMPCA,PMVK,POLDIP2,POLG,POLQ,POLR2H,POLR2J,POLR3D,POMP,POP5,POU2F2,PPCDC,PPM1G,PPP1CA,PPP1CB,PPP1CC,PPP1R11,PPP1R12C,PPP1R18,PPP2R2A,PPP2R5C,PPP2R5E,PPP3R1,PPP4C,PPP6R1,PRDM1,PRDX1,PRDX3,PRELID1,PRIM2,PRKAG1,PRKAR1A,PRKCD,PRKCH,PRKCSH,PRPF19,PRPF38A,PRPF40A,PRPF6,PRPS1,PRR13,PRRT3,PRUNE2,PSMA1,PSMA4,PSMA5,PSMA6,PSMB1,PSMB10,PSMB3,PSMB4,PSMB5,PSMB7,PSMB8,PSMC1,PSMC5,PSMD13,PSMD14,PSMD7,PSME1,PSMF1,PSRC1,PSTPIP1,PTBP3,PTP4A3,PTPN6,PTPN7,PTPRC,PTPRE,PTPRO,PTRH2,PTX3,PUS3,PXN,PYCARD,PYGL,PYROXD1,QPCT,RAB11A,RAB11B,RAB11FIP3,RAB14,RAB18,RAB20,RAB27A,RAB28,RAB2A,RAB32,RAB34,RAB35,RAB4B,RAB5A,RAB5B,RAB5C,RAB7A,RAB9A,RAC1,RAC2,RACK1,RAD23B,RAE1,RAP1GAP2,RARS,RASAL3,RASGRP4,RASSF2,RASSF5,RBM4,RBM42,RBM8A,RBMS1,REEP4,RER1,RETN,RFC1,RFK,RFXANK,RGS10,RGS19,RHOG,RIN3,RLIM,RNASE6,RNASEH1,RNASEH2A,RNASEH2C,RNASEK,RNF121,RNF175,RNF181,RNF185,RNF31,RNF6,ROCK1,RPF1,RPF2,RRP15,RTN1,RTN4,RTN4RL2,RUFY1,S100A4,S100A6,S100A9,S1PR4,SAMSN1,SAP30BP,SAR1B,SARNP,SART1,SASH3,SCAMP2,SCAMP3,SCNM1,SCO1,SCRN3,SDAD1,SDF2,SDHAF1,SDHAF2,SEC14L1,SEC24C,SEC31A,SEC61A1,SELPLG,SEMA4A,SERPING1,SF1,SF3A1,SF3B2,SF3B4,SF3B5,SFT2D1,SGTA,SH3BGRL,SH3BGRL3,SH3KBP1,SH3TC1,SHISA5,SHROOM1,SIGLEC12,SIGLEC14,SIGLEC5,SIGLEC7,SIPA1L3,SIRPB1,SIRPD,SIRT6,SKIL,SLA,SLC1A5,SLC25A11,SLC25A14,SLC25A24,SLC25A46,SLC2A9,SLC30A7,SLC30A9,SLC33A1,SLC38A10,SLC39A1,SLC39A11,SLC39A7,SLC39A9,SLC44A1,SLC7A7,SLX4IP,SMAP2,SMARCAD1,SMCO4,SMIM15,SMIM7,SMU1,SMUG1,SNAP23,SNF8,SNHG20,SNHG21,SNHG6,SNIP1,SNRNP27,SNRPB,SNRPB2,SNRPD2,SNRPE,SNUPN,SOCS3,SORT1,SPATA24,SPEN,SPI1,SPR,SPTLC1,SRBD1,SRF,SRFBP1,SRGN,SRP14,SRP19,SRP54,SRP72,SRPRA,SSBP1,ST3GAL2,STAB1,STARD3,STK38L,STK40,STOML2,STRAP,STX18,STXBP2,SUB1,SUMO1,SUN2,SUPT4H1,SYK,SYS1,SYTL3,SYVN1,SZRD1,TADA3,TAF10,TAF11,TAF1B,TAGAP,TAGLN2,TANK,TAPBPL,TAX1BP1,TBCC,TBRG1,TBXAS1,TCEAL8,TCF21,TCF25,TCP11L1,TEFM,TFE3,TFEB,TFG,TGFB1,TGIF1,TGIF2,TGM2,THAP1,THAP4,THAP5,THBD,THBS3,THOC2,THOC5,TIAL1,TIMM10,TIMM10B,TIMM23,TIMM8B,TIMMDC1,TINF2,TJP1,TLN1,TLR8,TM2D2,TM2D3,TM9SF1,TMA7,TMBIM1,TMBIM6,TMC8,TMED8,TMEM101,TMEM102,TMEM123,TMEM126B,TMEM161B-AS1,TMEM167A,TMEM167B,TMEM173,TMEM199,TMEM208,TMEM222,TMEM230,TMEM234,TMEM243,TMEM33,TMEM59,TMUB2,TMX2,TNFAIP8L2,TNFSF13,TNFSF13B,TNNI2,TOLLIP,TOMM5,TOR4A,TOX4,TP53I11,TP53RK,TPD52L2,TPM2,TPM3,TPM4,TPRG1L,TPRKB,TPRN,TRAF3IP3,TRAF7,TRAM1,TRAM2-AS1,TRAPPC1,TRAPPC3,TRAPPC4,TRIM5,TRIM8,TRIOBP,TRIP12,TRMT1,TSC22D4,TSG101,TSPAN14,TSPO,TSPOAP1,TSSC4,TTC1,TTC30B,TUBGCP2,TUSC2,TWF2,TXN2,TXNDC17,TXNDC9,TYROBP,UBA1,UBAC1,UBALD2,UBB,UBE2D3,UBE2E1,UBE2L3,UBE2L6,UBE2M,UBE2Q1,UBE2V1,UBE4A,UBL3,UBL5,UBOX5,UBXN2B,UBXN4,UCP2,UGDH,UGGT1,UHRF1,UMAD1,UNC45A,UQCRQ,URM1,USB1,USF1,USP1,USP30,UTP11,UTP18,UTP3,UVSSA,VASP,VAV1,VCAN,VIPAS39,VMA21,VPS16,VPS25,VPS26A,VPS45,VPS4B,VPS52,VRK3,VTA1,WAC,WARS,WAS,WASF2,WBP1L,WBP2,WDR45,WDR48,WDR55,WDR70,WDR83OS,WIPF1,WNT5A,WRAP73,WSB2,XRN2,YBX1,YES1,YIPF3,YIPF4,YKT6,YWHAB,YWHAE,YWHAZ,ZBTB45,ZBTB8OS,ZC3H10,ZC3H13,ZC3HC1,ZCCHC9,ZDHHC5,ZFC3H1,ZGPAT,ZKSCAN4,ZMAT2,ZMPSTE24,ZNF143,ZNF213,ZNF264,ZNF317,ZNF330,ZNF333,ZNF35,ZNF362,ZNF384,ZNF512,ZNF513,ZNF524,ZNF622,ZNF710,ZNF780A,ZNHIT1,ZRANB2,ZSCAN16-AS1
# CEBPB : ACTG1,ACTN4,ACVR1B,ADGRL2,ADIPOR1,ADIPOR2,AGPAT4,ALDH1A1,ANXA1,ARHGAP24,ARID5A,ASAH1,ASAP1,ATOH8,ATP6V0C,BIN1,BMPR1A,BNIP3L,BTG1,C15orf39,CARD16,CBX3,CCDC47,CCNI,CDKN1A,CEBPB,CHP1,CHSY1,CLEC4E,CMIP,CREB5,CSF3R,CSNK1G2,CTDSP1,CTNNA1,CXXC5,DDIT3,DDIT4,DEK,DOK3,DPYSL2,DUSP4,EED,EIF1,EIF4A1,EIF4EBP1,EIF5,ELF2,EMC2,EMD,EPB41L3,EPHB6,ESRRA,ETV6,FAM110B,FAM177A1,FAM49B,FBXO3,FGD5-AS1,FGR,FLOT1,FLOT2,FLVCR2,FMNL1,FOSB,FOSL2,FURIN,G0S2,GAPDH,GCA,GHITM,GLUL,GNB2,GPR85,GRB2,GYG1,HDAC4,HERPUD1,HMGN2,HNRNPA3,HPSE,IER2,IGFBP2,IRF2,ISG20,ITGB2,JAK1,JDP2,KDM6B,KIAA0040,LATS2,LMNB1,MBD6,MCL1,MIDN,MPP7,MSRB1,MTA2,NAA20,NACC2,NFATC1,NFE2L2,NFIL3,NUP50,OAZ1,OSM,PCED1B-AS1,PCGF5,PCTP,PDE4D,PDXK,PELI2,PER1,PFDN1,PFN1,PGK1,PGS1,PHC2,PHIP,PKM,PLBD1,PLPP5,PLPPR2,PNPLA2,PPM1B,PPP3CA,PRKAG2,PSTPIP1,PXK,QPCT,RAB1A,RAB2A,RAB6A,RAD21,RAP1A,RARA,RASGRP2,RB1,RBM15B,RBM27,RNASE2,RNF149,S100A8,S100A9,SEC23A,SEC62,SELL,SERPINA1,SEZ6,SFTPD,SH3KBP1,SKI,SLC16A3,SLC16A7,SLC25A15,SLC3A2,SLC8A1,SMARCD3,SMCR5,SNX14,SNX21,SOCS3,SSBP2,SSH2,STAU1,STX11,SULF2,SUMO2,TBL1X,TGFBI,THRB,TMED7,TMEM154,TMEM33,TMEM91,TNFRSF1A,TNFSF13B,TOP1,TP53TG1,TREM1,TRIB1,TRIM69,TRIOBP,TSPO,TXNDC12,UBQLN1,USP3,VAMP5,VGLL4,YWHAG,ZADH2,ZDHHC13,ZDHHC2,ZFP36L1,ZFP36L2,ZMIZ1,ZNF281,ZNF695
# NFKB1 : AATK,ABCA1,ABCA5,ABHD17C,ABL2,ACAT2,ACOT9,ACSL1,ACSL4,ACSL5,ACTN1,ACTR3,ACVR2A,ADA,ADCY9,ADORA2A,ADRA2B,ADRB2,ADSS,ADTRP,AEN,AGAP3,AGPAT4,AHCYL1,AKAP13,AKIRIN1,AKT3,ALCAM,AMPD3,ANKRD12,ANKRD28,ANKRD33B,ANTXR2,ANXA2,ANXA5,ANXA7,AP1S2,APOL6,AQP9,ARAP3,AREL1,ARF1,ARFGAP3,ARHGAP21,ARHGAP24,ARHGAP26-IT1,ARHGAP31,ARHGAP31-AS1,ARHGEF2,ARID3A,ARID5B,ARL5B,ARL8B,ASAP1,ASAP1-IT2,ASTL,ATF2,ATF5,ATF7IP2,ATP13A3,ATP1A1,ATP1B3,ATP2B1,ATP2C1,ATP6AP1,ATP6V1B2,ATP6V1C1,ATP8B2,ATXN1,AZIN1,B3GNT2,B3GNT5,B4GALT1,B4GALT5,BAALC,BACH1,BANP,BASP1,BATF,BAZ1A,BAZ2A,BBIP1,BCAR3,BCAT1,BCL2A1,BCL2L1,BCL2L14,BCL3,BCL6,BCOR,BEND3,BHLHE40,BID,BIRC2,BIRC3,BMP6,BMT2,BRAF,BSDC1,BTBD19,BTG1,BTG2,BTG3,C15orf48,C1orf21,C21orf62,C3,C6orf223,C9orf72,CA13,CAMK1G,CAMKK2,CANT1,CANX,CASP1,CCDC71L,CCDC93,CCL20,CCL23,CCL3,CCL3L1,CCL4,CCL4L2,CCL5,CCL7,CCNT2,CCR5,CCR7,CCRL2,CD109,CD200,CD40,CD44,CD58,CD82,CD83,CD9,CDC42BPG,CDC42EP2,CDC42EP3,CDC42EP4,CDC42SE1,CDCA4,CDK1,CDK14,CDKN1A,CDV3,CFLAR,CFLAR-AS1,CHD2,CHD4,CHEK1,CHMP2B,CHMP4B,CHST11,CHST15,CKB,CLCF1,CLEC12A,CLEC2D,CLEC4D,CLIC4,CLSTN3,CMTM6,CNKSR3,CORO1C,CPAMD8,CPEB4,CRADD,CREB5,CRISPLD2,CSF2,CSF3,CSGALNACT2,CSNK1A1,CSNK1E,CSRNP1,CTNNA1,CTNNB1,CTSL,CXCL1,CXCL2,CXCL3,CXCL5,CXCL8,CXCL9,CXXC5,CYB5R2,CYFIP1,CYTH1,DARS,DCP1A,DCTN4,DDX6,DDX60L,DEFB1,DENND4A,DHX34,DIAPH1,DICER1,DIDO1,DLGAP1-AS1,DLGAP1-AS2,DMXL2,DNAJB5,DNAJB6,DRAM1,DSE,DUSP2,DUSP3,DUSP5,DUSP6,DUSP8,DVL3,DYNLT3,DYRK3,E2F7,EBF1,EBI3,EBLN2,ECE1,EDEM1,EEPD1,EFR3A,EGR3,EHD1,EHD4,EHF,EIF2AK3,EIF4E,EIF4G1,EIF5,ELF4,ELL2,ELMO2,ELOVL7,EML3,EML4,EMP2,EMP3,ENDOV,EPAS1,EPB41L3,EPC1,EPC2,EPM2AIP1,ERBIN,ERCC6,EREG,ERGIC1,ERMN,ERN1,ESPL1,ESYT2,ETS1,ETV3,ETV3L,EZH2,F12,F3,FAM107B,FAM222A,FAM3C,FAM49A,FAM83G,FAS,FBRS,FEM1B,FEM1C,FERMT2,FGFR1,FGFRL1,FGR,FKRP,FLCN,FLNA,FLT1,FMNL3,FNBP1,FNDC3A,FNDC3B,FNIP2,FOSL1,FOSL2,FOXO3,FRMD6,FRMD7,FRY-AS1,FSTL3,FTH1,FXYD6,FZD4,FZR1,G0S2,GAB2,GABPB1,GADD45A,GAK,GALNT6,GAS2L3,GBE1,GBP2,GCH1,GCNT4,GCSAML,GEM,GGA1,GHRLOS,GJB2,GJC1,GK,GK5,GLA,GLIS3,GLS,GNA15,GNB1,GNG2,GNG4,GPATCH2L,GPR132,GPR137B,GPR157,GPR183,GPR35,GPR68,GPR84,GRAMD1A,GRASP,GSAP,GTF2IRD1,HAPLN3,HAS1,HDAC4,HDGF,HECW2,HES4,HEY1,HIPK2,HIPK3,HIVEP1,HIVEP2,HLA-A,HLA-F,HMBOX1,HMGCS1,HMGXB4,HNRNPC,HOMER1,HPCAL1,HRH1,HS3ST1,HS3ST3B1,HSP90B1,HYOU1,IBTK,ICAM1,IDO2,IER3,IFFO2,IFIH1,IFIT3,IFNGR2,IL10,IL10RB,IL1A,IL1B,IL1R1,IL1RN,IL20,IL23A,IL27,IL2RA,IL2RG,IL32,IL36G,IL36RN,IL4I1,IL4R,IL6,IL6R,IL6ST,IL7,IL7R,INHBA,INPP5A,INSM1,IPO5,IQCG,IQGAP2,ISG20,ITCH,ITGA1,ITGA3,ITGA5,ITGAV,ITGAX,ITGB8,ITK,ITPKB,ITPKB-IT1,ITPRIP,ITSN2,IVNS1ABP,JAG1,JAK1,JARID2,JMJD1C,JRKL,KANK1,KCNA3,KCNJ2,KCNN4,KCNQ1OT1,KDM4B,KDM5C,KDM6A,KDM6B,KIF1B,KIF21B,KIFC3,KLF10,KLF6,KLF7,KLHL21,KMO,KMT2C,KMT2E,KPNB1,KTN1,KYNU,LACC1,LAMB3,LAMP3,LAT,LCOR,LCP1,LCP2,LEMD2,LENG8,LEP,LGALS3,LHFPL2,LIMS1,LINC-PINT,LINC00158,LINC00211,LINC00309,LINC00520,LINC00622,LINC00656,LINC00852,LINC00862,LINC00884,LINC00910,LINC00926,LINC00937,LINC00963,LINC00971,LINC01093,LINC01136,LINC01160,LINC01215,LINC01262,LINC01270,LINC01353,LINC01366,LINC01465,LINC01476,LINC01588,LIPN,LITAF,LMNA,LMNB2,LONRF1,LONRF3,LOXHD1,LPCAT1,LRCH1,LRG1,LRP12,LRRC32,LRRC8C,LRRFIP2,LUCAT1,LYN,LYPD3,MAFF,MAFG,MALT1,MAML2,MAMLD1,MANF,MAP1LC3A,MAP2K1,MAP2K3,MAP3K8,MAP4K4,MAPK13,MAPK6,MAPK8,MAPKAPK2,MARCH3,MARCKS,MBD6,MCOLN2,MDGA1,MECP2,MED13,MELTF,MEPCE,MET,METRNL,MFGE8,MFHAS1,MFSD2A,MGAM,MGLL,MIAT,MICALL2,MIDN,MIR155HG,MIR3142HG,MIR3945HG,MLLT6,MMP1,MMP14,MMP19,MMP9,MN1,MOB3B,MOB3C,MPRIP,MREG,MROH1,MS4A7,MSANTD3,MSC,MSN,MT2A,MTRF1L,MXD1,MYH9,MYO1E,MYO1G,MYO9B,N4BP1,N4BP2L1,NAB1,NAB2,NABP1,NAMPT,NANOS3,NBPF10,NBPF11,NBPF14,NBPF15,NBPF20,NBPF9,NCF2,NCOR2,NCR3LG1,NCS1,NDRG1,NECTIN1,NECTIN2,NEDD4L,NEMP1,NEURL3,NFAT5,NFATC1,NFE2L2,NFKB1,NFKB2,NFKBIA,NFKBIB,NFKBID,NFKBIE,NFKBIZ,NINJ1,NIPAL4,NLRP3,NOTCH1,NOTCH2,NR3C1,NR4A3,NRIP3,NRP2,NSMAF,NUMB,NUP153,NUP188,NUP62,NUS1,NXT2,OAF,OAZ2,OGT,OPTN,ORAI2,OSBPL8,OSTM1,OTUD4,OTUD5,P2RX4,P2RX7,P4HA2,PABPC1,PABPC4,PACERR,PAFAH1B1,PALM2-AKAP2,PANK3,PANX1,PAPLN,PATL1,PBRM1,PCNX1,PDE4A,PDE4B,PDE4D,PDE4DIP,PDE7A,PDE8A,PDGFB,PDLIM5,PDLIM7,PEA15,PEAK1,PELI1,PER2,PFKFB3,PFKP,PHC1,PHF1,PHLDA1,PHLDA2,PHLDB1,PHTF1,PIAS1,PIGA,PIK3AP1,PIK3R5,PILRB,PIM1,PIM2,PIM3,PISD,PKNOX1,PLAGL2,PLAU,PLAUR,PLD1,PLEC,PLEK,PLEKHB2,PLEKHO2,PLIN4,PLK3,PLOD2,PLPP5,PLSCR1,PLXNA1,PLXNC1,PLXND1,PMAIP1,PNPLA1,PNPLA8,PNRC2,POGZ,POLG,POR,PPARD,PPFIA1,PPIF,PPP1R12A,PPP1R15B,PPP4R2,PPP6R3,PPTC7,PRDM1,PRDM8,PRKCH,PRLR,PRR7,PRRC2C,PSEN1,PSMA3-AS1,PSTPIP2,PTBP1,PTGES,PTGIR,PTGS2,PTK2B,PTPN1,PTPN14,PTPN2,PTPRE,PTPRJ,PTX3,PVR,PWWP2A,PXDC1,PXN,QKI,QSOX2,RAB10,RAB11FIP5,RAB12,RAB21,RAB35,RAB5A,RAB7A,RAB9A,RABEP1,RABGEF1,RACGAP1,RAD1,RAD21-AS1,RALGDS,RANBP9,RAP1B,RAP2C,RAP2C-AS1,RAPGEF1,RAPGEF2,RAPH1,RASGRP3,RASSF3,RASSF5,RB1CC1,RBBP8,RBM17,RBM47,RDX,REL,RELA,RELB,RELT,RETN,REV3L,RFX2,RGCC,RGL4,RGPD2,RHCG,RHEB,RHOF,RHOH,RHOU,RILPL2,RIN2,RNF144B,RNF145,RNF19B,RNF2,RNF24,ROCK1,ROCK2,RRAD,RREB1,RUFY3,RUNX1,RYBP,S100A2,SALL4,SAMSN1,SAR1A,SAT1,SATB1,SAV1,SBF2,SBNO2,SCAF4,SCARF1,SCN1B,SDC2,SDC4,SDCBP,SEC22B,SEC23B,SEC24A,SEMA3E,SEMA4A,SERF1A,SERPINB2,SERPINB9,SERPINB9P1,SESN2,SESTD1,SETD5,SFMBT2,SFR1,SGMS1,SGMS2,SGPP2,SH2D3A,SH3BP2,SHF,SIAH1,SIAH2,SIK2,SIK3,SIPA1L1,SIRPA,SIX5,SKA3,SKIL,SLAMF1,SLC11A2,SLC12A6,SLC14A2,SLC16A6,SLC1A3,SLC24A4,SLC25A13,SLC25A34,SLC2A6,SLC30A4,SLC30A7,SLC35F2,SLC38A5,SLC39A13,SLC39A8,SLC43A2,SLC43A3,SLC44A1,SLC45A4,SLC7A1,SLC8A2,SLC8B1,SLC9A8,SLFN5,SMAD2,SMAD3,SMAD7,SMG1,SMG7,SMG7-AS1,SMG9,SMIM3,SMOX,SMS,SNN,SNRK-AS1,SNX10,SOCS2,SOCS3,SOD2,SOX13,SP3,SPAG9,SPHK1,SPRED2,SPSB1,SPTBN5,SQSTM1,SRC,SRD5A1,SREBF2,SRP54,SRSF5,SSH1,ST18,ST20-AS1,ST3GAL1,ST3GAL2,STAG2,STARD10,STARD3NL,STARD4,STARD8,STAT1,STAT4,STAT5A,STK26,STK4,STX11,STX12,STX3,STX4,STX5,SUN2,SURF4,SUSD6,SVIL,SWAP70,SYNJ2,SYNPO2,TAB2,TAF4B,TANK,TAOK1,TAOK3,TAPBP,TBC1D22B,TBC1D30,TBK1,TCHH,TET2,TEX14,TFE3,TFEC,TGIF1,THAP2,THBS1,THUMPD3-AS1,TICAM1,TIFA,TJP2,TLE3,TLE4,TLR2,TM4SF19,TMBIM1,TMBIM6,TMEM106A,TMEM120A,TMEM120B,TMEM138,TMEM217,TMEM26,TMEM54,TMEM63B,TMEM88,TMEM8A,TMTC1,TNF,TNFAIP2,TNFAIP3,TNFAIP6,TNFAIP8,TNFRSF10B,TNFRSF14,TNFRSF18,TNFRSF1B,TNFRSF4,TNFRSF9,TNFSF14,TNFSF15,TNIK,TNIP1,TNIP2,TNIP3,TNKS2,TNS3,TOMM40L,TP53BP1,TP53BP2,TP53INP2,TPRA1,TRA2A,TRAF1,TRAF3,TRAF3IP2,TRAF3IP2-AS1,TRAF4,TRAK1,TRAPPC10,TRG-AS1,TRIM36,TRIP10,TSC22D2,TTC39A,TTL,TTYH2,TUBB6,UAP1,UBA1,UBAC2,UBALD2,UBASH3B,UBE2D1,UBE2E1,UBE2J1,UBE2Z,UBR4,UBXN7,UPB1,UPP1,URGCP,USP12,USP12-AS2,USP13,USP16,USP27X,USP27X-AS1,VASP,VAV2,VDR,VEZF1,VPS13A,VPS37C,WAPL,WDR1,WDR31,WIPF2,WTAP,XBP1,XIAP,YY1AP1,ZBTB10,ZBTB25,ZBTB5,ZBTB7A,ZC3H12A,ZC3H12C,ZC3H18,ZCCHC14,ZEB2,ZFHX3,ZFP36,ZFX,ZHX2,ZMIZ1,ZMIZ1-AS1,ZMIZ2,ZNF140,ZNF217,ZNF267,ZNF277,ZNF281,ZNF295-AS1,ZNF316,ZNF410,ZNF697,ZNFX1,ZSCAN5A,ZSWIM6
# RELA : ABCA1,ABCC11,ACOD1,ACSL1,ACSL5,ACTN1,ACVR2A,ADAMTSL4-AS1,ADARB1,ADD2,ADM,ADNP,ADPRH,AK4,AKIRIN1,AKT1,ALDOA,AMOTL1,AMPH,ANKRD12,ANKRD22,ANKRD28,ANPEP,ANTXR2,ANXA2,AP1G1,AP1S2,APOBEC3A,APOBR,APOL3,AQP9,ARAP2,ARCN1,ARHGAP26,ARHGAP29,ARHGEF10L,ARHGEF11,ARID3A,ARNTL2,ASGR2,ATP11C,ATP1A1,ATP2B4,ATP6V1C1,ATP8B2,B2M,B3GAT2,B3GNT2,B3GNT5,B4GALT5,BACH1,BARD1,BATF,BATF2,BAZ1A,BBIP1,BCAT1,BCL2A1,BCL3,BCL6,BIRC2,BPI,BTG3,BZW1,C11orf24,C16orf72,C19orf38,C2CD4B,C6orf223,C8orf76,CA13,CAHM,CANT1,CAPNS1,CARD19,CASP10,CASP4,CASP5,CBWD2,CBX2,CCDC71L,CCL2,CCL23,CCL3,CCL3L1,CCL4,CCL4L2,CCL5,CCL7,CCNA1,CCRL2,CD274,CD44,CD58,CD82,CD9,CDC42EP4,CDKN1A,CDKN2A,CDKN2B,CFLAR,CHEK1,CHIC2,CHMP4B,CITED4,CLEC4E,CLIC4,CNTLN,CPAMD8,CPM,CREBBP,CRELD2,CSF2,CSF3,CSRNP1,CTDNEP1,CUL1,CUL4B,CWC25,CXCL10,CXCL11,CYB5R2,CYBRD1,CYFIP1,CYP51A1,CYP7B1,DCTN2,DDX60L,DENND2D,DGAT2,DIDO1,DLGAP1-AS1,DLGAP1-AS2,DNAJA1,DNER,DNPEP,DOCK4,DPYSL3,DTX3L,DUSP2,E2F7,ECE1,EDEM1,EDN1,EGLN1,EHD1,EIF1B,EIF4E,EML4,EMP3,ENDOV,ENPP4,EPAS1,ERBIN,ESPL1,ETS2,F5,FAM124A,FAM126B,FAM214B,FAM3C,FAM71F2,FAR2,FAS,FASTKD5,FBRS,FBXO48,FCAR,FEM1C,FES,FFAR2,FLCN,FLNA,FLOT2,FLT1,FLVCR2,FOLR3,FOXK1,FOXP4,FPGS,FTSJ1,FURIN,FXR2,FYN,G6PD,GABRG2,GALNT6,GAN,GAPVD1,GARS,GATAD2A,GBF1,GBP1,GBP2,GBP3,GBP7,GCH1,GDPD5,GFM1,GHRLOS,GIGYF1,GK,GK3P,GLT1D1,GMPPB,GNA11,GNA15,GNE,GNG4,GNL2,GPD2,GPR63,GPR68,GPR84,GPX4,GRAMD1A,GRB10,GRINA,GSAP,GSK3A,H1FX,HELZ2,HEXB,HEY1,HIF1A,HIST1H4I,HIVEP1,HK3,HLA-B,HLX,HPSE,HSPA5,IFFO2,IFITM1,IFNGR2,IGF2BP3,IKBKB,IL15RA,IL17RA,IL1B,IL1RN,IL27,IL31RA,IL32,IL36G,IL36RN,IL4I1,IL4R,INO80C,INPP5A,INPP5K,IRF1,IRF7,ITGA5,ITGAM,ITGAX,ITGB1,ITGB8,ITPKC,IVNS1ABP,JAK1,JAK3,JCHAIN,KCNJ2,KDELR1,KIAA1191,KIF21B,KIF24,KLHL2,KLHL6,KMT2E,KREMEN1,KRT23,KRT7,KTN1,KYNU,LAT2,LATS2,LDHA,LDLRAD3,LENG9,LHFPL2,LIF,LILRA6,LIMK2,LINC00309,LINC00482,LINC00656,LINC00852,LINC00909,LINC00937,LINC01001,LINC01093,LINC01176,LINC01268,LINC01270,LINC01588,LINS1,LIPN,LMNA,LPAR2,LPGAT1,LRRC59,LRRC8B,LRRC8C,LRRFIP2,LYN,MAFF,MAK,MALT1,MAML2,MANBAL,MANF,MAP1LC3B2,MAP2K1,MAP3K5,MAP3K7CL,MAP3K8,MAPK8,MARCKS,MBNL3,MCFD2,MCTP1,MCTP2,MED1,MED13L,MED27,MED30,MEI1,MET,METTL9,MFHAS1,MGAT2,MIR155HG,MIR3945HG,MLLT1,MLLT6,MMP19,MN1,MPP6,MRPS18A,MSL1,MTF1,MUCL1,MX1,MYLPF,MYO1G,N4BP1,NAB1,NATD1,NBN,NBPF9,NCK2,NCOA1,NDE1,NDRG1,NECTIN2,NEDD1,NEDD4,NEIL3,NEURL3,NEXN,NF1,NFE2L2,NFKB1,NFKBIA,NFKBIB,NINJ1,NIPAL4,NIT1,NLRP3,NOCT,NPTN-IT1,NR1H2,NR3C1,NRF1,NRP2,NSMAF,NUP62,NXT2,OAS2,OPTN,OR2B11,ORAI2,OTULIN,OXSR1,P2RX4,P2RX7,P4HB,PACERR,PALM2-AKAP2,PCGF1,PCID2,PDLIM7,PDPK1,PELI1,PELO,PFDN1,PHF11,PHF5A,PHLDA2,PIAS3,PIGA,PIK3AP1,PIK3R5,PIM3,PIWIL4,PKNOX1,PLAGL2,PLAUR,PLEC,PLEK,PLEKHA8,PLEKHF2,PLEKHO2,PLGRKT,PLK3,PLSCR1,PLXND1,PMPCB,PNP,PNPLA1,POGK,POR,PORCN,POU2F1,PPCDC,PPM1K,PPP1R15A,PPP1R17,PPTC7,PRDM8,PRELID3B,PRKACA,PRLR,PRND,PROK2,PRPF3,PRRC2A,PSD,PSEN1,PSMD11,PSMD6,PSRC1,PTDSS1,PTGER2,PTGIR,PTGS2,PTK2,PTPN14,PTPRA,PTX3,QPCT,R3HDM4,RAB12,RAB13,RAB22A,RAB28,RAB5A,RAB5B,RAB8B,RABGEF1,RALGDS,RAP2C,RASSF5,RBM15B,RBM17,RBM4,RBMS1,RBP7,RELA,RELB,RFFL,RHOH,RHOU,RIN2,RIT1,RMI2,RNF145,RNF185,RNF19B,RNF213,RP2,RSPO3,RTN4,RUNX1,S100A2,SBF1,SBNO2,SCAF1,SCARF1,SCLT1,SCN1B,SDC4,SDE2,SEC22B,SEC23B,SEMA3E,SEMA4B,SERINC2,SESTD1,SETD5,SH2B2,SH3GL1,SIPA1L2,SIRPA,SKA3,SKIL,SLAMF1,SLC11A2,SLC12A6,SLC16A7,SLC25A37,SLC25A44,SLC26A2,SLC26A6,SLC2A5,SLC2A6,SLC35B3,SLC39A14,SLC39A8,SLC43A3,SLC45A4,SLC9A8,SLCO4A1,SMAD3,SMAD4,SMG5,SMG7,SMG9,SMIM3,SMPD1,SMPDL3B,SMURF1,SNRK,SNX20,SOCS2,SOD2,SPACA6,SPATA24,SPG21,SPHK1,SPPL2A,SPSB1,SPTA1,SRC,SRD5A1,SREBF2,SRGN,SRI,SRSF5,SSH1,SSTR2,ST8SIA4,STAP1,STARD3,STARD8,STAT1,STK26,STK3,STOM,STOML2,STRIP2,STT3A,STX1A,SUZ12,SVIL,SYK,SYS1,TAP1,TAPBP,TBC1D22B,TBC1D30,TBK1,TBX21,TDRD6,TDRD7,TES,TEX10,TEX14,TEX30,TICAM2,TIFA,TIGD6,TIMM17B,TIMM23,TJAP1,TLCD2,TLE4,TM4SF19,TMEM106A,TMEM158,TMEM185B,TMEM217,TMEM52B,TMEM54,TMEM71,TMEM8A,TMOD3,TMTC2,TMUB2,TNF,TNFAIP2,TNFAIP6,TNFAIP8,TNFRSF12A,TNFRSF1A,TNFRSF1B,TNFRSF8,TNFSF15,TNFSF8,TNIP1,TNIP2,TOM1,TOR4A,TP53BP2,TP53I11,TPCN2,TRAF1,TRAF3,TRAF3IP2,TRIB1,TRIM21,TRIM25,TRIM69,TRPS1,TRPV2,TSPAN31,TUBB6,UBASH3B,UBE2W,UBE2Z,UBLCP1,UBR4,UBXN7,UEVLD,UNC13D,UPP1,URM1,UXS1,VAMP4,VIM,VNN3,WDFY1,WIPF1,WSB2,WTAP,XPR1,XRN1,XRN2,YRDC,ZBED6,ZBTB17,ZBTB21,ZC3H12A,ZC3HAV1,ZCCHC14,ZDHHC5,ZFYVE1,ZNF140,ZNF213,ZNF333,ZNF492,ZNF667-AS1,ZNF697,ZNF74,ZNF778,ZNFX1,ZNRF2,ZPR1
# SMDA3 : AMPH,ARHGAP24,ARHGAP29,BPI,CCDC71L,DLGAP1-AS1,DLX4,FXYD6,GLIS3,ITGA1,LAMB3,LIF,LIPN,LONRF3,MCEMP1,NOCT,NRGN,P2RX1,PLCB1,PNPLA1,PPARG,SELL,SERPINB10,SERPINB2,SGCA,SMAD3,STARD4,SYT1,TMTC1,TPBG,VNN2
# IRF1 : AAGAB,ABCA6,ABCB7,ABCC11,ABCC2,ABHD13,ABHD17C,ACBD3,ACOD1,ACSL4,ACTN1,ACTR6,ADAM10,ADAMTSL4,ADAMTSL4-AS1,ADAR,ADM,ADPRHL1,AFF1,AFF2,AFTPH,AGPAT3,AHCTF1,AHNAK,AKAP13,AKR1C3,ALDOA,AMIGO2,AMMECR1L,ANKIB1,ANKRD22,ANKRD28,ANKRD40,ANKS1A,ANTXR2,ANXA1,AP5B1,APOBEC3B,APOBR,APOL1,APOL2,APOL3,APOL4,APOL6,APTX,ARAP1,ARAP2,ARCN1,ARFGEF1,ARFIP1,ARHGAP25,ARHGAP30,ARHGAP31,ARHGAP5,ARHGAP9,ARHGEF10L,ARHGEF11,ARHGEF3,ARHGEF7,ARID4B,ARID5A,ARID5B,ARL17A,ARL17B,ARL5B,ARNT2,ARNTL,ARNTL2,ASH1L,ASH1L-AS1,ATF3,ATF6,ATP11B,ATP11C,ATP6V1B2,ATP6V1D,ATR,ATXN7L3,AZI2,B3GAT2,B3GNT5,B4GALT4,BAALC-AS1,BACH1,BANK1,BARD1,BATF2,BAZ1A,BAZ2A,BAZ2B,BBIP1,BBX,BCL10,BCL2L1,BCL2L13,BCL2L14,BCL3,BCL6,BCOR,BLOC1S6,BMF,BMPR2,BPI,BRMS1L,BTBD11,BTN3A1,C21orf91,C2orf27A,C3AR1,C5,C5orf15,C5orf51,C5orf56,C6orf47,CAB39L,CALM1,CALM2,CAMK1D,CAMK1G,CAMSAP2,CAPG,CAPZA1,CARD16,CARD17,CARM1,CASP1,CASP10,CASP3,CASP4,CASP5,CASP8,CBX4,CCDC117,CCDC149,CCDC30,CCDC69,CCL7,CCL8,CCM2,CCNG1,CCNL1,CCNY,CCR5,CCR6,CCRL2,CCSER2,CD2,CD22,CD226,CD247,CD274,CD300E,CD69,CD80,CDC14A,CDC27,CDC42EP4,CDC42SE2,CDCP1,CDK13,CDK2,CDKN1A,CDYL,CEBPB,CEBPZOS,CEMIP,CENPC,CENPJ,CENPU,CEP162,CFLAR,CHD2,CHD8,CHI3L2,CHN2,CHORDC1,CHRFAM7A,CIITA,CKAP2L,CLEC2B,CLEC6A,CLIC4,CLIP4,CLMN,CMIP,CMPK2,CMTR1,CNOT6,COA6,COBLL1,COLEC12,CPNE5,CPQ,CREB5,CREBRF,CREBZF,CRIM1,CRK,CSF2RB,CSNK1D,CSNK1G1,CSTA,CSTF2T,CUBN,CUL1,CUL2,CUL4B,CUX1,CXCL10,CXCL11,CXCL13,CXCL9,CXorf21,CXXC5,CYP7B1,CYTH4,DAPP1,DBR1,DCAF10,DCAF6,DCP2,DCUN1D3,DDX39B,DDX58,DDX60,DDX60L,DEDD,DENND1B,DENND3,DENND4A,DENND4C,DGKA,DGKG,DHX15,DHX8,DIAPH3,DICER1,DIS3,DKK2,DLGAP1-AS1,DNAJA1,DNAJB1,DNAJB4,DNAJC11,DNAJC13,DNAJC14,DNER,DNM1,DOCK4,DOCK5,DPM1,DSE,DTNBP1,DTX3L,DUSP10,DUSP6,DYNLT1,DYRK1A,EDEM1,EDN1,EFR3A,EGLN1,EGR4,EIF2AK2,EIF4G1,EIF5AL1,ELF4,ELK3,ELK4,ELMO2,ELMSAN1,ELOVL6,EMILIN2,EML4,EMP1,ENAH,ENC1,ENTPD4,ENY2,EPC2,EPG5,EPN2,ERAP1,ERC1,ERO1B,ETV7,EXT1,EXTL2,EYA3,EZH2,FAM102B,FAM111A,FAM114A1,FAM118A,FAM122B,FAM124A,FAM126A,FAM13B,FAM149B1,FAM174B,FAM8A1,FAR1,FAS,FBXL3,FBXO11,FBXO30,FBXO34,FBXO6,FBXW4,FBXW7,FCER2,FCGR1B,FCGR2B,FCHSD2,FCRLA,FEM1C,FER,FFAR2,FGD2,FGD4,FHAD1,FLVCR2,FMNL2,FNBP1,FNDC3B,FOXN2,FPR2,FRRS1,FRS2,FRY,G3BP2,G6PD,GAB2,GABBR1,GABRG2,GADD45B,GADD45G,GALNT4,GAN,GATAD2B,GBF1,GBP1,GBP2,GBP3,GBP4,GBP5,GBP7,GCH1,GFI1B,GGCX,GGT5,GK,GK3P,GLUL,GNAI3,GNAQ,GNB4,GNG5,GNPNAT1,GPANK1,GPHN,GPN2,GPR65,GPR84,GPX3,GREM1,GRIPAP1,GSS,GTDC1,GTF2B,GTF2E1,GTF2F2,GTPBP1,H2AFY,HAPLN3,HCCS,HCFC1,HCP5,HDX,HECA,HEG1,HELB,HELZ2,HERC5,HERC6,HES1,HESX1,HGF,HHAT,HIF1A,HIPK1,HIRA,HIST1H1C,HIST1H1E,HIST1H2BD,HIST1H2BE,HIST1H4E,HIVEP1,HIVEP3,HLA-E,HLX,HOTAIRM1,HPGD,HPSE,HRH2,HS3ST3B1,HS6ST1,HSPA1B,HSPH1,ID2,IDO1,IDO2,IER2,IER3,IFFO2,IFI35,IFI44,IFI44L,IFI6,IFIH1,IFIT2,IFIT3,IFIT5,IFITM1,IGF2BP2,IGSF6,IKBKG,IKZF4,IL10RA,IL15,IL15RA,IL1B,IL1R2,IL1RAP,IL1RL2,IL1RN,IL27,IL2RA,IL31RA,IMPAD1,INHBA,INO80D,INPP1,INTS3,INTS6,INTS6L,IQSEC1,IRF1,IRF2,IRF2BPL,IRF5,IRF7,IRF8,IRF9,IRS2,ISG15,ITPRIP,ITSN2,IVNS1ABP,JAG1,JAK2,JAZF1,JCHAIN,JMJD1C,KAT2B,KAT6A,KAT6B,KAT7,KBTBD2,KCNA3,KCNC4,KCNJ2,KCNMB1,KCTD20,KHDRBS1,KIAA0040,KIAA1109,KIAA1841,KIF13A,KIF13B,KIF2A,KLC1,KLF10,KLHL2,KLHL6,KMT2A,KSR1,KYAT3,KYNU,LACTB,LAMTOR3,LAP3,LATS1,LATS2,LCOR,LCORL,LDLRAD3,LDLRAD4,LEMD3,LILRA6,LILRB2,LILRB3,LIMK2,LIN54,LINC00189,LINC00211,LINC00265,LINC00309,LINC00528,LINC00607,LINC00683,LINC00884,LINC00910,LINC00937,LINC01089,LINC01093,LINC01355,LINC01366,LINS1,LMNA,LMNB1,LMTK2,LNPEP,LONRF1,LONRF3,LRCH1,LRP10,LRP8,LRRC39,LRRK2,LTBR,LYST,MAATS1,MAD2L2,MAFK,MAML3,MAMLD1,MAP3K11,MAP3K5,MAP3K7CL,MAPKAPK2,MAPKAPK3,MARCH1,MARCO,MAST4,MAX,MBD5,MBNL1,MBNL2,MBOAT1,MBP,MCF2L2,MCM6,MED16,MEF2A,MEFV,MFN2,MFSD6,MICALL1,MICU1,MIOS,MIR3945HG,MIR4435-2HG,MKLN1,MKRN1,MLKL,MLLT3,MLLT6,MMGT1,MMP25-AS1,MOB3C,MR1,MRPS6,MS4A7,MSL1,MSRB1,MSTO1,MTMR1,MTRF1L,MUC1,MX1,MX2,MXRA7,MYBL1,MYD88,MYH9,MYLIP,MYO1E,MYO1G,MYOF,N4BP1,N4BP2,NAA35,NAA60,NABP1,NAF1,NAGK,NAMPT,NBN,NBPF10,NBPF14,NBR1,NCALD,NCF2,NCOA6,NDST1,NEAT1,NEBL,NECTIN2,NEK10,NEK7,NF1,NFE2L2,NIPAL4,NKD1,NKX3-1,NLRC5,NLRP3,NMI,NNT,NOD1,NPTN,NR2C2,NR4A3,NRBP1,NREP,NRGN,NT5C2,NUB1,NUBP1,NUFIP2,NUGGC,NUP153,NUP62,OAS2,OASL,ODF2L,OGFR,OR2B11,ORC5,OSBPL3,OSBPL9,OSER1,OTUD6B,P2RX7,PACSIN2,PALM2-AKAP2,PAN3,PANX1,PARP10,PARP12,PARP14,PARP8,PARP9,PATL1,PCED1B-AS1,PCNX1,PCNX4,PDE4B,PDE6H,PDE8B,PDHX,PDP1,PDXK,PEF1,PELI1,PFKP,PGK1,PGLYRP2,PHC2,PHTF1,PIAS1,PICALM,PIK3C3,PIK3CA,PIM1,PIP4K2B,PITPNM2,PLA1A,PLCG2,PLD2,PLEK,PLEKHA7,PLEKHN1,PLK2,PLSCR1,PLSCR4,PML,PMP22,PNPT1,POLD3,POLR3C,POMK,POU2F1,PPARD,PPARG,PPP1R10,PPP1R12A,PPP1R12C,PPP1R15A,PPP1R15B,PPP1R17,PPP1R18,PPP1R8,PPP4R2,PRDM1,PRDM2,PRDM4,PRKAG2,PRKD2,PRLR,PRNP,PRR16,PRR5L,PRRC2B,PSD3,PSD4,PSMA4,PSMB2,PSMB8,PSMB8-AS1,PSMB9,PSTPIP2,PTAFR,PTBP2,PTK2,PTPN6,PTPRA,PTPRE,PURA,PUS3,PXN,QSOX2,R3HCC1L,RAB13,RAB31,RAB39A,RAB5B,RAF1,RALB,RALGPS2,RANBP10,RAP2C,RAP2C-AS1,RARA,RASA2,RASA3,RASAL2,RASGEF1B,RASSF5,RBM33,RBM4,RBM43,RBMS1,RBMXL1,RDX,RELA,RELB,RERE,REV1,RFWD3,RFX7,RGMA,RGS17,RHBDF2,RHOH,RHOU,RIPK1,RIPK2,RLIM,RMRP,RNASE2,RNASEL,RNF111,RNF144B,RNF169,RNF19B,RNF213,RNF24,RNH1,ROCK1,ROCK2,RPH3A,RRAGC,RREB1,RSAD2,RSPO3,RSPRY1,RUNX1,S100A10,S1PR4,SAMD4A,SAMD9,SAMD9L,SAMHD1,SARS,SAT1,SBF2,SBNO2,SCAPER,SCARB2,SCARF1,SCIN,SCLT1,SCYL3,SDC2,SECTM1,SEMA4A,SEMA4D,SERPINB1,SERPINB9P1,SETD3,SETX,SH3BP2,SH3GLB1,SIK2,SIPA1L1,SIRPB2,SKIL,SLA,SLC11A2,SLC16A10,SLC16A3,SLC1A2,SLC1A3,SLC22A15,SLC25A28,SLC25A37,SLC2A9,SLC30A4,SLC38A2,SLC41A1,SLC43A2,SLC44A1,SLC5A3,SLC6A12,SLCO3A1,SLCO4A1,SMAD7,SMC5,SMG1,SMG7,SMIM11A,SMIM12,SMIM15,SNRNP27,SNTA1,SNTB1,SNTB2,SNX10,SNX12,SNX19,SNX20,SOCS1,SOS1,SP100,SP110,SP140,SPPL2A,SPTA1,SPTSSA,SRSF1,SSBP1,SSBP2,ST3GAL2,ST3GAL6,STAG1,STAMBPL1,STARD13,STAT1,STAT2,STAT3,STEAP3,STEAP4,STK17B,STK3,STK35,STOML1,STX12,STX17,STX3,SULT1B1,SVIL,SYNJ2,TAB2,TAB3,TAF1,TAF10,TAF1B,TAGAP,TANK,TAP1,TAP2,TARDBP,TBPL1,TBX21,TCF20,TCF4,TCF7L2,TDRD7,TERF2IP,TESK2,TET2,TFG,THAP5,THOC2,TICAM2,TIFA,TJAP1,TLCD2,TLK1,TLR1,TLR2,TLR4,TLR8,TMEM106A,TMEM131,TMEM173,TMEM43,TMEM54,TMEM64,TMOD3,TMPO,TMSB10,TMTC4,TNFAIP1,TNFAIP2,TNFAIP3,TNFAIP8,TNFRSF14,TNFSF10,TNFSF13,TNFSF13B,TNFSF8,TNKS,TNKS2,TNS3,TOB1,TOR1B,TP53BP2,TP53I11,TRAF3IP2-AS1,TRAF6,TRAFD1,TRAK1,TREM1,TRIM21,TRIM22,TRIM25,TRIM5,TRIM56,TRIM69,TRMT1,TRMT10B,TRMT6,TRPM2,TSC22D2,TSHZ1,TSPYL4,TTF2,TUBB2A,TUG1,TVP23B,TVP23C,TWF1,TXNDC15,TYMSOS,UBALD2,UBE2E1,UBE2H,UBR1,UBR4,UBR5,UBXN2A,UBXN7,UCHL5,UEVLD,UGGT2,UGP2,UHRF1BP1L,UIMC1,UNC45A,UROD,USP15,USP2,USP25,USP3,USP32,USP42,USP6NL,VAMP4,VAMP5,VAV1,VAV3,VCL,VCPIP1,VDR,VGLL4,VIM,VMP1,VNN3,VPS13B,VPS16,VPS18,VPS54,VPS8,WARS,WBP1L,WDFY1,WDFY3,WDR41,WDR43,WDR7,WIPF1,XAF1,XRN1,YME1L1,YTHDC2,ZADH2,ZBED5,ZBP1,ZBTB14,ZBTB18,ZBTB34,ZC3H12A,ZC3HAV1,ZCCHC14,ZDHHC5,ZER1,ZFC3H1,ZFP36,ZFP69B,ZFP90,ZMYM6,ZNF143,ZNF200,ZNF266,ZNF276,ZNF322,ZNF333,ZNF407,ZNF646,ZNF699,ZNF710,ZNF800,ZNF839,ZNFX1,ZNRF2,ZSWIM1
# STAT1 : ABCD1,ABR,ACAA2,ACBD5,ACSF3,ACSL1,ACSL4,ACTB,ADAM9,ADAR,ADGRE2,ADGRE5,ADPRH,AFF1,AFTPH,AGAP3,AGO2,AGO4,AGPAT2,AGPAT3,AIFM1,AIFM2,AIM2,AKT1S1,ALPK1,AMMECR1,AMZ2,ANKFY1,ANKRD17,ANKRD22,ANO6,ANXA1,ANXA4,ANXA6,AP5B1,APOBEC3C,APOBEC3D,APOL1,APOL2,APOL3,APOL4,APOL6,AQR,ARAP2,ARFGEF1,ARHGAP25,ARHGAP30,ARHGAP31,ARHGEF11,ARHGEF9,ARID4B,ARID5A,ARID5B,ARMC2,ARMCX3,ARNTL,ARPC1B,ASCL2,ASXL1,ATAD1,ATF5,ATG7,ATP11C,ATP13A1,ATP13A2,ATP2A2,ATP6V0E1,ATP6V1B2,ATXN1,ATXN7,ATXN7L3,B2M,B3GAT2,B4GALT4,BANK1,BATF2,BAZ1A,BAZ2A,BAZ2B,BBX,BCL2L1,BCL2L11,BCL2L14,BCOR,BLNK,BMF,BMP1,BMP4,BMPR2,BRWD3,BTBD10,BTN2A2,BTN3A3,C11orf58,C17orf53,C18orf25,C2,C2orf27A,C3orf38,C5orf56,C5orf66,C9orf72,CABP4,CADM1,CALM1,CALM2,CANX,CAPG,CARD6,CASK,CASP1,CASP3,CASP4,CASP5,CASP7,CBL,CCDC167,CCDC6,CCDC88B,CCL13,CCND3,CCNYL1,CCR1,CD207,CD276,CD38,CD47,CD59,CD80,CD86,CD93,CDC25B,CDC42SE1,CDC42SE2,CDCP1,CDK12,CDK17,CDK6,CDKL5,CEP162,CEP350,CEP57L1,CERS6,CFL1,CHMP5,CIITA,CLASP1,CLCN5,CLDN15,CLEC16A,CLIP4,CLK4,CLN5,CLTC,CMIP,CMPK2,CMTR1,CNDP2,CNP,CNPPD1,COLEC12,COLGALT2,COPB1,COPE,CORO1C,CPEB3,CPSF2,CREB1,CREBBP,CREBRF,CRISPLD2,CRLF3,CSF2RB,CSK,CSTF2T,CTNNA1,CTNND1,CTNS,CTSB,CTSL,CTSS,CTSZ,CUL1,CUL4B,CWC22,CXCL10,CXCL11,CXCL16,CXCL9,CXorf21,CXorf38,CYBB,CYTH1,CYTH4,DAPP1,DBR1,DCAF10,DCAF5,DCAF6,DCTN1,DDR1,DDX58,DDX60,DDX60L,DEDD,DENND1B,DENND2D,DENND3,DGKA,DGKH,DHX58,DHX8,DICER1,DIDO1,DLGAP4,DMTF1,DMXL2,DNAJC14,DNAJC5,DOCK10,DOCK4,DOCK7,DOK1,DPP9,DPYD,DPYSL2,DTNBP1,DTX2,DTX3L,DUSP10,DYNLT1,DYRK1A,DYSF,E2F3,EDIL3,EEPD1,EFR3A,EGLN2,EIF2AK2,EIF4E3,EIF4G1,EIF4G3,ELAVL4,ELF1,ELK4,ELMO2,ELMSAN1,EMILIN1,EML4,ENAH,ENPP4,ENTPD4,EP300,EPB41L3,EPDR1,EPS15,EPSTI1,ERAP1,ERF,ERMN,ERRFI1,ESR1,ETV6,ETV7,EVI2A,EXOC1,EYA3,EYS,FAM120B,FAM135A,FAM160B1,FAM193B,FAM204A,FAM225A,FAM3C,FAM49A,FAM49B,FAM53B,FAM91A1,FAR1,FAR2,FAS,FASTKD1,FBXO34,FBXO6,FCER2,FCGR1A,FCGR1B,FCGR2B,FCGR3A,FCHO2,FGR,FH,FHAD1,FKBP5,FLVCR2,FMN1,FMNL2,FMR1,FNDC3A,FRMD3,FRS2,FUBP3,FUT4,FYN,GAB2,GABBR1,GABPB1-AS1,GAS6,GBP1,GBP4,GBP5,GBP7,GCLM,GDI1,GFI1B,GGA1,GGT5,GIMAP2,GIMAP4,GLIPR2,GLS,GLUL,GNA13,GNAI3,GNB2,GNB4,GNG11,GNG2,GNS,GOLGA3,GOSR2,GPD2,GPM6A,GPR18,GPR68,GPX3,GRAMD1B,GRIPAP1,GRK5,GSAP,GSK3B,GTF2F1,GTF2I,GTPBP1,GTSF1,H6PD,HACE1,HAPLN3,HAT1,HAX1,HCP5,HELZ2,HERC5,HERC6,HIVEP3,HLA-A,HLA-B,HLA-E,HLX,HM13,HNRNPH2,HSD3B7,HSP90B1,HTATIP2,IDH1,IDO1,IDO2,IFFO2,IFI16,IFI35,IFI44,IFI44L,IFI6,IFIH1,IFIT2,IFIT3,IFIT5,IFITM1,IFITM3,IGF2BP2,IGF2BP3,IGF2R,IGSF8,IL10,IL12RB1,IL15,IL15RA,IL18BP,IL21R,IL27,IL2RA,IL4I1,IL6ST,INO80,INPPL1,INSM1,IPO9,IQSEC1,IRF1,IRF2,IRF2BPL,IRF5,IRF7,IRF8,IRF9,ISG15,ISG20,ITGAV,ITK,ITPKC,ITPR1,ITSN2,JAK2,JAK3,JCHAIN,JMJD1C,JUND,KAT2B,KBTBD2,KCNJ2,KCTD10,KCTD12,KDM6A,KHNYN,KIAA0232,KIAA1109,KIDINS220,KIF13B,KIF1B,KIF1C,KIFC3,KIZ,KLF10,KLF9,KLHDC10,KLHL6,KMT2C,KPNB1,KSR1,KTI12,LACTB,LAG3,LAMTOR5,LAP3,LARP4B,LCP1,LDLRAD4,LEMD3,LGALS3BP,LGALS8,LGALS9,LHFPL2,LILRA4,LILRB1,LILRB2,LILRB4,LINC00886,LINC00910,LINC00943,LINC01140,LINC01422,LINC01550,LMNA,LMNB1,LMO2,LNPEP,LPAR1,LPAR4,LPGAT1,LPP,LPXN,LRCH1,LRRC4,LRRC75A,LRRC8B,LRRK2,LSM14B,LY75,LYRM9,LZTR1,MAML3,MAN2A1,MAP3K11,MAP3K5,MAP4,MAPRE3,MARCH1,MAST4,MAT2B,MAX,MBD5,MCL1,MEF2A,MEF2D,MEFV,METRN,MFN2,MGAT1,MGAT5,MIA3,MICB,MICU1,MIEF1,MIER1,MIRLET7BHG,MITF,MLKL,MLLT6,MMAA,MMP25,MMP25-AS1,MOB3B,MOB3C,MOCS3,MOV10,MPEG1,MPP1,MR1,MRC1,MRFAP1L1,MSL1,MSN,MSR1,MSRA,MTMR1,MTMR10,MTMR2,MUC1,MX1,MX2,MXRA7,MYO1E,MYO6,MYOF,N4BP1,N4BP2,NADK,NAGK,NAPA,NBR1,NCF2,NCKAP1L,NCOA3,NCOA7,NCOR2,NDST1,NECTIN2,NEDD1,NEDD4,NEK7,NETO2,NF1,NFATC2IP,NFE2L1,NFE2L3,NIN,NINJ2,NIPBL,NKD1,NLK,NLRC5,NMI,NPAT,NPHP1,NR2C2,NR2F1,NRAS,NRN1,NRP1,NRP2,NT5C2,NT5C3A,NUB1,NUMB,NUP107,NUP153,NUP160,NUP205,NUP62,OAS2,OAS3,OGFR,OLFML2B,OSBPL3,OSGIN2,OTUD5,P2RY1,P2RY14,PABPC1,PAK2,PALLD,PAM,PAN2,PARP10,PARP11,PARP12,PARP14,PARP15,PARP8,PARP9,PBX3,PCGF5,PCNX3,PDE11A,PDIA3,PDIA4,PDP2,PDS5B,PEA15,PGRMC1,PHACTR1,PHACTR4,PHC2,PHC3,PHIP,PHKG1,PIAS1,PICALM,PIEZO1,PIK3C3,PIK3CA,PIKFYVE,PIM1,PINK1,PIP4K2A,PITRM1,PKN2,PLA1A,PLCL2,PLD1,PLEK,PLEKHA2,PLEKHA7,PLEKHA8,PLEKHG3,PLEKHM3,PLEKHO2,PLXDC2,PLXNC1,PLXND1,PML,PNPLA6,POC5,POGZ,POU5F2,PPCS,PPM1B,PPP1R12C,PPP1R18,PPP2R2A,PPP4R3A,PRDM1,PRDM2,PRDM4,PRICKLE4,PRKAG2,PRKCE,PRKCSH,PRKRIP1,PRLR,PRR14L,PRR5L,PSMA3,PSMA4,PSMA5,PSMB8,PSMB8-AS1,PSMB9,PSME2,PTAFR,PTMS,PTPN6,PTPRC,PTPRE,PTPRJ,PWP2,PXK,PYCARD,PYM1,QSOX1,RAB14,RAB20,RAB27A,RAB29,RAB35,RAB39A,RAB3IL1,RAB5A,RAB5C,RAD21,RALY,RASAL2,RASGEF1B,RASGRP4,RASSF2,RASSF4,RASSF5,RBM23,RBM47,RBMS1,RBMXL1,RC3H2,RCAN2,RCSD1,RFC1,RFFL,RFTN1,RFX5,RGL1,RGS5,RHEBL1,RHOH,RIPK1,RIPK2,RLIM,RNASEH2B,RNASEL,RNF111,RNF115,RNF13,RNF144B,RNF145,RNF19A,RNF2,RNF213,RNF217,RNF31,RNF4,ROCK1,ROCK2,RPTOR,RRBP1,RREB1,RSAD2,RSF1,RTCB,RTN4,RUFY3,RUSC1-AS1,S1PR2,SAMD4A,SAMD8,SAMD9,SAMD9L,SAMHD1,SAMSN1,SAP30BP,SAR1B,SARS,SASH1,SCAMP3,SCAPER,SCARB2,SCIMP,SCIN,SCLT1,SCN9A,SCYL3,SDC3,SDCCAG8,SEC23A,SEC24B,SEL1L,SELPLG,SEMA4A,SEMA4D,SERPING1,SERTAD1,SETX,SFT2D1,SH3BP1,SH3GLB1,SHISA5,SHOC2,SIGLEC1,SKAP2,SKIL,SLA,SLC12A5,SLC12A7,SLC15A3,SLC16A3,SLC1A3,SLC24A1,SLC25A25,SLC25A28,SLC39A11,SLC43A2,SLC4A7,SLC8A1,SLC9A3R1,SLF2,SLK,SLMAP,SMAD2,SMAP2,SMARCAL1,SMARCD3,SMCO4,SMG1,SMG7,SND1-IT1,SNX11,SNX12,SNX16,SNX20,SORT1,SOS1,SOX13,SP100,SP110,SP140L,SPATS2L,SPCS3,SPECC1L,SPHK1,SPIN1,SPRY2,SPTAN1,SPTLC1,SPTSSA,SRP54,SSBP3,ST3GAL2,ST3GAL5,STAC,STAG1,STAG2,STAM2,STARD13,STAT1,STAT2,STAT3,STAT6,STK3,STOM,STRN,STT3A,STX11,STX12,STX16,STX17,STX3,SVIL,SYPL1,SYTL1,TACSTD2,TADA3,TAF15,TAGAP,TAP1,TAP2,TAPBP,TAS2R14,TBC1D1,TBC1D14,TBC1D32,TBKBP1,TBL1XR1,TCF12,TCF20,TCF4,TDRD7,TEP1,TESK2,TFEC,TGFBI,TGM2,TGOLN2,THEMIS2,TIAM2,TICAM2,TIPARP,TIRAP,TJP2,TLK2,TLR2,TLR8,TMEM102,TMEM106A,TMEM109,TMEM123,TMEM140,TMEM236,TMEM30A,TMEM51,TMEM52B,TMEM62,TMEM63A,TMEM64,TMEM86A,TMOD3,TNFSF10,TNFSF13B,TNFSF15,TNKS,TNKS2,TNRC6C,TNS3,TOMM40L,TOP1,TOP2B,TOR1AIP1,TOR1B,TOR4A,TP53I11,TPM4,TRAF6,TRIM14,TRIM21,TRIM25,TRIM26,TRIM34,TRIM5,TRIM56,TRIM69,TRIOBP,TRIP12,TRPC4AP,TRPM2,TSC22D1,TSC22D4,TSPAN2,TTC7B,TTLL7,TWF1,TXNIP,UBA1,UBA6,UBALD2,UBD,UBE2L6,UBE3A,UBR1,UBR5,UBTF,UBXN4,UCP2,UEVLD,UNK,UQCRC2,URM1,USP15,USP3,USP32,USP38,USP42,USP6NL,USP9X,UTRN,UVRAG,VAMP5,VANGL1,VAV3,VCAM1,VCPIP1,VGLL4,VPS13C,VPS54,VPS8,WARS,WDFY1,WDR44,WDR47,WIPF1,WSB2,WWC2,WWP2,XAF1,XIAP,XPO6,XRN1,ZBP1,ZBTB26,ZBTB7B,ZC3H12D,ZEB2,ZFAT,ZFX,ZFYVE1,ZHX2,ZNF140,ZNF217,ZNF223,ZNF24,ZNF322,ZNF397,ZNF438,ZNF496,ZNF641,ZNF646,ZNF664,ZNF710,ZNF850,ZNFX1,ZNHIT1,ZNRF2,ZSCAN2
# HIF1A : ACVR1B,ARF6,ASPH,ATP11B,BCAT1,CPEB2,CRADD,DMXL2,ELF4,ERGIC1,ERP44,FZD4,GPAT4,HPSE,IL1R1,IQSEC1,JAK1,MORC3,NAMPT,NBN,PAG1,PATL1,PIK3AP1,PLXNC1,PNPLA8,PPP1R3B,PPP4R2,PTPN1,RAB5A,RABGEF1,RALA,RUNX1,SBNO2,SCYL2,SLC16A3,SMCHD1,SRF,SSH1,STK40,STRN3,TET2,TFE3,TICAM1,TPM4,TRIM21,USP32,VPS9D1

# RELB: AACS,ABTB2,ACHE,ACSL1,ACSL4,ACTN1,ACVR2A,ADA,ADGRE1,ADGRG1,ADNP2,ADORA2A,AGPAT4,AK4,ALAS1,AMIGO2,AMPD3,ANPEP,ANXA11,ANXA2,AP2B1,APBA3,APP,ARAP3,ARFGEF1,ARID3A,ARL13B,ARL8B,ATF2,ATP13A3,ATP1A1,ATP2C1,ATP6V1C1,ATXN1L,ATXN7L1,AZIN1,B3GNT2,BAALC,BASP1,BATF2,BAZ2A,BBIP1,BCL10,BCL2A1,BCL2L14,BCL3,BECN1,BEX1,BID,BIRC2,BPI,BTG1,BTG3,BZW1,C11orf80,C11orf96,C15orf39,C1orf198,C1orf21,C2orf49,CABLES2,CAMKK2,CAMTA1,CANT1,CARD16,CASP1,CASP5,CAST,CCDC107,CCDC71L,CCDC82,CCDC93,CCL2,CCL23,CCL4,CCL4L2,CCL5,CCL7,CCND2,CCNYL1,CCRL2,CCSER2,CD200,CD36,CD44,CD82,CDA,CDC42EP3,CDKN1A,CDKN2B,CELF1,CEP170B,CERK,CERS2,CFLAR,CHD4,CHERP,CHIC2,CHMP2B,CHST11,CHST2,CLCF1,CLEC4D,CLEC4E,CLGN,CLIC4,CLIP1,CLMN,CLSTN3,CMAS,CNKSR3,COCH,COG3,COX7A2L,CREB3,CREBRF,CSF2RB,CSF3,CSNK1G3,CTDSP1,CTNNB1,CUL4B,CXCL1,CXCL11,CXCL5,CYB5R2,CYFIP2,CYP7B1,DCUN1D3,DDX60L,DENND2D,DERL1,DGAT2,DIS3,DLEU7,DLG5,DLGAP1-AS1,DLL1,DMTF1,DNAAF1,DNAJA1,DNAJC1,DR1,DRAM1,DSE,DUSP18,DUSP2,DUSP3,DUSP8,EBF1,ECE1,EDNRB,EHD1,EHD4,EIF1,EIF1B,ELK3,EMP3,ENPP4,EPM2AIP1,ERF,ERRFI1,EYA3,F3,F5,FAM107B,FAM124A,FAM3C,FAM49A,FASTKD5,FBXL3,FBXO34,FKBP5,FLOT1,FMNL3,FOSL1,FOXC1,FTH1,FXYD6,FZD4,G0S2,GABRG2,GALK1,GATAD2A,GCLM,GGH,GK,GLRX3,GNA15,GNG5,GPBAR1,GPR35,GPR63,GPR84,GPSM1,GRINA,GUK1,GYPC,HECTD1,HEG1,HES4,HEY1,HIF1A,HIVEP1,HLA-F,HMGA1,HNRNPC,HNRNPH3,HOMER3,HOPX,HS3ST3B1,IER3,IFFO2,IFIH1,IFT43,IKBKB,IL15,IL1B,IL1R1,IL1RN,IL27,IL36G,IL36RN,IL4I1,IL4R,IL6,INHBA,INO80C,IRAK2,IRF2BPL,IRGQ,IST1,ITGA1,ITGA3,ITGB8,ITPRIP,JAK1,KAT6B,KCNJ2,KCTD21,KIAA0040,KIF13A,KLHL25,KRT17,KRT7,KYNU,LAMA2,LAMB3,LAMTOR3,LASP1,LCP2,LEP,LGALS3,LHFPL2,LIF,LIMS1,LINC00877,LINC01215,LINC01262,LINC01270,LINC01358,LITAF,LMNB2,LRP12,LRRFIP2,LYSMD2,MAK,MAML1,MAN1A2,MAP1LC3A,MAP2K3,MAP4K4,MARCH3,MARCKS,MARCKSL1,MCFD2,MCMBP,MCOLN3,MED1,MED10,MEG8,MET,MFHAS1,MICAL1,MICALL2,MIER3,MLLT1,MN1,MORF4L2,MPP3,MRPL52,MSANTD3,MSC,MSN,MTF1,MTMR14,MTMR3,MTSS1,MYH9,MYO1G,MYOM2,N4BP1,NAA50,NAMPT,NBEAL1,NCF1,NCR3LG1,NDUFV1,NEDD4,NEURL3,NFKB1,NFKB2,NFKBIA,NID1,NINJ1,NIPAL4,NRIP3,NRROS,NSMAF,NT5E,NUS1,OPTN,OR10G3,ORAI1,OSBP,OSGIN1,OSM,OTUD1,P2RX4,P4HA1,PAG1,PANX1,PCDH11X,PCGF1,PDE8A,PDLIM7,PEA15,PELI3,PELO,PEMT,PERP,PFKFB3,PGM3,PHF8,PHLDA1,PHLDA2,PIGA,PIK3CD,PIM3,PKNOX1,PLAC8,PLAGL2,PLD1,PLEC,PLEKHB2,PLEKHO2,PLPP1,PLPP3,PNPLA1,POU2F2,PPARD,PPP1CB,PPP1R18,PPP2CB,PPP2R5E,PRELID3B,PRPF3,PSD,PSMB1,PSMB2,PTGIR,PTX3,RAB11B,RAB12,RAB21,RAB22A,RAB5A,RAB8B,RAB9A,RALA,RALGDS,RALGPS2,RAP1B,RAP2C,RAP2C-AS1,RAPGEF2,RAPGEFL1,RASSF5,RDX,RELA,RELB,RFX2,RHBDD1,RHOC,RHOH,RHOU,RIN2,RIPK2,RMRP,RNF145,RNF152,RNF185,RNF19B,RNF24,RRAD,RRAGD,RUNX3,RYBP,SAMD4A,SAR1A,SAR1B,SAV1,SCN1B,SDC2,SDC4,SDCBP,SEC22B,SEC24A,SEMA3C,SEMA4F,SERPINA1,SERPINB8,SERPINB9,SESN2,SFR1,SGMS1,SGMS2,SGPP2,SH2D3A,SH3GLB1,SHB,SHC1,SIN3A,SIPA1L2,SLAIN2,SLAMF1,SLC11A2,SLC16A4,SLC26A2,SLC28A3,SLC2A14,SLC2A6,SLC39A8,SLC39A9,SLC41A2,SLC43A2,SLC43A3,SLC44A1,SLC4A2,SLCO3A1,SLCO4A1,SLPI,SMAD7,SMARCD2,SMURF1,SNN,SNX10,SNX20,SOCS2,SOCS3,SOD2,SPHK1,SPSB1,SRF,SRGN,SSPO,SSSCA1-AS1,SSTR2,ST3GAL1,STARD3NL,STARD4,STARD8,STEAP3,STK24,STK26,STK3,STK35,TAB2,TAGLN2,TANC2,TANK,TBC1D20,TBC1D30,TBX21,TET3,TEX10,TFPI,TGIF2,THAP12,THBS1,TICAM1,TIFA,TIMP3,TIPARP-AS1,TJP2,TLN1,TLR2,TMEM120A,TMEM159,TMEM185B,TMEM38B,TMEM63B,TNFAIP6,TNFAIP8,TNFRSF1B,TNFRSF9,TNFSF14,TNIP1,TNIP2,TNIP3,TOR4A,TPM4,TPST1,TRAF3IP2,TRAF6,TRAPPC3,TRIM36,TRIM72,TUBA1C,TWISTNB,TXN,TXNDC16,UAP1,UBA1,UBALD1,UBE2A,UBE2H,UBE2J1,UBE2Z,UBE3A,UPP1,USE1,UXS1,VCP,VDR,VIM,VPS37C,VWF,WDFY1,WDTC1,WIPF2,WNT5A,WTAP,XBP1,XRCC1,YBX1,YKT6,ZADH2,ZBTB17,ZBTB33,ZC3H12A,ZC3H12C,ZMIZ1,ZNF267,ZNF317,ZNF496,ZNF672,ZNF697,ZNF710,ZNF720,ZSWIM4
# HIVEP1_ext: AACS,AATK,ABCA1,ABCA5,ABHD17C,ABL2,ABLIM1,ABR,ACAT2,ACOT9,ACSL1,ACTN1,ACTR3,ACVR2A,ADA,ADAM17,ADAM19,ADM,ADORA2A,ADRB2,ADTRP,AGAP3,AGPAT4,AHCTF1,AK4,AKT3,ALCAM,AMPD3,ANKRD12,ANKRD28,ANKRD33B,ANPEP,ANTXR2,AQP9,ARAP1,ARAP2,ARAP3,ARFGAP3,ARHGAP26-IT1,ARHGAP31,ARHGAP31-AS1,ARHGEF10L,ARHGEF2,ARID3A,ARID5A,ARID5B,ARL5B,ARL8B,ARMC5,ARNTL2,ASAP1,ASH1L,ATF2,ATF5,ATF6,ATP13A3,ATP2A2,ATP6V1C1,AZIN1,BAALC-AS1,BACH1,BANP,BATF2,BAZ1A,BCAR3,BCAT1,BCL2,BCL2A1,BCL2L1,BCL3,BCL6,BCOR,BEND3,BID,BIRC2,BIRC3,BMP6,BPI,BRPF3,BTBD19,BTG3,BZW1,C1orf61,C1QTNF1,C21orf62,C6orf223,C6orf62,CA13,CAMK1G,CAMKK2,CANT1,CASZ1,CCDC154,CCDC71L,CCDC93,CCL20,CCL3,CCL3L1,CCL4,CCL4L2,CCL5,CCM2L,CCNL1,CCR7,CD109,CD274,CD300E,CD44,CD58,CD82,CD93,CDC42EP3,CDK1,CDK12,CDK14,CDKL5,CDKN1A,CDKN2B,CELF1,CFLAR,CFLAR-AS1,CHD4,CHEK1,CHMP4B,CHST15,CLEC2D,CLIC4,CLIP1,CNOT4,CPAMD8,CPEB4,CRADD,CRISPLD2,CRY1,CSF1,CSF2,CSF3,CSGALNACT2,CSRNP1,CTNNAL1,CUX1,CXCL1,CXCL2,CXCL3,CXCL5,CXCL8,CXXC5,CYB561A3,CYB5R2,CYFIP1,CYTH1,DDX3X,DDX58,DDX60L,DENND2D,DENND4A,DENND5A,DGKH,DIAPH1,DICER1,DLEU2,DLGAP1-AS1,DLGAP4,DMWD,DMXL2,DNAJB5,DRAM1,DUSP10,DUSP3,DUSP5,DUSP6,DYRK1A,DYRK3,E2F3,E2F6,E2F7,EBF1,EBI3,ECE1,EGOT,EGR3,EHD1,EID3,EIF1B,EIF1B-AS1,EIF4G2,ELL2,ELOVL7,EMILIN2,EML3,EMP1,EPC2,EPG5,ERBIN,ERGIC1,ETS1,ETV3,ETV6,EYA3,EZH2,F3,FAM107B,FAM124A,FAM83G,FBRS,FERMT2,FFAR2,FGF13,FGR,FLCN,FLNA,FLT1,FMN1,FMNL3,FNDC3A,FNDC3B,FNIP2,FOSL1,FOXP4,FRMD6,FRY-AS1,FSCN1,FSTL3,FURIN,FUT6,FYN,GABPB1,GALNT6,GBE1,GBP2,GCH1,GHRLOS,GK5,GLIS3,GLS,GNG2,GOLGA8A,GPATCH2L,GPR132,GPR137B,GPR157,GPR68,GPR84,GPRC5A,GRAMD1A,GSAP,HDGF,HIP1R,HIPK2,HIPK3,HIVEP1,HIVEP2,HMGCS1,HMGXB4,HNRNPC,HOPX,HRH1,HS3ST1,HS3ST3B1,HUWE1,HYMAI,HYOU1,ICAM1,ICAM4,IER3,IFFO2,IFIH1,IFNLR1,IGF2R,IL1A,IL1B,IL1R1,IL1RN,IL20,IL23A,IL2RA,IL32,IL36G,IL36RN,IL4R,IL6,IL7,IL7R,INF2,INHBA,INPP5A,INTS6,IRAK3,IRF1,ITGA5,ITGB8,ITPRIP,IVNS1ABP,JAG1,JAK1,JRKL,KANK1,KBTBD2,KCNA3,KCNJ2,KCNN4,KDM4B,KDM5C,KDM6A,KDM6B,KDM7A,KIDINS220,KIF13A,KIF1B,KIF21B,KIFC3,KLF6,KLF7,KLHL21,KLHL6,KMT2E,KREMEN1,KYNU,LACC1,LACTB2-AS1,LAMB3,LAMP3,LAT,LCOR,LCP2,LFNG,LGALS3,LHFPL2,LIMK2,LINC-PINT,LINC00158,LINC00189,LINC00309,LINC00346,LINC00622,LINC00909,LINC00910,LINC00937,LINC01215,LINC01268,LINC01465,LINC01476,LINC01588,LITAF,LMNB2,LONRF3,LRCH3,LRP12,LUCAT1,LYN,MAFF,MAFG,MALT1,MAML2,MAMLD1,MAP2K3,MAP3K5,MAP4K3,MAP4K4,MAPK13,MAPK6,MARCH3,MARCKS,MB21D2,MCOLN2,MCTP1,MDGA1,MECP2,MED13,MEFV,MELTF,MFHAS1,MFSD2A,MGAM,MGLL,MGRN1,MICALL1,MICALL2,MIR155HG,MIR3142HG,MIR3945HG,MLLT6,MMP14,MOB3B,MOB3C,MORF4L2,MREG,MSANTD3,MSC,MTMR3,MYBPC3,MYH9,MYO1E,MYO1G,MYO9B,N4BP1,NAB1,NAB2,NABP1,NBN,NBPF10,NBPF14,NBPF20,NBPF9,NCALD,NCOA5,NCOR2,NCR3LG1,NCS1,NDRG1,NECTIN2,NEDD4L,NEK7,NEMP1,NEURL3,NEUROD1,NFAT5,NFATC1,NFE2L2,NFKB1,NFKB2,NFKBIA,NFKBIB,NFKBIE,NFKBIZ,NFYA,NINJ1,NIPAL4,NKD1,NLRP3,NOTCH1,NOTCH2,NPTN-IT1,NR3C1,NR4A3,NRIP1,NRIP3,NRP2,NSMAF,NSUN6,NT5E,NUMB,NUP188,OAZ2,OGT,OPTN,OSBP,OSBPL8,OSGIN1,P2RX4,PAG1,PALM2-AKAP2,PAM,PANK3,PANX1,PAPSS2,PARP14,PATL1,PBX4,PCDH11X,PCF11,PCGF3,PCNX1,PDE4A,PDE4B,PDE4DIP,PDGFB,PDLIM7,PEAK1,PELI1,PHEX,PHF1,PHLDB1,PIK3AP1,PIK3R5,PIM1,PIM3,PISD,PLAGL2,PLAUR,PLCG2,PLD1,PLEC,PLEK,PLEKHG2,PLIN4,PLXNC1,PLXND1,PMEPA1,PNPLA1,POLG,PPARD,PPARG,PPFIA1,PPP1R12A,PPP1R13B,PPP1R15A,PPP1R15B,PPP1R18,PRDM11,PRDM8,PRKAG2,PRLR,PRRC2C,PSEN1,PSTPIP2,PTGIR,PTGS2,PTK2B,PTPN1,PTPRE,PTPRJ,PTX3,QKI,R3HCC1L,RAB10,RAB12,RAB3IP,RAB5A,RABGEF1,RALGAPA1,RALGAPA2,RALGDS,RAP1B,RAP2C,RAPGEF1,RAPGEF2,RAPH1,RARG,RASA2,RASA3,RASGRP3,RBBP8,RBM47,RBMS1,RDX,REL,RELA,RELB,RELL1,RFFL,RFTN1,RGCC,RGPD2,RHCG,RHOH,RHOQ,RIMKLB,RIN2,RND1,RNF144B,RNF19A,RNF19B,RNF213,RNF24,ROCK1,ROCK2,RPGR,RRAD,RYBP,SAMD4B,SAMSN1,SATB1,SAV1,SBNO2,SCAF4,SCARF1,SCN1B,SDC2,SDC4,SEC22B,SEC24A,SEMA3E,SERPINB2,SERPINB9,SESTD1,SETD5,SFI1,SFR1,SGIP1,SGMS2,SGPP2,SH2D3A,SIK2,SIPA1L1,SIRPA,SKIL,SLAMF1,SLAMF7,SLC12A6,SLC16A3,SLC16A6,SLC22A4,SLC24A4,SLC25A13,SLC2A6,SLC30A4,SLC30A7,SLC35F2,SLC39A8,SLC41A2,SLC43A2,SLC43A3,SLC6A12,SLC6A6,SLC9A8,SLCO3A1,SLCO4A1,SMAD7,SMG1,SMG7-AS1,SMG9,SMOX,SMS,SOD2,SPACA6,SPATA13,SPHK1,SPIN4,SPRED2,SPSB1,SQSTM1,SRC,SS18,SSTR2,ST3GAL1,ST3GAL2,ST8SIA4,STAG2,STARD8,STAT4,STAT5A,STEAP1B,STK26,STK3,STX3,SUMO4,SUSD6,SVIL,SYNJ2,SYNPO2,TAB2,TANK,TAPBP,TBC1D23,TBC1D30,TBC1D8,TET2,TEX14,THUMPD3-AS1,TICAM1,TIFA,TJP2,TLCD2,TLE4,TLR2,TM4SF1,TM4SF19,TMEM106A,TMEM217,TMEM265,TMEM45A,TMEM54,TMEM63B,TNF,TNFAIP1,TNFAIP2,TNFAIP6,TNFAIP8,TNFRSF14,TNFRSF18,TNFRSF1B,TNFRSF4,TNFRSF9,TNFSF14,TNIK,TNIP1,TNIP2,TNIP3,TNS3,TP53BP1,TP53BP2,TPRA1,TRA2A,TRAF1,TRAF3,TRAF3IP2,TRAK1,TREML4,TRIM36,TRIP10,TSC22D2,TTL,TUBB1,TUBB6,TUG1,TWISTNB,TXNRD1,UBA1,UBALD2,UBE2E1,UBE2J1,UBE2Z,UBR4,UBR5,URGCP,USF3,USP12,USP12-AS2,USP24,USP34,USP7,USP9X,VASH2,VASP,VEGFA,VIM,VPS13A,WDR31,WIPF2,WTAP,WWC2,XIAP,XRN1,ZBTB1,ZBTB17,ZBTB5,ZC3H12A,ZC3H12C,ZCCHC14,ZEB2,ZFAT,ZFX,ZFYVE1,ZHX2,ZMIZ1,ZMIZ1-AS1,ZNF267,ZNF316,ZNF638,ZNF674,ZNF697,ZNF800,ZNFX1,ZNRF2,ZSCAN5A,ZSWIM4,ZSWIM6
# STAT3: ADAM10,ADAMTSL4-AS1,APOL6,ARID4B,ASH1L-AS1,ATP11B,B3GNT5,BAZ2A,CTBP2,DENND3,DNAJA1,ENTPD4,ERRFI1,FUT4,GPX3,HCP5,HIPK3,HSPA13,IGFBP2,IPPK,IRF1,IRF7,KBTBD2,KDM2A,KMT2B,KMT2E,MAP3K11,MCU,MLKL,MMP25-AS1,N4BP2,NAMPT,NATD1,NRDC,NUFIP2,OR10G3,PARP14,PCNX1,PICALM,PJA2,PLEKHG2,PPARD,PPP1R12C,PRNP,PTPN7,RELA,RLIM,RNF213,RUNX1,RXRA,SBNO2,SEMA4A,SLCO3A1,SMCHD1,SSBP2,STAT2,STAT3,STK3,TBX21,TDRD7,TJAP1,TNFAIP1,TP53BP1,TRAK2,TRIP12,TXLNA,UBA6,UBR1,USP32,USP6NL,VCPIP1,WDFY3,XAF1,YY1AP1,ZBTB25,ZNF222,ZNF692
# STAT5A: ADA,ADAMTSL3,ADRB2,AMPD3,AREL1,ATG7,ATP2C1,BAALC-AS1,BCL2,BCL2L14,BHLHE40,BOD1,CD40,CFB,CLCF1,CLDN5,CLIC2,CLIC4,CREBRF,CRK,DGKH,DNAJB11,DTX2,E2F7,EDEM1,FAM160A1,FAS,FOXP1,GPR132,ICAM1,IL6,INF2,IPO9-AS1,IRF8,KLF6,LINC01270,LPP,LRRC8B,MAFF,MAN2C1,MAP3K5,MAP3K8,MAST4,MFHAS1,MICALL2,MIR155HG,MMP2,MOB3C,MUC1,NEURL3,NRP2,OSBPL3,PANX1,PAPSS2,PDE4DIP,PDGFB,PIM2,PITPNB,PNISR,PNRC1,PSD,PTTG1,RANGAP1,RAPGEF2,RASGRP3,RDX,RHOH,RNF19A,RTN4RL2,SAV1,SCN1B,SH2D3A,SHC1,SIRPA,SLC12A7,SLC15A3,SNX9,SOCS2,SPRED2,SRD5A1,ST3GAL6,STARD8,STAT5A,STAT5B,TIPRL,TNFRSF9,TRAK1,TRIP10,TTC23,UBE2D1,ZC3HAV1,ZFYVE1,ZNF189,ZNF24,ZNF655

# PROGENy gene lists:
# PROGENy JAK/STAT: OAS1,HERC6,OAS3,PLSCR1,DDX60,TRIM21,SP110,DDX58,RSAD2,STAT1,IFI44,IFIT1,IFI44L,OAS2,SP100,IFI35,MX2,SAMHD1,ISG15,DDX60L,IFI6,NMI,USP18,IFIT3,IFI16,SAMD9,EIF2AK2,XAF1,MX1,BST2,IFIH1,SLC25A28,C19orf66,RTP4,TDRD7,IFI27,IFIT5,IFITM1,UBE2L6,PML,BATF2,EPSTI1,APOL2,OASL,LAP3,SLC15A3,IRF7,DHX58,HELZ2,ADAR,PARP14,TMEM62,TRIM14,SECTM1,RNF114,APOL6,ETV7,TAP2,IFIT2,PARP9,IRF9,PHF11,TAP1,TRIM5,PSMB8,TRIM38,APOL4,PARP12,SAMD9L,APOBEC3G,CMTR1,HERC5,STAT2,OGFR,GBP1,MYD88,TRIM22,TMEM140,RNF213,CMPK2,LAMP3,NLRC5,TRANK1,ZBP1,PRKD2,APOL1,PSMB9,GBP4,CYP2J2,TRIM25,PSMB8-AS1,ZC3HAV1,ZNFX1,DTX3L,CNP,SP140L,APOL3,CXCL11,PLSCR2,CASP4,IL22RA1,TRIM69,TLR3,PLEKHA4,CXCL9,WARS,GMPR,SERPING1,DUOX2,TNFSF13B,JAK2,CASP7,AIM2,CASP10,MLKL,ANKFY1,IL18BP,TRAFD1,BLZF1,BCL2L14,THEMIS2,PARP10,LGALS3BP,OR52K3P,ERAP2,N4BP1,IL15,STARD5,TRIM26,ZCCHC2,RASGRP3,NUB1,IFITM3,HLA-F,HRASLS2,HLA-E,SOCS1,ZNF107,ISG20,TNFSF10,IDO1,LY6E,REC8,FAM46A,CX3CL1,DDO,LGALS9,IFITM2,NAPA,LGMN,TRIM56,GTPBP1,ERV3-2,PPM1K,C5orf56,B2M,DNPEP,APOBEC3F,PNPT1,GCA,XRN1,LYPD5,HDX,GBP6,PDCD1LG2,CD274,CD47,TLK2,MOV10,RARRES3,CASP1,IL12RB2,IRF2,TYMP,FAM122C,GBP5,IRF1,SLFN12,EDNRA,CTNNBL1,NOD2,NOD1,CBR3,SIDT1,BTC,STX17,GNB4,GBP1P1,CHMP5,ART3,UNC93B1,IL15RA,BISPR,HAPLN3,RBCK1,RBM11,FCGR1B,BTN3A3,SP140,LGALS17A,CEACAM1,IRF8,PDZD2,CALCOCO2,RIPK1,C1S,USP42,ADPRHL2,CASP8,TLE4,PCGF5,ZFYVE26,GBP2,FBXO6,GBP7,RABGAP1L,WDR25,MSRB1,DUOXA2,UBD,TRIM40,FZD5,MAK,LOC101927027,NT5C3A,ODF3B,C2,TEX29,SPATS2L,ZBTB42,GSDMD,CD38,VAMP5,ACE2,RUBCN,IL7,CNDP2,C4orf33,SCARB2,EHD4,PSMB10,NDUFA9,MCUB,MS4A6E,ACKR4,C1R,SERPINB9P1,TMEM106A,CXorf21,BTN3A1,SHISA5,CXorf38,PHACTR4,CTSS,SNX6,HLA-DOB,HLA-C,MOB3C,UBQLNL,CFH,RNF31,ETV6,CCL8,STK3,LOC101928809,RIPK2,PROCR,PSME2,GIMAP2,RNF19B,SERPINB1,HSH2D,CD74,STOML1,ABHD16A,NUDCD1,ASPHD2,MMAA,UBA6,CCL18,STAP1,CXCL10,HLA-DRA,HLA-G,PSMA4,NAT8B,ALPK1,PSMA3,TICAM1,USF1,ZNF620,PIK3AP1,DNAJA1,BAZ2A,ATP10A,CPEB3,FMR1,DCP1A,PSME1,BTN3A2,ANGPTL1,LMO2,RP2,OR51L1,CTRL,DYNLT1,BCL2L13,GTF2B,LGALS8,USP30-AS1,AZI2,STEAP4,PARP8,LYN,SBNO2,MED25,OTUD4,FRMD4B,CIITA,RFX5,GNG5,PGLYRP4,C6orf62,SCLT1,PSMA2,LRRC3,TAPBPL,MDGA1,HLA-B,INTS12,SCO2,VSIG10L,STAT3,CAND2,PI4K2B,SOCS3,KIAA0040,ARHGAP20,CD164,MYCBP2,IL1RN,SCYL3,NBN,NECTIN2,UBFD1,HLA-J,GBP3,CCDC68,BAK1,HCP5,MAP2,PATL1,MB21D1,CLDN23,CSNK1G1,ANKRD62,C5orf15,FLT3LG,TTC39B,MAB21L2,SLC6A14,SRD5A3-AS1,FBXO39,ANKRD22,VPS9D1,SPSB1,SLC38A5,ACSM5,CCR8,RBM43,SPTLC2,APOBEC3D,UBA7,CCND3,STX11,CLEC2B,C12orf74,CREM,CARD11,LRRTM2,TTC38,DAPP1,HELB,GPR180,IFNLR1,ERAP1,IFI30,KRT24,MAD2L1BP,JADE2,IL36RN,DYRK4,LRCH2,MMP25-AS1,TOP1P2,FLVCR2,NNMT,ZRSR2,LOC100506274,C1QB,CHRNA6,ELF1,RNASE10,CTSO
# PROGENy NFkB: TNFAIP3,CXCL2,CCL20,NFKBIA,CXCL8,RELB,NFKB1,TNFAIP2,ICAM1,IL6,NFKBIE,NFKB2,CXCL3,SLC2A6,NINJ1,WTAP,BCL2A1,BIRC3,CD83,TICAM1,CXCL1,RIPK2,INHBA,CSF2,TNIP1,DRAM1,TRAF1,IRAK2,TLR2,TNF,MSC,IKBKE,GCH1,CXCL5,IL1A,BID,PTX3,EHD1,CFLAR,SDC4,TNIP2,TNFSF15,IFNGR2,MAP3K8,IL1RN,CYLD,IL23A,IL36G,CCL2,KYNU,STAT5A,B4GALT1,IL1B,TNIP3,IFNAR2,HIVEP1,NKX3-1,CFB,TRIP10,CSF3,IRF1,TNFRSF9,CXCL6,NFKBIZ,CD40,OLR1,EFNA1,C11orf96,LTB,EBI3,IL7R,UBE2Z,SOCS3,ELOVL7,TIFA,REL,IL32,NFKBIB,STX11,KMO,ANKLE2,NECAP2,TNFAIP8,TNFAIP6,PPP4R4,RNF19B,HIVEP2,SOD2,ANKRD33B,SGPP2,BIRC2,CXCL11,SLC11A2,PRDM1,SLC41A2,TMEM217,IL15,CD80,MARCKS,IFIH1,EDN1,CLIC4,RELA,DENND5A,PLAU,N4BP1,TP53BP2,DTX2,DAPP1,TAPBP,CLDN1,IL15RA,AMPD3,RND1,CD44,SERPINA3,MRGPRX3,MMP10,CD70,GPR132,SERPINB9,RHCG,PLAGL2,UXS1,ZC3H12C,CLUHP3,C3,VCAM1,NOCT,C1QTNF1,PANX1,BMP2,PPP3CC,C8orf4,PSTPIP2,DUSP16,TYMP,JUNB,DNAAF1,RAP2C,TAP1,APBA3,IER3,HLA-F,CHST2,GSAP,HMGN2P46,LINC01465,CSF1,CD82,NFE2L3,ABTB2,LYN,IRGQ,NBN,POU2F2,MMP12,PIM3,ITGB8,SQRDL,GBP3,TANK,CCL5,SLC39A8,PTGIR,CCL7,CREB3,ICOSLG,GBP2,IRAK3,SQSTM1,WNT5A,SOCS1,RNF144B,ETS1,TNC,MMP9,IL12B,PLAUR,PPP1R15A,LINC01215,IL18R1,CLEC4E,CD69,DDX58,LACTB,CYP7B1,UPB1,CDC42EP2,SRC,ARNTL2,SUSD6,OPTN,RHBDF2,LRP12,C15orf48,TNFAIP1,SBNO2,PLA2G4C,NAMPT,PDE4B,PML,ZMIZ2,CD274,PMAIP1,PTGS2,PDGFB,TSLP,CTHRC1,LOC440934,SRGN,CXCL10,GBP1,LOC100130476,C10orf55,STAT4,SLC25A37,FEZ1,IL6ST,ARL5B,BTG3,IL4I1,PARP12,USP12,LINC00158,DAXX,NIPAL4,HS3ST3B1,TRIM47,S100A12,CASP5,TMEM106A,SERPINB2,IRF7,P2RX4,CLEC2D,G0S2,CYP27B1,BDKRB1,PLK3,FMNL3,LIF,NECTIN2,USP43,SLC31A2,CASP7,LINC00936,TRAF3,CHMP4B,HCK,GRAMD1A,ADPRH,LHFP,PELI1,B4GALT5,CCL4,CCL8,PATL1,BTN2A1,MYO10,GP1BA,ISG20,CDKN2B,IGFL1,NR4A3,NEURL3,GPR84,ATP2B1,B2M,SCARF1,RILPL2,LOC100288675,LOC101929709,FBRS,ZBTB10,RAP1B,FCAR,PNPLA1,MSANTD3,MMP13,SLAMF7,JAK3,MT1M,FFAR2,PDPN,DHX58,CANT1,LOC399900,CTSS,RCAN1,GRINA,CH25H,PDGFRL,IFIT3,SIX5,CD48,KCNJ2,MAFF,CCL19,LOC101928554,KDM6B,PTAFR,SAMSN1,LAMB3,STARD5,DCUN1D3,FNDC3B,GPR176,EREG,KRT6B,TRAF2,TRIM36,ZNF697,VNN1,CD47,KLF6,BTN2A2,MREG,MESDC1,NFKBID,NCOA7,CCNL1,AKR1B1,RHBDL2,ZBTB17,PPARD,SLC25A28,DNAJB5,CSRNP1,RARRES1,PTGER4,GRAMD3,GBP5,CYTH1,GNA15,PIK3R5,HELZ2,KLF9,ZXDA,IRF2,PTPN1,PLA1A,ATF5,SYNGR3,PIK3AP1,MMP14,MTF1,CHST7,IDO1,SDCBP,BAZ1A,EGOT,TNFSF9,ARID5A,ETS2,RNF207,SLC7A2,CCR7,GBP4,ZNF267,NOD2,APOL3,SLC1A2,ZBTB46,CCL18,PARP14,CEBPD,APOL6,CD58,SLC12A7,IFIT2,GATA6,MMP19,PPP1R18,SAV1,CFLAR-AS1,PI3,OAS3,SLC9A8,LGALSL,RRAD,IL10RA,SELE,CASP4
# PROGENy MAPK: DUSP6,SPRED2,SPRY2,ETV5,EPHA2,PHLDA2,FOSL1,SPRY4,DUSP4,SLCO4A1,ETV4,PHLDA1,ETV1,SLC20A1,DUSP5,STX1A,ELK3,NT5E,PHC2,SPRED1,TCOF1,BZW1,CCND1,AREG,UBA2,TNFRSF10A,PWP2,PNP,TRIB1,GPR3,HMGA2,POLR3G,YRDC,DCBLD2,SLC35F2,LRP8,AEN,TOP1,SDCCAG3,CDCA4,FERMT1,SLC4A7,GRPEL1,XRCC3,E2F3,PRNP,LDLR,MCL1,WDR3,CBFB,MCMBP,MIER2,KRR1,NUP50,GATAD2A,BRIX1,UBASH3B,DUSP7,RNF138,KMT5A,STK17A,SMURF2,HMGA1,MPP6,CAPN15,TEX10,TAF1A,KCTD5,MMP1,QSOX2,GTPBP4,SH2B3,RNF126,LIF,TSR1,RRP1,NIP7,MCFD2,SERPINB8,SLC25A32,EGR1,SPRY4-IT1,DDX21,CCDC86,USP36,RASSF8,PUM3,QTRT2,CD3EAP,TNFRSF12A,PITPNC1,EIF4E,PDSS1,NPAS2,PLK3,MRGBP,UBIAD1,TEX30,NUP153,MYEOV,BLM,SH3TC2,MAK16,SLC25A15,IER2,PVR,ODC1,POLR1E,TOMM34,RFWD3,ITGA6,LONRF3,SEH1L,TRMT6,ANO1,CRY1,NLE1,UBE2S,GART,PPAT,NOC3L,IPO7,MAGOHB,NDUFAF4,SLC45A3,WDR74,LINC00973,URB2,PFDN2,FIP1L1,FJX1,FOXD1,DSCC1,TGFA,PTPN12,SNRPB,CDC45,NAA15,RFC3,WDR43,CDC42EP2,IER3,ORC6,IRAK1,SLC35G2,MYO19,RPP25,MCM10,CRLF3,ECT2,MYC,GNL2,ZNF696,MAFF,HPCAL1,ESF1,TAGLN3,CHD1,UPP1,STEAP1,SLC25A19,MAD2L1,BYSL,PPRC1,AMD1,NOM1,WNK4,ERRFI1,NCL,POLR2D,LYAR,CEP72,KBTBD2,FAM84B,CCNE1,EIF5A2,PLAUR,DHODH,ENTPD6,TFAM,CHAF1A,ARPC5L,ISG20L2,ZNF473,TRA2B,GJB3,PLEK2,MRTO4,NOLC1,NSUN2,NOL6,S100A16,SMTN,EIF2B2,FBXO45,TFB1M,TMEM185B,SOWAHC,HBEGF,SSFA2,KIF23,CDK5R1,DAZAP1,ITPR3,NUP54,TCERG1,NOP16,DIAPH1,GEMIN4,PPM1G,UCHL5,ZNF215,ENOPH1,RAD18,AXIN1,UTP18,UBE3C,DNAJC2,WWTR1,HES1,PRMT3,PNO1,ATAD2,HAUS6,HEATR1,HSPA14,PSMC3IP,SMAGP,PUS1,STK10,HJURP,EGR3,KCNN4,TXLNG,PXN,CHSY1,ERCC6L,HAS3,WDR12,TLR4,SURF2,HNRNPAB,NTSR1,TOMM40,ARHGEF2,IFRD2,MLX,PAK1IP1,SNRPD1,SYNCRIP,DCAF13,RRP1B,TIPIN,TRIP13,CDC25A,RRM2,RIF1,KPNA4,RBM15,POLR3D,NOP14,ACOT7,RAN
# PROGENy Hypoxia: PDK1,ANKRD37,FAM162A,BNIP3L,ENO2,ANKZF1,NDRG1,INSIG2,FUT11,KDM3A,PFKFB4,GBE1,EGLN1,PGK1,HILPDA,PLOD2,RLF,ALDOC,ERO1A,BHLHE40,ALDOA,BNIP3,GPI,PLOD1,SLC2A1,KDM4B,ZNF395,MIR210HG,P4HA1,PPFIA4,EGLN3,VLDLR,VEGFA,EFNA3,PFKFB3,MXI1,WDR45B,WSB1,LDHA,AK4,NCKIPSD,PGM1,SLC25A36,YEATS2,ANGPTL4,TMEM45A,NGLY1,VKORC1,SEC61G,ENO1,C8orf58,LOC154761,C4orf47,MAPK7,PDK3,CA9,ZNF292,GAPDH,PFKP,NRN1,SPAG4,RIOK3,P4HA2,C4orf3,ZNF160,WDR54,MPI,CLK3,PFKL,ADM,NOL3,DARS,PKM,PGAM1,RAB20,SAP30,LOX,RBPJ,CRKL,NARF,FBXO42,RORA,GYS1,RNF24,ZNF654,PRSS53,FAM13A,PPP2R5B,PAM,ZBTB25,CCNG2,LOXL2,KDM4C,BCKDK,TBC1D9B,TPI1,STC2,SLC35E1,CNOT8,HK1,PLEKHA2,RIMKLA,RSBN1,BEND5,RAB11FIP5,MAP2K1,XPNPEP1,GOLGA8A,FAM13A-AS1,DDX41,FEM1C,BCKDHA,PPP1R13L,ITPR1,SH3D21,SYDE1,MIF,RRAGD,RNASET2,GPR146,JMJD6,WDR60,CXCR4,PPP1R3B,ANG,GUK1,ALKBH5,RNASE4,LDHC,S100A10,ARFGEF3,ARID3A,PIGA,SFXN3,RASSF7,KDM5B,PCAT6,SLC2A3,DPCD,UBE2O,CEP250,NFIL3,PPP1R3C,SMAP1,LNPK,DPYSL4,KCTD11,BHLHE41,SERTAD2,CLK1,ZMYM2,SMYD2,EFEMP2,RNMT,CRYBB2P1,TMSB10,PLXNA3,HK2,ITGA5,PPME1,HCG4,FGF11,SLC6A8,CSRP2,PLIN2,FKBP15,LRP1,OBSL1,NDUFA4L2,NDNF,FAM210A,CECR5,PNRC1,PAPD7,KISS1R,ARRDC3,INHA,METTL21B,LGALS1,IPMK,JAM2,SEMA4B,LRP2BP,CIART,COL6A2,IGFBP3,PIAS2,SERGEF,RYBP,TNIP1,UPK1A,SCAND2P,QSOX1,TMEM74B,BBX,MGEA5,FLNB,C3orf58,TPD52,PRKAA2,SNAPC3,ING1,CRABP2,UNC119,MOB3A,FOSL2,PRELID2,NPEPPS,KDM7A,BTG1,CDK19,STC1,PLAGL1,KLF7,MFSD10,PHF21A,LOC100630923,TET1,CABIN1,SNTA1,P4HB,LGALS8,MZF1,HDAC3,RNF113A,BMT2,NBAS,PNMA2,NUP58,RRAGA,MTFP1,B3GNT4,MKNK2,OXSR1,LONP1,EHD2,SLC6A6,ARHGAP45


