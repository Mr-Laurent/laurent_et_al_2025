# app_Macs_scatter
setwd("~/laurent_et_al_2025/Viz_Apps/")

#
library(shiny)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(dplyr)
library(glue)
library(ComplexHeatmap)
library(Matrix)


load("MacsEnr_z_mymod_150_metadata_goodlog.rd") 
load("Macs_dgCmc_16nov22_counts.rd")      # dgcmtxcounts    matrix of gene counts by metacells
hm_zmod<-as.matrix(metadf_z[,c(2:3,(length(colnames(metadf_z))-150):length(colnames(metadf_z)))])
colnames(hm_zmod)<-gsub("scaled","",colnames(hm_zmod))
hm_zmod<-as.data.frame(hm_zmod)
hm_zmod$metacell<-rownames(hm_zmod)

load("genelist_16nov22_Macs.rd")  # genelist, colnames should be "Mod" & "Gens"
genelist<-gnlscenics
mnpcounts<-dgcmtxcounts



pal=rev(rainbow(15))[4:15]
linebreaks <- function(n){HTML(strrep(br(), n))}

findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}


#library(edgeR)
#data.cpm<-cpm(mnpcounts, normalized.lib.sizes=TRUE, log=F)   # If log=TRUE : calcul is detailled here: https://support.bioconductor.org/p/107719/
# I'M USING 1e05 INSTEAD OF 1e06 AS MEAN METACELL READ COUNT IS ~ 100 000
data.cpm<-1e+05*sweep(mnpcounts,2,colSums(mnpcounts),`/`)  
l.data.cpm<-log2(1+data.cpm)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  navbarPage("Myelo_enrich",
             tabPanel("Classic",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(4, 
                               wellPanel(
                                 h4("V2 vs V3 modules score distribution"),
                                 selectInput("var", label = "Choose a module to display", choices=colnames(hm_zmod),selected="My_Mod_1"),
                                 #sliderInput(inputId = "dotSize", label = "Set point size", value=0.3, min=0, max=1,step=0.1),
                                 #linebreaks(4),
                                 textInput("geneask", h4("Looking for a gene ?"),value = "no gene selected"),
                                 verbatimTextOutput("vgen"),
                                 linebreaks(3),
                                 h4("ScatterPlot parameters"),
                                 #https://stackoverflow.com/questions/36906265/how-to-color-sliderbar-sliderinput
                                 tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #6e6e6e}")),
                                 tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: #bdbdbd}")),
                                 tags$style(HTML(".js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {background: #bdbdbd}")),
                                 tags$style(HTML(".js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {background: #6e6e6e}")),             
                                 sliderInput(inputId = "Uhei", label = "Height", value=500, min=0, max=1200,step=100),
                                 sliderInput(inputId = "Uwid", label = "Width", value=500, min=0, max=1200,step=100),
                                 sliderInput(inputId = "Upts", label = "Dot size", value=1, min=0, max=5,step=0.1),
                                 selectInput("um_x", label = "Module as UMAP x axis", choices=colnames(hm_zmod),selected="SigMac_Mono2"),
                                 selectInput("um_y", label = "Module as UMAP y axis", choices=colnames(hm_zmod),selected="SigMac_Mono1"),
                                 # textInput("um_x", h4("Module as UMAP x axis"),value = "14"),
                                 # textInput("um_y", h4("Module as UMAP y axis"),value = "16")
                               )
                        ),
                        mainPanel(
                          #plotOutput("vln"),
                          textOutput("text"),
                          verbatimTextOutput("verb"),
                          tags$hr(),
                          tabsetPanel(type = "tabs",
                                      # tabPanel("ScatterPlot", plotOutput("sct")),
                                      tabPanel("Custom ScatterPlot", plotOutput("sct2"))
                          )
                          
                        )
                      )
             ),
             tabPanel("z-scored",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Z-scaled module score by V2 or V3 cells, displayed as mean in metacells"),
                                 selectInput("varz", label = "Choose a module to display", choices=colnames(hm_zmod),selected="My_Mod_1"),
                                 textInput("um_xz", h4("Module as UMAP x axis"),value = "SigMac_Mono2"),
                                 textInput("um_yz", h4("Module as UMAP y axis"),value = "SigMac_Mono1"),
                                 h4("ScatterPlot parameters"),
                                 #https://stackoverflow.com/questions/36906265/how-to-color-sliderbar-sliderinput
                                 tags$style(HTML(".js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {background: #bdbdbd}")),
                                 tags$style(HTML(".js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {background: #bdbdbd}")),
                                 tags$style(HTML(".js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {background: #6e6e6e}")),
                                 sliderInput(inputId = "Uheiz", label = "Height", value=400, min=0, max=1200,step=100),
                                 sliderInput(inputId = "Uwidz", label = "Width", value=700, min=0, max=1200,step=100),
                                 sliderInput(inputId = "Uptsz", label = "Dot size", value=1.5, min=0, max=5,step=0.1),
                                 selectInput("typlotz", label = "Plot design", choices=c("Dark_rnbw","Dark_viri","Light_rnbw_alpha"),selected="Light_rnbw_alpha")
                               )
                        ),
                        mainPanel(
                          textOutput("textz"),
                          verbatimTextOutput("verbz"),
                          tags$hr(),
                          plotOutput("sctzvall"),
                          tags$hr(),
                          plotOutput("sctzv2")
                        )
                      )
                      
             ),
             tabPanel("Genes",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(4, 
                               wellPanel(
                                 h4("Personnal gene or modules on metacells"),
                                 h5("Chose a single gene, multiple genes for a module score, or one of the following parameters: 
                                    nCount_RNA, nFeature_RNA, V2vV3ratio, PBMCratio, V1prop, V2prop, V3prop, cellnbr, Nantes_prop, Powrie_prop, Sinai_prop, clustannot, clustSC"),
                                 textInput("select3", h4("Gene, Module or Param selection:"),value = "nCount_RNA"),
                                 selectInput("um_x3", label = "Module as UMAP x axis", choices=colnames(hm_zmod),selected="SigMac_Mono2"),
                                 selectInput("um_y3", label = "Module as UMAP y axis", choices=colnames(hm_zmod),selected="SigMac_Mono1"),
                                 # textInput("um_x3", h4("Module as UMAP x axis"),value = "My_Mod_21"),
                                 # textInput("um_y3", h4("Module as UMAP y axis"),value = "My_Mod_39"),
                                 h4("ScatterPlot parameters"),
                                 #https://stackoverflow.com/questions/36906265/how-to-color-sliderbar-sliderinput
                                 tags$style(HTML(".js-irs-7 .irs-single, .js-irs-7 .irs-bar-edge, .js-irs-7 .irs-bar {background: #6e6e6e}")),             
                                 sliderInput(inputId = "Uhei3", label = "Height", value=700, min=0, max=1200,step=100),
                                 sliderInput(inputId = "Uwid3", label = "Width", value=700, min=0, max=1200,step=100),
                                 sliderInput(inputId = "Upts3", label = "Dot size", value=1.5, min=0, max=5,step=0.1),
                                 selectInput("typlot3", label = "Plot design", choices=c("Dark_rnbw","Dark_viri"),selected="Dark_viri"),
                                 selectInput("rawlog3", label = "One gene expression", choices=c("Log cpm 1e05","Cpm 1e05 on MC exp raw","Log","Raw"),selected="Log cpm 1e05")
                               )
                        ),
                        mainPanel(
                          #tags$a(href="https://mr-laurent.shinyapps.io/MyeloEnrich_v1_scatterfind", "Click here to have the working personnalized scatterplot"),
                          plotOutput("sct3")
                          
                        )
                      )
                      
             )
             
             
  ))

hm_env = new.env()

server <- function(input, output) {
  output$vln <- renderPlot({
    ggplot(minidf, aes_string(x = "chem", y = input$var)) +
      geom_violin(draw_quantiles = c(.25, .5, .75),aes(fill=chem))+theme_minimal()+
      geom_jitter(shape=16, position=position_jitter(0.2),size=input$dotSize)+geom_boxplot(width=rev(0.2), fill="white")
  })
  output$text <- renderText({paste0(input$var," genes: ")})
  output$verb <- renderText({genelist$Gens[which(genelist$Mod==input$var)] })
  output$vgen <- renderText({
    if(input$geneask%in%unlist(strsplit(genelist$Gens,split=","))){
      print(paste0("Module ",findgene(strsplit(genelist$Gens,split=","),input$geneask)))
    }else{print("Not found")}
  })
  #Can't have reactive sizes for the plotting, so get them in a reactive function
  u_hei<-function(){return(input$Uhei)}
  u_wid<-function(){return(input$Uwid)}
  
  output$sct2 <- renderPlot({
    ggplot(hm_zmod, aes(x=!!sym(paste0(input$um_x)) , y=!!sym(paste0(input$um_y)) , color=hm_zmod[,input$var])) +
      geom_point(size=input$Upts)+
      ggtitle(paste("projection of ",input$var," score in Macs metacells",sep=""))+
      scale_color_gradientn(colours = pal,name="score")+theme_dark()
  },height=u_hei, width=u_wid)
  
  # Tab 2 : z-scores
  
  output$textz <-renderText({paste0(input$varz," genes: ")})
  output$verbz <- renderText({genelist$Gens[which(genelist$Mod==varz)] })
  u_heiz<-function(){return(input$Uheiz)}
  u_widz<-function(){return(input$Uwidz)}
  u_widz2<-function(){return(2*input$Uwidz)}
  
  #v2sc_Mod_136
  
  output$sctzvall <- renderPlot({
    pvall<-ggplot(metadf_z, aes(x=!!sym(paste0("My_Mod_",input$um_xz,"score")) , y=!!sym(paste0("My_Mod_",input$um_yz,"score")) , 
                                color=(metadf_z$`V2prop`*metadf_z[,paste0("v2sc_Mod_",input$varz)]%>%(function(x) { x[is.na(x)] <- 0; return(x) }))+(metadf_z$`V3prop`*metadf_z[,paste0("v3sc_Mod_",input$varz)]%>%(function(x) { x[is.na(x)] <- 0; return(x) }))   )) +
      geom_point(size=input$Uptsz)+
      ggtitle(paste("projection of My_Mod_",input$varz," V3cells mean of z-score in Macs metacells ",sep=""))
    
    if(input$typlotz=="Dark_rnbw"){
      pvallb<-pvall+scale_color_gradientn(colours = pal)+theme_dark()+ labs(color = "score")
    }else if(input$typlotz=="Dark_viri"){
      pvallb<-pvall+scale_color_viridis()+theme_dark()+ labs(color = "score")
    }else if(input$typlotz=="Light_rnbw_alpha"){
      pvallb<-ggplot(metadf_z, aes(x=!!sym(paste0("My_Mod_",input$um_xz,"score")) , y=!!sym(paste0("My_Mod_",input$um_yz,"score")) , 
                                   color=(metadf_z$`V2prop`*metadf_z[,paste0("v2sc_Mod_",input$varz)]%>%(function(x) { x[is.na(x)] <- 0; return(x) }))+(metadf_z$`V3prop`*metadf_z[,paste0("v3sc_Mod_",input$varz)]%>%(function(x) { x[is.na(x)] <- 0; return(x) }))  )) +
        geom_point(size=input$Uptsz)+
        ggtitle(paste("projection of My_Mod_",input$varz," z-score corrected in Macs metacells",sep=""))+
        scale_color_gradientn(colours = pal)+ labs(color = "score",alpha="V3 prop")  
    }
    pvallb
    
  },height=u_heiz, width=u_widz)
  
  output$sctzv2 <- renderPlot({
    pv2<-ggplot(metadf_z, aes(x=!!sym(paste0("My_Mod_",input$um_xz,"score")) , y=!!sym(paste0("My_Mod_",input$um_yz,"score")) , color=metadf_z[,paste0("v2sc_Mod_",input$varz)])) +
      geom_point(size=input$Uptsz)+theme(legend.text = "score")+
      ggtitle(paste("projection of My_Mod_",input$varz," V2cells mean of z-score in Macs metacells",sep=""))
    
    if(input$typlotz=="Dark_rnbw"){
      pv2b<-pv2+scale_color_gradientn(colours = pal)+theme_dark()+ labs(color = "score")
    }else if(input$typlotz=="Dark_viri"){
      pv2b<-pv2+scale_color_viridis()+theme_dark()+ labs(color = "score")
    }else if(input$typlotz=="Light_rnbw_alpha"){
      pv2b<-ggplot(metadf_z, aes(x=!!sym(paste0("My_Mod_",input$um_xz,"score")) , y=!!sym(paste0("My_Mod_",input$um_yz,"score")) , color=metadf_z[,paste0("v2sc_Mod_",input$varz)],alpha=metadf_z[,"V2prop"])) +
        geom_point(size=input$Uptsz)+
        ggtitle(paste("projection of My_Mod_",input$varz," V2cells mean of z-score in Macs metacells",sep=""))+
        scale_color_gradientn(colours = pal)+ labs(color = "score",alpha="V2 prop") 
    }
    
    pv3<-ggplot(metadf_z, aes(x=!!sym(paste0("My_Mod_",input$um_xz,"score")) , y=!!sym(paste0("My_Mod_",input$um_yz,"score")) , color=metadf_z[,paste0("v3sc_Mod_",input$varz)])) +
      geom_point(size=input$Uptsz)+
      ggtitle(paste("projection of My_Mod_",input$varz," V3cells mean of z-score in Macs metacells ",sep=""))
    
    if(input$typlotz=="Dark_rnbw"){
      pv3b<-pv3+scale_color_gradientn(colours = pal)+theme_dark()+ labs(color = "score")
    }else if(input$typlotz=="Dark_viri"){
      pv3b<-pv3+scale_color_viridis()+theme_dark()+ labs(color = "score")
    }else if(input$typlotz=="Light_rnbw_alpha"){
      pv3b<-ggplot(metadf_z, aes(x=!!sym(paste0("My_Mod_",input$um_xz,"score")) , y=!!sym(paste0("My_Mod_",input$um_yz,"score")) , color=metadf_z[,paste0("v3sc_Mod_",input$varz)],alpha=metadf_z[,"V3prop"])) +
        geom_point(size=input$Uptsz)+
        ggtitle(paste("projection of My_Mod_",input$varz," V3cells mean of z-score in Macs metacells",sep=""))+
        scale_color_gradientn(colours = pal)+ labs(color = "score",alpha="V3 prop")  
    }
    pv2b+pv3b
  },height=u_heiz, width=u_widz2)
  
  
  
  # Tab 3 : z-scores heatmap
  
  
  hm_heiz<-function(){return(input$HMheiz)}
  hm_widz<-function(){return(input$HMwidz)}
  
  output$hmzall <- renderPlot({
    # Conditional removal of modules :
    if(input$selectrm==""){hm_zmod_rm<-hm_zmod}else{
      hm_zmod_rm<-hm_zmod[,-eval(parse(text=paste0("c(",input$selectrm,")")))] }
    #Compute the heatmap figure
    set.seed(42)
    hm_env$ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                             row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),
                             column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = input$HMkc,row_km = input$HMkr ) )
    hm_env$ht_pos = ht_pos_on_device(hm_env$ht)
  }) #,height=hm_heiz, width=hm_widz)
  
  output$hmzoom = renderPlot({
    if(is.null(input$ht_brush)) {
      grid.newpage()
      grid.text("No region is selected.", 0.5, 0.5)
    } else {
      lt = ComplexHeatmap:::get_pos_from_brush(input$ht_brush)
      pos1 = lt[[1]]
      pos2 = lt[[2]]
      
      ht = hm_env$ht
      pos = selectArea(ht, mark = FALSE, pos1 = pos1, pos2 = pos2, 
                       verbose = FALSE, ht_pos = hm_env$ht_pos)
      
      row_index = unlist(pos[1, "row_index"])
      column_index = unlist(pos[1, "column_index"])
      m = ht@ht_list[[1]]@matrix
      set.seed(42)
      ht_select = Heatmap(m[row_index, column_index, drop = FALSE],
                          col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
                          show_heatmap_legend = FALSE,
                          cluster_rows = FALSE, cluster_columns = FALSE)
      draw(ht_select)
      
    }
  })
  
  output$ht_click_content = renderText({
    if(is.null(input$ht_click)) {
      "Not selected."
    } else {
      pos1 = ComplexHeatmap:::get_pos_from_click(input$ht_click)
      lt = ComplexHeatmap:::get_pos_from_brush(input$ht_brush)
      posm1 = lt[[1]]
      posm2 = lt[[2]]
      posmod <- selectArea(ht, mark = FALSE, pos1 = posm1, pos2 = posm2, 
                           verbose = FALSE, ht_pos = hm_env$ht_pos)
      
      ht = hm_env$ht
      pos = selectPosition(ht, mark = FALSE, pos = pos1, 
                           verbose = FALSE, ht_pos = hm_env$ht_pos)
      
      row_index = pos[1, "row_index"]
      column_index = pos[1, "column_index"]
      row_mod=unlist(posmod[1,"row_index"])
      test1= input$ht_brush$coords_css$xmin
      test2= input$ht_brush$coords_css$ymin
      m = ht@ht_list[[1]]@matrix
      v = m[row_index, column_index]
      glue('Selected modules: {paste(row_mod, collapse=",")}',
           "xmin: {test1}",
           "ymin: {test2}",
           #"row index: {row_index}",
           #"column index: {column_index}",
           "value: {v}", .sep = "\n")
      
    }
  })
  
  
  
  u_hei3<-function(){return(input$Uhei3)}
  u_wid3<-function(){return(input$Uwid3)}
  output$sct3 <- renderPlot({
    ilist<-unlist(strsplit(input$select3,split = ","))
    auth<-c("nCount_RNA","nFeature_RNA","clustSC","V2vV3ratio","PBMCratio","V1prop","V2prop","V3prop","cellnbr","Nantes_prop","Powrie_prop","Sinai_prop")
    hm_zmod$currentscore<-NA
    if(length(ilist)<2){
      if(ilist%in%auth){ #Avant il était avant le else if <2 mais ça marchait plsu sur les version de R/shiny récentes donc je le rentre dedans pour pas de soucis
        eval(parse(text=paste0("hm_zmod$currentscore[match(hm_zmod$metacell,rownames(hm_zmod))]<-hm_zmod$",ilist)))
        sct3p<-ggplot(hm_zmod, aes(x=!!sym(paste0(input$um_x3)) , y=!!sym(paste0(input$um_y3)) , color=hm_zmod[,"currentscore"])) +
          geom_point(size=input$Upts3)+
          ggtitle(paste("projection of ",ilist," in Macs metacells",sep=""))
      }else{
        if(input$rawlog3=="Cpm 1e05 on MC exp raw"){
          # Take the "raw" counts
          eval(parse(text=paste0("hm_zmod$currentscore[match(hm_zmod$metacell,colnames(mnpcounts))]<-data.cpm[ilist,,drop=F]")))
        }else if(input$rawlog3=="Log cpm 1e05"){
          eval(parse(text=paste0("hm_zmod$currentscore[match(hm_zmod$metacell,colnames(mnpcounts))]<-l.data.cpm[ilist,,drop=F]")))
        }else if(input$rawlog3=="Raw"){
          eval(parse(text=paste0("hm_zmod$currentscore[match(hm_zmod$metacell,colnames(mnpcounts))]<-t(mnpcounts[ilist,])")))
        }else if(input$rawlog3=="Log"){
          eval(parse(text=paste0("hm_zmod$currentscore[match(hm_zmod$metacell,colnames(mnpcounts))]<-log(1+t(mnpcounts[ilist,]))")))
        }
        sct3p<-ggplot(hm_zmod, aes(x=!!sym(paste0(input$um_x3)) , y=!!sym(paste0(input$um_y3)) , color=hm_zmod[,"currentscore"])) +
          geom_point(size=input$Upts3)+
          ggtitle(paste("projection of ",ilist," expression in Macs metacells",sep=""))
      }
    }else{
      values<-Matrix::colSums(log(1+mnpcounts[ilist,,drop=F]))/Matrix::colSums(log(1+mnpcounts))
      eval(parse(text=paste0("hm_zmod$currentscore[match(hm_zmod$metacell,names(values))]<-values")))
      sct3p<-ggplot(hm_zmod, aes(x=!!sym(paste0(input$um_x3)) , y=!!sym(paste0(input$um_y3)) , color=hm_zmod[,"currentscore"])) +
        geom_point(size=input$Upts3)+
        ggtitle(paste("projection of the current module score in Macs metacells",sep=""))
    }
    
    if(input$typlot3=="Dark_rnbw"){
      sct3pb<-sct3p+scale_color_gradientn(colours = pal)+theme_dark()+ labs(color = "score")
    }else if(input$typlot3=="Dark_viri"){
      sct3pb<-sct3p+scale_color_viridis()+theme_dark()+ labs(color = "score")
    }
    
    sct3pb
    
  },height=u_hei3, width=u_wid3)
  
}

shinyApp(ui = ui, server = server)



