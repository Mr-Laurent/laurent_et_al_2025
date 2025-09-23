#### app_XeniumSlideViewer_v1_2.R : added annotation selection, change cluster order etc 
# 1_3 : Switch from h5ad objects to rd objects for sharing the app
# 1_3_1 : clean code + comments
# 1_3_2 : new page for subclusters, can also take additive gene expression ("=GENE+GENE+GENE")
# 1_3_3 : new page for subclusters from an annotation file ( + fix color scale in signature / genes) => only in winodw 3 as of 30/04: need to implement it for windows 1 & 2
# 1_3_4 : new secret use for page 3 (log norm '_XXX' gene expression, '#XXXXXX' color of subtype) + Lineages buttons for easy selection in page 3 + coord_fix option

# THOMAS : 
googlepath="~/Viz_Apps/"


setwd(paste0(googlepath))

library(Matrix)
library(ggplot2)
library(reticulate)
library(dplyr)
library(scales)
library(viridis)
library(shiny)
`%ni%` <- Negate(`%in%`)

options(shiny.maxRequestSize = 500 * 1024^2)



###----------------------------------------------------------### 
#### Function to transform list of genes to vector readable #### 
#### by the program, with exclusion of unknown genes        ####
###----------------------------------------------------------### 
gene_input2vect_t<-function(genelist,dgcmtx){
  gnshow<-unlist(strsplit(genelist,split=","))
  gnshow<-gsub(" ","",gnshow)
  notfound<-setdiff(gnshow,colnames(dgcmtx))
  print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  return(setdiff(gnshow,notfound) )
}



pal_hiro<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0","#74BCD4","#548FAC","#396794","#20466D")


ui <- fluidPage(
  navbarPage("Xenium slide viewer v1.3.4",
             tabPanel("Gene and signatures",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 fileInput(inputId = "file", "Upload Xenium dgcmtx_XXX.rd files", accept = ".rd"),
                                 textInput("geneask", h4("Gene / Signature\nto plot"),placeholder = "no gene selected"),
                                 verbatimTextOutput("notf_gen"),
                                 actionButton("click", "Show slide"),
                                 sliderInput(inputId = "bx_hei", label = "Height", value=800, min=0, max=2000,step=50),
                                 sliderInput(inputId = "bx_wid", label = "Width", value=1000, min=0, max=2000,step=50),
                                 sliderInput(inputId = "bx_siz", label = "Dot size", value=0.2, min=0, max=3,step=0.05),
                                 sliderInput(inputId = "bx_interval", label = "Color interval", value=c(-2,2), min=-5, max=5,step=0.5),
                                 checkboxInput("select_scale","Scale on subset", value=TRUE),
                                 downloadButton('export_plot'),
                                 verbatimTextOutput("obsm_names"),
                                 uiOutput("chk_bx"),
                                 actionButton("select_all", "All"),
                                 actionButton("select_none", "None")
                               )
                        ),
                        mainPanel(
                          
                          plotOutput("slide_plt", height = 800, width = 1000),
                          
                        )
                      )
                      
             ),
             tabPanel("Same but w/ sublineages\n(Jake annotations)",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 #                     # fileInput(inputId = "file", "Upload Xenium dgcmtx_XXX.rd files", accept = ".rd"),
                                 textInput("geneask2", h4("Gene / Signature\nto plot"),placeholder = "no gene selected"),
                                 verbatimTextOutput("notf_gen2"),
                                 actionButton("click2", "Show slide"),
                                 sliderInput(inputId = "bx_hei2", label = "Height", value=800, min=0, max=2000,step=50),
                                 sliderInput(inputId = "bx_wid2", label = "Width", value=1000, min=0, max=2000,step=50),
                                 sliderInput(inputId = "bx_siz2", label = "Dot size", value=0.2, min=0, max=3,step=0.05),
                                 sliderInput(inputId = "bx_interval2", label = "Color interval", value=c(-2,2), min=-5, max=5,step=0.5),
                                 checkboxInput("select_scale2","Scale on subset", value=TRUE),
                                 downloadButton('export_plot2'),
                                 verbatimTextOutput("obsm_names2"),
                                 uiOutput("chk_bx2"),
                                 actionButton("select_all2", "All"),
                                 actionButton("select_none2", "None")
                               )
                        ),
                        mainPanel(
                          
                          plotOutput("slide_plt2", height = 800, width = 1000),
                          
                        )
                      )
                      
             ),
             tabPanel("Same but w/ sublineages\n(my annotations)",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 fileInput(inputId = "file_annot", "Upload subtype annotation file", accept = ".rd"),
                                 div(
                                   style = "display: flex; align-items: center;",
                                   textInput("geneask3", h4("Gene / Signature\nto plot"),placeholder = "no gene selected"),
                                   actionLink("helpLink", "(?)", style = "margin-left: 5px; font-size: 20px;")
                                 ),
                                 verbatimTextOutput("notf_gen3"),
                                 actionButton("click3", "Show slide"),
                                 sliderInput(inputId = "bx_hei3", label = "Height", value=800, min=0, max=2000,step=50),
                                 sliderInput(inputId = "bx_wid3", label = "Width", value=1000, min=0, max=2000,step=50),
                                 checkboxInput("select_coordfix","Fix coordinates", value=TRUE),
                                 sliderInput(inputId = "bx_siz3", label = "Dot size", value=0.2, min=0, max=3,step=0.05),
                                 sliderInput(inputId = "bx_interval3", label = "Color interval", value=c(-2,2), min=-5, max=15,step=0.5),
                                 # checkboxInput("select_scale3","Scale on subset", value=TRUE),
                                 actionButton("toggle_scale3", "Scaled on subset.\nClick to scale on all"),
                                 downloadButton('export_plot3'),
                                 verbatimTextOutput("obsm_names3"),
                                 uiOutput("chk_bx3"),
                                 actionButton("select_all3", "All"),
                                 actionButton("select_none3", "None"),
                                 uiOutput("dynamic_buttons3")
                               )
                        ),
                        mainPanel(
                          
                          plotOutput("slide_plt3", height = 800, width = 1000),
                          
                        )
                      )
                      
             )
             
  ))

hm_env = new.env()

server <- function(input, output, session) {
  # Reactive object (count matrix + metadata dataframe) initiation
  adata <- reactiveVal(NULL)
  meta_d <- reactiveVal(NULL)
  meta_ann <- reactiveVal(NULL)
  meta_names <- reactiveVal(NULL)
  grouplin_names <- reactiveVal(NULL)
  
  observeEvent(input$file_annot, {
    req(input$file_annot)
    tryCatch({  # Load the manual annotations of subtypes 
      load(input$file_annot$datapath, envir = hm_env)
      
      loaded_medata <- get("meta_annot_perso", envir = hm_env) 
      meta_ann(loaded_medata)
      showNotification("Annotations loaded successfully!", type = "message")
      
      group_names3 <- unique(meta_ann()$my_annots)
      output$chk_bx3 <- renderUI({
        checkboxGroupInput("groups3", "Select Groups:",
                           choices = sort(group_names3),
                           selected = group_names3)
      })
    }, error = function(e) {
      showNotification("Error loading file: Make sure it's a valid count object .rd", type = "error")
      meta_ann(NULL)
    })
    tryCatch({  # Add the lineage buttons for easier selections
      rd_file_path_3 <- paste0(googlepath,"./../Grouped_objects/Xenium/meta_annot_",gsub("meta_annot_([^_]*)_.*$","\\1",input$file_annot$name), "_names.rd")
      load(rd_file_path_3, envir = hm_env)
      print("Catch2done")
      df_re4 <- get("df_re4", envir = hm_env)
      meta_names(df_re4)  # Store the loaded meta names in reactiveVal
      print(unique(meta_names()$group)[!is.na(unique(meta_names()$group))])
      grplin<-unique(meta_names()$group)[!is.na(unique(meta_names()$group))]
      grouplin_names(grplin)
    }, error = function(e) {
      showNotification("Meta names not found.", type = "warning")
      meta_names(NULL)
      grouplin_names(NULL)
    })
  })
  
  output$dynamic_buttons3 <- renderUI({  # Create additionnal buttons to None and All, with lineage groups defined in the file meta_annot_names_XXX.rd
    req(grouplin_names())
    tagList(
      lapply(grouplin_names(), function(name) {
        actionButton(inputId = paste0("btn_", name), label = name)
      })
    )
  })
  
  
  observeEvent(input$file, {
    req(input$file)  # Ensure a file is selected
    rd_file_path <- paste0(googlepath,"./../Grouped_objects/Xenium/meta_annot_",gsub("dgcmtx_(.*).rd$","\\1",input$file$name), ".rd")
    
    tryCatch({ # Try to load a count matrix
      load(input$file$datapath, envir = hm_env)
      print(input$file$datapath)
      loaded_data <- get("dgcmtxcounts", envir = hm_env) 
      adata(loaded_data)  # Store the loaded data in reactiveVal
      print(dim(loaded_data))
      print('is ok')
      print(dim(adata() ))
      showNotification("Anndata loaded successfully!", type = "message")
      print((rd_file_path))
      if (file.exists(rd_file_path)) { # Try to load the metadata table associated to the count matrix
        load(rd_file_path, envir = hm_env)
        meta_annot <- get("meta_annot", envir = hm_env)
        meta_annot$cell_id<-rownames(meta_annot)
        meta_annot$currentscore<-NA
        meta_d(meta_annot)  # Store the loaded meta data in reactiveVal
        print(unique(meta_d()$Ann_lin_AndJak))
        showNotification("Meta object file loaded successfully!", type = "message")
        group_names <- unique(meta_d()$Ann_lin_AndJak)
        output$chk_bx <- renderUI({
          checkboxGroupInput("groups", "Select Groups:",
                             choices = group_names,
                             selected = group_names)
        })
        group_names2 <- unique(meta_d()$Annotations)
        output$chk_bx2 <- renderUI({
          checkboxGroupInput("groups2", "Select Groups:",
                             choices = group_names2,
                             selected = group_names2)
        })
        
      } else {
        showNotification("Meta object file not found.", type = "warning")
        meta_d(NULL)  # Reset meta_object if the file is not found
      }
      
    }, error = function(e) {
      showNotification("Error loading file: Make sure it's a valid count object .rd", type = "error")
      adata(NULL)  # Reset adata if an error occurs
      meta_d(NULL)
    })
  })
  
  # DEBUG: Check that the object is loaded by displaying its size
  output$obsm_names <- renderPrint({
    req(adata())  # Ensure data is loaded before accessing
    print(paste0(dim(adata())[1]," cells\n",dim(adata())[2]," probes"))
  })
  output$obsm_names2 <- renderPrint({
    req(adata())  # Ensure data is loaded before accessing
    print(paste0(dim(adata())[1]," cells\n",dim(adata())[2]," probes"))
  })
  output$obsm_names3 <- renderPrint({
    req(adata())  # Ensure data is loaded before accessing
    print(paste0(dim(adata())[1]," cells\n",dim(adata())[2]," probes"))
  })
  # All vs None buttons to check/uncheck all
  observeEvent(input$select_all, {
    if (!is.null(meta_d())) {
      group_names <- unique(meta_d()$Ann_lin_AndJak)
      updateCheckboxGroupInput(session, "groups", selected = group_names)
    }
  })
  observeEvent(input$select_none, {
    updateCheckboxGroupInput(session, "groups", selected = character(0))
  })
  observeEvent(input$select_all2, {
    if (!is.null(meta_d())) {
      group_names2 <- unique(meta_d()$Annotations)
      updateCheckboxGroupInput(session, "groups2", selected = group_names2)
    }
  })
  observeEvent(input$select_none2, {
    updateCheckboxGroupInput(session, "groups2", selected = character(0))
  })
  observeEvent(input$select_all3, {
    if (!is.null(meta_ann())) {
      group_names3 <- unique(meta_ann()$my_annots)
      updateCheckboxGroupInput(session, "groups3", selected = group_names3)
    }
  })
  observeEvent(input$select_none3, {
    updateCheckboxGroupInput(session, "groups3", selected = character(0))
  })
  observe({   # For each additionnal lineage button, allows to check the associated annotations when clicked
    lapply(grouplin_names(), function(name) {
      observeEvent(input[[paste0("btn_", name)]], {
        group_names3 <- unique(meta_ann()$my_annots)
        updateCheckboxGroupInput(session, "groups3", selected = group_names3[group_names3%in%meta_names()$my_annots[which(meta_names()$group==name)]] )
      })
    })
  })
  
  observeEvent(input$helpLink, {
    showModal(modalDialog(
      title = "Help Information",
      HTML("<h3>How to use the search bar:</h3><br/>
           <b>- '#000000'</b> <br/>
           Enter an Hex code to color the selected subtypes<br/>
           <br/>
           <b>- 'CD14'</b><br/>
           Enter a single gene name to show this gene expression in the selected subtypes<br/>
           <br/>
           <b>- '_CD14'</b><br/>
           Enter a single gene preceded by a '_' to show the log2 1+x normalized expression (scale factor: 10^4)<br/>
           <br/>
           <b>- 'GZMB,PRF1,GNLY,IFNG,NKG7,GZMA,CXCR3'</b><br/>
           Enter multiple genes separated by commas to compute an enrichment score<br/>
           <br/>
           <b>- '=EPCAM+AGR3+GPX2+MMP7'</b><br/>
           Enter multiple genes with '=' and separated by '+' to get an additive expression score<br/><br/><br/>
           "),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  # Output text to inform if genes put by user are not found
  output$notf_gen = renderText({
    gnshow<-unlist(strsplit(input$geneask,split=","))
    gnshow<-gsub(" ","",gnshow)
    notfound<-setdiff(gnshow,colnames(adata() ))
    print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  })
  output$notf_gen2 = renderText({
    gnshow<-unlist(strsplit(input$geneask2,split=","))
    gnshow<-gsub(" ","",gnshow)
    notfound<-setdiff(gnshow,colnames(adata() ))
    print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  })
  output$notf_gen3 = renderText({
    gnshow<-unlist(strsplit(input$geneask3,split=","))
    gnshow<-gsub(" ","",gnshow)
    notfound<-setdiff(gnshow,colnames(adata() ))
    print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  })
  
  
  # Initiation of plot parameters
  bxhei<-function(){return(input$bx_hei)}
  bxwid<-function(){return(input$bx_wid)}
  bxsiz<-function(){return(input$bx_siz)}
  interval_min<-function(){return(input$bx_interval[1])}
  interval_max<-function(){return(input$bx_interval[2])}
  bxhei2<-function(){return(input$bx_hei2)}
  bxwid2<-function(){return(input$bx_wid2)}
  bxsiz2<-function(){return(input$bx_siz2)}
  interval_min2<-function(){return(input$bx_interval2[1])}
  interval_max2<-function(){return(input$bx_interval2[2])}
  bxhei3<-function(){return(input$bx_hei3)}
  bxwid3<-function(){return(input$bx_wid3)}
  bxsiz3<-function(){return(input$bx_siz3)}
  interval_min3<-function(){return(input$bx_interval3[1])}
  interval_max3<-function(){return(input$bx_interval3[2])}
  
  output$slide_plt <- renderPlot({
    plot_data()
  },height=bxhei, width=bxwid)
  output$slide_plt2 <- renderPlot({
    plot_data2()
  },height=bxhei2, width=bxwid2)
  output$slide_plt3 <- renderPlot({
    plot_data3()
  },height=bxhei3, width=bxwid3)
  
  plot_data<- eventReactive(input$click,{ # T prevent print at first lauch
    req(meta_d())  # Ensure data is loaded before accessing
    req(input$groups)  # Ensure groups are selected
    select_lin <- input$groups
    
    isolate({
      allsum<-Matrix::rowSums(adata())
      print(input$geneask)
      # if(grep("^=",input$geneask)==1){
      if(grepl('^=', input$geneask) ){
        print("y a un =")
        gn2sum<-unlist(strsplit(gsub("^=","",input$geneask),split="\\+"))
        legd<-paste0("Additive expression")
        
        meta_data<-meta_d()
        # Scale on the subsetted cells or not
        if(input$select_scale==T){
          df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$Ann_lin_AndJak%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F]))[which(meta_data$Ann_lin_AndJak%in%select_lin)])
          meta_data$currentscore<-df_inter$score[match(meta_data$cell_id,df_inter$cell_name)]
        }else{
          meta_data$currentscore[which(meta_data$Ann_lin_AndJak%in%select_lin)]<-Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F])[which(meta_data$Ann_lin_AndJak%in%select_lin)]
        }
        
        brplt_title<-paste0("projection of current score in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
        
        plot_sc<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(size = 1*input$bx_siz) +
          scale_color_gradientn(colours = rev(pal_hiro),na.value ="#EDEDED", limits=c(input$bx_interval[1],input$bx_interval[2]), oob=squish)+
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        if(input$select_scale==T){
          # plot_score <- plot_sc+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
          pal_hiro_min<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0")
          plot_score <- plot_sc+scale_color_gradientn(colors = rev(pal_hiro_min),na.value ="#EDEDED" )
        }else{
          pal_hiro_min<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0")
          # plot_score <- plot_sc+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max(meta_data$currentscore) ) )
          plot_score <- plot_sc+scale_color_gradientn(colors = rev(pal_hiro_min),na.value ="#EDEDED", limits = c(0, max(meta_data$currentscore) ) )
        }
      }else{
        
        
        ilist<-gene_input2vect_t(input$geneask,adata())
        print(ilist)
        if(length(ilist)==1){   # If only 1 gene is selected, display its flat expression
          print("ilist est de 1")
          legd<-paste0(ilist," \nExpression")
          meta_data<-meta_d()  # Need to get out of the reactive variable or else I will have "Error in <-: invalid (NULL) left side of assignment" 
          # eval(parse(text=paste0("meta_data$currentscore<-as.vector(as.data.frame(adata())[,which(colnames(adata())==ilist)])")))                              ############### CHANGE IT : matrix subset before dataframe ??
          eval(parse(text=paste0("meta_data$currentscore<-as.vector(adata()[,which(colnames(adata())==ilist)])")))  
          max_cur_sco<-max(meta_data$currentscore)
          meta_data$currentscore[which(meta_data$Ann_lin_AndJak%ni%select_lin)]<-NA
          brplt_title<-paste0("projection of ",ilist," expression in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
          
          plot_sc<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
            geom_point(data = subset(meta_data, is.na(currentscore) ), aes(x = x, y = y), 
                       color = "#EDEDED", size = 0.2*input$bx_siz) +
            geom_point(data = subset(meta_data, currentscore == 0), aes(x = x, y = y), 
                       color = "#D6D6D6", size = 0.5*input$bx_siz) +
            geom_point(data = subset(meta_data, currentscore > 0), 
                       aes(color = currentscore), size = 1*input$bx_siz)  +
            theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            labs(title = brplt_title,
                 color = legd )+
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank())
          if(input$select_scale==T){
            plot_score <- plot_sc+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
          }else{
            plot_score <- plot_sc+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
          }
          
          
        }else if(length(ilist)>=2){  # If multiple genes are selected, compute and plot its score of enrichment
          print("ilist is more than 1")
          meta_data<-meta_d()
          # Scale on the subsetted cells or not
          if(input$select_scale==T){
            df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$Ann_lin_AndJak%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$Ann_lin_AndJak%in%select_lin)])
            df_inter$scale_score<-scale(df_inter$score)
            meta_data$currentscore<-df_inter$scale_score[match(meta_data$cell_id,df_inter$cell_name)]
          }else{
            meta_data$currentscore[which(meta_data$Ann_lin_AndJak%in%select_lin)]<-scale(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$Ann_lin_AndJak%in%select_lin)]
          }
          
          legd<-"Signature\nenrichment"
          brplt_title<-paste0("projection of current score in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
          
          plot_score<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
            geom_point(size = 1*input$bx_siz) +
            scale_color_gradientn(colours = rev(pal_hiro),na.value ="#EDEDED", limits=c(input$bx_interval[1],input$bx_interval[2]), oob=squish)+
            theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            labs(title = brplt_title,
                 color = legd )+
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank())
          
        }
      }
      print(plot_score)
      vals$plt1<-plot_score
      
      
    })
  })
  
  vals <- reactiveValues(plt1=NULL)
  
  output$export_plot = downloadHandler(  # Save the plot as a PDF
    filename = function() {paste0("plots_",Sys.Date(),".pdf")},
    content = function(file) {
      pdf(file, onefile = TRUE, width = 0.016*input$bx_wid, height = 0.016*input$bx_hei)
      print(vals$plt1)
      
      dev.off()
    })
  
  plot_data2<- eventReactive(input$click2,{ # T prevent print at first lauch
    req(meta_d())  # Ensure data is loaded before accessing
    req(input$groups2)  # Ensure groups are selected
    select_lin <- input$groups2
    
    isolate({
      allsum<-Matrix::rowSums(adata())
      ilist<-gene_input2vect_t(input$geneask2,adata())
      print(ilist)
      if(length(ilist)==1){   # If only 1 gene is selected, display its flat expression
        print("ilist est de 1")
        legd<-paste0(ilist," \nExpression")
        meta_data<-meta_d()  # Need to get out of the reactive variable or else I will have "Error in <-: invalid (NULL) left side of assignment" 
        # eval(parse(text=paste0("meta_data$currentscore<-as.vector(as.data.frame(adata())[,which(colnames(adata())==ilist)])")))                              ############### CHANGE IT : matrix subset before dataframe ??
        eval(parse(text=paste0("meta_data$currentscore<-as.vector(adata()[,which(colnames(adata())==ilist)])")))  
        max_cur_sco<-max(meta_data$currentscore)
        meta_data$currentscore[which(meta_data$Annotations%ni%select_lin)]<-NA
        brplt_title<-paste0("projection of ",ilist," expression in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
        
        plot_sc2<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(data = subset(meta_data, is.na(currentscore) ), aes(x = x, y = y), 
                     color = "#EDEDED", size = 0.2*input$bx_siz2) +
          geom_point(data = subset(meta_data, currentscore == 0), aes(x = x, y = y), 
                     color = "#D6D6D6", size = 0.5*input$bx_siz2) +
          geom_point(data = subset(meta_data, currentscore > 0), 
                     aes(color = currentscore), size = 1*input$bx_siz2)  +
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        if(input$select_scale2==T){
          plot_score2 <- plot_sc2+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
        }else{
          plot_score2 <- plot_sc2+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
        }
        
        
        
      }else if(length(ilist)>=2){  # If multiple genes are selected, compute and plot its score of enrichment
        print("ilist is more than 1")
        meta_data<-meta_d()
        # Scale on the subsetted cells or not
        if(input$select_scale2==T){
          df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$Annotations%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$Annotations%in%select_lin)])
          df_inter$scale_score<-scale(df_inter$score)
          meta_data$currentscore<-df_inter$scale_score[match(meta_data$cell_id,df_inter$cell_name)]
        }else{
          meta_data$currentscore[which(meta_data$Annotations%in%select_lin)]<-scale(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$Annotations%in%select_lin)]
        }
        
        legd<-"Signature\nenrichment"
        brplt_title<-paste0("projection of current score in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
        
        plot_score2<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(size = 1*input$bx_siz2) +
          scale_color_gradientn(colours = rev(pal_hiro),na.value ="#EDEDED", limits=c(input$bx_interval2[1],input$bx_interval2[2]), oob=squish)+
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        
      }
      print(plot_score2)
      vals$plt2<-plot_score2
      
      
    })
  })
  
  vals <- reactiveValues(plt2=NULL)
  
  output$export_plot2 = downloadHandler(  # Save the plot as a PDF
    filename = function() {paste0("plots_",Sys.Date(),".pdf")},
    content = function(file) {
      pdf(file, onefile = TRUE, width = 0.016*input$bx_wid2, height = 0.016*input$bx_hei2)
      print(vals$plt2)
      
      dev.off()
    })
  
  
  ##################
  #### Window 3 ####
  ##################
  
  leg_sc3_state <- reactiveVal("Scaled on subset.\nClick to scale on all")
  observeEvent(input$toggle_scale3, {
    # Cycle through the three states
    new_state3 <- switch(leg_sc3_state(),
                         "Scaled on subset.\nClick to scale on all" = "Scaled on all.\nClick to scale manually",
                         "Scaled on all.\nClick to scale manually" = "Scaled manually.\nClick to scale on subset",
                         "Scaled manually.\nClick to scale on subset" = "Scaled on subset.\nClick to scale on all") 
    leg_sc3_state(new_state3) 
    # Update button label based on state
    new_label3 <- new_state3
    updateActionButton(session, "toggle_scale3", label = new_label3)
  })
  
  plot_data3<- eventReactive(input$click3,{ # T prevent print at first lauch
    req(meta_d())  # Ensure data is loaded before accessing
    req(meta_ann())
    req(input$groups3)  # Ensure groups are selected
    select_lin <- input$groups3
    
    isolate({
      allsum<-Matrix::rowSums(adata())
      ilist<-gene_input2vect_t(input$geneask3,adata())
      print(ilist)
      if(grepl('^#', input$geneask3)){
        
        print("color")
        meta_data<-meta_ann()  # Need to get out of the reactive variable or else I will have "Error in <-: invalid (NULL) left side of assignment" 
        meta_da<-meta_d()
        meta_data$cell_id<-rownames(meta_data)
        meta_data$x<-meta_da$x[match(rownames(meta_data),rownames(meta_da) )]
        meta_data$y<-meta_da$y[match(rownames(meta_data),rownames(meta_da) )]
        meta_data$currentscore[which(meta_data$my_annots%ni%select_lin)]<-FALSE
        meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-TRUE
        brplt_title<-paste0("projection of selected celltypes in ", gsub("dgcmtx_(.*).rd$","\\1",input$file_annot$name) ," slide")
        
        plot_sc3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(data = subset(meta_data, currentscore==FALSE ), aes(x = x, y = y), 
                     color = "#EDEDED", size = 0.2*input$bx_siz3) +
          geom_point(data = subset(meta_data, currentscore ==TRUE), aes(x = x, y = y),  
                     color = as.character(input$geneask3), size = 1*input$bx_siz3)  +
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),legend.position = "none")
        
        if(input$select_coordfix==T){
          plot_sc3 <- plot_sc3+coord_fixed()
        }
        plot_score3 <- plot_sc3
        
      }else if(grepl('^=', input$geneask3) ){
        print("y a un =")
        gn2sum<-unlist(strsplit(gsub("^=","",input$geneask3),split="\\+"))
        legd<-paste0("Additive expression")
        
        meta_data<-meta_ann()  
        meta_da<-meta_d()
        meta_data$cell_id<-rownames(meta_data)
        meta_data$x<-meta_da$x[match(rownames(meta_data),rownames(meta_da) )]
        meta_data$y<-meta_da$y[match(rownames(meta_data),rownames(meta_da) )]
        # Scale on the subsetted cells or not
        # if(input$select_scale==T){
        #   df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$my_annots%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F]))[which(meta_data$my_annots%in%select_lin)])
        #   meta_data$currentscore<-df_inter$score[match(meta_data$cell_id,df_inter$cell_name)]
        # }else{
        #   meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F])[which(meta_data$my_annots%in%select_lin)]
        # }
        
        if (leg_sc3_state() == "Scaled on subset.\nClick to scale on all") {
          meta_data$currentscore<-NA
          df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$my_annots%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F]))[which(meta_data$my_annots%in%select_lin)])
          meta_data$currentscore<-df_inter$score[match(meta_data$cell_id,df_inter$cell_name)]
        } else if (leg_sc3_state() == "Scaled on all.\nClick to scale manually") { # If other selection: score is computed on all, but legend scale value is either 0 to max of total score or manually set 
          meta_data$currentscore<-NA
          meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F])[which(meta_data$my_annots%in%select_lin)]
        }else if (leg_sc3_state() == "Scaled manually.\nClick to scale on subset") { # If other selection: score is computed on all, but legend scale value is either 0 to max of total score or manually set 
          meta_data$currentscore<-NA
          meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F])[which(meta_data$my_annots%in%select_lin)]
        }
        
        max_cur_sco<-max( Matrix::rowSums(adata()[,which(colnames(adata())%in%gn2sum),drop=F]) )
        
        brplt_title<-paste0("projection of current score in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
        
        plot_sc3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(size = 1*input$bx_siz) +
          # scale_color_gradientn(colours = rev(pal_hiro),na.value ="#EDEDED", limits=c(input$bx_interval3[1],input$bx_interval3[2]), oob=squish)+
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        
        if(input$select_coordfix==T){
          plot_sc3 <- plot_sc3+coord_fixed()
        }
        # if(input$select_scale==T){
        #   plot_score <- plot_sc+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
        # }else{
        #   
        #   plot_score <- plot_sc+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max(meta_data$currentscore) ) )
        # }
        pal_hiro_min<-c("#000005","#230058","#5F0C6E","#B12362","#E56157","#ED894D","#FECF75","#DEDEDE")
        if (leg_sc3_state() == "Scaled on subset.\nClick to scale on all") {
          # plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
          plot_score3 <- plot_sc3+scale_color_gradientn(colors = rev(pal_hiro_min), na.value ="#EDEDED" )
        } else if (leg_sc3_state() == "Scaled on all.\nClick to scale manually") {
          # plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
          plot_score3 <- plot_sc3+scale_color_gradientn(colors = rev(pal_hiro_min), na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
        } else if (leg_sc3_state() == "Scaled manually.\nClick to scale on subset"){
          # plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits= c(input$bx_interval3[1],input$bx_interval3[2]), oob = scales::squish )
          plot_score3 <- plot_sc3+scale_color_gradientn(colors = rev(pal_hiro_min), na.value ="#EDEDED", limits= c(input$bx_interval3[1],input$bx_interval3[2]), oob = scales::squish )
        }
        
        
      }else if(grepl('^_', input$geneask3) ){
        print("y a un _")
        gn2look<-unlist(strsplit(gsub("^_","",input$geneask3),split="\\+"))
        legd<-paste0("log2 1+x \nnormalized expression\n(scale factor: 10^4")
        
        meta_data<-meta_ann()  
        meta_da<-meta_d()
        meta_data$cell_id<-rownames(meta_data)
        meta_data$x<-meta_da$x[match(rownames(meta_data),rownames(meta_da) )]
        meta_data$y<-meta_da$y[match(rownames(meta_data),rownames(meta_da) )]
        #  log2(1+(10000*count_matrix[,"CD8A"]/Matrix::rowSums(count_matrix) ))   Normalize a bit like NormalizeData() in Seurat
        eval(parse(text=paste0("meta_data$currentscore<-as.vector( log2(1+(10000*adata()[,which(colnames(adata())==gn2look)]/allsum)) )")))  
        max_cur_sco<-max(meta_data$currentscore)
        meta_data$currentscore[which(meta_data$my_annots%ni%select_lin)]<-NA
        brplt_title<-paste0("projection of ",ilist," normalized expression in ", gsub("dgcmtx_(.*).rd$","\\1",input$file_annot$name) ," slide")
        
        plot_sc3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(data = subset(meta_data, is.na(currentscore) ), aes(x = x, y = y), 
                     color = "#EDEDED", size = 0.2*input$bx_siz3) +
          geom_point(data = subset(meta_data, currentscore == 0), aes(x = x, y = y), 
                     color = "#D6D6D6", size = 0.5*input$bx_siz3) +
          geom_point(data = subset(meta_data, currentscore > 0), 
                     aes(color = currentscore), size = 1*input$bx_siz3)  +
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        
        if(input$select_coordfix==T){
          plot_sc3 <- plot_sc3+coord_fixed()
        }
        
        if (leg_sc3_state() == "Scaled on subset.\nClick to scale on all") {
          plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
        } else if (leg_sc3_state() == "Scaled on all.\nClick to scale manually") {
          plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
        } else if (leg_sc3_state() == "Scaled manually.\nClick to scale on subset"){
          plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits= c(input$bx_interval3[1],input$bx_interval3[2]), oob = scales::squish  )
          
        }
        
        
      }else if(length(ilist)==1){   # If only 1 gene is selected, display its flat expression
        print("ilist est de 1")
        legd<-paste0(ilist," \nExpression")
        meta_data<-meta_ann()  # Need to get out of the reactive variable or else I will have "Error in <-: invalid (NULL) left side of assignment" 
        meta_da<-meta_d()
        meta_data$cell_id<-rownames(meta_data)
        meta_data$x<-meta_da$x[match(rownames(meta_data),rownames(meta_da) )]
        meta_data$y<-meta_da$y[match(rownames(meta_data),rownames(meta_da) )]
        # eval(parse(text=paste0("meta_data$currentscore<-as.vector(as.data.frame(adata())[,which(colnames(adata())==ilist)])")))                              ############### CHANGE IT : matrix subset before dataframe ??
        eval(parse(text=paste0("meta_data$currentscore<-as.vector(adata()[,which(colnames(adata())==ilist)])")))  
        max_cur_sco<-max(meta_data$currentscore)
        meta_data$currentscore[which(meta_data$my_annots%ni%select_lin)]<-NA
        brplt_title<-paste0("projection of ",ilist," expression in ", gsub("dgcmtx_(.*).rd$","\\1",input$file_annot$name) ," slide")
        
        plot_sc3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          geom_point(data = subset(meta_data, is.na(currentscore) ), aes(x = x, y = y), 
                     color = "#EDEDED", size = 0.2*input$bx_siz3) +
          geom_point(data = subset(meta_data, currentscore == 0), aes(x = x, y = y), 
                     color = "#D6D6D6", size = 0.5*input$bx_siz3) +
          geom_point(data = subset(meta_data, currentscore > 0), 
                     aes(color = currentscore), size = 1*input$bx_siz3)  +
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        
        if(input$select_coordfix==T){
          plot_sc3 <- plot_sc3+coord_fixed()
        }
        
        # if(input$select_scale3==T){
        #   plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
        # }else{
        #   plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
        # }
        
        if (leg_sc3_state() == "Scaled on subset.\nClick to scale on all") {
          plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED" )
        } else if (leg_sc3_state() == "Scaled on all.\nClick to scale manually") {
          plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
        } else if (leg_sc3_state() == "Scaled manually.\nClick to scale on subset"){
          plot_score3 <- plot_sc3+scale_color_viridis_c(option = "magma", direction = -1,na.value ="#EDEDED", limits= c(input$bx_interval3[1],input$bx_interval3[2]), oob = scales::squish  )
          
        }
        
        
        
      }else if(length(ilist)>=2){  # If multiple genes are selected, compute and plot its score of enrichment
        print("ilist is more than 1")
        meta_data<-meta_ann()
        meta_da<-meta_d()
        meta_data$cell_id<-rownames(meta_data)
        meta_data$x<-meta_da$x[match(rownames(meta_data),rownames(meta_da) )]
        meta_data$y<-meta_da$y[match(rownames(meta_data),rownames(meta_da) )]
        # Scale on the subsetted cells or not
        # if(input$select_scale3==T){
        #   df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$my_annots%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$my_annots%in%select_lin)])
        #   df_inter$scale_score<-scale(df_inter$score)
        #   meta_data$currentscore<-df_inter$scale_score[match(meta_data$cell_id,df_inter$cell_name)]
        # }else{
        #   meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-scale(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$my_annots%in%select_lin)]
        # }
        # 
        
        if (leg_sc3_state() == "Scaled on subset.\nClick to scale on all") {
          meta_data$currentscore<-NA
          df_inter<-data.frame("cell_name"=meta_data$cell_id[which(meta_data$my_annots%in%select_lin)],"score"=(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$my_annots%in%select_lin)])
          df_inter$scale_score<-scale(df_inter$score)
          meta_data$currentscore<-df_inter$scale_score[match(meta_data$cell_id,df_inter$cell_name)]
        } else if (leg_sc3_state() == "Scaled on all.\nClick to scale manually") {
          meta_data$currentscore<-NA
          meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-scale(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$my_annots%in%select_lin)]
        } else if (leg_sc3_state() == "Scaled manually.\nClick to scale on subset"){
          meta_data$currentscore<-NA
          meta_data$currentscore[which(meta_data$my_annots%in%select_lin)]<-scale(Matrix::rowSums(adata()[,which(colnames(adata())%in%ilist),drop=F])/allsum)[which(meta_data$my_annots%in%select_lin)]
          
        }
        
        
        legd<-"Signature\nenrichment"
        brplt_title<-paste0("projection of current score in ", gsub("dgcmtx_(.*).rd$","\\1",input$file$name) ," slide")
        
        plot_score3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
          # geom_point(size = 1*input$bx_siz3) +
          
          geom_point(data = subset(meta_data, is.na(currentscore) ), aes(x = x, y = y), 
                     color = "#EDEDED", size = 0.2*input$bx_siz3) +
          geom_point(data = subset(meta_data, !is.na(currentscore)), 
                     aes(color = currentscore), size = 1*input$bx_siz3)   + 
          
          scale_color_gradientn(colours = rev(pal_hiro), limits=c(input$bx_interval3[1],input$bx_interval3[2]), oob=squish)+ #na.value ="#EDEDED",
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = legd )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        if(input$select_coordfix==T){
          plot_score3 <- plot_score3+coord_fixed()
        }
        
      }
      
      print(plot_score3)
      vals$plt3<-plot_score3
      
      
    })
  })
  
  vals <- reactiveValues(plt3=NULL)
  
  output$export_plot = downloadHandler(  # Save the plot as a PDF
    filename = function() {paste0("plots_",Sys.Date(),".pdf")},
    content = function(file) {
      pdf(file, onefile = TRUE, width = 0.016*input$bx_wid3, height = 0.016*input$bx_hei3)
      print(vals$plt3)
      
      dev.off()
    })
  
}


shinyApp(ui = ui, server = server)













# Infl_Mono1 sig remaining genes: GATA2,CD80,IL7R,IDO1,SPP1,EDN1,MET,CCL5,CSF3,CD274,SLAMF1,STEAP4,IL3RA




