library(cowplot)
library(dplyr)
library(gginnards)
library(ggplot2)
library(ggpp)
library(ggpubr)
library(grid)
library(gridExtra)
library(gtable)
library(patchwork)
library(reshape2)
library(rstatix)
library(tidyverse)


# Load the metadata and count matrix from Ngollo cohort
load("./Grouped_objects/RNAseq_REMIND/minipheno_AllezNgollo_postop_goodRemRec_sig.rd")
load("./Grouped_objects/genelist_16nov22_Macs.rd")

miniphenoM0I<-minipheno[which(minipheno$`location:ch1`%in%c("M0I")),]
miniphenoCTRL<-minipheno[which(minipheno$`location:ch1`%in%c("Ctrl")),]


# signatures <- c("data_Siggrp_I","data_Siggrp_II","data_Siggrp_III","data_Siggrp_IV","data_Siggrp_V",
#                 "data_Siggrp_VI","data_Siggrp_VII","data_Siggrp_VIII","data_Siggrp_IX","data_Siggrp_X",   
#                 "data_Siggrp_XI","data_Siggrp_XII","data_Siggrp_XIII","data_Siggrp_XIV","data_Siggrp_XV",  
#                 "data_Siggrp_XVI","data_Siggrp_XVII")
signatures <- c("data_Siggrp_VIII","data_Siggrp_XV")
levels(signatures)<-signatures





##---------------------------------------##
##--##  Camembert of Rutgeerts prop  ##--##
##---------------------------------------##

camembert_rut<-function(data_set,cat_prop,col_rut,col_rut_txt,plot_title){
  
  mltab<-melt(table(data_set$EnrichSupInf,data_set$new_rut) )
  perc_df <-mltab[which(mltab$Var2==cat_prop),c(1,3)]%>% mutate(perc = `value` / sum(`value`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  # Reorder row order + levels to always have the same order
  perc_df <- perc_df[match(names(col_rut_txt), perc_df$Var1), ]
  perc_df$Var1 <- factor(perc_df$Var1, levels = levels(perc_df$Var1)[match(names(col_rut_txt), levels(perc_df$Var1))])
  return(ggplot(perc_df, aes(x="", y=value, fill=Var1))+
           coord_polar("y", start=0)+
           geom_bar(stat="identity", size=rel(0.8), color="black") + 
           scale_fill_manual(values=col_rut ) +  theme_void()+
           theme(axis.text.x=element_blank(),
                 plot.title=element_text(size=rel(1.4), face="bold"),
                 legend.title = element_text(size=rel(1.4)),
                 legend.text = element_text(size=rel(1.2))) +
           geom_text(aes(label = paste0(labels,"\nn=",value)),color=col_rut_txt,size=rel(4),
                     position = position_stack(vjust = 0.5))+guides(fill = guide_legend(title = "Status")) +
           ggtitle(plot_title)
  )
}

##---------------------------------------##
##--##   Table of rutgeerts counts   ##--##
##---------------------------------------##

count_rutgeerts<-function(data_set,fnt_size=24){
  
  g2 <- tableGrob(table(data_set$EnrichSupInf,data_set$new_rut),theme=ttheme_minimal(base_size = fnt_size) )
  g2 <- gtable_add_grob(g2,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = nrow(g2), l = 2, r = ncol(g2))
  g2 <- gtable_add_grob(g2,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 1, l = 2, r = ncol(g2) )
  g2$widths <- 1.2*g2$widths
  plt2<-ggplot() +
    annotation_custom(grob = g2, xmin = 0.9, xmax = 0, ymin = -Inf, ymax = Inf) +
    #ggtitle(paste0("Proportions of\n",val1," x ",val2,"\ndouble positive samples"))+ 
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5,vjust = -rel(6)))
  
  return(plt2)
  
}


##---------------------------------------##



all_sig_comp_CAMrutCT<-function(sub_set){
  n <- length(signatures)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      x <- signatures[i]
      y <- signatures[j]
      lgd1=paste0(gsub("data_Siggrp_","",x))
      lgd2=paste0(gsub("data_Siggrp_","",y))
      print(paste0(lgd1," ",lgd2))
      
      # Distrib DP vs Others in each Rec rutgeerts score (i2a, i2b, i3, i4)
      eval(parse(text=paste0("data_set=minipheno",sub_set)))
      data_set$EnrichSupInf <- with(data_set, ifelse(eval(parse(text=paste0(x)))>= max(eval(parse(text=paste0("miniphenoCTRL$",x))))
                                                     & eval(parse(text=paste0(y)))>= max(eval(parse(text=paste0("miniphenoCTRL$",y)))) , "Double Pos","Others"))
      # To avoid cases where one condition is 0% of data, force the labels so melt() will create a 0 for them instead of deleting it
      data_set$EnrichSupInf<-factor(x=data_set$EnrichSupInf, levels = c("Others","Double Pos"))
      data_set$new_rut<-data_set$`rutgeerts:ch1`
      data_set<-data_set[-which(data_set$new_rut%in%c("0","1")),]
      levels(data_set$EnrichSupInf)<-c("Others","Double Pos")
      
      cat_prop="i2a"
      col_rut<-c("Others"="#E4C784","Double Pos"="#d19a02")
      col_rut_txt<-c("Others"="black","Double Pos"="black")
      i2a_dist<-camembert_rut(data_set=data_set,cat_prop=cat_prop,col_rut=col_rut,col_rut_txt=col_rut_txt,
                              plot_title=paste0(cat_prop," distribution in ",lgd1," & ",lgd2) )
      cat_prop="i2b"
      col_rut<-c("Others"="#D3A46C","Double Pos"="#ba7000")
      col_rut_txt<-c("Others"="black","Double Pos"="black")
      i2b_dist<-camembert_rut(data_set=data_set,cat_prop=cat_prop,col_rut=col_rut,col_rut_txt=col_rut_txt,
                              plot_title=paste0(cat_prop," distribution in ",lgd1," & ",lgd2) )
      cat_prop="i3"
      col_rut<-c("Others"="#C77B53","Double Pos"="#b34b02")
      col_rut_txt<-c("Others"="white","Double Pos"="white")
      i3_dist<-camembert_rut(data_set=data_set,cat_prop=cat_prop,col_rut=col_rut,col_rut_txt=col_rut_txt,
                             plot_title=paste0(cat_prop," distribution in ",lgd1," & ",lgd2) )
      cat_prop="i4"
      col_rut<-c("Others"="#A2594B","Double Pos"="#801c01")
      col_rut_txt<-c("Others"="white","Double Pos"="white")
      i4_dist<-camembert_rut(data_set=data_set,cat_prop=cat_prop,col_rut=col_rut,col_rut_txt=col_rut_txt,
                             plot_title=paste0(cat_prop," distribution in ",lgd1," & ",lgd2) )
      cnt_2ab34<-count_rutgeerts(data_set=data_set,fnt_size=rel(12))
      
      # Distrib DP vs Others in Mild (i2a+i2b) vs Severe (i3+i4)
      eval(parse(text=paste0("data_set=minipheno",sub_set)))
      data_set$EnrichSupInf <- with(data_set,ifelse(eval(parse(text=paste0(x)))>= max(eval(parse(text=paste0("miniphenoCTRL$",x))))
                                                    & eval(parse(text=paste0(y)))>= max(eval(parse(text=paste0("miniphenoCTRL$",y)))), "Double Pos","Others"))
      # To avoid cases where one condition is 0% of data, force the labels so melt() will create a 0 for them instead of deleting it
      data_set$EnrichSupInf<-factor(x=data_set$EnrichSupInf, levels = c("Others","Double Pos"))
      data_set$new_rut<-data_set$`rutgeerts:ch1`
      data_set<-data_set[-which(data_set$new_rut%in%c("0","1")),]
      data_set$new_rut[which(data_set$new_rut%in%c("i2a","i2b"))]<-"i2"
      data_set$new_rut[which(data_set$new_rut%in%c("i3","i4"))]<-"3_4"
      levels(data_set$EnrichSupInf)<-c("Others","Double Pos")
      
      cat_prop="i2"
      col_rut<-c("Others"="#DDB87A","Double Pos"="#C28019")
      col_rut_txt<-c("Others"="black","Double Pos"="black")
      mild_dist<-camembert_rut(data_set=data_set,cat_prop=cat_prop,col_rut=col_rut,col_rut_txt=col_rut_txt,
                               plot_title=paste0(cat_prop," distribution in ",lgd1," & ",lgd2) )
      cat_prop="3_4"
      col_rut<-c("Others"="#B76D4F","Double Pos"="#9D360D")
      col_rut_txt<-c("Others"="white","Double Pos"="white")
      severe_dist<-camembert_rut(data_set=data_set,cat_prop=cat_prop,col_rut=col_rut,col_rut_txt=col_rut_txt,
                                 plot_title=paste0(cat_prop," distribution in ",lgd1," & ",lgd2) )
      cnt_mi_sev<-count_rutgeerts(data_set=data_set,fnt_size=rel(12))
      
      
      # print((((i2a_dist/i3_dist)|(i2b_dist/i4_dist))/cnt_2ab34)|(mild_dist/severe_dist/cnt_mi_sev))
      design_perso="
      AACCFF
      AACCFF
      BBDDGG
      BBDDGG
      #EE#HH
      #EE#HH
      "
      
      print(i2a_dist+i3_dist+i2b_dist+i4_dist+cnt_2ab34+mild_dist+severe_dist+cnt_mi_sev+plot_layout(design = design_perso) )
      
    }}}


pdf(paste0("./Figures_print/FigS5Ha_M0I_DP_CamRut_tmCT.pdf"),width =9,height = 7)
all_sig_comp_CAMrutCT("M0I")
dev.off()

