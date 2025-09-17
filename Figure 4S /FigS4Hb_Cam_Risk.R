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

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")

load("./Grouped_objects/r_dsgn_InfUninf_cd68norm_Momacsscores.rd")
r_dsgn_Inf<-r_dsgn[which(r_dsgn$status%in%c("Inflamed")),]
r_dsgn_Uninf<-r_dsgn[which(r_dsgn$status%in%c("Uninflamed")),]



##--------------------------------------##
##--##  Camembert of Response prop  ##--##
##--------------------------------------##

camembert_rep<-function(data_set,cat_prop,col_rep,col_rep_txt,plot_title){
  
  mltab<-melt(table(data_set$EnrichSupInf,data_set$response) )
  perc_df <-mltab[which(mltab$Var2==cat_prop),c(1,3)]%>% mutate(perc = `value` / sum(`value`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  # Reorder row order + levels to always have the same order
  perc_df <- perc_df[match(names(col_rep_txt), perc_df$Var1), ]
  perc_df$Var1 <- factor(perc_df$Var1, levels = levels(perc_df$Var1)[match(names(col_rep_txt), levels(perc_df$Var1))])
  return(ggplot(perc_df, aes(x="", y=value, fill=Var1))+
           coord_polar("y", start=0)+
           geom_bar(stat="identity", size=rel(0.8), color="black") + 
           scale_fill_manual(values=col_rep ) +  theme_void()+
           theme(axis.text.x=element_blank(),
                 plot.title=element_text(size=rel(1.4), face="bold"),
                 legend.title = element_text(size=rel(1.4)),
                 legend.text = element_text(size=rel(1.2))) +
           geom_text(aes(label = labels),color=col_rep_txt,size=rel(4),
                     position = position_stack(vjust = 0.5))+guides(fill = guide_legend(title = "Status")) +
           ggtitle(plot_title)
  )
}

camembert_rep_dp<-function(data_set,cat_prop,col_rep,col_rep_txt,plot_title){
  
  mltab<-melt(table(data_set$EnrichSupInf,data_set$response) )
  perc_df <-mltab[which(mltab$Var1==cat_prop),c(2,3)]%>% mutate(perc = `value` / sum(`value`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  # Reorder row order + levels to always have the same order
  perc_df <- perc_df[match(names(col_rep_txt), perc_df$Var2), ]
  perc_df$Var2 <- factor(perc_df$Var2, levels = levels(perc_df$Var2)[match(names(col_rep_txt), levels(perc_df$Var2))])
  return(ggplot(perc_df, aes(x="", y=value, fill=Var2))+
           coord_polar("y", start=0)+
           geom_bar(stat="identity", size=rel(0.8), color="black") + 
           scale_fill_manual(values=col_rep ) +  theme_void()+
           theme(axis.text.x=element_blank(),
                 plot.title=element_text(size=rel(1.4), face="bold"),
                 legend.title = element_text(size=rel(1.4)),
                 legend.text = element_text(size=rel(1.2))) +
           geom_text(aes(label = paste0(labels,"\nn=",value)),color=col_rep_txt,size=rel(4),
                     position = position_stack(vjust = 0.5))+guides(fill = guide_legend(title = "Status")) +
           ggtitle(plot_title)
  )
}


##--------------------------------------##
##--##   Table of response counts   ##--##
##--------------------------------------##

count_response<-function(data_set,fnt_size=24){
  
  g2 <- tableGrob(table(data_set$EnrichSupInf,data_set$response),theme=ttheme_minimal(base_size = fnt_size) )
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


r_dsgn_r<-r_dsgn[r_dsgn$response%in%c("NR","R"),]
r_dsgn_r_Inf<-r_dsgn_r[which(r_dsgn_r$status%in%c("Inflamed")),]
r_dsgn_r_Uninf<-r_dsgn_r[which(r_dsgn_r$status%in%c("Uninflamed")),]
signatures<-c("Siggrp_VIII","Siggrp_XV")


all_sig_comp_CAMrepCT<-function(sub_set){
  n <- length(signatures)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      x <- signatures[i]
      y <- signatures[j]
      lgd1=paste0(gsub("Siggrp_","",x))
      lgd2=paste0(gsub("Siggrp_","",y))
      print(paste0(lgd1," ",lgd2))
      
      # Distrib DP vs Others in each Rec rutgeerts score (i2a, i2b, i3, i4)
      eval(parse(text=paste0("data_set=r_dsgn_r",sub_set)))
      data_set$EnrichSupInf <- with(data_set, ifelse(eval(parse(text=paste0("dataCD68n_",x)))>= 0
                                                     & eval(parse(text=paste0("dataCD68n_",y)))>= 0 , "Double Pos","Others"))
      
      data_set$EnrichSupInf <- with(data_set, ifelse(eval(parse(text=paste0("dataCD68n_",x)))>= max(eval(parse(text=paste0("r_dsgn_r_Uninf$dataCD68n_",x))))
                                                     & eval(parse(text=paste0("dataCD68n_",y)))>= max(eval(parse(text=paste0("r_dsgn_r_Uninf$dataCD68n_",y)))) , "Double Pos","Others"))
      # To avoid cases where one condition is 0% of data, force the labels so melt() will create a 0 for them instead of deleting it
      data_set$EnrichSupInf<-factor(x=data_set$EnrichSupInf, levels = c("Others","Double Pos"))
      levels(data_set$EnrichSupInf)<-c("Others","Double Pos")
      
      cat_prop="R"
      col_rep<-c("Others"="#8FD0EC","Double Pos"="#20A1DA")
      col_rep_txt<-c("Others"="black","Double Pos"="black")
      r_dist<-camembert_rep(data_set=data_set,cat_prop=cat_prop,col_rep=col_rep,col_rep_txt=col_rep_txt,
                            plot_title=paste0(lgd1," & ",lgd2,"\ndistribution in ",cat_prop) )
      cat_prop="NR"
      col_rep<-c("Others"="#DD938F","Double Pos"="#BB2820")
      col_rep_txt<-c("Others"="black","Double Pos"="white")
      nr_dist<-camembert_rep(data_set=data_set,cat_prop=cat_prop,col_rep=col_rep,col_rep_txt=col_rep_txt,
                             plot_title=paste0(lgd1," & ",lgd2,"\ndistribution in ",cat_prop) )
      cnt_r_nr<-count_response(data_set=data_set,fnt_size=rel(12))
      
      
      cat_prop="Double Pos"
      col_rep<-c("R"="#20A1DA","NR"="#BB2820")
      col_rep_txt<-c("R"="black","NR"="white")
      dp_g_dist<-camembert_rep_dp(data_set=data_set,cat_prop=cat_prop,col_rep=col_rep,col_rep_txt=col_rep_txt,
                                  plot_title=paste0(cat_prop," distribution in\n",lgd1," & ",lgd2) )
      cat_prop="Others"
      col_rep<-c("R"="#8FD0EC","NR"="#DD938F")
      col_rep_txt<-c("R"="black","NR"="black")
      ot_g_dist<-camembert_rep_dp(data_set=data_set,cat_prop=cat_prop,col_rep=col_rep,col_rep_txt=col_rep_txt,
                                  plot_title=paste0(cat_prop," distribution in\n",lgd1," & ",lgd2) )
      cnt_mi_sev<-count_response(data_set=data_set,fnt_size=rel(12))
      
      
      design_perso="
      AABBDD
      AABBDD
      #CC#EE
      #CC#EE
      ####FF
      ####FF
      "
      
      print(r_dist+nr_dist+cnt_r_nr+dp_g_dist+ot_g_dist+cnt_mi_sev+plot_layout(design = design_perso) )
      
    }}}

pdf(paste0("./Figure 4S/FigS4Hb_DP_CamResp_tmCT.pdf"),width =9,height = 7)
all_sig_comp_CAMrepCT("_Inf")
dev.off()

