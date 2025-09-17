library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
load("./Grouped_objects/r_dsgn_InfUninf_cd68norm_Momacsscores.rd")
r_dsgn_Inf<-r_dsgn[which(r_dsgn$status%in%c("Inflamed")),]
r_dsgn_Uninf<-r_dsgn[which(r_dsgn$status%in%c("Uninflamed")),]



r_dsgn_r<-r_dsgn[r_dsgn$response%in%c("NR","R"),]
r_dsgn_r_Inf<-r_dsgn_r[which(r_dsgn_r$status%in%c("Inflamed")),]
r_dsgn_r_Uninf<-r_dsgn_r[which(r_dsgn_r$status%in%c("Uninflamed")),]

##------------------------------##
## Boxplots w/ Signature scores ##
##------------------------------##

bx_rep<-function(data_set,sub_set,val2,col_rep,lgd2,plot_title,comp_rep,stat){
  eval(parse(text=paste0("min_y_value<-min(data_set$",val2,")" )))
  eval(parse(text=paste0("max_y_value<-max(data_set$",val2,")" )))
  #val_2<-paste0("data",sub_set,"_",val2)
  val_2<-val2
  plot<-ggplot(data_set, aes(x =response, y =!!sym(val_2), color =response, fill=response)) +
    geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
    scale_x_discrete(limits = names(col_rep) )+
    geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
    ggprism::theme_prism(base_size = 20)+
    labs(y= sig_rec[lgd2]) +
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.line = element_line(size = rel(1.6)),
          axis.text = element_text(face="bold",size=rel(1.5), color="black"),
          axis.ticks = element_line(size = rel(1.6)),
          axis.ticks.length =unit(16,"pt"),
          axis.title = element_text(face="bold",size=rel(1.8), color="black"),
          axis.text.x = element_text(margin = margin(t = 8)),
          axis.text.y = element_text(margin = margin(r = 8)),
          panel.grid.minor = element_blank(),panel.grid.major = element_blank() )+
      scale_color_manual(values=col_rep, name="Response")+
      scale_fill_manual(values=col_rep, name="Response")+
    ggtitle(plot_title)+
    coord_cartesian(ylim = c(min_y_value, max_y_value+(0.1*(max_y_value-min_y_value) )))#Get 10% more space above the plot so the stat text is not cut
  
  if(stat==T){
    #It does a wilcoxon test to compare means between 2 groups
    return(plot+stat_compare_means(comparisons = comp_rep,size = rel(10),vjust=-0.2,
                                   bracket.size = rel(2)) )
  }else{return(plot)}
}


##---------------------------------------##

#### D6) Signature Boxplot R vs NR ####

col_rep<-c("R"="#20A1DA","NR"="#BB2820" )
comp_rep=list(c("R","NR"))

signatures <- c("dataCD68n_Siggrp_I","dataCD68n_Siggrp_II","dataCD68n_Siggrp_III","dataCD68n_Siggrp_IV","dataCD68n_Siggrp_V",
                "dataCD68n_Siggrp_VI","dataCD68n_Siggrp_VII","dataCD68n_Siggrp_VIII","dataCD68n_Siggrp_IX","dataCD68n_Siggrp_X",   
                "dataCD68n_Siggrp_XI","dataCD68n_Siggrp_XII","dataCD68n_Siggrp_XIII","dataCD68n_Siggrp_XIV","dataCD68n_Siggrp_XV",  
                "dataCD68n_Siggrp_XVI","dataCD68n_Siggrp_XVII")
levels(signatures)<-signatures

all_sig_comp_BR<-function(sub_set){
  data_set=r_dsgn_r
  #signatures<-colnames(r_dsgn_r)[grep(paste0("dataCD68n_Siggrp"),colnames(r_dsgn_r))]
  n <- length(signatures)
  for (j in 1:n) {
    y <- signatures[j]
    lgd2=paste0("Signature ",gsub(paste0("dataCD68n_Siggrp_"),"",y))
    print(lgd2)
    print(bx_rep(data_set =data_set,sub_set=sub_set,val2=y,lgd2=lgd2,plot_title=sub_set,
                 col_rep=col_rep,comp_rep=comp_rep,stat=T))
  }
}

sig_rec<-c("Signature I"="I_Transition A","Signature II"="II_OXPHOS","Signature III"="Program III","Signature IV"="IV_Transition B","Signature V"="V_Macrophage","Signature VI"="VI_Transition C",
           "Signature VII"="VII_Mono A","Signature VIII"="VIII_Mono B","Signature IX"="IX_Interferon","Signature X"="Program X","Signature XI"="Program XI","Signature XII"="XII_Immunoreg/repair",
           "Signature XIII"="XIII_Inflammatory A","Signature XIV"="XIV_Inflammatory B","Signature XV"="XV_Inflammatory Mono1","Signature XVI"="XVI_Stress","Signature XVII"="XVII_Hypoxia")

pdf(paste0("./Figure 4S/FigS4Ib_Bxplt_Risk.pdf"),width =6, height = 7.5)
all_sig_comp_BR("_Inf")
dev.off()


