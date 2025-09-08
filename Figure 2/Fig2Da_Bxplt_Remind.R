library(Rfast)
library(ggplot2)
library(ggprism)
library(ggpubr)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
# Load the metadata and count matrix from Ngollo cohort
load("Grouped_objects/minipheno_AllezNgollo_postop_goodRemRec_sig.rd")
load("Grouped_objects/genelist_16nov22_Macs.rd")

miniphenoM0I<-minipheno[which(minipheno$`location:ch1`%in%c("M0I")),]
signatures <- c("data_Siggrp_XV")
levels(signatures)<-signatures



##---------------------------------------##
## Boxplots w/ Rutgeerts scores (CTM0I)  ##
##---------------------------------------##

bx_rutgeerts<-function(data_set,val2,col_rut,lgd2,plot_title,comp_rut,stat){
  plot<-ggplot(data_set, aes(x =new_rut, y =!!sym(val2), color =new_rut, fill=new_rut)) +
    geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
    scale_x_discrete(limits = names(col_rut) )+
    geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
    ggprism::theme_prism(base_size = 25)+
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30))+
    labs(y= lgd2) +
    scale_color_manual(values=col_rut, name="Rutgeerts\nendocopic\nrecurrence\nscore")+
    scale_fill_manual(values=col_rut, name="Rutgeerts\nendocopic\nrecurrence\nscore")+
    ggtitle(plot_title)
  if(stat==T){
    #It does a wilcoxon test to compare means between 2 groups
    return(plot+stat_compare_means(comparisons = comp_rut,size = 6,
                                   bracket.size = rel(2)) )
  }else{return(plot)}
}  


##---------------------------------------##


col_rut<-c("i2"="#ba7000",
           "i3_4"="#801c01")
comp_rut=list(c("i2","i3_4"))


all_sig_comp_BR<-function(sub_set){
  eval(parse(text=paste0("data_set=minipheno",sub_set)))
  data_set$new_rut<-data_set$`rutgeerts:ch1`
  data_set<-data_set[-which(data_set$new_rut%in%c("0","1")),]
  data_set$new_rut[which(data_set$new_rut%in%c("i2a","i2b"))]<-"i2"
  data_set$new_rut[which(data_set$new_rut%in%c("i3","i4"))]<-"i3_4"
  n <- length(signatures)
  for (j in 1:n) {
    y <- signatures[j]
    # lgd2=paste0("Signature ",gsub("data_Siggrp_","",y))
    lgd2="XV_Infl._III"
    print(paste0(lgd2))
    print(bx_rutgeerts(data_set =data_set ,val2=y,lgd2=lgd2,plot_title=sub_set,
                       col_rut=col_rut,comp_rut=comp_rut,stat=T))
  }
}



pdf("./Figure 2/Fig2Da_XV_Bxplt.pdf",width =6,height = 7.5)
all_sig_comp_BR("M0I")
dev.off()



