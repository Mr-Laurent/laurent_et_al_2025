# Figures Upa cytok.
library(tidyr)
library(tidyverse)
library(ggh4x)
library(ggplot2)
library(ggpubr)
library(Rfast)
library(edgeR)
library(ggrepel)
library(factoextra)
library(viridis)


load("./Grouped_objects/DGEseq_43_44_45_46/data_logcpm_mac43to46.rd")
load("./Grouped_objects/DGEseq_43_44_45_46/metadata_mac43to46.rd") 

genes_lst<-"CXCL10,SEMA4A,LTB,EBI3,CCL19,CD274,INHBA,SPP1,CCL2,CXCL5,CLU,CXCL9,CXCL1,CSF1,CXCL11,IL6,TNFSF14,CSF3,IL1A,TNFSF15,CSF2,IL1RN,CCL5,CCL22,CXCL8,TNF,CLCF1,CCL3,CCL4,CCL20,IL23A,IL1B,CCL4L2,CCL7,EREG,ADM,VEGFA,PDGFB,IL10,TNFSF12,HBEGF,APLP2,THBS1,SEMA4D,CXCL14,IL15,OSM,TNFSF13,TNFSF13B,AREG,CXCL16,PLAU,TGFB1,TNFSF10,NRG1,RNASET2,LRPAP1,IL16,IL18,PDGFC,VEGFB,CXCL12,IGF1,CCL18,APOE"
genes_lst<-c(genes_lst,"IL12B","IL6")

tl_data<-as.data.frame(as.matrix(t(l.data.cpm)))

m_43to46<-meta_43to46[,c(3,6,7)]
rownames(m_43to46)<-m_43to46$Samples
m_43to46<-m_43to46[,-1]
merged_df <- merge(m_43to46, tl_data, by = "row.names")   # don't keep LPS/TNF/GM detail as it might interfer with genes
rownames(merged_df) <- merged_df$Row.names
merged_df <- merged_df[, -1] 

## MAC44-5 and MAC45-3 low counts : exclude
ord_x1<-c("MAC43-1","MAC44-1","MAC45-1","MAC46-1","MAC43-2","MAC44-2","MAC45-2","MAC46-2",
          "MAC43-3","MAC44-3","MAC46-3","MAC43-4","MAC44-4","MAC45-4","MAC46-4","MAC43-5","MAC45-5","MAC46-5")
to_remv<-c("MAC45-3","MAC44-5")

ord_x3=c("1","2","3","4","5")

Upts2=1.4
y_txt_size=3


`%ni%` <- Negate(`%in%`)
unlist(strsplit(genes_lst,split=","))[which(unlist(strsplit(genes_lst,split=","))%ni%colnames(tl_data))]
# CXCL14 not found
genes_list<-unlist(strsplit(genes_lst,split=","))[which(unlist(strsplit(genes_lst,split=","))%in%colnames(tl_data))]


mean_sem<-function(x,digits= 3,na.rm=FALSE){
  if(na.rm==TRUE) {x<-x[!is.na(x)]}
  y <- mean(x)
  sem_val<-round(sd(x)/sqrt(length(x)),digits)
  ymin<-y-sem_val
  ymax<-y+sem_val
  return(c(y = y, ymin = ymin, ymax = ymax))
}

genes_list_sel=c("CXCL9","CXCL10","CXCL11","CSF1","IL23A","CCL7","TNF","SPP1","INHBA","IL1A","IL1B","IL6")
pdf(file="./Figures_print/FigS10F_selectCytok_Bulk43to46.pdf",width = 7, height = 9)

for(y_name in genes_list_sel){
  eval(parse(text=paste0("merged_df$i_y2<-merged_df[,y_name]")))  
  
  dfAtt<-plyr::ddply(merged_df,~condition,function(x) mean_sem(x$i_y2,na.rm=T))
  colnames(dfAtt)[2]<-"i_y2"
  
  brbx_2<-ggplot(dfAtt, aes(x = condition, y = i_y2, group=i_y2)) + 
    geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
    geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
    geom_point(size=rel(7.2))+theme_minimal()+ggtitle(label=y_name,)+
    xlab(label = "Upadacitinib (µM)")+
    theme(axis.line = element_line(size = rel(3.2)),
          title = element_text(size=rel(3)),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(face="bold",size=rel(2.5), color="black"),
          axis.ticks = element_line(size = rel(3.5)),
          axis.ticks.length =unit(16,"pt"),
          axis.title = element_text(face="bold",size=rel(0.7), color="black"),
          axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
          axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),axis.title.y = element_blank() )+
    scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1") )
  
  print(brbx_2)
}
dev.off()
  
# pdf(file="G./Figure 10S/FigS10F_allCytok_Bulk43to46.pdf",width = 7, height = 8)
# 
# for(y_name in genes_list){
#   eval(parse(text=paste0("merged_df$i_y2<-merged_df[,y_name]")))  
#   
#   dfAtt<-plyr::ddply(merged_df,~condition,function(x) mean_sem(x$i_y2,na.rm=T))
#   colnames(dfAtt)[2]<-"i_y2"
#   
#   brbx_2<-ggplot(dfAtt, aes(x = condition, y = i_y2, group=i_y2)) + 
#     geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
#     geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
#     geom_point(size=rel(7.2))+theme_minimal()+ggtitle(label=y_name)+
#     xlab(label = "Upadacitinib (µM)")+
#     theme(axis.line = element_line(size = rel(3.2)),
#           title = element_text(size=rel(2)),
#           axis.text = element_text(face="bold",size=rel(2.5), color="black"),
#           axis.ticks = element_line(size = rel(3.5)),
#           axis.ticks.length =unit(16,"pt"),
#           axis.title = element_text(face="bold",size=rel(1.6), color="black"),
#           axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
#           axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),axis.title.y = element_blank() )+
#     scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1") )
#   
#   print(brbx_2)
#   
# }
# dev.off()
