library(xlsx)
library(ggplot2)
library(ggpubr)
library(Rfast)
library(edgeR)
library(ggrepel)
library(factoextra)
library(DESeq2)


`%ni%` <- Negate(`%in%`)
# Compute SEM values for the error bars
mean_sem<-function(x,digits= 3,na.rm=FALSE){
  if(na.rm==TRUE) {x<-x[!is.na(x)]}
  y <- mean(x)
  sem_val<-round(sd(x)/sqrt(length(x)),digits)
  ymin<-y-sem_val
  ymax<-y+sem_val
  return(c(y = y, ymin = ymin, ymax = ymax))
}


# Table S18: Data points for Figure 4E & S10G 
ccl7_data<-read.xlsx("./Figure 10S/FigS10G_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 1:6)
rownames(ccl7_data)<-ccl7_data$NA.
ccl7_data$NA.<-NULL
colnames(ccl7_data)<-1:5

opn_data<-read.xlsx("./Figure 10S/FigS10G_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 8:13)
rownames(opn_data)<-opn_data$NA.
opn_data$NA.<-NULL
colnames(opn_data)<-1:5

il1a_data<-read.xlsx("./Figure 10S/FigS10G_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 15:20)
rownames(il1a_data)<-il1a_data$NA.
il1a_data$NA.<-NULL
colnames(il1a_data)<-1:5

il6_data<-read.xlsx("./Figure 10S/FigS10G_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 22:27)
rownames(il6_data)<-il6_data$NA.
il6_data$NA.<-NULL
colnames(il6_data)<-1:5



pdf(file="./Figures_print/FigS10G_invitro_gene.pdf",width = 6, height = 9)
for(i in c("ccl7","opn","il1a","il6")){
  eval(parse(text=paste0("gene_data <- ",i,"_data")))
  
  dfAtt<-plyr::ddply(reshape2::melt(gene_data),~variable,function(x) mean_sem(x$value,na.rm=T))
  colnames(dfAtt)[2]<-"i_y2"  # so I can merge stats on meta_43to46 with dfAtt
  hline_layer<- if (min(range(dfAtt$ymin)) < 0 && max(range(dfAtt$ymax)) > 0) {
    geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))}else{NULL}
  brbx_2<-ggplot(dfAtt, aes(x = variable, y = i_y2, group=i_y2)) + 
    hline_layer+
    geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
    geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
    geom_point(size=rel(7.2))+theme_minimal()+
    xlab(label = "Upadacitinib (?M)")+ #coord_cartesian(ylim = c(minval2, (maxval2)))+ 
    theme(axis.line = element_line(size = rel(3.2)),
          axis.text = element_text(face="bold",size=rel(2.5), color="black"),
          axis.ticks = element_line(size = rel(3.5)),
          axis.ticks.length =unit(16,"pt"),
          axis.title = element_text(face="bold",size=rel(3.2), color="black"),
          axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
          axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),axis.title.y = element_blank() )+
    scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1") )
  
  m_g_d<-reshape2::melt(gene_data)
  colnames(m_g_d)[2]<-"i_y2"
  brbx_3 <- brbx_2 +stat_compare_means(data=m_g_d, comparisons = list(c("1","5")),size = 9, tip.length = 0,
                                       label.y= max(dfAtt$ymax)+0.1*(max(dfAtt$ymax)-min(dfAtt$ymin)),
                                       bracket.size = rel(1.5),paired = FALSE,method = "wilcox.test")+  #method.args = list(alternative = "greater")
    ggtitle(paste0(i," Wilcox test, \nnot paired, bilateral 1!=5\n"))+ expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )
  print(brbx_3)
}
dev.off()




# Table S19: Data points for Figure 4F & S10G 
ccl7_data<-read.xlsx("./Figure 10S/FigS10G_exvivo_datapoints.xlsx","Feuil1", rowIndex = 2:6, colIndex = 1:7)
rownames(ccl7_data)<-ccl7_data$NA.
ccl7_data$NA.<-NULL
colnames(ccl7_data)<-1:6

opn_data<-read.xlsx("./Figure 10S/FigS10G_exvivo_datapoints.xlsx","Feuil1", rowIndex = 2:6, colIndex = 9:15)
rownames(opn_data)<-opn_data$NA.
opn_data$NA.<-NULL
colnames(opn_data)<-1:6

il1a_data<-read.xlsx("./Figure 10S/FigS10G_exvivo_datapoints.xlsx","Feuil1", rowIndex = 2:6, colIndex = 17:23)
rownames(il1a_data)<-il1a_data$NA.
il1a_data$NA.<-NULL
colnames(il1a_data)<-1:6

il6_data<-read.xlsx("./Figure 10S/FigS10G_exvivo_datapoints.xlsx","Feuil1", rowIndex = 2:6, colIndex = 25:31)
rownames(il6_data)<-il6_data$NA.
il6_data$NA.<-NULL
colnames(il6_data)<-1:6



pdf(file="./Figures_print/FigS10G_exvivo_gene.pdf",width = 6, height = 9)
for(i in c("ccl7","opn","il1a","il6")){
  eval(parse(text=paste0("gene_data <- ",i,"_data")))
  
  dfAtt<-plyr::ddply(reshape2::melt(gene_data),~variable,function(x) mean_sem(x$value,na.rm=T))
  colnames(dfAtt)[2]<-"i_y2"  # so I can merge stats on meta_43to46 with dfAtt
  hline_layer<- if (min(range(dfAtt$ymin)) < 0 && max(range(dfAtt$ymax)) > 0) {
    geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))}else{NULL}
  brbx_2<-ggplot(dfAtt, aes(x = variable, y = i_y2, group=i_y2)) + 
    hline_layer+
    geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
    geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
    geom_point(size=rel(7.2))+theme_minimal()+
    xlab(label = "Upadacitinib (?M)")+ #coord_cartesian(ylim = c(minval2, (maxval2)))+ 
    theme(axis.line = element_line(size = rel(3.2)),
          axis.text = element_text(face="bold",size=rel(2.5), color="black"),
          axis.ticks = element_line(size = rel(3.5)),
          axis.ticks.length =unit(16,"pt"),
          axis.title = element_text(face="bold",size=rel(3.2), color="black"),
          axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
          axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),axis.title.y = element_blank() )+
    scale_x_discrete(limits=c("1","2","3","4","5","6"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1","6"="10") )
  
  m_g_d<-reshape2::melt(gene_data)
  colnames(m_g_d)[2]<-"i_y2"
  brbx_3 <- brbx_2 +stat_compare_means(data=m_g_d, comparisons = list(c("1","5")),size = 9, tip.length = 0,
                                       label.y= max(dfAtt$ymax)+0.1*(max(dfAtt$ymax)-min(dfAtt$ymin)),
                                       bracket.size = rel(1.5),paired = FALSE,method = "wilcox.test")+  #method.args = list(alternative = "greater")
    ggtitle(paste0(i," Wilcox test, \nnot paired, bilateral 1!=5\n"))+ expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )
  print(brbx_3)
}
dev.off()
