#From Analysis_DGE43to46.R

# EFS samples 
# All 4 samples are cultured in GMCSF + LPS + TNF + IFNg
# 
# Conditions:
# 1: No Upadacitinib
# 2: 0.001 mg/mL  Upa
# 3: 0.01 mg/mL
# 4: 0.1 mg/mL
# 5: 1 mg/mL
# 


library(ggplot2)
library(ggpubr)
library(Rfast)
library(edgeR)
library(ggrepel)
library(factoextra)



load("./Grouped_objects/DGEseq_43_44_45_46/data_counts_mac43to46.rd") 
load("./Grouped_objects/DGEseq_43_44_45_46/data_logcpm_mac43to46.rd")
load("./Grouped_objects/Modsig_MoMacsEnr.rd")
colorpalette<-c("1"="#DDBBBB","2"="#CE9988","3"="#D07766","4"="#C15544","5"="#BB4030")
`%ni%` <- Negate(`%in%`)
obj<-readRDS("./Grouped_objects/DGEseq_43_44_45_46/MM250_DGEseq_DDS.rds")
meta_43to46<-data.frame(Wells=obj$wells,Barcodes=obj$barcodes,Samples=obj$samples,RIN=obj$RIN,UMI=obj$UMI)


allsum<-Rfast::colsums(l.data.cpm)  
to_remv<-c("MAC45-3","MAC44-5")  # too few reads to have trustable outputs  

#Measured values, for stats
wilcox.test(c(1.22191103,0.20585640,0.69181546,0.37733470), c(-1.61697427,0.13499555,-0.43523240,-1.20075527), paired = TRUE, alternative = "greater") # we test if value in 1 are greater than 5
wilcox.test(c(1.22191103,0.69181546,0.37733470), c(-1.61697427,-0.43523240,-1.20075527), paired = TRUE, alternative = "greater") # we test if value in 1 are greater than 5
# unilat: p-value = 0.0625    Si que 3 values (pas MAC44-5): 0.125
# bilat:  p-value = 0.125     Si que 3 values              : 0.25
wilcox.test(c(1.22191103,0.69181546,0.37733470), c(-1.61697427,-0.43523240,-1.20075527), paired = F, alternative = "greater")
# unilat: p-value = 0.0625    Si que 3 values (pas MAC44-5): 0.05
wilcox.test(c(1.22191103,0.20585640,0.69181546,0.37733470), c(-1.61697427,NA,-0.43523240,-1.20075527), paired = F, alternative = "greater")
# unilat: p-value = 0.02857
# bilat : p-value = 0.05714

mean_sem<-function(x,digits= 3,na.rm=FALSE){
  if(na.rm==TRUE) {x<-x[!is.na(x)]}
  y <- mean(x)
  sem_val<-round(sd(x)/sqrt(length(x)),digits)
  ymin<-y-sem_val
  ymax<-y+sem_val
  return(c(y = y, ymin = ymin, ymax = ymax))
}






#Get min and max score values across all score to homogeneize the plots:
minval=0
maxval=0
for(mods in names(Modsig)[-length(names(Modsig))]){
  eval(parse(text=paste0("gnshow_y<-unlist(strsplit(Modsig$",mods,",split = ','))")))
  notfound_y<-setdiff(gnshow_y,rownames(l.data.cpm))
  gnshow_y2<-setdiff(gnshow_y,notfound_y)
  
  currentscore<-(Rfast::colsums(l.data.cpm[gnshow_y2,])/allsum )
  names(currentscore)<-colnames(l.data.cpm)
  currentscore<-currentscore[setdiff(names(currentscore),to_remv)]
  crt_sc<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))
  names(crt_sc)<-names(currentscore)
  if(min(crt_sc)<minval){minval<-min(crt_sc)}
  if(max(crt_sc)>minval){maxval<-max(crt_sc)}
}  


pdf(file="./Figures_print/FigS10D_Bulk43to46_Sig.pdf",width = 5.5, height = 10)
for(i in names(Modsig)[1:length(names(Modsig))-1]){

    print(i)

    eval(parse(text=paste0("gnshow_y<-unlist(strsplit(Modsig$",i,",split = ','))")))
    notfound_y<-setdiff(gnshow_y,rownames(l.data.cpm))
    gnshow_y2<-setdiff(gnshow_y,notfound_y)
    y_name=i

    currentscore<-(Rfast::colsums(l.data.cpm[gnshow_y2,])/allsum )
    names(currentscore)<-colnames(l.data.cpm)
    currentscore<-currentscore[setdiff(names(currentscore),to_remv)]
    crt_sc<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))
    names(crt_sc)<-names(currentscore)
    
    meta_43to46$i_y2<-crt_sc[match(meta_43to46$Samples,names(crt_sc))]
    meta_43to46$condition<-as.character(gsub("^.*\\-","",meta_43to46$Samples))
    meta_43to46$batch<-as.character(gsub("^(.*)\\-.*","\\1",meta_43to46$Samples))
    pairtest_meta_43to46<-meta_43to46[which(!is.na(meta_43to46$i_y2)),]
    pairtest_meta_43to46<-pairtest_meta_43to46[which(pairtest_meta_43to46$Samples%ni%c("MAC44-1","MAC44-5")),]  # MAC44-5 is missing so for paired wilcoxon test between condition 1 and 5, I need to remove MAC44-1
    
    dfAtt<-plyr::ddply(meta_43to46,~condition,function(x) mean_sem(x$i_y2,na.rm=T))
    colnames(dfAtt)[2]<-"i_y2"  # so I can merge stats on meta_43to46 with dfAtt
    
    brbx_2<-ggplot(dfAtt, aes(x = condition, y = i_y2, group=i_y2)) + 
      geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
      geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
      geom_point(size=rel(7.2))+theme_minimal()+
      xlab(label = "Upadacitinib (µM)")+coord_cartesian(ylim = c(minval, (maxval+ 0.5)))+ 
      theme(axis.line = element_line(size = rel(3.2)),
            axis.text = element_text(face="bold",size=rel(2.5), color="black"),
            axis.ticks = element_line(size = rel(3.5)),
            axis.ticks.length =unit(16,"pt"),
            axis.title = element_text(face="bold",size=rel(3.2), color="black"),
            axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
            axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),axis.title.y = element_blank() )+
      scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1") )
    
    if (min(range(pairtest_meta_43to46$i_y2)) < 0 && max(range(pairtest_meta_43to46$i_y2)) > 0) {
      brbx_3 <- brbx_2 + geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))+
        stat_compare_means(data=meta_43to46, comparisons = list(c("1","5")),size = 9,bracket.size = rel(1.5),label.y = 1.6,
                           paired = FALSE,method = "wilcox.test",method.args = list(alternative = "greater"))+
        ggtitle(paste0(i," Wilcox test, \nnot paired (only MAC44-5 removed), unilateral 1>5"))
      brbx_4 <- brbx_2 + geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))+
        stat_compare_means(data=pairtest_meta_43to46, comparisons = list(c("1","5")),size = 9,bracket.size = rel(1.5),label.y = 1.6,paired = TRUE,method = "wilcox.test",method.args = list(alternative = "greater"))+
        ggtitle(paste0(i," Wilcox test, \nPaired (only MAC44-5 removed), unilateral 1>5")) 
      brbx_5 <- brbx_2 + geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))+
        stat_compare_means(data=meta_43to46, comparisons = list(c("1","5")),size = 9,bracket.size = rel(1.5),label.y = 1.6,paired = FALSE,method = "wilcox.test")+
        ggtitle(paste0(i," Wilcox test, \nnot paired (only MAC44-5 removed), bilateral 1!=5")) }
    print(brbx_5)
    print(brbx_3)
    print(brbx_4)
    }
dev.off()



minval2=0
maxval2=0
for(mods in names(Modsig)[c(14,15)]){
  eval(parse(text=paste0("gnshow_y<-unlist(strsplit(Modsig$",mods,",split = ','))")))
  notfound_y<-setdiff(gnshow_y,rownames(l.data.cpm))
  gnshow_y2<-setdiff(gnshow_y,notfound_y)
  
  currentscore<-(Rfast::colsums(l.data.cpm[gnshow_y2,])/allsum )
  names(currentscore)<-colnames(l.data.cpm)
  currentscore<-currentscore[setdiff(names(currentscore),to_remv)]
  crt_sc<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))
  names(crt_sc)<-names(currentscore)
  if(min(crt_sc)<minval2){minval2<-min(crt_sc)}
  if(max(crt_sc)>minval2){maxval2<-max(crt_sc)}
}  


pdf(file="./Figures_print/FigS10D_Bulk43to46_Sig_grand.pdf",width = 7, height = 9.5)
for(i in names(Modsig)[c(14,15)]){
  
  print(i)
  
  eval(parse(text=paste0("gnshow_y<-unlist(strsplit(Modsig$",i,",split = ','))")))
  notfound_y<-setdiff(gnshow_y,rownames(l.data.cpm))
  gnshow_y2<-setdiff(gnshow_y,notfound_y)
  y_name=i
  
  currentscore<-(Rfast::colsums(l.data.cpm[gnshow_y2,])/allsum )
  names(currentscore)<-colnames(l.data.cpm)
  currentscore<-currentscore[setdiff(names(currentscore),to_remv)]
  crt_sc<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))
  names(crt_sc)<-names(currentscore)
  
  meta_43to46$i_y2<-crt_sc[match(meta_43to46$Samples,names(crt_sc))]
  meta_43to46$condition<-as.character(gsub("^.*\\-","",meta_43to46$Samples))
  meta_43to46$batch<-as.character(gsub("^(.*)\\-.*","\\1",meta_43to46$Samples))
  pairtest_meta_43to46<-meta_43to46[which(!is.na(meta_43to46$i_y2)),]
  pairtest_meta_43to46<-pairtest_meta_43to46[which(pairtest_meta_43to46$Samples%ni%c("MAC44-1","MAC44-5")),]  # MAC44-5 is missing so for paired wilcoxon test between condition 1 and 5, I need to remove MAC44-1
  
  dfAtt<-plyr::ddply(meta_43to46,~condition,function(x) mean_sem(x$i_y2,na.rm=T))
  colnames(dfAtt)[2]<-"i_y2"  # so I can merge stats on meta_43to46 with dfAtt
    hline_layer<- if (min(range(pairtest_meta_43to46$i_y2)) < 0 && max(range(pairtest_meta_43to46$i_y2)) > 0) {
    geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))}else{NULL}
    
  brbx_2<-ggplot(dfAtt, aes(x = condition, y = i_y2, group=i_y2)) + 
    hline_layer+
    geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
    geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
    geom_point(size=rel(7.2))+theme_minimal()+
    xlab(label = "Upadacitinib (µM)")+coord_cartesian(ylim = c(minval2, (maxval2)))+ 
    theme(axis.line = element_line(size = rel(3.2)),
          axis.text = element_text(face="bold",size=rel(2.5), color="black"),
          axis.ticks = element_line(size = rel(3.5)),
          axis.ticks.length =unit(16,"pt"),
          axis.title = element_text(face="bold",size=rel(3.2), color="black"),
          axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
          axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),axis.title.y = element_blank() )+
    scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1") )
  
    brbx_3 <- brbx_2 +stat_compare_means(data=meta_43to46, comparisons = list(c("1","5")),size = 9, 
                                         label.y= max(dfAtt$ymax)+0.1*(max(dfAtt$ymax)-min(dfAtt$ymin)),
                         bracket.size = rel(1.5),paired = FALSE,method = "wilcox.test",method.args = list(alternative = "greater"))+
      ggtitle(paste0(i," Wilcox test, \nnot paired (only MAC44-5 removed), unilateral 1>5\n"))+ expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )
    print(brbx_3)
  
}
dev.off()

