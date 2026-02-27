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
library(DESeq2)
library(xlsx)
library(rstatix)



load("./Grouped_objects/DGEseq_43_44_45_46/data_counts_mac43to46.rd") 
load("./Grouped_objects/DGEseq_43_44_45_46/data_logcpm_mac43to46.rd")
load("./Grouped_objects/Modsig_MoMacsEnr.rd")
colorpalette<-c("1"="#DDBBBB","2"="#CE9988","3"="#D07766","4"="#C15544","5"="#BB4030")
`%ni%` <- Negate(`%in%`)
obj<-readRDS("./Grouped_objects/DGEseq_43_44_45_46/MM250_DGEseq_DDS.rds")
meta_43to46<-data.frame(Wells=obj$wells,Barcodes=obj$barcodes,Samples=obj$samples,RIN=obj$RIN,UMI=obj$UMI)

allsum<-Rfast::colsums(l.data.cpm)  
to_remv<-c("MAC45-3","MAC44-5")  # too few reads to have trustable outputs  

# Compute SEM values for the error bars
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


# pdf(file="./Figures_print/FigS9C_Bulk43to46_Sig.pdf",width = 5.5, height = 10)
# for(i in names(Modsig)[1:length(names(Modsig))-1]){
# 
#     print(i)
# 
#     eval(parse(text=paste0("gnshow_y<-unlist(strsplit(Modsig$",i,",split = ','))")))
#     notfound_y<-setdiff(gnshow_y,rownames(l.data.cpm))
#     gnshow_y2<-setdiff(gnshow_y,notfound_y)
#     y_name=i
# 
#     currentscore<-(Rfast::colsums(l.data.cpm[gnshow_y2,])/allsum )
#     names(currentscore)<-colnames(l.data.cpm)
#     currentscore<-currentscore[setdiff(names(currentscore),to_remv)]
#     crt_sc<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))
#     names(crt_sc)<-names(currentscore)
#     
#     meta_43to46$i_y2<-crt_sc[match(meta_43to46$Samples,names(crt_sc))]
#     meta_43to46$condition<-as.character(gsub("^.*\\-","",meta_43to46$Samples))
#     meta_43to46$batch<-as.character(gsub("^(.*)\\-.*","\\1",meta_43to46$Samples))
#     pairtest_meta_43to46<-meta_43to46[which(!is.na(meta_43to46$i_y2)),]
#     pairtest_meta_43to46<-pairtest_meta_43to46[which(pairtest_meta_43to46$Samples%ni%c("MAC44-1","MAC44-5")),]  # MAC44-5 is missing so for paired wilcoxon test between condition 1 and 5, I need to remove MAC44-1
#     
#     dfAtt<-plyr::ddply(meta_43to46,~condition,function(x) mean_sem(x$i_y2,na.rm=T))
#     colnames(dfAtt)[2]<-"i_y2"  # so I can merge stats on meta_43to46 with dfAtt
#     
#     brbx_2<-ggplot(dfAtt, aes(x = condition, y = i_y2, group=i_y2)) + 
#       geom_line(aes(group=1),size=rel(2) ) +  #only 1 group to show so need aes(group=1) 
#       geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,size=rel(0.9)) +
#       geom_point(size=rel(7.2))+theme_minimal()+
#       xlab(label = "Upadacitinib (?M)")+coord_cartesian(ylim = c(minval, (maxval+ 0.5)))+ 
#       theme(axis.line = element_line(size = rel(3.2)),
#             axis.text = element_text(face="bold",size=rel(2.5), color="black"),
#             axis.ticks = element_line(size = rel(3.5)),
#             axis.ticks.length =unit(16,"pt"),
#             axis.title = element_text(face="bold",size=rel(3.2), color="black"),
#             axis.text.x = element_text(margin = margin(t = 8),angle = 90,vjust=0.5,hjust=1),
#             axis.title.x = element_text(margin = margin(t = 26)),panel.grid.minor = element_blank(),
#             panel.grid.major = element_blank(),axis.title.y = element_blank() )+
#       scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1"="0","2"="0.001","3"="0.01","4"="0.1","5"="1") )
#     
#     if (min(range(pairtest_meta_43to46$i_y2)) < 0 && max(range(pairtest_meta_43to46$i_y2)) > 0) {
#       # Testing unilateral 1 (0?M) > 5 (1?M)
#       brbx_3 <- brbx_2 + geom_hline(yintercept = 0, linetype="dashed", color = "#B3B3B3",size=rel(1.5))+
#         stat_compare_means(data=meta_43to46, comparisons = list(c("1","5")),size = 9,bracket.size = rel(1.5),label.y = 1.6,
#                            paired = FALSE,method = "wilcox.test",method.args = list(alternative = "greater"))+
#         ggtitle(paste0(i," Wilcox test, \nnot paired (only MAC44-5 removed), unilateral 1>5"))
#     print(brbx_3)
#     }
# }
# dev.off()


# Same but only for Fig 4F with Sig_XIV and Sig_XV
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


pdf(file="./Figures_print/Fig4E_Bulk43to46_Sig_grand.pdf",width = 4, height = 9)
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
    xlab(label = "Upadacitinib (?M)")+coord_cartesian(ylim = c(minval2, (maxval2)))+ 
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


###--- Take values from excel files ---####


# Table S18: Data points for Figure 4F & S10G
cxcl10_data<-read.xlsx("./Figure 4/Fig4E_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 1:6)
rownames(cxcl10_data)<-cxcl10_data$NA.
cxcl10_data$NA.<-NULL
colnames(cxcl10_data)<-1:5

il23_data<-read.xlsx("./Figure 4/Fig4E_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 8:13)
rownames(il23_data)<-il23_data$NA.
il23_data$NA.<-NULL
colnames(il23_data)<-1:5

ccl2_data<-read.xlsx("./Figure 4/Fig4E_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 15:20)
rownames(ccl2_data)<-ccl2_data$NA.
ccl2_data$NA.<-NULL
colnames(ccl2_data)<-1:5

mcsf_data<-read.xlsx("./Figure 4/Fig4E_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 22:27)
rownames(mcsf_data)<-mcsf_data$NA.
mcsf_data$NA.<-NULL
colnames(mcsf_data)<-1:5

il1b_data<-read.xlsx("./Figure 4/Fig4E_datapoints.xlsx","Feuil1", rowIndex = 2:7, colIndex = 29:34)
rownames(il1b_data)<-il1b_data$NA.
il1b_data$NA.<-NULL
colnames(il1b_data)<-1:5

# Comme CD150 je peux faire Nemenyi car appari?, mais jsp si on veut "comparer" au ex-vivo, lui est incomplet docn besoin de Wilcox 1>5 pour stat 

# data_set<-reshape2::melt(as.matrix(cxcl10_data))
# colnames(data_set)<-c("sample","variable","value")
# Fried_rstatix<-data_set %>%rstatix::friedman_test(value ~ variable |sample)
# pwc_label <- bquote(paste("pwc: ", bold("Friedman's test")) )  #pwc= PairWise Comparison
# test_label <- get_test_label(Fried_rstatix, detailed = TRUE)
# combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))
# Fried_rstatix
# library(PMCMRplus)
# # Pairwise comparisons using Nemenyi-Wilcoxon-Wilcox all-pairs test for a two-way balanced complete block design
# nemen_test<-PMCMRplus::frdAllPairsNemenyiTest(value ~ variable | sample, data=data_set)
# # encore plus stats que Wilcox 1>5 unpaired


pdf(file="./Figures_print/Fig4F_Bulk43to46_Gene_grand.pdf",width = 4, height = 9)
for(i in c("cxcl10","il23","ccl2","mcsf","il1b")){
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
                                     bracket.size = rel(1.5),paired = FALSE,method = "wilcox.test",method.args = list(alternative = "greater"))+
  ggtitle(paste0(i," Wilcox test, \nnot paired, unilateral 1>5\n"))+ expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )
print(brbx_3)
}
dev.off()





CD150_data<-read.xlsx("./Figure 4/Fig4E_datapoints.xlsx","Feuil1", rowIndex =11:15, colIndex = 8:13)
rownames(CD150_data)<-CD150_data$NA.
CD150_data$NA.<-NULL
colnames(CD150_data)<-1:5

data_set<-reshape2::melt(as.matrix(CD150_data))
colnames(data_set)<-c("sample","variable","value")
## Here all conditios use the same samples, so even if non parametric, it's paired, so use Friedman's not KW :
Fried_rstatix<-data_set %>%rstatix::friedman_test(value ~ variable |sample)
pwc_label <- bquote(paste("pwc: ", bold("Friedman's test")) )  #pwc= PairWise Comparison
test_label <- get_test_label(Fried_rstatix, detailed = TRUE)
combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))

Fried_rstatix
## Fried is stat, perform paired Wilcoxon signed-rank test ?  (Dunn is not for paired)
# It's an adapted test but we only have 4 samples so we will never get anything lower than p=0.125, it's useless to do a stat test with n=4, even in "one-sided" it's 0.0625
# Also recommended: Nemenyi and Conover

# wlcx_results <- data_set %>%
# wilcox_test(value ~ variable, paired = TRUE, p.adjust.method = "holm", alternative= "greater")
# wlcx_results$padj_round<-round(wlcx_results$p.adj,4)

## Testing Nemenyi:
# https://sites.google.com/site/rgraphiques/4--stat/comparaison-de-moyennes-avec-r/les-tests-post-hoc-sous-r
library(PMCMRplus)
# Pairwise comparisons using Nemenyi-Wilcoxon-Wilcox all-pairs test for a two-way balanced complete block design
nemen_test<-PMCMRplus::frdAllPairsNemenyiTest(value ~ variable | sample, data=data_set)

# Stat of condition 1 (0uM) vs 5 (1uM) --->>>
stat_0v1<-nemen_test$p.value["5", "1"]

#make a bar :
y_max <- max(data_set$value, na.rm = TRUE)
y_bar <- y_max * 1.1
y_text <- y_max * 1.2


i="CD150"
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
brbx_3 <- brbx_2 +
  scale_y_continuous(expand = c(0, 0),n.breaks = 5)+
  coord_cartesian(ylim = c(0, 40000),clip = 'off')+
  geom_segment( aes(x = 1, xend = 5, y = y_bar, yend = y_bar),
    inherit.aes = FALSE ,size=1.5 ) +
  geom_segment( aes(x = 1, xend = 1, y = y_bar, yend = y_bar * 0.97),
    inherit.aes = FALSE ,size=1.5  ) +
  geom_segment( aes(x = 5, xend = 5, y = y_bar, yend = y_bar * 0.97),
    inherit.aes = FALSE ,size=1.5  ) +
  annotate("text",x = 3,y = y_text,size=8.5 ,
    label = paste0("p = ", signif(stat_0v1, 3))  )+
  ggtitle(paste0(i," Friedman -> Nemenyi post-hoc test, paired data"))+ expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )



pdf(file="./Figures_print/Fig4E_CD150MFI_nemenyi.pdf",width = 4, height = 9)
print(brbx_3)
dev.off()


pdf(file="./Figures_print/Fig4E_CD150MFI.pdf",width = 4, height = 9)
i="CD150"
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
                                       bracket.size = rel(1.5),paired = FALSE,method = "wilcox.test",method.args = list(alternative = "greater"))+
    ggtitle(paste0(i," Wilcox test, \nnot paired, unilateral 1>5\n"))+ expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )+ 
    expand_limits(y = max(dfAtt$ymax)+0.3*(max(dfAtt$ymax)-min(dfAtt$ymin)) )+
    scale_y_continuous(expand = c(0, 0),n.breaks = 5)+
    coord_cartesian(ylim = c(0, 40000),clip = 'off')
  print(brbx_3)

dev.off()

