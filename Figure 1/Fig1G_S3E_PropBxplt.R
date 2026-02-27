library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(tidyverse)
library(ggpubr)
library(rstatix)
library(openxlsx)


# dfinfo has chemistry, sample, disease, location, etc informations
dfinfo<-read.csv2(file="./Grouped_objects/sample_annots.csv",sep=",",header=T)
colnames(dfinfo)[1]<-"sample"
listobj<-c("df3hcco","df3hcil","df3CDuni","df3CD","df3UC")

load("./Grouped_objects/MacsTot_z_mymod_150_metadata_goodlog_v3.rd") # metadf_z     dataframe with module scores and scaled scores by metacells
df_zmod<-metadf_z[,154:length(colnames(metadf_z))]
load("./Grouped_objects/ht_Total_Macs.rd")
load("./Grouped_objects/df_mcall_Macs_select2.rd")
metadata<-merge(df[,1:3], dfinfo[,c("sample","Disease","tissue","status","biotherapy_status","CHEMISTRY","location","subset")], by="sample")

###----------------------------------------------### 
#### Function to put call clusters on metacells ####
###----------------------------------------------### 
ht_clustering<-function(ht,df_zmod){
  metacells_kclust<- column_order(ht) 
  for (i in 1:length(metacells_kclust)){
    if (i == 1) {
      clu <- t(t(rownames(df_zmod)[metacells_kclust[[i]]]))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("Metacell", "Cluster")
    } else {
      clu <- t(t(rownames(df_zmod)[metacells_kclust[[i]]]))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }
  out<-as.data.frame(out)
  out<-out[match(rownames(df_zmod),out$Metacell),]
  df_zmod$cluster<-out$Cluster  # So we don't have the 3 clusters of unknown
  df_zmod$cluster<-factor(df_zmod$cluster, levels = paste0("cluster", 1:40))
  setClass(Class="Myobj1",
           representation(
             out="data.frame",
             df_zmod="data.frame"
           ))
  return(new("Myobj1",out=out,df_zmod=df_zmod))
}


###----------------------------------------------### 
#### Function to get the samples with > X cells ####
###----------------------------------------------### 

TL_samples_gut<-function(nbrcel,df_celtyp,obj_celtype,t_dfn){
  
  df2<-df_celtyp[which(df_celtyp$mc>=0),]
  ## put the corresponding cluster to each metacell
  df2$Cluster<-obj_celtype@out$Cluster[match(df2$mc2,obj_celtype@out$Metacell)]
  df2$Group<-NA
  for(i in t_dfn$Grp_cell){
    df2$Group[which(df2$Cluster%in%paste0('cluster',unlist(strsplit(t_dfn$Grp_cls[which(t_dfn$Grp_cell==i)],split=','))) )]<-i
  }
  
  df3<-cbind(metadata[match(df2$names, metadata$names),c(1:2,4:10)],df2[,c("Cluster","Group")])
  
  df3<-df3[-which(df3$sample=="n14"),]
  df3$status[which(df3$status=="Involved?")]<-"Involved"
  df3_nbr<-df3[which(df3$sample%in%names(which(table(df3$sample)>=nbrcel))),] 
  smpls<- names(which(table(df3$sample)>=nbrcel)) 
  ## Cut each subset distribution :
  df3hc<-df3_nbr[which(df3_nbr$Disease=="CTRL"&df3_nbr$biotherapy_status!="After"),]                            #
  df3UCuni<-df3_nbr[which(df3_nbr$Disease=="UC"&df3_nbr$status=="Uninvolved"&df3_nbr$biotherapy_status!="After"),]  #
  #All Uninvolved UC are from Ileon : take them as Healthy controls
  df3hc_ucu<-rbind(df3hc,df3UCuni)
  df3hcil<-df3hc_ucu[which(df3hc_ucu$tissue=="ILEUM"),]
  df3hcco<-df3hc_ucu[which(df3hc_ucu$tissue=="COLON"),]
  df3UC<-df3_nbr[which(df3_nbr$Disease=="UC"&df3_nbr$status=="Involved"&df3_nbr$biotherapy_status!="After"),]       #
  df3CDuni<-df3_nbr[which(df3_nbr$Disease=="CD"&df3_nbr$status=="Uninvolved"&df3_nbr$biotherapy_status!="After"),]  #
  df3CD<-df3_nbr[which(df3_nbr$Disease=="CD"&df3_nbr$status=="Involved"&df3_nbr$biotherapy_status!="After"),]       #
  
  clust_ord<-levels(as.factor(df2$Group))
  
  
  setClass(Class="Myobj2",
           representation(
             smpls="character",
             df3_nbr="data.frame",
             df3hcco="data.frame",
             df3hcil="data.frame",
             df3CDuni="data.frame",
             df3CD="data.frame",
             df3UC="data.frame",
             clust_ord="character"
           ))
  return(new("Myobj2",smpls=smpls,df3_nbr=df3_nbr,df3hcco=df3hcco,
             df3hcil= df3hcil,df3CDuni=df3CDuni,df3CD=df3CD,df3UC=df3UC,clust_ord=clust_ord))
  
}


###--------------------------------------------### 


# Do the stats for grouped clusters : 
grp_mac<-c("Mono1", 
           "Mono3", "Mono4", "Mono5", 
           "Macro6","Macro7", "Macro8","Macro9")
grp_cls<-c("38,37,31,26,32,34,35,39,40", 
           "36,22", "20,21,3","4,11,5,7,6,29,28",
           "1,33,25", "30,23,27,24", "15,2,14,9,10,8", "12,16,13,19,17,18")
mac_dfn<-data.frame(Grp_cell=grp_mac,Grp_cls=grp_cls)
celtypannot=c('Mono1'='#874037',
              'Mono3'='#f77774','Mono4'='#ffa3a1','Mono5'='#B4DEC6',  
              'Macro6'='#c967eb','Macro7'='#5fc2ed','Macro8'='#4F88B9','Macro9'='#5C538B')

# Samples ordered by category:
df3hcil_o1<-c("214","226","216","218","217","235","227","220") # in Mo/Macs:"215" "225" "234" "236" HC are "After" treatment     "220" has 28 cells
df3CDuni_o1<-c("68","192","189","159","195","208","129","186","180","135")   # in Mo/Macs:"129" has 28 cells
df3CD_o1<-c("190","196","193","158","209","69","n7","n1","n3","n17","n47","n49","n6","n33","n37","128","138","187","181","n13")   # in Mo/Macs: n14 is removed, n49 has 2 cells, n6 has 23 cells
df3UC_o1<-c("238","200","198","202")
df3hcco_o1<-c("237","219","221","204","206")

Macs_myc<-ht_clustering(ht,df_zmod)
Macs_df15<-TL_samples_gut(15,df,Macs_myc,mac_dfn)
#Force order of subsets in the barplot:
Macs_df15@clust_ord<-names(celtypannot)
levels(Macs_df15@clust_ord)<-names(celtypannot)

# Get the MoMacs distribution by sample:
prop_by_pat_CD<-as.data.frame.matrix(prop.table(table(Macs_df15@df3CD$sample,Macs_df15@df3CD$Group),margin = 1) )
prop_by_pat_CD$category<-"Inf."
prop_by_pat_CDuni<-as.data.frame.matrix(prop.table(table(Macs_df15@df3CDuni$sample,Macs_df15@df3CDuni$Group),margin = 1))
prop_by_pat_CDuni$category<-"Uninf."
prop_by_pat_hcil<-as.data.frame.matrix(prop.table(table(Macs_df15@df3hcil$sample,Macs_df15@df3hcil$Group),margin = 1) )
prop_by_pat_hcil$category<-"Ctrls"

# bind_rows allow to put "NA" if column was not present in the dataframe
combined_df <- dplyr::bind_rows(prop_by_pat_hcil, prop_by_pat_CDuni, prop_by_pat_CD)
combined_df[is.na(combined_df)]<-0
combined_df$patient<-rownames(combined_df)

long_df <- gather(combined_df, key = "Subtype", value = "Proportion", -patient, -category)


# Facetwrap helps having freedom of color, with separation by subtype
# Use a factor to force the order of sybtypes
momacpalette=c('Mono1'='#874037',
               'Mono3'='#f77774','Mono4'='#ffa3a1','Mono5'='#77c799',  
               'Macro6'='#c967eb','Macro7'='#41b5e8','Macro8'='#4F88B9','Macro9'='#5C538B')
long_df$Subtype<-factor(long_df$Subtype, levels=names(celtypannot) )


# test without patient 69
long_df2<-long_df[-which(long_df$patient=="69"),]

kw_results<-long_df2 %>%
  group_by(Subtype) %>%
  rstatix::kruskal_test(Proportion ~ category)
dunn_results<-long_df2 %>%
  group_by(Subtype) %>%
  rstatix::dunn_test(Proportion ~ category, p.adjust.method = "holm")

# Keep Dunn results in KW signif subtypes only: 
Sig_subtypes <- kw_results %>%
  filter(p < 0.05) %>%
  pull(Subtype)

dunn_results_filtered <- dunn_results %>%
  filter(Subtype %in% Sig_subtypes)

dunn_results_filtered <- dunn_results_filtered %>% add_y_position()
dunn_results_filtered$y.position<-rep(c(0.9,0.95,1),nrow(dunn_results_filtered)/3)  # better if I manually tell them where to be
dunn_results_filtered$padj_round<-round(dunn_results_filtered$p.adj,4)


Sig_subtypes <- kw_results %>% pull(Subtype)
dunn_results_all <- dunn_results %>% filter(Subtype %in% Sig_subtypes)
dunn_results_all <- dunn_results_all %>% add_y_position()
dunn_results_all$y.position<-rep(c(0.9,0.95,1),nrow(dunn_results_all)/3) 
dunn_results_all$padj_round<-round(dunn_results_all$p.adj,4)
  
pdf("./Figures_print/FigS3E_PropBxpt.pdf",width = 8,height = 5)
ggplot(long_df2, aes(x = category, y = Proportion, fill = Subtype,color = Subtype)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA,size=rel(0.6) )+ # I remove outlier dots as I have jitter, it will avoid misunderstanding
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", fatten = 2) + # Add black median line
  facet_wrap(~ Subtype,ncol = 1,) +
  theme_minimal() +
  theme(strip.text = element_blank(), axis.text.y = element_text(size=rel(1.3)), axis.title.y = element_blank()  )+    #,panel.border = element_rect(colour = "black", fill=NA, size=1) )+  # To remove subplot titles, needs to be put after theme_minimal
  labs(title = "Proportions of Subtypes by Group (only patients with >15 MoMacs, no pat.69)",
       x = "Category",
       y = "Proportion",
       fill = "Subtype") +
  scale_fill_manual(values =celtypannot)+scale_x_discrete(limits=c("Inf.","Uninf.","Ctrls"))+
  scale_color_manual(values =celtypannot)+
  geom_jitter(height = 0,size=rel(1), width = 0.25,alpha=0.5,aes(color = Subtype))+
  stat_pvalue_manual(dunn_results_filtered, label = "p.adj.signif", tip.length = 0,hide.ns = TRUE,coord.flip = TRUE,size=rel(8.5),bracket.size = rel(1)) +
  coord_flip()

ggplot(long_df2, aes(x = category, y = Proportion, fill = Subtype,color = Subtype)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA,size=rel(0.6) )+ # I remove outlier dots as I have jitter, it will avoid misunderstanding
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", fatten = 2) + # Add black median line
  facet_wrap(~ Subtype,ncol = 1,) +
  theme_minimal() +
  theme(strip.text = element_blank(), axis.text.y = element_text(size=rel(1.3)), axis.title.y = element_blank()  )+    #,panel.border = element_rect(colour = "black", fill=NA, size=1) )+  # To remove subplot titles, needs to be put after theme_minimal
  labs(title = "Proportions of Subtypes by Group (only patients with >15 MoMacs, no pat.69)",
       x = "Category",
       y = "Proportion",
       fill = "Subtype") +
  scale_fill_manual(values =celtypannot)+scale_x_discrete(limits=c("Inf.","Uninf.","Ctrls"))+
  scale_color_manual(values =celtypannot)+
  geom_jitter(height = 0,size=rel(1), width = 0.25,alpha=0.5,aes(color = Subtype))+
  stat_pvalue_manual(dunn_results_filtered, label = "padj_round", tip.length = 0,hide.ns = TRUE,coord.flip = TRUE,size=rel(2.5),bracket.size = rel(1)) +
  coord_flip()

# Show all p-val even the non signif ones:
ggplot(long_df2, aes(x = category, y = Proportion, fill = Subtype,color = Subtype)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA,size=rel(0.6) )+ # I remove outlier dots as I have jitter, it will avoid misunderstanding
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", fatten = 2) + # Add black median line
  facet_wrap(~ Subtype,ncol = 1,) +
  theme_minimal() +
  theme(strip.text = element_blank(), axis.text.y = element_text(size=rel(1.3)), axis.title.y = element_blank()  )+    #,panel.border = element_rect(colour = "black", fill=NA, size=1) )+  # To remove subplot titles, needs to be put after theme_minimal
  labs(title = "Proportions of Subtypes by Group (only patients with >15 MoMacs, no pat.69)",
       x = "Category",
       y = "Proportion",
       fill = "Subtype") +
  scale_fill_manual(values =celtypannot)+scale_x_discrete(limits=c("Inf.","Uninf.","Ctrls"))+
  scale_color_manual(values =celtypannot)+
  geom_jitter(height = 0,size=rel(1), width = 0.25,alpha=0.5,aes(color = Subtype))+
  stat_pvalue_manual(dunn_results_all, label = "padj_round", tip.length = 0,hide.ns = FALSE,coord.flip = TRUE,size=rel(2.5),bracket.size = rel(1)) +
  coord_flip()

dev.off()  

# write.csv2(long_df2,file="Gaelle_Fig1G_datapoints.csv",quote=F,row.names=F)


old.cluster.ids<-c("Mono1","Mono2","Mono3","Mono4","Mono5",
                   "Macro6","Macro7","Macro8","Macro9")    

new.cluster.ids <- c("IFIM","IFN_Mono","IRM","Early_Mono","Late_Mono",
                     "IFN_Macro","Infl_Macro","FOLR2neg_Macro","FOLR2pos_Macro")

kw_results$Subtype<- plyr::mapvalues(x = kw_results$Subtype, from = old.cluster.ids, to = new.cluster.ids)
dunn_results_filtered$Subtype<- plyr::mapvalues(x = dunn_results_filtered$Subtype, from = old.cluster.ids, to = new.cluster.ids)
dunn_results_all$Subtype<- plyr::mapvalues(x = dunn_results_all$Subtype, from = old.cluster.ids, to = new.cluster.ids)

wb <- createWorkbook()
addWorksheet(wb, "Dunn_Prop")
writeData(wb, "Dunn_Prop", kw_results[,c(1,3:6)], startRow = 2)
writeData(wb, "Dunn_Prop", dunn_results_filtered[,c(1,3:10)], startRow = nrow(kw_results)+4)
saveWorkbook(wb, "./Figure 1/Fig1G_S3E_PropBxpt_stats.xlsx", overwrite = TRUE)







# Prism-like figure :
colnames(dunn_results_filtered)[1]<-"category" 
m1_df_pvals<-dunn_results_filtered[which(dunn_results_filtered$category=="IFIM"),]
m1_long_df2<-long_df2[which(long_df2$Subtype=="Mono1"),]
m3_long_df2<-long_df2[which(long_df2$Subtype=="Mono3"),]
m4_long_df2<-long_df2[which(long_df2$Subtype=="Mono4"),]
m5_df_pvals<-dunn_results_filtered[which(dunn_results_filtered$category=="Late_Mono"),]
m5_long_df2<-long_df2[which(long_df2$Subtype=="Mono5"),]
m6_long_df2<-long_df2[which(long_df2$Subtype=="Macro6"),]
m7_long_df2<-long_df2[which(long_df2$Subtype=="Macro7"),]
m8_long_df2<-long_df2[which(long_df2$Subtype=="Macro8"),]
m8_df_pvals<-dunn_results_filtered[which(dunn_results_filtered$category=="FOLR2neg_Macro"),]
m9_long_df2<-long_df2[which(long_df2$Subtype=="Macro9"),]
m9_df_pvals<-dunn_results_filtered[which(dunn_results_filtered$category=="FOLR2pos_Macro"),]

palette<-c("#EEEEEE","#999999","#000000")
# m1_df_pvals$y.position<-df_pvals$y.position-0.125   # lower the stat bar

  
prismplot<-function(long_df2){
  ggplot(long_df2, aes(x =category, y =Proportion, fill=category))+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
  scale_x_discrete(limits = c("Ctrls","Uninf.", "Inf."), labels=c("Ctrls.","Uninfl.","Infl."))+
  scale_color_manual(values = palette)+
  scale_y_continuous(n.breaks = 3) +
  ggprism::theme_prism(base_size = 16)+
  labs(y= "Monocyte frequency\nin Macrophage subset")+
  stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = category), fatten = 4)+ # Visualize the median on top, allows to adapt color
  scale_color_manual(values = c("black", "white", "black")) +
  geom_jitter(shape=21, position=position_jitter(width = 0.24),size=5, fill="white",stroke=1.8)+
  theme(legend.position = "none",axis.title.x =element_blank(),
        axis.text.y = element_text(size=rel(2.75),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(1)),
        axis.ticks.length = unit(15,"points"),
        plot.margin = unit(c(40, 4, 4, 4),"points"))+   #add margin above plot so asterisks are not cropped
  coord_cartesian(ylim = c(0, 1),clip = 'off')#+  #cut the plot at 1
# stat_compare_means(label.y = max(df4$Freq)+0.5,size = 7) # show K-W value
}

for(i in c("m1_df_pvals","m5_df_pvals","m8_df_pvals","m9_df_pvals")){
  eval(parse(text=paste0(i,"$padj_round<-round(",i,"$p.adj,4)" )))
}


m1_plot_prism<-prismplot(m1_long_df2)+
  stat_pvalue_manual(m1_df_pvals, label = "p.adj.signif",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(25))
m1_plot_prism2<-prismplot(m1_long_df2)+
  stat_pvalue_manual(m1_df_pvals, label = "padj_round",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(10))
m3_plot_prism<-prismplot(m3_long_df2)
m4_plot_prism<-prismplot(m4_long_df2)
m5_plot_prism<-prismplot(m5_long_df2)+
  stat_pvalue_manual(m5_df_pvals, label = "p.adj.signif",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(25))
m5_plot_prism2<-prismplot(m5_long_df2)+
  stat_pvalue_manual(m5_df_pvals, label = "padj_round",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(10))
m6_plot_prism<-prismplot(m6_long_df2)
m7_plot_prism<-prismplot(m7_long_df2)
m8_plot_prism<-prismplot(m8_long_df2)+
  stat_pvalue_manual(m8_df_pvals, label = "p.adj.signif",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(25))
m8_plot_prism2<-prismplot(m8_long_df2)+
  stat_pvalue_manual(m8_df_pvals, label = "padj_round",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(10))
m9_plot_prism<-prismplot(m9_long_df2)+
  stat_pvalue_manual(m9_df_pvals, label = "p.adj.signif",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(25))
m9_plot_prism2<-prismplot(m9_long_df2)+
  stat_pvalue_manual(m9_df_pvals, label = "padj_round",tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(10))



pdf("./Figures_print/Fig1G_PropBxpt_stats.pdf",width = 5,height = 7.5)
print(m1_plot_prism+ labs(y= "IFIM frequency\nin Macrophage subset"))
print(m1_plot_prism2+ labs(y= "IFIM frequency\nin Macrophage subset"))
print(m3_plot_prism+ labs(y= "IRM frequency\nin Macrophage subset"))
print(m4_plot_prism+ labs(y= "Early Mono frequency\nin Macrophage subset"))
print(m5_plot_prism+ labs(y= "Late Mono frequency\nin Macrophage subset"))
print(m5_plot_prism2+ labs(y= "Late Mono frequency\nin Macrophage subset"))
print(m6_plot_prism+ labs(y= "IFN Mac frequency\nin Macrophage subset"))
print(m7_plot_prism+ labs(y= "Infl Mac frequency\nin Macrophage subset"))
print(m8_plot_prism+ labs(y= "FOLR2- Mac frequency\nin Macrophage subset"))
print(m8_plot_prism2+ labs(y= "FOLR2- Mac frequency\nin Macrophage subset"))
print(m9_plot_prism+ labs(y= "FOLR2+ Mac frequency\nin Macrophage subset"))
print(m9_plot_prism2+ labs(y= "FOLR2+ Mac frequency\nin Macrophage subset"))
dev.off()
