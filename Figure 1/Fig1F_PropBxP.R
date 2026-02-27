library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(ggpubr)
library(rstatix)
library(openxlsx)

dfinfo<-read.csv2(file="./Grouped_objects/sample_annots.csv",sep=",",header=T) # dfinfo has sample metadata: chemistry, sample, disease, location, etc informations
colnames(dfinfo)[1]<-"sample"
listobj<-c("df3hcco","df3hcil","df3CDuni","df3CD","df3UC") # The 5 conditions: healthy control (colon), healthy control (ileum), uninflamed ileum, inflamed ileum, inflamed colon 
# Loading data of the total cells dataset 
load("./Grouped_objects/MacsTot_z_mymod_150_metadata_goodlog_v3.rd") # metadf_z dataframe with module scores and scaled scores by metacells
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


###------------------------------------### 
#### Function for the proportion plot ####
###------------------------------------### 

TL_prop_plot_gut<-function(celtypannot,df3){
  
  # V4: To accomodate with empty spaces when I give my own order including patients with less cells than expected :
  design_perso=paste0(c(rep("A",length(unique(df3hcil_o1))),rep("E",length(unique(df3CDuni_o1))),rep("F",length(unique(df3CD_o1))),"\n",
                        rep("B",length(unique(df3hcil_o1))),rep("E",length(unique(df3CDuni_o1))),rep("F",length(unique(df3CD_o1))),"\n",
                        rep("C",length(unique(df3hcil_o1))),rep("E",length(unique(df3CDuni_o1))),rep("F",length(unique(df3CD_o1))),"\n",
                        rep("D",length(unique(df3hcil_o1))),rep("E",length(unique(df3CDuni_o1))),rep("F",length(unique(df3CD_o1)))),collapse = "")
  
  
  noyaxe<-theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  
  for(i in 1:length(listobj)){
    df4<-eval(parse(text=paste0("df3@",listobj[i]) ))  
    df4$Group<-factor(df4$Group,levels=df3@clust_ord)
    eval(parse(text=paste0("pt",i,"<-ggplot(df4,aes(x= sample,fill= Group)) +  geom_bar(position = 'fill')+
  scale_fill_manual(values=celtypannot )+  
  theme_classic()+labs( x = ' ', y = 'Proportion',fill = 'Celltype')+ 
                         scale_x_discrete(limits =",listobj[i],"_o1)+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))")))
    # Define df for metadata rows
    dfannot<-data.frame(patient=levels(as.factor(df4$sample)),origin=df4[match(levels(as.factor(df4$sample)), df4$sample), "location"],
                        chem=df4[match(levels(as.factor(df4$sample)), df4$sample), "CHEMISTRY"],
                        disease=paste0(df4[match(levels(as.factor(df4$sample)), df4$sample), "Disease"],df4[match(levels(as.factor(df4$sample)), df4$sample), "status"]),
                        tissue=df4[match(levels(as.factor(df4$sample)), df4$sample), "tissue"] )
    dfannot2 <- dfannot %>%mutate(disease = if_else(disease == "CDInvolved", "CD",disease))
    dfannot <- dfannot2 %>%mutate(disease = if_else(disease == "CTRLCONTROL", "HC", disease))
    dfannot2 <- dfannot %>%mutate(disease = if_else(disease == "CDUninvolved", "CDuninv",disease))
    dfannot <- dfannot2 %>%mutate(disease = if_else(disease == "UCInvolved", "UC", disease))
    dfannot2 <- dfannot %>%mutate(disease = if_else(disease == "UCUninvolved", "UCuninv",disease))
    dfannot <- dfannot2 %>%mutate(tissue = if_else(tissue == "ILEUM", "ileum", tissue))
    dfannot2 <- dfannot %>%mutate(tissue = if_else(tissue == "COLON", "colon",tissue))
    dfannot<-dfannot2
    #disease
    eval(parse(text=paste0("pt",i,"_2 <- ggplot(dfannot, aes(x=patient, y=1, fill = disease)) +
    geom_tile() + guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.4))+
    theme_void() + scale_x_discrete(limits =",listobj[i],"_o1)+
    scale_fill_manual(values=c(HC='#abcbd6',CD='#c22d4f', CDuninv='#c97185',UC='#314d73', UCuninv='#5877a1'))  ")))
    #tissue
    eval(parse(text=paste0("pt",i,"_3 <- ggplot(dfannot, aes(x=patient, y=1, fill = tissue)) +
    geom_tile() + guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.4))+
    theme_void() + scale_x_discrete(limits =",listobj[i],"_o1)+
    scale_fill_manual(values=c(ileum='#b3aca4', colon='#5e5b57'))  ")))
    #origin
    eval(parse(text=paste0("pt",i,"_4 <- ggplot(dfannot, aes(x=patient, y=1, fill = origin)) +
    geom_tile() + guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.4))+
    theme_void() + scale_x_discrete(limits =",listobj[i],"_o1)+
    scale_fill_manual(values=c(mtsinai='#17b6ff', nantes='#cf8013'))   ")))
  }
  pt1_1<-pt1+guides(fill= "none")+ggtitle("HC colon")
  pt1_2<-pt1_2+guides(fill= "none")
  pt1_3<-pt1_3+guides(fill= "none")
  pt1_4<-pt1_4+guides(fill= "none")
  plot1<-((pt1_1 / pt1_2) / pt1_3) / pt1_4 + 
    plot_layout(widths = c(2, 1), nrow=4, heights = c(6,1.2,1.2,1.2))
  pt2_1<-pt2+guides(fill= "none")+ggtitle("HC ileum") # +noyaxe   pour Mono/Mac ratio j'enleve les colons donc besoin annot ici
  pt2_2<-pt2_2+guides(fill= "none")
  pt2_3<-pt2_3+guides(fill= "none")
  pt2_4<-pt2_4+guides(fill= "none")
  plot2<-((pt2_1 / pt2_2) / pt2_3) / pt2_4 + 
    plot_layout(widths = c(2, 1), nrow=4, heights = c(6,1.2,1.2,1.2))
  pt3_1<-pt3+guides(fill= "none")+ggtitle("CDuninv")+noyaxe
  pt3_2<-pt3_2+guides(fill= "none")
  pt3_3<-pt3_3+guides(fill= "none")
  pt3_4<-pt3_4+guides(fill= "none")
  plot3<-((pt3_1 / pt3_2) / pt3_3) / pt3_4 + 
    plot_layout(widths = c(2, 1), nrow=4, heights = c(6,1.2,1.2,1.2))
  pt4_1<-pt4+ggtitle("CD")+noyaxe   #+guides(fill= "none") pour Mono/Mac ratio j'enleve les UC donc besoin l?gende ici
  pt4_2<-pt4_2#+guides(fill= "none") 
  pt4_3<-pt4_3#+guides(fill= "none")
  pt4_4<-pt4_4#+guides(fill= "none")
  plot4<-((pt4_1 / pt4_2) / pt4_3) / pt4_4 + 
    plot_layout(widths = c(2, 1), nrow=4, heights = c(6,1.2,1.2,1.2))
  pt5_1<-pt5+ggtitle("UC")+noyaxe
  pt5_2<-pt5_2
  pt5_3<-pt5_3
  pt5_4<-pt5_4
  plot5<-((pt5_1 / pt5_2) / pt5_3) / pt5_4 + 
    plot_layout(widths = c(2, 1), nrow=4, heights = c(6,1.2,1.2,1.2))
  
  return(list(pt1_1=pt1_1,pt1_2=pt1_2,pt1_3=pt1_3,pt1_4=pt1_4,plot1=plot1,
              pt2_1=pt2_1,pt2_2=pt2_2,pt2_3=pt2_3,pt2_4=pt2_4,plot2=plot2,
              pt3_1=pt3_1,pt3_2=pt3_2,pt3_3=pt3_3,pt3_4=pt3_4,plot3=plot3,
              pt4_1=pt4_1,pt4_2=pt4_2,pt4_3=pt4_3,pt4_4=pt4_4,plot4=plot4,
              pt5_1=pt5_1,pt5_2=pt5_2,pt5_3=pt5_3,pt5_4=pt5_4,plot5=plot5,
              design_perso=design_perso))
}


###--------------------------------------------### 
#### Function for the SMALLER proportion plot ####
###--------------------------------------------### 

TL_prop_plot_gut_smol<-function(celtypannot,df3){
  
  # V4: To accomodate with empty spaces when I give my own order including patients with less cells than expected :
  design_perso=paste0(c(rep("A",length(unique(df3hcil_o1))),rep("B",length(unique(df3CDuni_o1))),rep("C",length(unique(df3CD_o1)))),collapse = "")
  
  noyaxe<-theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  for(i in 1:length(listobj)){
    df4<-eval(parse(text=paste0("df3@",listobj[i]) ))  
    df4$Group<-factor(df4$Group,levels=df3@clust_ord)
    if(listobj[i]=="df3hcil"){df4$title="Ctrls"}else if(listobj[i]=="df3CDuni"){df4$title="Uninfl."}else if(listobj[i]=="df3CD"){df4$title="Inflamed"}else{df4$title="not found"}
    eval(parse(text=paste0("pt",i,"<-ggplot(df4,aes(x= sample,fill= Group)) +  geom_bar(position = 'fill')+
  scale_fill_manual(values=celtypannot )+  
  theme_classic()+labs( x = ' ', y = 'Proportion',fill = 'Celltype')+ 
                         scale_x_discrete(limits =",listobj[i],"_o1,labels=c('n1'='GIM7','n3'='GIM8','n7'='GIM23','n13'='GIM31','n17'='GIM33','n33'='GIM35','n37'='GIM36','n47'='GIM38'))+
                           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))")))

  }
  

  pt2_1<-pt2+guides(fill= "none")+#ggtitle("HC ileum")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#ABCBD6",0.5), color = alpha("#ABCBD6",0.5)),  
          strip.text = element_text(size = rel(2.2) )) 

  pt3_1<-pt3+guides(fill= "none")+noyaxe+#ggtitle("CDuninv")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#C97185",0.5), color = alpha("#C97185",0.5)),  
          strip.text = element_text(size = rel(2.2) )) 

  pt4_1<-pt4+noyaxe+ #ggtitle("CD")+ 
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#C22D4F",0.5), color = alpha("#C22D4F",0.5)),  
          strip.text = element_text(size = rel(2.2) )) 

  
  return(list(pt2_1=pt2_1,
              pt3_1=pt3_1,
              pt4_1=pt4_1,
              design_perso=design_perso))
}


###------------------------------------###


# Assign cluster to metacell + get module scores:
Macs_myc<-ht_clustering(ht,df_zmod)
# Define which clusters are Mono, which are Macro:
grp_mac<-c("Monocytes","Macrophages")
grp_cls<-c("38,37,31,26,32,34,35,39,40,36,22,20,21,3,4,11,5,7,6,29,28",
           "1,33,25,30,23,27,24,15,2,14,9,10,8,12,16,13,19,17,18")
celtypannot=c('Monocytes'='#CF6A11','Macrophages'='#4080D9')
mac_dfn2<-data.frame(Grp_cell=grp_mac,Grp_cls=grp_cls)

# Group cells by condition with metadata, threshold for "good" sample: n=30 cells min
Macs_df3<-TL_samples_gut(30,df,Macs_myc,mac_dfn2)


# Sort samples by conditions based on the Macs %:
ordsimple<-function(data=NULL){
  protab<-prop.table(table(data$Group,data$sample),margin=2)
  ord2<-colnames(protab)[order(protab["Macrophages",])]
  return(ord2)}
###------------------------------------### 

df3hcil_o1<-ordsimple(data=Macs_df3@df3hcil)
df3hcco_o1<-ordsimple(data=Macs_df3@df3hcco)
df3UC_o1<-ordsimple(data=Macs_df3@df3UC)
df3CDuni_o1<-ordsimple(data=Macs_df3@df3CDuni)
df3CD_o1<-ordsimple(data=Macs_df3@df3CD)

#Proportion plot
plot<-TL_prop_plot_gut(celtypannot,Macs_df3)
plot_smol<-TL_prop_plot_gut_smol(celtypannot,Macs_df3)



# Figure 1E: proportions as boxplots:

df4hcil<-as.data.frame(prop.table(table(Macs_df3@df3hcil$sample,Macs_df3@df3hcil$Group),margin = 1))
df4hcil$category<-"HC ileum"
df4CD<-as.data.frame(prop.table(table(Macs_df3@df3CD$sample,Macs_df3@df3CD$Group),margin = 1))
df4CD$category<-"CD"
df4CDuni<-as.data.frame(prop.table(table(Macs_df3@df3CDuni$sample,Macs_df3@df3CDuni$Group),margin = 1))
df4CDuni$category<-"CD uninv."

df4<-rbind(rbind(df4hcil,df4CD),df4CDuni)
df4<-df4[-which(df4$Var2=="Macrophages"),]


my_comparisons <- list( c("HC ileum", "CD uninv."), c("CD uninv.", "CD"), c("HC ileum", "CD") )
df4$category<- factor(df4$category, levels= c("HC ileum","CD uninv.", "CD") )
palette<-c("#EEEEEE","#999999","#000000")


# Stats:
# K-W + Dunn post hoc test
KW_rstatix<-df4 %>%rstatix::kruskal_test(Freq ~ category)
dunn_rstatix<-df4 %>%rstatix::dunn_test(Freq ~ category, p.adjust.method = "holm")
df_pvals <- dunn_rstatix %>% add_xy_position() %>% mutate(category="HC ileum") # creates a "category" so results can be read in my ggpubr/ggplot2 way of plotting


df_pvals$y.position<-df_pvals$y.position-0.125   # lower the stat bar
df_pvals$padj_round<-round(df_pvals$p.adj,4)

plot_prism<-function(p_lab=p_lab,psiz=25){
  print(ggplot(df4, aes(x =category, y =Freq, fill=category, p_lab=p_lab))+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
  scale_x_discrete(limits = c("HC ileum","CD uninv.", "CD"), labels=c("Ctrls.","Uninfl.","Infl."))+
  scale_color_manual(values = palette)+
  scale_y_continuous(n.breaks = 3) +
  ggprism::theme_prism(base_size = 16)+
  labs(y= "Monocyte frequency\nin Macrophage subset")+
  stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = category), fatten = 4)+ # Visualize the median on top, allows to adapt color
  scale_color_manual(values = c("black", "black", "white")) +
  geom_jitter(shape=21, position=position_jitter(width = 0.24),size=5, fill="white",stroke=1.8)+
  stat_pvalue_manual(df_pvals, label = p_lab,tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(psiz))+
  theme(legend.position = "none",axis.title.x =element_blank(),
        axis.text.y = element_text(size=rel(2.75),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(1)),
        axis.ticks.length = unit(15,"points"),
        plot.margin = unit(c(40, 4, 4, 4),"points"))+   #add margin above plot so asterisks are not cropped
  coord_cartesian(ylim = c(0, 1),clip = 'off') )#+  #cut the plot at 1
  # stat_compare_means(label.y = max(df4$Freq)+0.5,size = 7) # show K-W value

}
  
  

pdf("./Figures_print/Fig1F_BoxP.pdf",width = 5,height = 7.5)
plot_prism(p_lab="padj_round",psiz=14)
plot_prism(p_lab="p.adj.signif")
dev.off()

# Save the stats values
wb <- createWorkbook()
addWorksheet(wb, "KW_Dunn")
xl_row <- 1
writeData(wb, "KW_Dunn", paste0("Kruskal-Wallis:"), startRow = xl_row, colNames = FALSE)
xl_row <- xl_row + 1
writeData(wb, "KW_Dunn", KW_rstatix[,2:5], startRow = xl_row)
xl_row <- xl_row + nrow(KW_rstatix) + 1
writeData(wb, "KW_Dunn", paste0("Dunn post-hoc test:"), startRow = xl_row, colNames = FALSE)
xl_row <- xl_row + 1
writeData(wb, "KW_Dunn", df_pvals[,c(2:9)], startRow = xl_row)
saveWorkbook(wb, "./Figure 1/Fig1F_BoxP_stats.xlsx", overwrite = TRUE)

