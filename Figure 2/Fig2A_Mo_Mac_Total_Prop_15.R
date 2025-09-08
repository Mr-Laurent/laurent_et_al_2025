library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(ggpubr)
library(rstatix)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
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


###--------------------------------------------------### 
#### Function to get the samples with > X cells (v2)####
###--------------------------------------------------### 


TL_samples_gut<-function(nbrcel,df_celtyp,obj_celtype,t_dfn,remove_ct=NULL){
  
  df2<-df_celtyp[which(df_celtyp$mc>=0),]
  ## put the corresponding cluster to each metacell
  df2$Cluster<-obj_celtype@out$Cluster[match(df2$mc2,obj_celtype@out$Metacell)]
  df2$Group<-NA
  for(i in t_dfn$Grp_cell){
    df2$Group[which(df2$Cluster%in%paste0('cluster',unlist(strsplit(t_dfn$Grp_cls[which(t_dfn$Grp_cell==i)],split=','))) )]<-i
  }
  # If there are celltypes to remove:
  if(!is.null(remove_ct)){
    df2<-df2[-which(df2$Group%in%remove_ct),] 
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


# Annotation in the total cells cohort:
grp_mac<-c("Mono1",
           "Mono3","Mono4","Mono5",
           "Macro6","Macro7","Macro8","Macro9")
grp_cls<-c("38,37,31,26,32,34,35,39,40", 
           "36,22", "20,21,3","4,11,5,7,6,29,28",
           "1,33,25", "30,23,27,24", "15,2,14,9,10,8", "12,16,13,19,17,18")
mac_dfn<-data.frame(Grp_cell=grp_mac,Grp_cls=grp_cls)


Mono_celtypannot=c('Mono1'='#874037','Mono3'='#f77774','Mono4'='#ffa3a1','Mono5'='#B4DEC6')
Macs_celtypannot=c('Macro6'='#c967eb','Macro7'='#5fc2ed','Macro8'='#4F88B9','Macro9'='#5C538B')
Mono_names=c("Mono1","Mono3","Mono4","Mono5")
Macs_names=c("Macro6","Macro7","Macro8","Macro9")

# Order by condition were computed on the total cell composition
df3hcil_o1<-c("214","226","216","218","217","235","227","220") # in Mo/Macs:"215" "225" "234" "236" HC are "After" treatment     "220" has 28 cells
df3CDuni_o1<-c("68","192","189","159","195","208","129","186","180","135")   # in Mo/Macs:"129" has 28 cells
df3CD_o1<-c("190","196","193","158","209","69","n7","n1","n3","n17","n47","n49","n6","n33","n37","128","138","187","181","n13")   # in Mo/Macs: n14 is removed, n49 has 2 cells, n6 has 23 cells
df3UC_o1<-c("238","200","198","202")
df3hcco_o1<-c("237","219","221","204","206")



###--------------------------------------------### 
#### Function for the SMALLER proportion plot ####
###--------------------------------------------### 

TL_prop_plot_gut_ann<-function(celtypannot,df3){
  
  # To accomodate with empty spaces when I give my own order including patients with less cells than expected :
  design_perso=paste0(c(rep("A",length(unique(df3hcco_o1))),rep("B",length(unique(df3hcil_o1))),rep("C",length(unique(df3CDuni_o1))),rep("D",length(unique(df3CD_o1))),rep("E",length(unique(df3UC_o1))) ),collapse = "")
  
  noyaxe<-theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  for(i in 1:length(listobj)){
    df4<-eval(parse(text=paste0("df3@",listobj[i]) ))  
    df4$Group<-factor(df4$Group,levels=df3@clust_ord)
    if(listobj[i]=="df3hcil"){df4$title="Ctrls\nileum"}else if(listobj[i]=="df3hcco"){df4$title="Ctrls\ncolon"}
    else if(listobj[i]=="df3CDuni"){df4$title="Uninfl.\nileum"}else if(listobj[i]=="df3UC"){df4$title="Inflam.\ncolon"}
    else if(listobj[i]=="df3CD"){df4$title="Inflamed\nileum"}else{df4$title="not found"}
    eval(parse(text=paste0("pt",i,"<-ggplot(df4,aes(x= sample,fill= Group)) +  geom_bar(position = 'fill')+
  scale_fill_manual(values=celtypannot )+  
  theme_classic()+labs( x = ' ', y = 'Proportion',fill = 'Celltype')+ 
                         scale_x_discrete(limits =",listobj[i],"_o1,labels=c('n1'='GIM7','n3'='GIM8','n6'='GIM21','n7'='GIM23','n13'='GIM31','n17'='GIM33','n33'='GIM35','n37'='GIM36','n47'='GIM38','n49'='GIM39'))+
                           theme(axis.text.y = element_text(size=rel(1.8)),axis.text.x = element_text(size=rel(1.5),angle = 90, hjust = 1, vjust = 0.5))")))
    
  }
  
  pt1_1<-pt1+guides(fill= "none")+#ggtitle("HC colon")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#5877a1",0.5), color = alpha("#5877a1",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  
  pt2_1<-pt2+guides(fill= "none")+noyaxe+#ggtitle("HC ileum")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#ABCBD6",0.5), color = alpha("#ABCBD6",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  
  pt3_1<-pt3+guides(fill= "none")+noyaxe+#ggtitle("CDuninv")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#C97185",0.5), color = alpha("#C97185",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  
  pt4_1<-pt4+guides(fill= "none")+noyaxe+ #ggtitle("CD")+ 
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#C22D4F",0.5), color = alpha("#C22D4F",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  pt5_1<-pt5+noyaxe+ #ggtitle("CD")+ 
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#314d73",0.5), color = alpha("#314d73",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  
  return(list(pt1_1=pt1_1,
              pt2_1=pt2_1,
              pt3_1=pt3_1,
              pt4_1=pt4_1,
              pt5_1=pt5_1,
              design_perso=design_perso))
}


###--------------------------------------------### 
#### Function for the SMALLER proportion plot ####
###--------------------------------------------### 

TL_prop_plot_gut_ann_smol<-function(celtypannot,df3){
  
  # To accomodate with empty spaces when I give my own order including patients with less cells than expected :
  design_perso=paste0(c(rep("A",length(unique(df3hcil_o1))),rep("B",length(unique(df3CDuni_o1))),rep("C",length(unique(df3CD_o1))) ),collapse = "")
  
  noyaxe<-theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  for(i in 1:length(listobj)){
    df4<-eval(parse(text=paste0("df3@",listobj[i]) ))  
    df4$Group<-factor(df4$Group,levels=df3@clust_ord)
    if(listobj[i]=="df3hcil"){df4$title="Ctrls\nileum"}else if(listobj[i]=="df3hcco"){df4$title="Ctrls\ncolon"}
    else if(listobj[i]=="df3CDuni"){df4$title="Uninfl.\nileum"}else if(listobj[i]=="df3UC"){df4$title="Inflam.\ncolon"}
    else if(listobj[i]=="df3CD"){df4$title="Inflamed\nileum"}else{df4$title="not found"}
    eval(parse(text=paste0("pt",i,"<-ggplot(df4,aes(x= sample,fill= Group)) +  geom_bar(position = 'fill')+
  scale_fill_manual(values=celtypannot )+  
  theme_classic()+labs( x = ' ', y = 'Proportion',fill = 'Celltype')+ 
                         scale_x_discrete(limits =",listobj[i],"_o1,labels=c('n1'='GIM7','n3'='GIM8','n6'='GIM21','n7'='GIM23','n13'='GIM31','n17'='GIM33','n33'='GIM35','n37'='GIM36','n47'='GIM38','n49'='GIM39'))+ 
                           theme(axis.text.y = element_text(size=rel(1.8)),axis.text.x = element_text(size=rel(1.5),angle = 90, hjust = 1, vjust = 0.5))")))
    
  }
  ## Remove proportions from Colon samples, not needed in main figures:
  
  # pt1_1<-pt1+guides(fill= "none")+#ggtitle("HC colon")+
  #   facet_grid(. ~ title) +
  #   theme(strip.background = element_rect(fill = alpha("#5877a1",0.5), color = alpha("#5877a1",0.5)),  
  #         strip.text = element_text(size = rel(1.3) )) 
  
  pt2_1<-pt2+guides(fill= "none")+#ggtitle("HC ileum")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#ABCBD6",0.5), color = alpha("#ABCBD6",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  
  pt3_1<-pt3+guides(fill= "none")+noyaxe+#ggtitle("CDuninv")+
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#C97185",0.5), color = alpha("#C97185",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  
  pt4_1<-pt4+guides(fill= "none")+noyaxe+ #ggtitle("CD")+ 
    facet_grid(. ~ title) +
    theme(strip.background = element_rect(fill = alpha("#C22D4F",0.5), color = alpha("#C22D4F",0.5)),  
          strip.text = element_text(size = rel(1.3) )) 
  # pt5_1<-pt5+noyaxe+ #ggtitle("CD")+ 
  #   facet_grid(. ~ title) +
  #   theme(strip.background = element_rect(fill = alpha("#314d73",0.5), color = alpha("#314d73",0.5)),  
  #         strip.text = element_text(size = rel(1.3) )) 
  
  return(list(
              pt2_1=pt2_1,
              pt3_1=pt3_1,
              pt4_1=pt4_1,
              
              design_perso=design_perso))
}


###------------------------------------###
# 30 cells doesn't work : too few cells to see plots in HC 


Mono_myc<-ht_clustering(ht,df_zmod)

Mono_df15<-TL_samples_gut(15,df,Mono_myc,mac_dfn, remove_ct=Macs_names)
#Force order of subsets in the barplot:
Mono_df15@clust_ord<-Mono_names
levels(Mono_df15@clust_ord)<-Mono_names
Mono_plot15<-TL_prop_plot_gut_ann(Mono_celtypannot,Mono_df15)
Mono_s_plot15<-TL_prop_plot_gut_ann_smol(Mono_celtypannot,Mono_df15)

Macs_df15<-TL_samples_gut(15,df,Mono_myc,mac_dfn, remove_ct=Mono_names)
#Force order of subsets in the barplot:
Macs_df15@clust_ord<-Macs_names
levels(Macs_df15@clust_ord)<-Macs_names
Macs_plot15<-TL_prop_plot_gut_ann(Macs_celtypannot,Macs_df15)
Macs_s_plot15<-TL_prop_plot_gut_ann_smol(Macs_celtypannot,Macs_df15)


pdf("./Figure 2/Fig2A_MoMacTot_Prop_15.pdf",width = 9,height = 6)
Mono_s_plot15$pt2_1+Mono_s_plot15$pt3_1+Mono_s_plot15$pt4_1+guides(fill= "none")+plot_layout(design = Mono_s_plot15$design_perso)
Macs_s_plot15$pt2_1+Macs_s_plot15$pt3_1+Macs_s_plot15$pt4_1+guides(fill= "none")+plot_layout(design = Macs_s_plot15$design_perso)
Mono_plot15$pt1_1+Mono_plot15$pt2_1+Mono_plot15$pt3_1+Mono_plot15$pt4_1+Mono_plot15$pt5_1+plot_layout(design = Mono_plot15$design_perso)
Macs_plot15$pt1_1+Macs_plot15$pt2_1+Macs_plot15$pt3_1+Macs_plot15$pt4_1+Macs_plot15$pt5_1+plot_layout(design = Macs_plot15$design_perso)
dev.off()



