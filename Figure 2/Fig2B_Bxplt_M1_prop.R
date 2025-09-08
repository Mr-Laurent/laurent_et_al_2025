library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
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


###------------------------------------### 

# Annotation in the total cells cohort:
grp_mac<-c("Mono1",
           "Mono3","Mono4","Mono5",
           "Macro6","Macro7","Macro8","Macro9")
grp_cls<-c("38,37,31,26,32,34,35,39,40", 
           "36,22", "20,21,3","4,11,5,7,6,29,28",
           "1,33,25", "30,23,27,24", "15,2,14,9,10,8", "12,16,13,19,17,18")
mac_dfn<-data.frame(Grp_cell=grp_mac,Grp_cls=grp_cls)
Macs_names=c("Macro6","Macro7","Macro8","Macro9")

Mono_myc<-ht_clustering(ht,df_zmod)
Mono_df15<-TL_samples_gut(15,df,Mono_myc,mac_dfn, remove_ct=Macs_names)

# Keep only ratio in Infl Ileum, and only Mono1/Mono3
ratiom1m3_ibd<-data.frame(prop_m1=prop.table(table(Mono_df15@df3_nbr$sample,Mono_df15@df3_nbr$Group)[unique(Mono_df15@df3_nbr$sample[which(Mono_df15@df3_nbr$status=="Involved"&Mono_df15@df3_nbr$Disease=="CD")]),c("Mono1","Mono3")], margin=1)[,1],origin="nan_sin")

                                                        
ratiom1m3_ibd<-ratiom1m3_ibd[-which(rownames(ratiom1m3_ibd)=="158"),] # it has 0 m1, 0 m3


# Load info from Buckley et al 
df = read.csv("./Grouped_objects/momacs_cell_metadata.csv")
annots = read.csv("./Grouped_objects/Buckley_momacs_annot_241217.csv")
df = df[df$metacell_name != "Outliers",]
df$annot = unlist(lapply(X = df$metacell_name, FUN = function(x) {return(annots$annotation[annots$metacell == x])} ))
df = df[df$annot != "undefined" ,]   # Remove the "undefined" annotation


# Select only the patients of interest:
df_CD = df[df$Disease == "CD",]
df_CD_Infl = df_CD[df_CD$Inflammation == "Inflamed",]
df_CD_Infl_Pre = df_CD_Infl[df_CD_Infl$Treatment == "Pre",]
df_CD_Infl_Pre_Ileum = df_CD_Infl_Pre[df_CD_Infl_Pre$Site == "Terminal_Ileum",]

table(rowSums(table(df_CD_Infl_Pre_Ileum$sample_id,df_CD_Infl_Pre_Ileum$annot)[,5:9])>15)
ratiom1m3_tau<-data.frame(prop_m1=prop.table(table(df_CD_Infl_Pre_Ileum$sample_id,df_CD_Infl_Pre_Ileum$annot)[-11,c(5,7)],margin=1)[,1],origin="taurus") #CID006558-1 has less than 15 cells, remove it



ratiom1m3<-rbind(ratiom1m3_ibd,ratiom1m3_tau)

col_m1m3<-c("taurus"="#E29453","nan_sin"="#A83584")
comp_m1m3=list(c("nan_sin","taurus"))


##  /!\  9 of the 29 points are patients with 10 Mono1+Mono3 or less
# 
# ggplot(ratiom1m3, aes(x =origin, y =prop_m1, color =origin, fill=origin)) +
#   geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
#   scale_x_discrete(limits = names(col_m1m3) )+
#   geom_jitter(shape=16, position=position_jitter(0.2),size=4,alpha=0.8)+
#   ggprism::theme_prism(base_size = 20)+
#   theme(legend.position = "none",axis.title.x =element_blank())+
#   scale_color_manual(values=col_m1m3)+
#   scale_fill_manual(values=col_m1m3)+
#   stat_compare_means(comparisons = comp_m1m3,size = 6, bracket.size = rel(1.5),method = "wilcox",paired=F)+ 
#   scale_y_continuous(labels = scales::percent)


# Removing samples with <10 cells available for the ratios :
tibd<-table(Mono_df15@df3_nbr$sample,Mono_df15@df3_nbr$Group)[unique(Mono_df15@df3_nbr$sample[which(Mono_df15@df3_nbr$status=="Involved"&Mono_df15@df3_nbr$Disease=="CD")]),c("Mono1","Mono3")]
gd_tibd<-tibd[which(rowSums(tibd)>9),]
gd_ratiom1m3_ibd<-data.frame(prop_m1=prop.table(gd_tibd, margin=1)[,1],origin="nan_sin")

ttau<-table(df_CD_Infl_Pre_Ileum$sample_id,df_CD_Infl_Pre_Ileum$annot)[-11,c(5,7)]
gd_ttau<-ttau[which(rowSums(ttau)>9),]
gd_ratiom1m3_tau<-data.frame(prop_m1=prop.table(gd_ttau,margin=1)[,1],origin="taurus") 

gd_ratiom1m3<-rbind(gd_ratiom1m3_ibd,gd_ratiom1m3_tau)

pdf("./Figure 2/Fig2B_Bxplt_M1prop.pdf",width =5,height = 7.5)
ggplot(gd_ratiom1m3, aes(x =origin, y =prop_m1, color =origin, fill=origin)) +
  geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
  scale_x_discrete(limits = names(col_m1m3) )+
  geom_jitter(shape=16, position=position_jitter(0.2),size=4,alpha=0.8)+
  ggprism::theme_prism(base_size = 20)+
  theme(legend.position = "none",axis.title.x =element_blank())+
  scale_color_manual(values=col_m1m3)+
  scale_fill_manual(values=col_m1m3)+
  stat_compare_means(comparisons = comp_m1m3,size = 6, bracket.size = rel(1.5),method = "wilcox",paired=F)+ 
  scale_y_continuous(labels = scales::percent)
dev.off()
