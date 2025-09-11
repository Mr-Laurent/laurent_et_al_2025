library(ComplexHeatmap)
library(ggplot2)
library(Rfast)
library(ggrepel)
library(ggpubr)
library(scales)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
dfinfo<-read.csv2(file="./Grouped_objects/sample_annots.csv",sep=",",header=T)
colnames(dfinfo)[1]<-"sample"
listobj<-c("df3hcco","df3hcil","df3CDuni","df3CD","df3UC")

load("./Grouped_objects/MacsTot_z_mymod_150_metadata_goodlog_v3.rd") # metadf_z     dataframe with module scores and scaled scores by metacells
load("./Grouped_objects/genelist_16nov22_Macs.rd")
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
#### Function to split genes in signature lists ####
###----------------------------------------------###

automate_strsplit <- function(df, col, indices) {
  unlist(strsplit(paste0(df[[col]][indices], collapse = ","), split = ","))
}

###----------------------------------------------### 
#### Function to order with hierarc. clustering ####
###----------------------------------------------###

ordering<-function(data=NULL){
  res2<-cor(as.matrix(t(data)), method = c("pearson"))
  res2_dist=parallelDist::parDist(res2,method = 'euclidean')
  res2_clust=hclust(res2_dist,method = 'complete')
  return(res2_clust$order)}

###----------------------------------------------### 


#Define which modules are used in each programs
grp_I<-c(43,128,70,38)
grp_II<-c(8,13)
grp_III<-c(131,133,95,40,66)
grp_IV<-c(99,25,33)
grp_V<-c(55,6,50,21,31)
grp_VI<-c(16,81,58,119)
grp_VII<-c(143,27,12,14)
grp_VIII<-c(19,71,123,18,144)
grp_IX<-c(28,83,48,103,49,141,53,126)
grp_X<-c(140,72,94,106)
grp_XI<-c(4,10,73)
grp_XII<-c(93,111,90,107)
grp_XIII<-c(11,7,65,74,102)
grp_XIV<-c(20,52,54,34,113,108,37,60,142,121,116,87,120)
grp_XV<-c(46,77,125,2,3,44,56)
grp_XVI<-c(39)
grp_XVII<-c(84,88)
Modused<-c()
Sig_list<-list()
for(i in ls(pattern = "^grp_")){
  eval(parse(text=paste0("Sig_list$Sig",gsub("_ids$","",i),"<- automate_strsplit(genelist, 'Gens', ",i,")" )))
  Modused<-append(Modused,paste0("Mods ",paste0(eval(as.name(i)),collapse = ",")  ) )
}

Modsig<-data.frame(signatures=names(Sig_list),modused=Modused)


df_zmod<-metadf_z[,154:length(colnames(metadf_z))]
#Attribute cluster number to each metacell
Macs_myc<-ht_clustering(ht,df_zmod)



load("./Grouped_objects/MacsfromMNP_select2_dgCmc_counts.rd")

for(Mod_score in Modsig$signatures){
  geneshow<-eval(parse(text=paste0("Sig_list$",Mod_score)))
  gnshow<-unlist(strsplit(geneshow,split=","))
  gnshow<-gsub(" ","",gnshow)
  ilist<-intersect(gnshow,rownames(dgcmtxcounts))
  currentscore<-Matrix::colSums(log(1+dgcmtxcounts[ilist,,drop=F]))/Matrix::colSums(log(1+dgcmtxcounts))
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
}
df_zmod<-metadf_z[,(ncol(metadf_z)-16):(ncol(metadf_z))]
df_zmod$cluster<-Macs_myc@out$Cluster[match(rownames(df_zmod),Macs_myc@out$Metacell)]
df_zmod$cluster<-factor(df_zmod$cluster)

meddf <- aggregate(df_zmod[,1:(length(colnames(df_zmod))-1)], by = list(df_zmod$cluster), median)
minimelt<-reshape::melt(meddf)  
#Reorder the plot manually based on clustering + similarities with Fig 1B:
clnmord3<-paste0("cluster",c(38,37,39,31,40,26,32,34,35,36,22,3,20,21,7,6,29,28,5,11,4,25,33,1,30,27,24,23,15,2,14,9,10,8,13,12,16,19,17,18) )
ord_sig_perso2<-c("Siggrp_X","Siggrp_XI","Siggrp_XII","Siggrp_I","Siggrp_III","Siggrp_IV","Siggrp_VI",
                  "Siggrp_II","Siggrp_XVI","Siggrp_XVII","Siggrp_XV","Siggrp_XIV","Siggrp_XIII","Siggrp_IX","Siggrp_VIII","Siggrp_VII","Siggrp_V")



plot_reorderselectedrow_selectedcol<-ggplot(minimelt, aes(x = variable, y = Group.1, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish, name = "Median Score") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "right")+ geom_hline(yintercept = 0.5 + 0:40, colour = "black", size = 0.05) 


pdf("./Figure 2S/FigS2A_HM.pdf",height =10,width = 6.8)

plot_reorderselectedrow_selectedcol+
  scale_x_discrete(limits = ord_sig_perso2,position = "top", labels=gsub("Siggrp_","",ord_sig_perso2))+
  scale_y_discrete(limits = rev(clnmord3),                   labels=gsub("cluster","",rev(clnmord3)) )+
  theme(axis.text.y =element_text(size = rel(1.5),face="bold"),
        axis.text.x =element_text(hjust=0.5,angle = 0,size = rel(1.25),face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()  )+
  ggtitle("Relative enrichment of signatures by clusters \n(Macs from total samples, manual reordering)")


dev.off()
