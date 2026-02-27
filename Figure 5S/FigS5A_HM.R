# Refaire Fig S5A : peut etre ht momacs plutot que myelo ? Sur Buckley
library(ComplexHeatmap)
library(ggplot2)
library(Rfast)
library(ggrepel)
library(ggpubr)
library(scales)



dfinfo<-read.csv2(file="./Grouped_objects/Buckley_momacs_annot_241217.csv",sep=",",header=T)


object="100mods_momacs_ThomasBuckley"
load(paste0("./Grouped_objects/Taurus/genelist_",object,".rd"))
load(paste0("./Grouped_objects/Taurus/metadata_",object,".rd"))
load(paste0("./Grouped_objects/Taurus/counts_",object,".rd"))
load(paste0("./Grouped_objects/Taurus/ht_",object,".rd"))

load("./Grouped_objects/Modsig_MoMacsEnr.rd")

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




df_zmod<-metadf_z[,154:length(colnames(metadf_z))]
#Attribute cluster number to each metacell
Macs_myc<-ht_clustering(ht,df_zmod)


ds<-as.matrix(dgcmtxcounts)
allsum<-Rfast::colsums(log(1+ds))

for(Mod_score in names(Modsig)[1:17] ){
  print(Mod_score)
  geneshow<-unlist(Modsig[Mod_score])
  print(setdiff(geneshow,rownames(ds)) ) # Only keeps genes detected in the count matrix
  ilist<-intersect(geneshow,rownames(ds))
  if(length(ilist)<2){currentscore<-log(1+ds[ilist,])/allsum 
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }else{currentscore<-(Rfast::colsums(log(1+ds[ilist,,drop=F]))/allsum )
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }
}




df_zmod<-metadf_z[,(ncol(metadf_z)-16):(ncol(metadf_z))]
df_zmod$cluster<-Macs_myc@out$Cluster[match(rownames(df_zmod),Macs_myc@out$Metacell)]
df_zmod$cluster<-factor(df_zmod$cluster)

meddf <- aggregate(df_zmod[,1:(length(colnames(df_zmod))-1)], by = list(df_zmod$cluster), median)
minimelt<-reshape::melt(meddf)  
#Reorder the plot manually based on clustering + similarities with Fig 1B:
clnmord3<-paste0("cluster",c(39,40,28,29,31,37,30,18,16,15,35,38,34,36,32,33,14,3,4,7,8,5,27,13,20,12,6,21,9,17,11,26,10,1,23,24,22,2,25) )  #19,
 
ord_sig_perso2<-c("Sig_X","Sig_XI","Sig_XII","Sig_I","Sig_III","Sig_IV","Sig_VI",
                  "Sig_II","Sig_XVI","Sig_XVII","Sig_XV","Sig_XIV","Sig_XIII","Sig_IX","Sig_VIII","Sig_VII","Sig_V")



plot_reorderselectedrow_selectedcol<-ggplot(minimelt, aes(x = variable, y = Group.1, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish, name = "Median Score") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "right")+ geom_hline(yintercept = 0.5 + 0:40, colour = "black", size = 0.05) 


pdf("./Figures_print/FigS5A_HM.pdf",height =10,width = 6.8)
plot_reorderselectedrow_selectedcol+
  scale_x_discrete(limits = ord_sig_perso2,position = "top", labels=gsub("Sig_","",ord_sig_perso2))+
  scale_y_discrete(limits = (clnmord3),                   labels=gsub("cluster","",(clnmord3)) )+
  theme(axis.text.y =element_text(size = rel(1.5),face="bold"),
        axis.text.x =element_text(hjust=0.5,angle = 0,size = rel(1.25),face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()  )+coord_flip()+
  ggtitle("Relative enrichment of signatures by clusters \n(Macs from total samples, manual reordering)")
dev.off()






