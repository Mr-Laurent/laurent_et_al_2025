# Scatter plot of metacells distributed by Program VII (MonoA) and Program VIII (MonoB) with each subtype highlighted
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
# dfinfo has chemistry, sample, disease, location, etc informations
dfinfo<-read.csv2(file="./Grouped_objects/sample_annots.csv",sep=",",header=T)
colnames(dfinfo)[1]<-"sample"
listobj<-c("df3hcco","df3hcil","df3CDuni","df3CD","df3UC")
load("./Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v3.rd") # metadf_z     dataframe with module scores and scaled scores by metacells
df_zmod<-metadf_z[,154:length(colnames(metadf_z))]
load("./Grouped_objects/ht_Enrich_Macs.rd")
load("./Grouped_objects/df_mcall_50clust_16nov22_Macs.rd")

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

# Annotations of MoMacs in the total cohort:
mac_dfn<-data.frame(Grp_cell=c("Mono-like 1", "Mono-like 2", "Mono-like 3", "Mono-like 4", "Mono-like 5", 
                              "Macrophage-like 6",  "Macrophage-like 7", "Macrophage-like 8","Macrophage-like 9")
                    ,Grp_cls=c("2,40,6,34,35,36,37,38,39,33",   "3,9,11,25,4,13", "7,8,5,10",  "24,26",  "27,28",
                               "12,14", "1,16,18", "19,20,23,29,30,32", "15,17,21,22,31") )



Momacs_myc<-ht_clustering(ht,df_zmod)


# df is the object with the info Cell ID to metacell
df2<-df[which(df$mc>=0),]
df2$mc2<-paste0("mc",df2$mc)

# put the corresponding cluster to each metacell
df2$Cluster<-Momacs_myc@out$Cluster[match(df2$mc2,Momacs_myc@out$Metacell)]
df2$Group<-NA

for(i in mac_dfn$Grp_cell){
  df2$Group[which(df2$Cluster%in%paste0('cluster',unlist(strsplit(mac_dfn$Grp_cls[which(mac_dfn$Grp_cell==i)],split=','))) )]<-i
  
}


pdf("./Figure 2S/Paper_FigS2H_Sctplt_MoMac.pdf", width = 6, height = 6)
for (i in 1:length(mac_dfn$Grp_cell)) {
  current_group <- mac_dfn$Grp_cell[i]
  current_df <- metadf_z[as.vector(df2$mc2[df2$Group == current_group]), ]
  other_df <- metadf_z[as.vector(df2$mc2[df2$Group != current_group]), ]
  print( ggplot(other_df, aes(x =  SigMac_Mono2 , y =  SigMac_Mono1 )) + geom_point( color = "grey")+
           geom_point(data = current_df, aes(x = SigMac_Mono2 , y =  SigMac_Mono1), color = "red") + 
           ggtitle(paste("Scatter Plot for Group", current_group))+
           theme_bw()+
           theme(axis.line = element_line(colour = "black",size = 0.6)) ) 
}
dev.off()
