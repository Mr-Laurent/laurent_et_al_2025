library(Matrix)
library(ComplexHeatmap)
library(reshape)
library(ggplot2)
library(scales)
library(Rfast)


load("./Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v3.rd") # metadf_z     dataframe with module scores and scaled scores by metacells
hm_zmod<-metadf_z[,154:303]
load("./Grouped_objects/df_mcall_50clust_16nov22_Macs.rd")  # df contains the metadatas: cell ID, metacell ID, sample of origin, cluster and lab of origin
load("./Grouped_objects/Macs_dgCmc_16nov22_counts.rd")      # dgcmtxcounts is the dgCMatrix (condensed) with counts by metacells
load("./Grouped_objects/genelist_16nov22_Macs.rd")          # genelist has the gene lists associated to each of the 150 modules 
load("./Grouped_objects/ht_Enrich_Macs.rd")                 # ht has the cluster + metacells order defined with ComplexHeatmap 
colnames(hm_zmod)<-gsub("scaled","",colnames(hm_zmod))


###-------------------------------------------------------------### 
#### Function to assign vector of splitted genes from genelist #### 
###-------------------------------------------------------------### 
automate_strsplit <- function(df, col, indices) {
  unlist(strsplit(paste0(df[[col]][indices], collapse = ","), split = ","))
}
###-------------------------------------------------------------### 

# Our programs are a combination of modules selected as follows:
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


Signatures<-c()
Modused<-c()
for(i in ls(pattern = "^grp_")){
  eval(parse(text=paste0("Sig",gsub("_ids$","",i),"<- automate_strsplit(genelist, 'Gens', ",i,")" )))
  Signatures<-append(Signatures,paste0("Sig",gsub("_ids$","",i)) )
  Modused<-append(Modused,paste0("Mods ",paste0(eval(as.name(i)),collapse = ",")  ) )
}

Modsig<-data.frame(signatures=Signatures,modused=Modused)


metadf_z[,304:(303+length(colnames(metadf_z)))]<-NULL

ds<-as.matrix(dgcmtxcounts) 
allsum<-Rfast::colsums(log(1+ds)) 

## Compute program enrichment scores by metacells
for(Mod_score in Modsig$signatures){
  print(Mod_score)
  geneshow<-eval(parse(text=Mod_score))
  ilist<-intersect(geneshow,rownames(ds))
  if(length(ilist)<2){currentscore<-(log(1+ds[intersect(ilist,rownames(ds)),])/allsum )
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }else{currentscore<-(Rfast::colsums(log(1+ds[ilist,]))/allsum )
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }}
df_zmod<-metadf_z[,304:(303+length(Modsig$signatures))]


## Attribute each metacell to their cluster
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

df_zmod$cluster<-out$Cluster 
df_zmod$cluster<-factor(df_zmod$cluster, levels = paste0("cluster", 1:length(metacells_kclust)))

# The final heatmap show the median enrichment score of each program in each of the 40 metacell clusters
meddf <- aggregate(df_zmod[,1:(length(colnames(df_zmod))-1)], by = list(df_zmod$cluster), median)
minimelt<-reshape::melt(meddf)   

ordclu<-paste0("cluster",c(1:40))

# Heatmap function 
plot_reorderselectedrow_selectedcol<-ggplot(minimelt, aes(x = variable, y = Group.1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish, name = "Median Score") + # oob=squish : out of bound keep the limit color
  scale_y_discrete(limits = ordclu) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")+ geom_hline(yintercept = 0.5 + 0:40, colour = "black", size = 0.05) # adds lines between clusters


clnmord2<-paste0("cluster",c(6,40,2,33,39,38,37,34,35,36,3,9,11,25,4,13,7,8,5,10,24,26,27,28,12,14,16,18,1,19,20,23,29,30,32,15,17,21,22,31) )
ord_sig_perso2<-c("Siggrp_X","Siggrp_XI","Siggrp_XII","Siggrp_I","Siggrp_III","Siggrp_IV","Siggrp_VI",
                  "Siggrp_II","Siggrp_XVI","Siggrp_XVII","Siggrp_XV","Siggrp_XIV","Siggrp_XIII","Siggrp_IX","Siggrp_VIII","Siggrp_VII","Siggrp_V")



# Print the figure
lab_alt1<-c("Program X","Program XI","XII   Immunoreg./Repair","I      Transition A",
            "Program III","IV    Transition B","VI    Transition C","II     OXPHOS",
            "XVI  HSF","XVII Hypoxia","XV   IFIM program","XIV  Inflammatory B",
            "XIII  Inflammatory A","IX    IFN","VIII  MonoB","VII   MonoA","V     Macrophage")
pdf("./Figures_print/Fig1B_HM.pdf",height =6.5,width = 8)
plot_reorderselectedrow_selectedcol+ coord_flip()+
  scale_x_discrete(limits = ord_sig_perso2,position = "top", labels=lab_alt1 )+
  scale_y_discrete(limits = (clnmord2),                   labels=gsub("cluster","",clnmord2) )+
  theme(axis.text.y =element_text(size = rel(1.5),face="bold", hjust=0.5),
        axis.text.x =element_text(hjust=1,vjust=0.4,angle = 90,size = rel(0.75),face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()  )+ 
  ggtitle("Relative enrichment of signatures by clusters \n(Macs from enriched samples, manual reordering)")
dev.off()



