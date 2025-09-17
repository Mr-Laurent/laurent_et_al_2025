googlepath="G:/Mon Drive/UG_metacells/Xenium/"
setwd(paste0(googlepath))

library(Matrix)
library(ggplot2)
library(reticulate)
library(dplyr)
library(scales)
library(viridis)
library(shiny)
library(patchwork)
library(ggpubr)
library(ggrastr)
`%ni%` <- Negate(`%in%`)

#------------------------------------------------------#

makeplot <- function(meta_annot_perso,smalldot_size=0.1,bigdot_size=0.5){
  plot<-ggplot(meta_annot_perso, aes(x = x, y = y, color = my_annots2)) +
    rasterise(geom_point(data = subset(meta_annot_perso, my_annots2 == "others"), aes(size = smalldot_size))) +
    geom_point(data = subset(meta_annot_perso, my_annots2 != "others"), aes(size = bigdot_size )) +
    scale_color_manual(values = c(lin_palette3),breaks=levels(lin_palette3))+
    guides(color = guide_legend(override.aes = list(size = 4)))+
    theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(title = brplt_title,
         color = "Lineage" )+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_size_identity()+coord_fixed()
  return(plot)
}      

#------------------------------------------------------#






for(spl in c("GIM23_NonInfl","GIM23_InfROI1","GIM23_InfROI2","GIM33_InfROI1","GIM33_InfROI2","GIM38_NonInfl","GIM38_InfROI1","GIM38_InfROI2")){

  
  load(file=paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("G:/Mon Drive/UG_metacells/Xenium/SubsetsAppObj/meta_annot_june5_",spl,".rd"))
  load(paste0(googlepath,"/SubsetsAppObj/meta_annot_",spl, ".rd"))
  set.seed(42)
  
  
  
  
  
  # unique(meta_annot_perso$my_annots)
  old.cluster.ids<-c("Submucosa Myofibro","FIB","B_PC","Epith","CD3D TILC","Endo",
                     "SM RERGL","DC1","Plasma Cell","Macs","Infl Fib","CD4 Tcell",
  "CD8 Tcell","DC2_3","GClike B cell","TILC","Mem. B cell","Sm. Muscle",
  "MNP_N","CD8 eff.","Glial","Lymphat","Naive B cell","unannotated",
  "Mast","addGoblet","CD4 Tph","Mono","Neutro","aDC",
  "iaDC","Mono5","Mono1","Mono3","addNKlike","pDC",
  "Momacs_transition","addEosinoBaso")     

  new.cluster.ids <- c("Fibroblasts","Fibroblasts","B cells / PC","Epithelial cells","T cells / ILC","Endothelial cells",
                       "Fibroblasts","MNP / Neutro","B cells / PC","MNP / Neutro","Fibroblasts","T cells / ILC",
                       "T cells / ILC","MNP / Neutro","B cells / PC","T cells / ILC","B cells / PC","Smooth Muscle",
                       "MNP / Neutro","T cells / ILC","others","others","B cells / PC","others",
                       "others","Epithelial cells","T cells / ILC","MNP / Neutro","MNP / Neutro","MNP / Neutro",
                       "MNP / Neutro","MNP / Neutro","MNP / Neutro","MNP / Neutro","T cells / ILC","MNP / Neutro",
                       "MNP / Neutro","others")
  
  
  
  # [1] "Fibroblasts"       "B cells / PC"      "Epithelial cells"  "T cells / ILC"     "Endothelial cells" "MNP / Neutro"     
  # [7] "Smooth Muscle"     "others"       
  meta_annot_perso$my_annots2<- plyr::mapvalues(x = meta_annot_perso$my_annots, from = old.cluster.ids, to = new.cluster.ids)
 
  
  ct_names <- names(table(meta_annot_perso$my_annots2))
  # lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  # lin_palette3["B cells / PC"] <-"#FFCC56"
  # lin_palette3["Epithelial cells"] <- "#60BEB0"
  # lin_palette3["MNP / Neutro"] <- "#E5755F"
  # lin_palette3["Endothelial cells"] <- "#2E7B6F" 
  # lin_palette3["T cells / ILC"] <- "#5B85DE"
  # lin_palette3["Fibroblasts"] <- "#E590BC"
  # lin_palette3["Smooth Muscle"] <- "#E5408C"
  # 
  # 
  # colorpal2<-c("CD4 T cells"="#1A6C92","MAIT"="#7873B4","CD8 T cells"="#74A8DB","NK proliferating"="#2E7B6F","NK cells"="#60BEB0","CD56bright NK"="#A6D7D6",
  #              "CD16 Monocytes"="#A9AB35","CD14 Monocytes"="#CACC65","cDC"="#F3A7BB","pDC"="#E18399","Naive B cells"="#EE8866","Memory B cells"="#EABC7E","others"="#BCBCBB")
  # 
  # lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  # lin_palette3["B cells / PC"] <-"#74A8DB"
  # lin_palette3["Epithelial cells"] <- "#A9AB35"
  # lin_palette3["MNP / Neutro"] <- "#E5755F"
  # lin_palette3["Endothelial cells"] <- "#2E7B6F" 
  # lin_palette3["T cells / ILC"] <- "#3B55AE"
  # lin_palette3["Fibroblasts"] <- "#F3A7BB"
  # lin_palette3["Smooth Muscle"] <- "#E1608C"
  
  lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  lin_palette3["B cells / PC"] <-"#FFCC56"
  lin_palette3["Epithelial cells"] <- "#893B45"
  lin_palette3["MNP / Neutro"] <- "#E5753F"
  lin_palette3["Endothelial cells"] <- "#60BEB0" 
  lin_palette3["T cells / ILC"] <- "#3B55AE"
  lin_palette3["Fibroblasts"] <- "#E590BC"
  lin_palette3["Smooth Muscle"] <- "#E5408C"
  
  
  levels(lin_palette3)<-c("MNP / Neutro","B cells / PC","T cells / ILC","Fibroblasts","Smooth Muscle",             
                          "Endothelial cells","Epithelial cells","others" )           
                           
                          
                 
  
  pdf(file=paste0("G:/Mon Drive/UG_metacells/Figures paper Aout/Figure 5S/FigS5A_Lin_",spl,".pdf"),width = 13,height = 10)
  brplt_title<-paste0(spl) 
  print( makeplot(meta_annot_perso) )
  
  # meta_annot_perso$my_annots2 <-meta_annot_perso$my_annots
  # meta_annot_perso$my_annots2 <- ifelse(meta_annot_perso$my_annots2 %in% c("CD4 Tph","Plasma Cell","GClike B cell","Naive B cell"), meta_annot_perso$my_annots2, "others")
  # ct_names <- names(table(meta_annot_perso$my_annots2))
  # lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  # lin_palette3["Naive B cell"] <-"#DDCC77"
  # lin_palette3["GClike B cell"] <- "#5FB2E7"
  # lin_palette3["Plasma Cell"] <- "#E5855F"
  # lin_palette3["CD4 Tph"] <- "#26547C"
  # levels(lin_palette3)<-c("Naive B cell","GClike B cell","Plasma Cell","CD4 Tph","others")
  # print( makeplot(meta_annot_perso) )
  # 

  
  
  
  dev.off()
}