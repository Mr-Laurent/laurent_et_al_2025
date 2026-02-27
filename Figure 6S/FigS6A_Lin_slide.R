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



# Plot the main lineages together on each slides
for(spl in c("GIM23_NonInfl","GIM23_InfROI1","GIM23_InfROI2","GIM33_InfROI1","GIM33_InfROI2","GIM38_NonInfl","GIM38_InfROI1","GIM38_InfROI2")){

  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
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
                          
  pdf(file=paste0("./Figures_print/FigS6A_Lin_",spl,".pdf"),width = 13,height = 10)
  brplt_title<-paste0(spl) 
  print( makeplot(meta_annot_perso) )
  
  dev.off()
}