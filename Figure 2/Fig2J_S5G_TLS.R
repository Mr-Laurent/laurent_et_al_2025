library(ggplot2)
library(ggpubr)
library(ggrastr)
setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")


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






for(spl in c("GIM23_InfROI2","GIM33_InfROI1","GIM38_InfROI1","GIM38_InfROI2")){
  

  
  
load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))

set.seed(42)

  
  
  

# > ct_names 
# [1] "Act_DC"         "addEosinoBaso"  "addGoblet"      "addNKlike"      "Bcell"          "CD4 Tcell"      "CD8 Tcell"      "DC1"            "doublet"       
# [10] "Endo"           "Epith"          "Fibro"          "Glial"          "Infl. Fib"      "Infl. FRC-like" "Lymphat"        "Macs"           "Mast"          
# [19] "MNP_N"          "Mono"           "Mono1"          "Neutro"         "pDC"            "PlasmaC"        "Smooth Muscle"  "TILC"           "unannotated"   


meta_annot_perso$my_annots2 <-meta_annot_perso$my_annots
meta_annot_perso$my_annots2 <- ifelse(meta_annot_perso$my_annots2 %in% c("addNKlike","CD4 Tcell","CD4 Tph","CD8 eff.","CD8 Tcell",
                                                                         "GClike B cell","Mem. B cell","Plasma Cell","Naive B cell"), meta_annot_perso$my_annots2, "others")

ct_names <- names(table(meta_annot_perso$my_annots2))
lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
lin_palette3["Naive B cell"] <-"#FFCC56"
lin_palette3["GClike B cell"] <- "#DDA493"
lin_palette3["Mem. B cell"] <- "#AB7453"
lin_palette3["Plasma Cell"] <- "#E5756F"
lin_palette3["addNKlike"] <- "#2B9542" 
lin_palette3["CD8 Tcell"] <- "#699EEB"
lin_palette3["CD8 eff."] <- "#1B458E"
lin_palette3["CD4 Tcell"] <- "#E590BC"
lin_palette3["CD4 Tph"] <- "#E5408C"



levels(lin_palette3)<-names(lin_palette3)

pdf(file=paste0("./Figure 2/Fig2J_TLS_",spl,".pdf"),width = 13,height = 10)
brplt_title<-paste0(spl) 
print( makeplot(meta_annot_perso) )

meta_annot_perso$my_annots2 <-meta_annot_perso$my_annots
meta_annot_perso$my_annots2 <- ifelse(meta_annot_perso$my_annots2 %in% c("CD4 Tph","Plasma Cell","GClike B cell","Naive B cell"), meta_annot_perso$my_annots2, "others")
ct_names <- names(table(meta_annot_perso$my_annots2))
lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
lin_palette3["Naive B cell"] <-"#DDCC77"
lin_palette3["GClike B cell"] <- "#5FB2E7"
lin_palette3["Plasma Cell"] <- "#E5855F"
lin_palette3["CD4 Tph"] <- "#26547C"
levels(lin_palette3)<-c("Naive B cell","GClike B cell","Plasma Cell","CD4 Tph","others")
print( makeplot(meta_annot_perso) )

for(ct in c("CD4 Tph","Plasma Cell","GClike B cell","Naive B cell")){
  meta_annot_perso$my_annots2 <-meta_annot_perso$my_annots
  meta_annot_perso$my_annots2 <- ifelse(meta_annot_perso$my_annots2 %in% ct , meta_annot_perso$my_annots2, "others")
  ct_names <- names(table(meta_annot_perso$my_annots2))
  
  if(spl=="GIM23_InfROI2"){ #Zoom selection by slide
    sub_map<-meta_annot_perso[which(meta_annot_perso$x>24000&meta_annot_perso$x<25500&
                                      meta_annot_perso$y>11000&meta_annot_perso$y<12500),]
  }else if(spl=="GIM33_InfROI1"){ #Zoom selection by slide
    sub_map<-meta_annot_perso[which(meta_annot_perso$x>2500&meta_annot_perso$x<4250&
                                      meta_annot_perso$y>14500&meta_annot_perso$y<16500),]
  }else if(spl=="GIM38_InfROI2"){ #Zoom selection by slide
    sub_map<-meta_annot_perso[which(meta_annot_perso$x>25000&meta_annot_perso$x<27000&
                                      meta_annot_perso$y>3000&meta_annot_perso$y<5000),]
  }else{sub_map<-meta_annot_perso}
  
  lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  lin_palette3["Naive B cell"] <-"#FFD166"
  lin_palette3["Plasma Cell"] <- "#E5855F"
  lin_palette3["GClike B cell"] <- "#5FB2E7"
  lin_palette3["CD4 Tph"] <- "#26547C"
  levels(lin_palette3)<-c(ct,"others")
  print( makeplot(sub_map,smalldot_size=1,bigdot_size=3) )
}



dev.off()
}