library(Matrix)
library(ggplot2)
library(reticulate)
library(dplyr)
library(scales)
library(patchwork)
library(ggpubr)
library(ggrastr)
`%ni%` <- Negate(`%in%`)


setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")


#Figures
for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  
  ### Slides panels:
  # 1) Visualize Mono vs Macs
  # 2) Visualize Mono1, Mono3, Mono5
  # 3) Visualize Mono1, Mono3, Mono5, and other non annotated "mono"-like and Mono in transition
  
  # 1)
  ct=c("Macs","Mono","Mono1","Mono3","Mono5","Momacs_transition")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["Macs"] <- "#3F7FD9"
  lin_palette3["Mono"] <- "#CF6910"
  lin_palette3["Mono1"] <- "#CF6910"
  lin_palette3["Mono3"] <- "#CF6910"
  lin_palette3["Mono5"] <- "#CF6910"
  lin_palette3["Momacs_transition"] <- "#CF6910"
  meta_annot_perso$current_ct<-NA
  meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot%in%ct)]<-meta_annot_perso$merge_annot[which(meta_annot_perso$merge_annot%in%ct)]
  
  
  pdf(file=paste0("./Figure 5S/FigS5B_Xen_Lin_myannot_Mono_",spl,"_v5bscaleonall.pdf"),width = 13,height = 10)
  levels(lin_palette3)<-names(lin_palette3)
  brplt_title<-paste0(spl," -- projection of Mono vs Macs" )
  print(ggplot(meta_annot_perso, aes(x = x, y = y, color = merge_annot)) +
          rasterise(geom_point(data = subset(meta_annot_perso, is.na(current_ct) ), aes(x = x, y = y),
                               color = "#EDEDED", size = 0.2)) +
          geom_point(data = subset(meta_annot_perso, current_ct %in% names(lin_palette3)),
                     aes(x=x,y=y), color =  lin_palette3[subset(meta_annot_perso, current_ct %in% names(lin_palette3))$current_ct], size = 0.8) +
          scale_color_manual(values = lin_palette3)+
          guides(color = guide_legend(override.aes = list(size = 4)))+
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = "Lineage" )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())+ coord_fixed() ) 
  # 2)
  ct=c("Mono1","Mono3","Mono5")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["Mono1"] <- "#874037"
  lin_palette3["Mono3"] <- "#F77774"
  lin_palette3["Mono5"] <- "#50AE86"
  levels(lin_palette3)<-names(lin_palette3)
  brplt_title<-paste0(spl," -- projection of main subtypes in Mono" )
  print(ggplot(meta_annot_perso, aes(x = x, y = y, color = merge_annot)) +
          rasterise(geom_point(data = subset(meta_annot_perso, is.na(current_ct) ), aes(x = x, y = y),
                               color = "#EDEDED", size = 0.2)) +
          geom_point(data = subset(meta_annot_perso, current_ct %in% names(lin_palette3)),
                     aes(x=x,y=y), color =  lin_palette3[subset(meta_annot_perso, current_ct %in% names(lin_palette3))$current_ct], size = 0.8) +
          scale_color_manual(values = lin_palette3)+
          guides(color = guide_legend(override.aes = list(size = 4)))+
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = "Lineage" )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())+ coord_fixed() )
  # 3)
  ct=c("Mono","Mono1","Mono3","Mono5","Momacs_transition")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["Mono"] <- "#CF6910"
  lin_palette3["Mono1"] <- "#874037"
  lin_palette3["Mono3"] <- "#F77774"
  lin_palette3["Mono5"] <- "#50AE86"
  lin_palette3["Momacs_transition"] <- "#458875"
  print(ggplot(meta_annot_perso, aes(x = x, y = y, color = merge_annot)) +
          rasterise(geom_point(data = subset(meta_annot_perso, is.na(current_ct) ), aes(x = x, y = y),
                               color = "#EDEDED", size = 0.2)) +
          geom_point(data = subset(meta_annot_perso, current_ct %in% names(lin_palette3)),
                     aes(x=x,y=y), color =  lin_palette3[subset(meta_annot_perso, current_ct %in% names(lin_palette3))$current_ct], size = 0.8) +
          scale_color_manual(values = lin_palette3)+
          guides(color = guide_legend(override.aes = list(size = 4)))+
          theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          labs(title = brplt_title,
               color = "Lineage" )+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())+ coord_fixed() )
  
  
  dev.off()
}

