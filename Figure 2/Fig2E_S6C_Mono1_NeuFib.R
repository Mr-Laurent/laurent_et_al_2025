library(Matrix)
library(ggplot2)
library(reticulate)
library(dplyr)
library(scales)
library(patchwork)
library(ggpubr)
library(ggrastr)
`%ni%` <- Negate(`%in%`)


pdf(file=paste0("./Figures_print/Fig2E_S6C_Mono1_Neutro_v5bscaleonall.pdf"),width = 13,height = 10)

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  
  brplt_title<-paste0(spl," -- projection of Mono1 vs Neutro" )
  
  ct=c("Mono1","Neutro")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["Neutro"] <- "#7484EC"
  lin_palette3["Mono1"] <- "#874037"
  levels(lin_palette3)<-names(lin_palette3)
  meta_annot_perso$current_ct<-NA
  meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot%in%ct)]<-meta_annot_perso$merge_annot[which(meta_annot_perso$merge_annot%in%ct)]
  
  
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
  
}
dev.off()


pdf(file=paste0("./Figures_print/Fig2E_S6C_Mono1_InflFibro_v5bscaleonall.pdf"),width = 13,height = 10)

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  
  brplt_title<-paste0(spl," -- projection of Mono1 vs Infl Fib" )
  
  ct=c("Mono1","Infl Fib")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["Infl Fib"] <- "#C967EA"
  lin_palette3["Mono1"] <- "#874037"
  levels(lin_palette3)<-names(lin_palette3)
  meta_annot_perso$current_ct<-NA
  meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot%in%ct)]<-meta_annot_perso$merge_annot[which(meta_annot_perso$merge_annot%in%ct)]
  
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
  
}
dev.off()

pdf(file=paste0("./Figures_print/Fig2E_S6C_iaDC_then_EffCD8_v5bscaleonall.pdf"),width = 13,height = 10)

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  
  brplt_title<-paste0(spl," -- projection of iaDC" )
  
  ct=c("iaDC")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["iaDC"] <- "#642624"

  levels(lin_palette3)<-names(lin_palette3)
  meta_annot_perso$current_ct<-NA
  meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot%in%ct)]<-meta_annot_perso$merge_annot[which(meta_annot_perso$merge_annot%in%ct)]
  
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
  brplt_title<-paste0(spl," -- projection of Effector T CD8" )
  
  ct=c("CD8 eff.")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["CD8 eff."] <- "#133E8A"
  
  levels(lin_palette3)<-names(lin_palette3)
  meta_annot_perso$current_ct<-NA
  meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot%in%ct)]<-meta_annot_perso$merge_annot[which(meta_annot_perso$merge_annot%in%ct)]
  
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
  
}
dev.off()



pdf(file=paste0("./Figures_print/FigS6F_NK_v5bscaleonall.pdf"),width = 13,height = 10)

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  
  brplt_title<-paste0(spl," -- projection of NK cells" )
  
  ct=c("addNKlike")
  lin_palette3 <- setNames(rep("#EDEDED", length(ct)), ct)
  lin_palette3["addNKlike"] <- "#20B138"
  
  levels(lin_palette3)<-names(lin_palette3)
  meta_annot_perso$current_ct<-NA
  meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot%in%ct)]<-meta_annot_perso$merge_annot[which(meta_annot_perso$merge_annot%in%ct)]
  
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
  
}
dev.off()

