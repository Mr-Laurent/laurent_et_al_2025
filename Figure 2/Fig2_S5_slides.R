# Sub annotations, based on single lineages across samples

#####################################
#######  Subtypes 
#####################################

library(Matrix)
library(ggplot2)
library(reticulate)
library(dplyr)
library(scales)
library(patchwork)
library(ggpubr)
library(ggrastr)
`%ni%` <- Negate(`%in%`)


setwd("G:/Mon Drive/UG_metacells/laurent_et_al_2025/")

grp_dmrc<-function(sel_lin){
  spl="GIM23_InfROI1"
  load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_my_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
  bigmap<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
  big_dmrc<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%bigmap$cell),]
  rownames(big_dmrc)<-rownames(bigmap)
  
  
  for(spl in c("GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
    print(spl)
    load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
    load(paste0("./Grouped_objects/Xenium/meta_annot_my_",spl,".rd"))
    load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
    set.seed(42)
    
    sub_ma_mnpn<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
    sub_dmrc_mnpn<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%sub_ma_mnpn$cell),colnames(big_dmrc)] #colnames to ensure the gene order is the same before merging
    rownames(sub_dmrc_mnpn)<-rownames(sub_ma_mnpn)
    
    bigmap<-rbind(bigmap,sub_ma_mnpn)
    big_dmrc<-rbind(big_dmrc,sub_dmrc_mnpn)
  }  
  return( list(bigmap=bigmap, big_dmrc=big_dmrc))
}

### MNP:
sel_lin=c("MNPNeutro","MNPNeutroInfl")
MNPN_obj<-grp_dmrc(sel_lin)
bigmap<-MNPN_obj$bigmap
big_dmrc<-MNPN_obj$big_dmrc


sig_Macs<-unlist(strsplit("FOLR2,IGF1,KCNMA1,C1QA,AXL,VSIG4,DNASE1L3,GATM,SIGLEC1",split=','))
sig_Mono<-unlist(strsplit("CLEC4E,RETN,FCN1,LILRA5,MARCO,MCEMP1,FCAR,CD93,ANPEP,APOBEC3A,CCR1,VCAN,CD300E",split=','))
sig_Neutro<-unlist(strsplit("CXCR2,MME",split=','))
sig_pDC<-unlist(strsplit("GZMB,CXCR3,LILRA4,TCL1A",split=','))
sig_aDC<-unlist(strsplit("CCR7,LAD1,FSCN1,LAMP3",split=','))
sig_iaDC<-unlist(strsplit("CXCL9,CXCL10,C15orf48,SRC",split=','))
sig_DC1<-unlist(strsplit("MYLK,CADM1,PPY",split=','))
sig_DC2_3<-unlist(strsplit("CD1C,FCER1A,CLEC10A",split=','))

allsum<-Matrix::rowSums(big_dmrc) 
bigmap$cell_id<-rownames(bigmap)
for(i in c("sig_Macs","sig_Mono","sig_Neutro",
           "sig_pDC","sig_aDC","sig_iaDC","sig_DC1","sig_DC2_3")){
  
  eval(parse(text=paste0("ilist<-",i)))

  bigmap$currentscore<-NA
  # Add all gene counts from the signature
  df_inter<-data.frame("cell_name"=bigmap$cell_id,
                       "score"=(Matrix::rowSums(big_dmrc[,ilist])/allsum) )
  df_inter$scale_score<-scale(df_inter$score)
  eval(parse(text=paste0("bigmap$",i,"<-as.vector(df_inter$scale_score)")))
}
dim(bigmap) # 206467 cells in c("MNPNeutro","MNPNeutroInfl")



for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){

  # First annotation round:

  bigmap$pos_Macs<-ifelse(bigmap$sig_Macs>=0.5,TRUE,FALSE)
  bigmap$pos_Mono<-ifelse(bigmap$sig_Mono>=0.5,TRUE,FALSE)
  bigmap$pos_Neutro<-ifelse(bigmap$sig_Neutro>=0.5,TRUE,FALSE)
  bigmap$pos_pDC<-ifelse(bigmap$sig_pDC>=0.5,TRUE,FALSE)
  bigmap$pos_aDC<-ifelse(bigmap$sig_aDC>=0.5,TRUE,FALSE)
  bigmap$pos_iaDC<-ifelse(bigmap$sig_iaDC>=0.5,TRUE,FALSE)
  bigmap$pos_DC1<-ifelse(bigmap$sig_DC1>=0.5,TRUE,FALSE)
  bigmap$pos_DC2_3<-ifelse(bigmap$sig_DC2_3>=0.5,TRUE,FALSE)
  
  
  load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_my_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
  meta_annot_perso$x<-meta_annot$x
  meta_annot_perso$y<-meta_annot$y
  meta_annot_perso$currentscore<-NA
  match3<-bigmap$sig_Macs[match(rownames(meta_annot_perso), bigmap$cell_id)]
  meta_annot_perso$currentscore[!is.na(match3)]<-match3[!is.na(match3)]
  meta_annot_perso$my_annots_lin<-meta_annot_perso$my_annots
  
  pal_hiro<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0","#74BCD4","#548FAC","#396794","#20466D")
  

  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"] <- "MNP_N"
  # Condition: put "Neutro" in MNPN_ann for the cells that are MNP_N, not annotated yet, and positive to "pos_neutro" (but match names of the object slide with all lineages with the MNPN object with all slides)
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Neutro"])] <-"Neutro"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Macs"]& 
                                        !bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Mono"])] <-"Macs"  #Pos to Macs but not Neutro nor Mono
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Mono"]& 
                                        !bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Macs"])] <-"Mono"  #Pos to Mono but not Neutro nor Macs
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Macs"]& 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Mono"] )] <-"Momacs_transition"  # Pos to both Mono and Macs but not Neutro

  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_pDC"])] <-"pDC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_aDC"])] <-"aDC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="aDC" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_iaDC"])] <-"iaDC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_DC1"])] <-"DC1"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_DC2_3"])] <-"DC2_3"

  # Second threshold lower to annotate more cells
  bigmap$pos_Neutro<-ifelse(bigmap$sig_Neutro>=0.1,TRUE,FALSE)
  bigmap$pos_Macs<-ifelse(bigmap$sig_Macs>=0.1,TRUE,FALSE)
  bigmap$pos_Mono<-ifelse(bigmap$sig_Mono>=0.1,TRUE,FALSE)
  bigmap$pos_pDC<-ifelse(bigmap$sig_pDC>=0.1,TRUE,FALSE)
  bigmap$pos_aDC<-ifelse(bigmap$sig_aDC>=0.1,TRUE,FALSE)
  bigmap$pos_iaDC<-ifelse(bigmap$sig_iaDC>=0.1,TRUE,FALSE)
  bigmap$pos_DC1<-ifelse(bigmap$sig_DC1>=0.1,TRUE,FALSE)
  bigmap$pos_DC2_3<-ifelse(bigmap$sig_DC2_3>=0.1,TRUE,FALSE)
  
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Neutro"])] <-"Neutro"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Macs"]& 
                                        !bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Mono"])] <-"Macs" 
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Mono"]& 
                                        !bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Macs"])] <-"Mono" 
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Macs"]& 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Mono"] )] <-"Momacs_transition"

  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_pDC"])] <-"pDC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_aDC"])] <-"aDC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="aDC" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_iaDC"])] <-"iaDC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_DC1"])] <-"DC1"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]=="MNP_N" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_DC2_3"])] <-"DC2_3"

  
  print(table(meta_annot_perso$MNPN_ann))
  print(round(prop.table(table(meta_annot_perso$MNPN_ann))*100,2))
  
  sel_lin2<-c("iaDC","aDC","DC2_3","DC1","Macs","MNP_N","Mono","Neutro","pDC","Momacs_transition")
  print(paste0(spl,": ",sum(table(meta_annot_perso$MNPN_ann[which(meta_annot_perso$MNPN_ann!="MNP_N")]))," annotated cells for ",sum(table(meta_annot_perso$MNPN_ann)),
               " MNP/N cells (",round(100*sum(table(meta_annot_perso$MNPN_ann[which(meta_annot_perso$MNPN_ann!="MNP_N")]))/sum(table(meta_annot_perso$MNPN_ann)),1),")%"))
  
  meta_annot_perso$merge_annot<-meta_annot_perso$my_annots
  meta_annot_perso$merge_annot[which(!is.na(meta_annot_perso$MNPN_ann))]<-meta_annot_perso$MNPN_ann[which(!is.na(meta_annot_perso$MNPN_ann))]
  table(meta_annot_perso$merge_annot)
  
  
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  save(meta_annot_perso,file=paste0("./Grouped_objects/Xenium/meta_annot_june1_",spl,".rd"))  
} # End of MNPN loop
        
  
  
  
  
######## Mono subsets


sel_lin=c("Mono","Momacs_transition")
spl="GIM23_InfROI1"
load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_june1_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
bigmap<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
big_dmrc<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%bigmap$cell),]
rownames(big_dmrc)<-rownames(bigmap)


for(spl in c("GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  print(spl)
  load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_june1_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
  set.seed(42)
  
  sub_ma_mono<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
  sub_dmrc_mono<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%sub_ma_mono$cell),colnames(big_dmrc)] #colnames to ensure the gene order is the same before merging
  rownames(sub_dmrc_mono)<-rownames(sub_ma_mono)
  
  bigmap<-rbind(bigmap,sub_ma_mono)
  big_dmrc<-rbind(big_dmrc,sub_dmrc_mono)
}  




sig_InflMono1<-unlist(strsplit("GATA2,CD80,IL7R,IDO1,SPP1,EDN1,MET,CCL5,CSF3,CD274,SLAMF1,STEAP4,IL3RA",split=','))
sig_TransiC<-unlist(strsplit("PLD4,CLEC10A,SLAMF8,CSF1R,HLA-DRB1,PLA2G7,GPR183",split=','))
sig_ImmReg<-unlist(strsplit("IL10,TFPI,CDK1,FSTL3,DUSP2,PRDM1,THAP2",split=','))

allsum<-Matrix::rowSums(big_dmrc) 
bigmap$cell_id<-rownames(bigmap)
for(i in c("sig_InflMono1","sig_TransiC","sig_ImmReg")){
  eval(parse(text=paste0("ilist<-",i)))
  bigmap$currentscore<-NA
  df_inter<-data.frame("cell_name"=bigmap$cell_id,
                       "score"=(Matrix::rowSums(big_dmrc[,ilist])/allsum) )
  df_inter$scale_score<-scale(df_inter$score)
  eval(parse(text=paste0("bigmap$",i,"<-as.vector(df_inter$scale_score)")))
}
dim(bigmap) # 47440 cells in c("Mono","Momacs_transition")


for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){

  bigmap$pos_InflMono1<-ifelse(bigmap$sig_InflMono1>=0.5,TRUE,FALSE)  
  bigmap$pos_ImmReg<-ifelse(bigmap$sig_ImmReg>=0.5,TRUE,FALSE)
  bigmap$pos_TransiC<-ifelse(bigmap$sig_TransiC>=0.5,TRUE,FALSE)

  
  load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_june1_",spl,".rd"))
  load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
  meta_annot_perso$x<-meta_annot$x
  meta_annot_perso$y<-meta_annot$y
  meta_annot_perso$currentscore<-NA
  match3<-bigmap$sig_Macs[match(rownames(meta_annot_perso), bigmap$cell_id)]
  meta_annot_perso$currentscore[!is.na(match3)]<-match3[!is.na(match3)]
  meta_annot_perso$my_annots_lin<-meta_annot_perso$my_annots
  
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]%in%sel_lin &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_InflMono1"])] <-"Mono1"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]%in%sel_lin &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_ImmReg"])] <-"Mono3"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]%in%sel_lin &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_TransiC"])] <-"Mono5"

  
  table(meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"my_annots"], meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"])
  
  bigmap$pos_InflMono1<-ifelse(bigmap$sig_InflMono1>=0.1,TRUE,FALSE)  
  bigmap$pos_TransiC<-ifelse(bigmap$sig_TransiC>=0.1,TRUE,FALSE)
  bigmap$pos_ImmReg<-ifelse(bigmap$sig_ImmReg>=0.1,TRUE,FALSE)
  
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]%in%sel_lin &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_InflMono1"])] <-"Mono1"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]%in%sel_lin &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_ImmReg"])] <-"Mono3"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "MNPN_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"MNPN_ann"]%in%sel_lin &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_TransiC"])] <-"Mono5"
  
  print(table(meta_annot_perso$MNPN_ann))
  print(round(prop.table(table(meta_annot_perso$MNPN_ann))*100,2))
  

  sel_lin2<-c("iaDC","aDC","DC2_3","DC1","Macs","MNP_N","Mono","Neutro","pDC","Momacs_transition","Mono1","Mono3","Mono5")
  print(paste0(spl,": ",sum(table(meta_annot_perso$MNPN_ann[which(meta_annot_perso$MNPN_ann!="MNP_N")]))," annotated cells for ",sum(table(meta_annot_perso$MNPN_ann)),
               " MNP/N cells (",round(100*sum(table(meta_annot_perso$MNPN_ann[which(meta_annot_perso$MNPN_ann!="MNP_N")]))/sum(table(meta_annot_perso$MNPN_ann)),1),")%"))
  
  meta_annot_perso$merge_annot<-meta_annot_perso$my_annots
  meta_annot_perso$merge_annot[which(!is.na(meta_annot_perso$MNPN_ann))]<-meta_annot_perso$MNPN_ann[which(!is.na(meta_annot_perso$MNPN_ann))]
  table(meta_annot_perso$merge_annot)
  
  
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  save(meta_annot_perso,file=paste0("./SubsetsAppObj/meta_annot_june2_",spl,".rd"))  
  

} # End of Mono loop




  


###########################

grp_dmrc<-function(sel_lin){
  spl="GIM23_InfROI1"
  load(paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_my_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_",spl, ".rd"))
  bigmap<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
  big_dmrc<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%bigmap$cell),]
  rownames(big_dmrc)<-rownames(bigmap)
  
  for(spl in c("GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
    print(spl)
    load(paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
    load(paste0("./SubsetsAppObj/meta_annot_my_",spl,".rd"))
    load(paste0("./SubsetsAppObj/meta_annot_",spl, ".rd"))
    set.seed(42)
    
    sub_ma_mnpn<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
    sub_dmrc_mnpn<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%sub_ma_mnpn$cell),colnames(big_dmrc)] #colnames to ensure the gene order is the same before merging
    rownames(sub_dmrc_mnpn)<-rownames(sub_ma_mnpn)
    
    bigmap<-rbind(bigmap,sub_ma_mnpn)
    big_dmrc<-rbind(big_dmrc,sub_dmrc_mnpn)
  }  
  return( list(bigmap=bigmap, big_dmrc=big_dmrc))
}


### Tcells:
# table(meta_annot_perso$my_annots_lin)
sel_lin=c("TILC")
TILC_obj<-grp_dmrc(sel_lin)
bigmap<-TILC_obj$bigmap
big_dmrc<-TILC_obj$big_dmrc

allsum<-Matrix::rowSums(big_dmrc) 
bigmap$cell_id<-rownames(bigmap)


bigmap$sig_n_CD8A<-as.vector( log2(1+(10000*big_dmrc[,"CD8A"]/allsum)) ) 
bigmap$sig_n_CD4<-as.vector( log2(1+(10000*big_dmrc[,"CD4"]/allsum)) ) 
bigmap$sig_n_CD3D<-as.vector( log2(1+(10000*big_dmrc[,"CD3D"]/allsum)) ) 
bigmap$sig_n_CXCL13<-as.vector( log2(1+(10000*big_dmrc[,"CXCL13"]/allsum)) ) 


sig_CD8eff<-unlist(strsplit("GZMB,PRF1,GNLY,IFNG,NKG7,GZMA,CXCR3",split=','))
# sig_CD4treg<-unlist(strsplit("IL2RA,TIGIT,CTLA4,FOXP3",split=','))


allsum<-Matrix::rowSums(big_dmrc) 
bigmap$cell_id<-rownames(bigmap)
for(i in c("sig_CD8eff")){
  eval(parse(text=paste0("ilist<-",i)))
  bigmap$currentscore<-NA
  # Add all gene counts from the signature
  df_inter<-data.frame("cell_name"=bigmap$cell_id,"score"=(Matrix::rowSums(big_dmrc[,ilist])/allsum) )
  df_inter$scale_score<-scale(df_inter$score)
  eval(parse(text=paste0("bigmap$",i,"<-as.vector(df_inter$scale_score)")))
}


dim(bigmap) # 308535 cells in c("TILC")
pal_hiro<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0","#74BCD4","#548FAC","#396794","#20466D")

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  # First annotation round:
  bigmap$pos_n_CD8A<-ifelse(bigmap$sig_n_CD8A>=3,TRUE,FALSE)
  bigmap$pos_n_CD4<-ifelse(bigmap$sig_n_CD4>=3,TRUE,FALSE)
  bigmap$pos_n_CD3D<-ifelse(bigmap$sig_n_CD3D>=1,TRUE,FALSE)
  bigmap$pos_n_CXCL13<-ifelse(bigmap$sig_n_CXCL13>=2,TRUE,FALSE)
  bigmap$pos_CD8eff<-ifelse(bigmap$sig_CD8eff>=0.5,TRUE,FALSE)
  
  load(paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_june2_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_",spl, ".rd"))
  
  meta_annot_perso$currentscore<-NA
  

  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"] <- "TILC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="TILC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CD8A"])] <-"CD8 Tcell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="TILC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CD4"])] <-"CD4 Tcell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="CD8 Tcell" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_CD8eff"])] <-"CD8 eff."
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="CD4 Tcell" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CXCL13"])] <-"CD4 Tph"
  
  
  
  # print(table(meta_annot_perso$TILC_ann))
  # print(round(prop.table(table(meta_annot_perso$TILC_ann))*100,2))

  #2nd annotation round
  bigmap$pos_n_CD8A<-ifelse(bigmap$sig_n_CD8A>=2,TRUE,FALSE)
  bigmap$pos_n_CD4<-ifelse(bigmap$sig_n_CD4>=2,TRUE,FALSE)
  bigmap$pos_n_CXCL13<-ifelse(bigmap$sig_n_CXCL13>=1.5,TRUE,FALSE)  # don't lower too uch as it's very high in Tph and expressed lower in Treg
  bigmap$pos_CD8eff<-ifelse(bigmap$sig_CD8eff>=0.1,TRUE,FALSE)
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="TILC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CD8A"])] <-"CD8 Tcell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="TILC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CD4"])] <-"CD4 Tcell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="CD8 Tcell" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_CD8eff"])] <-"CD8 eff."
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="CD4 Tcell" & 
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CXCL13"])] <-"CD4 Tph"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "TILC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"TILC_ann"]=="TILC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_CD3D"])] <-"CD3D TILC"
  
  print(table(meta_annot_perso$TILC_ann))
  print(round(prop.table(table(meta_annot_perso$TILC_ann))*100,2))
  
  sel_lin2<-c("CD4 Tcell","CD4 Tph","CD8 eff.","CD8 Tcell","TILC","CD3D TILC")
  print(paste0(spl,": ",sum(table(meta_annot_perso$TILC_ann[which(meta_annot_perso$TILC_ann!="TILC")]))," annotated cells for ",sum(table(meta_annot_perso$TILC_ann)),
               " T/ILC cells (",round(100*sum(table(meta_annot_perso$TILC_ann[which(meta_annot_perso$TILC_ann!="TILC")]))/sum(table(meta_annot_perso$TILC_ann)),1),")%"))
  
  
  
  meta_annot_perso$merge_annot<-meta_annot_perso$my_annots
  # meta_annot_perso$merge_annot[which(meta_annot_perso$my_annots%in%sel_lin)]<-"MNP_N"
  meta_annot_perso$merge_annot[which(!is.na(meta_annot_perso$TILC_ann))]<-meta_annot_perso$TILC_ann[which(!is.na(meta_annot_perso$TILC_ann))]
  table(meta_annot_perso$merge_annot)
  
  objec="T_ILC"
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  
  ct_names <- names(table(meta_annot_perso$my_annots))
  lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  lin_palette3["TILC"] <- "#86756f"
  lin_palette3["CD4 Tcell"] <- "#F77794"
  lin_palette3["CD4 Tph"] <- "#52A38C"
  lin_palette3["CD8 eff."] <- "#447AB1"
  lin_palette3["CD8 Tcell"] <- "#8Bb4dD"
  lin_palette3["CD3D TILC"] <- "#e99a14"

  save(meta_annot_perso,file=paste0("./SubsetsAppObj/meta_annot_june3_",spl,".rd"))
  
  
} # End of TILC loop


###########################




### Bcells /  Plasma cells:
# table(meta_annot_perso$my_annots_lin)
sel_lin=c("Bcell","PlasmaC")
B_PC_obj<-grp_dmrc(sel_lin)
bigmap<-B_PC_obj$bigmap
big_dmrc<-B_PC_obj$big_dmrc

allsum<-Matrix::rowSums(big_dmrc) 
bigmap$cell_id<-rownames(bigmap)

bigmap$sig_n_TNFRSF13B<-as.vector( log2(1+(10000*big_dmrc[,"TNFRSF13B"]/allsum)) ) 

sig_PC<-unlist(strsplit("MZB1,TNFRSF17,DERL3,TENT5C",split=','))
sig_NaiveB<-unlist(strsplit("SELL,C1orf162,TCL1A",split=','))
sig_GClike<-unlist(strsplit("BCL6,RGS16,SPI1,CD9",split=','))

for(i in c("sig_PC","sig_NaiveB","sig_GClike")){
  eval(parse(text=paste0("ilist<-",i)))
  bigmap$currentscore<-NA
  # Add all gene counts from the signature
  df_inter<-data.frame("cell_name"=bigmap$cell_id,"score"=(Matrix::rowSums(big_dmrc[,ilist])/allsum) )
  df_inter$scale_score<-scale(df_inter$score)
  eval(parse(text=paste0("bigmap$",i,"<-as.vector(df_inter$scale_score)")))
}


dim(bigmap) # 280089 cells in c("B_PC")
pal_hiro<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0","#74BCD4","#548FAC","#396794","#20466D")

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  # First annotation round:
  bigmap$pos_n_TNFRSF13B<-ifelse(bigmap$sig_n_TNFRSF13B>=7,TRUE,FALSE)
  bigmap$pos_PC<-ifelse(bigmap$sig_PC>=0.5,TRUE,FALSE)
  bigmap$pos_NaiveB<-ifelse(bigmap$sig_NaiveB>=0.5,TRUE,FALSE)
  bigmap$pos_GClike<-ifelse(bigmap$sig_GClike>=0.5,TRUE,FALSE)
 
  
  load(paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_june3_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_",spl, ".rd"))
  
  meta_annot_perso$currentscore<-NA
  
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"] <- "B_PC"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==T &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"Plasma Cell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==T  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"Naive B cell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==T &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"GClike B cell"  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==T &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"Mem. B cell" 
  
  # print(table(meta_annot_perso$B_PC_ann))
  # print(round(prop.table(table(meta_annot_perso$B_PC_ann))*100,2))
  
  #2nd annotation round
  bigmap$pos_n_TNFRSF13B<-ifelse(bigmap$sig_n_TNFRSF13B>=5,TRUE,FALSE)
  bigmap$pos_PC<-ifelse(bigmap$sig_PC>=0.2,TRUE,FALSE)
  bigmap$pos_NaiveB<-ifelse(bigmap$sig_NaiveB>=0.2,TRUE,FALSE)
  bigmap$pos_GClike<-ifelse(bigmap$sig_GClike>=0.2,TRUE,FALSE)
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==T &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"Plasma Cell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==T  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"Naive B cell"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==T &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==F )] <-"GClike B cell"  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "B_PC_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"B_PC_ann"]=="B_PC" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_PC"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_NaiveB"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_GClike"]==F &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_TNFRSF13B"]==T )] <-"Mem. B cell" 
  
  print(table(meta_annot_perso$B_PC_ann))
  print(round(prop.table(table(meta_annot_perso$B_PC_ann))*100,2))
  
  sel_lin2<-c("B_PC","Plasma Cell","Naive B cell","GClike B cell","Mem. B cell")
  print(paste0(spl,": ",sum(table(meta_annot_perso$B_PC_ann[which(meta_annot_perso$B_PC_ann!="B_PC")]))," annotated cells for ",sum(table(meta_annot_perso$B_PC_ann)),
               " B/PC cells (",round(100*sum(table(meta_annot_perso$B_PC_ann[which(meta_annot_perso$B_PC_ann!="B_PC")]))/sum(table(meta_annot_perso$B_PC_ann)),1),")%"))
  
  
  
  meta_annot_perso$merge_annot<-meta_annot_perso$my_annots
  meta_annot_perso$merge_annot[which(!is.na(meta_annot_perso$B_PC_ann))]<-meta_annot_perso$B_PC_ann[which(!is.na(meta_annot_perso$B_PC_ann))]
  table(meta_annot_perso$merge_annot)
  
  objec="B_PC"
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  
  ct_names <- names(table(meta_annot_perso$my_annots))
  lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  lin_palette3["B_PC"] <- "#b6958f"
  lin_palette3["Plasma Cell"] <- "#E56157"
  lin_palette3["Naive B cell"] <- "#FFA01C"
  lin_palette3["GClike B cell"] <- "#CC9F94"
  lin_palette3["Mem. B cell"] <- "#4273B7"

  
  
  
  # 
  # pdf(file=paste0("./Figures/Xen_Lin_myannot_B_PC_",spl,"_v3scaleonall.pdf"),width = 13,height = 10)
  # levels(lin_palette3)<-names(lin_palette3)
  # brplt_title<-paste0(spl," -- projection of main cell types in ", objec )
  # print(ggplot(meta_annot_perso, aes(x = x, y = y, color = merge_annot)) +
  #         rasterise(geom_point(size = 0.1)) +
  #         scale_color_manual(values = lin_palette3)+
  #         guides(color = guide_legend(override.aes = list(size = 4)))+
  #         theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #         labs(title = brplt_title,
  #              color = "Lineage" )+
  #         theme(axis.text = element_blank(),
  #               axis.ticks = element_blank(),
  #               axis.title = element_blank()) )
  # for(ct in c("B_PC","Plasma Cell","Naive B cell","GClike B cell","Mem. B cell")){
  #   meta_annot_perso$current_ct<-NA
  #   meta_annot_perso$current_ct[which(meta_annot_perso$my_annots%in%sel_lin2)]<-"B_PC"
  #   meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot==ct)]<-ct
  #   pt1<-ggplot(meta_annot_perso, aes(x = x, y = y, color = current_ct)) +
  #     rasterise(geom_point(data = subset(meta_annot_perso, is.na(current_ct) ), aes(x = x, y = y),
  #                          color = "#EDEDED", size = 0.2)) +
  #     rasterise(geom_point(data = subset(meta_annot_perso, current_ct == "B_PC"), aes(x = x, y = y),
  #                          color = "#D0D0D0", size = 0.5)) +
  #     geom_point(data = subset(meta_annot_perso, current_ct == ct), aes(x=x,y=y), color =  lin_palette3[ct], size = 0.8)  +
  #     theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #     labs(title = paste0(spl," -- projection of ",ct," in ", objec ),
  #          color = "Lineage" )+
  #     theme(axis.text = element_blank(),
  #           axis.ticks = element_blank(),
  #           axis.title = element_blank())
  #   print(pt1)
  # }
  # 
  # dev.off()
  save(meta_annot_perso,file=paste0("./SubsetsAppObj/meta_annot_june4_",spl,".rd"))
  
  
} # End of B_PC loop

###########################




### Fibroblasts:
# table(meta_annot_perso$my_annots_lin)
sel_lin=c("Fibro","addPeri")
FIB_obj<-grp_dmrc(sel_lin)
bigmap<-FIB_obj$bigmap
big_dmrc<-FIB_obj$big_dmrc

allsum<-Matrix::rowSums(big_dmrc) 
bigmap$cell_id<-rownames(bigmap)


bigmap$sig_n_SFRP2<-as.vector( log2(1+(10000*big_dmrc[,"SFRP2"]/allsum)) ) 

sig_SM<-unlist(strsplit("MYH11,ACTA2,ACTG2,DES",split=','))
sig_SM_RERGL<-unlist(strsplit("RERGL,NTN4,MEF2C",split=','))
sig_Infl_Fib<-unlist(strsplit("CHI3L1,CAV1,FKBP11,GLIPR1,PDPN,PCOLCE,NPDC1,RAMP2,CXCL13,CCL3,CFB,CD40",split=','))
# sig_Infl_FRC<-unlist(strsplit("BASP1,MYLK,OGN,DPT,SFRP2,FSTL3,PTGDS,CCL19",split=','))


for(i in c("sig_SM","sig_Infl_Fib","sig_SM_RERGL")){
  eval(parse(text=paste0("ilist<-",i)))
  bigmap$currentscore<-NA
  # Add all gene counts from the signature
  df_inter<-data.frame("cell_name"=bigmap$cell_id,"score"=(Matrix::rowSums(big_dmrc[,ilist])/allsum) )
  df_inter$scale_score<-scale(df_inter$score)
  eval(parse(text=paste0("bigmap$",i,"<-as.vector(df_inter$scale_score)")))
}


dim(bigmap) # 475739 cells in c("FIB")
pal_hiro<-c("#E56157","#ED894D","#F6A95E","#FECF75","#FEE6B9","#ABDCE0","#74BCD4","#548FAC","#396794","#20466D")

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  # First annotation round:
  bigmap$pos_n_SFRP2<-ifelse(bigmap$sig_n_SFRP2>=2,TRUE,FALSE)
  
  bigmap$pos_SM<-ifelse(bigmap$sig_SM>=0.5,TRUE,FALSE)
  bigmap$pos_Infl_Fib<-ifelse(bigmap$sig_Infl_Fib>=0.5,TRUE,FALSE)
  bigmap$pos_SM_RERGL<-ifelse(bigmap$sig_SM_RERGL>=0.5,TRUE,FALSE)
  
  
  load(paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_june4_",spl,".rd"))
  load(paste0("./SubsetsAppObj/meta_annot_",spl, ".rd"))
  
  meta_annot_perso$currentscore<-NA
  
  
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"] <- "FIB"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==T &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM_RERGL"]==F  &
                                        bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==F  )] <-"Sm. Muscle"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==F &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM_RERGL"]==T  &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==F  )] <-"SM RERGL"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==F &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM_RERGL"]==F  &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==T  )] <-"Infl Fib"
  meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
                   "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==F &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM_RERGL"]==F  &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==F  &
                                       bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_n_SFRP2"]==T  )] <-"Submucosa Myofibro"
  
  # print(table(meta_annot_perso$FIB_ann))
  # print(round(prop.table(table(meta_annot_perso$FIB_ann))*100,2))
  
  # #2nd annotation round
  # bigmap$pos_SM<-ifelse(bigmap$sig_SM>=0.2,TRUE,FALSE)
  # bigmap$pos_Infl_Fib<-ifelse(bigmap$sig_Infl_Fib>=0.5,TRUE,FALSE)
  # bigmap$pos_Infl_FRC<-ifelse(bigmap$sig_Infl_FRC>=0.5,TRUE,FALSE)
  # 
  # meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
  #                  "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==T &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==F  &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_FRC"]==F  )] <-"Sm. Muscle"
  # meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
  #                  "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==F &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==T  &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_FRC"]==F  )] <-"Infl Fib"
  # meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],
  #                  "FIB_ann"][which( meta_annot_perso[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"FIB_ann"]=="FIB" &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_SM"]==F &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_Fib"]==F  &
  #                                      bigmap[rownames(meta_annot_perso)[meta_annot_perso$my_annots%in%sel_lin],"pos_Infl_FRC"]==T  )] <-"Infl FRC"
  # 
  print(table(meta_annot_perso$FIB_ann))
  print(round(prop.table(table(meta_annot_perso$FIB_ann))*100,2))
  
  sel_lin2<-c("FIB","Sm. Muscle","Infl Fib","SM RERGL","Submucosa Myofibro")
  print(paste0(spl,": ",sum(table(meta_annot_perso$FIB_ann[which(meta_annot_perso$FIB_ann!="FIB")]))," annotated cells for ",sum(table(meta_annot_perso$FIB_ann)),
               " Fibro cells (",round(100*sum(table(meta_annot_perso$FIB_ann[which(meta_annot_perso$FIB_ann!="FIB")]))/sum(table(meta_annot_perso$FIB_ann)),1),")%"))
  
  
  
  meta_annot_perso$merge_annot<-meta_annot_perso$my_annots
  meta_annot_perso$merge_annot[which(!is.na(meta_annot_perso$FIB_ann))]<-meta_annot_perso$FIB_ann[which(!is.na(meta_annot_perso$FIB_ann))]
  table(meta_annot_perso$merge_annot)
  
  objec="FIB"
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  
  ct_names <- names(table(meta_annot_perso$my_annots))
  lin_palette3 <- setNames(rep("#EDEDED", length(ct_names)), ct_names)
  lin_palette3["FIB"] <- "#F77774"
  lin_palette3["Sm. Muscle"] <- "#197F4F"
  lin_palette3["SM RERGL"] <- "#194FCF"
  lin_palette3["Infl Fib"] <- "#864036"
  lin_palette3["Submucosa Myofibro"] <- "#C967EA"
  
  
  
  
  
  # pdf(file=paste0("./Figures/Xen_Lin_myannot_FIB_",spl,"_v3scaleonall.pdf"),width = 13,height = 10)
  # levels(lin_palette3)<-names(lin_palette3)
  # brplt_title<-paste0(spl," -- projection of main cell types in ", objec )
  # print(ggplot(meta_annot_perso, aes(x = x, y = y, color = merge_annot)) +
  #         rasterise(geom_point(size = 0.1)) +
  #         scale_color_manual(values = lin_palette3)+
  #         guides(color = guide_legend(override.aes = list(size = 4)))+
  #         theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #         labs(title = brplt_title,
  #              color = "Lineage" )+
  #         theme(axis.text = element_blank(),
  #               axis.ticks = element_blank(),
  #               axis.title = element_blank()) )
  # for(ct in sel_lin2){
  #   meta_annot_perso$current_ct<-NA
  #   meta_annot_perso$current_ct[which(meta_annot_perso$my_annots%in%sel_lin2)]<-"FIB"
  #   meta_annot_perso$current_ct[which(meta_annot_perso$merge_annot==ct)]<-ct
  #   pt1<-ggplot(meta_annot_perso, aes(x = x, y = y, color = current_ct)) +
  #     rasterise(geom_point(data = subset(meta_annot_perso, is.na(current_ct) ), aes(x = x, y = y),
  #                          color = "#EDEDED", size = 0.2)) +
  #     rasterise(geom_point(data = subset(meta_annot_perso, current_ct == "FIB"), aes(x = x, y = y),
  #                          color = "#D0D0D0", size = 0.5)) +
  #     geom_point(data = subset(meta_annot_perso, current_ct == ct), aes(x=x,y=y), color =  lin_palette3[ct], size = 0.8)  +
  #     theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #     labs(title = paste0(spl," -- projection of ",ct," in ", objec ),
  #          color = "Lineage" )+
  #     theme(axis.text = element_blank(),
  #           axis.ticks = element_blank(),
  #           axis.title = element_blank())
  #   print(pt1)
  # }
  # 
  # dev.off()
  save(meta_annot_perso,file=paste0("./SubsetsAppObj/meta_annot_june5_",spl,".rd"))
  
  
} # End of FIB loop





################################
library(dplyr)
spl="GIM23_InfROI1"

load(paste0("./SubsetsAppObj/meta_annot_june5_",spl,".rd"))

df_re4<-data.frame(my_annots=names(table(meta_annot_perso$my_annots)),group=NA)
# df_re4$group[df_re4$my_annots%in%c("Act_DC","DC1","Macs","MNP_N","Mono","Mono1","Neutro","pDC")]<-"MNP_Neutro_pDC"
df_re4$group[df_re4$my_annots%in%c("iaDC","aDC","DC2_3","DC1","Macs","MNP_N","Mono","Mono1","Mono3","Mono5","Momacs_transition","Neutro","pDC")]<-"MNP_Neutro_pDC"
df_re4$group[df_re4$my_annots%in%c("addNKlike","CD3D TILC","CD4 Tcell","CD4 Tph","CD8 eff.","CD8 Tcell","TILC")]<-"T_ILC_NK"
df_re4$group[df_re4$my_annots%in%c("B_PC","GClike B cell","Mem. B cell","Naive B cell","Plasma Cell")]<-"B_PC"
# df_re4$group[df_re4$my_annots%in%c("FIB","Infl Fib","Infl FRC","Sm. Muscle")]<-"Stromal"
df_re4$group[df_re4$my_annots%in%c("FIB","Sm. Muscle","Infl Fib","SM RERGL","Submucosa Myofibro")]<-"Stromal"
names(table(df_re4$group))
save(df_re4,file="./SubsetsAppObj/meta_annot_june5_names.rd")

