# Sub annotations, based on single lineages across samples
# Hypothesis that batch correction not needed for probe levels ?

#####################################
#######  Subtypes annotation
#####################################

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")

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


grp_dmrc<-function(sel_lin){
  spl="GIM23_InfROI1"
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_my_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
  bigmap<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
  big_dmrc<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%bigmap$cell),]
  rownames(big_dmrc)<-rownames(bigmap)
  
  
  for(spl in c("GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
    print(spl)
    load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
    load(file=paste0("./Grouped_objects/meta_annot_my_",spl,".rd"))
    load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
    set.seed(42)
    
    sub_ma_mnpn<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
    sub_dmrc_mnpn<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%sub_ma_mnpn$cell),colnames(big_dmrc)] #colnames to ensure the gene order is the same before merging
    rownames(sub_dmrc_mnpn)<-rownames(sub_ma_mnpn)
    
    bigmap<-rbind(bigmap,sub_ma_mnpn)
    big_dmrc<-rbind(big_dmrc,sub_dmrc_mnpn)
  }  
  return( list(bigmap=bigmap, big_dmrc=big_dmrc))
}


##------------------------------##
##         MNP + Neutro         ##
##------------------------------##

sel_lin=c("MNPNeutro","MNPNeutroInfl")
MNPN_obj<-grp_dmrc(sel_lin)
bigmap<-MNPN_obj$bigmap
big_dmrc<-MNPN_obj$big_dmrc


for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  print(spl)
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  print(dim(dgcmtx_raw_counts))

}

sig_Macs<-unlist(strsplit("FOLR2,IGF1,KCNMA1,C1QA,AXL,VSIG4,DNASE1L3,GATM,SIGLEC1",split=','))
sig_Mono<-unlist(strsplit("CLEC4E,RETN,FCN1,LILRA5,MARCO,MCEMP1,FCAR,CD93,ANPEP,APOBEC3A,CCR1,VCAN,CD300E",split=','))
sig_Neutro<-unlist(strsplit("CXCR2,MME",split=','))
sig_pDC<-unlist(strsplit("GZMB,CXCR3,LILRA4,TCL1A",split=','))
sig_aDC<-unlist(strsplit("CCR7,LAD1,FSCN1,LAMP3",split=','))
sig_iaDC<-unlist(strsplit("CXCL9,CXCL10,C15orf48,SRC",split=',')) # permet de se restreindre à ceux proches de l'ulcération
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
  
  
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_my_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
  meta_annot_perso$x<-meta_annot$x
  meta_annot_perso$y<-meta_annot$y
  meta_annot_perso$currentscore<-NA
  match3<-bigmap$sig_Macs[match(rownames(meta_annot_perso), bigmap$cell_id)]
  meta_annot_perso$currentscore[!is.na(match3)]<-match3[!is.na(match3)]
  meta_annot_perso$my_annots_lin<-meta_annot_perso$my_annots
  
  
  # Annotation attribution based on score positivity:
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
  

  # 2nd pass with lower threshold to annotate more cells:
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
  save(meta_annot_perso,file=paste0("G:/Mon Drive/UG_metacells/Xenium/SubsetsAppObj/meta_annot_june1_",spl,".rd"))  
} # End of MNPN loop





##------------------------------##
##          Monocytes           ##
##------------------------------##
rm(list=ls())

sel_lin=c("Mono","Momacs_transition")
spl="GIM23_InfROI1"
load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
load(file=paste0("./Grouped_objects/meta_annot_june1_",spl,".rd"))
load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
bigmap<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
big_dmrc<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%bigmap$cell),]
rownames(big_dmrc)<-rownames(bigmap)


for(spl in c("GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  print(spl)
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_june1_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
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
  
  
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_june1_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
  meta_annot_perso$x<-meta_annot$x
  meta_annot_perso$y<-meta_annot$y
  meta_annot_perso$currentscore<-NA
  match3<-bigmap$sig_Macs[match(rownames(meta_annot_perso), bigmap$cell_id)]
  meta_annot_perso$currentscore[!is.na(match3)]<-match3[!is.na(match3)]
  meta_annot_perso$my_annots_lin<-meta_annot_perso$my_annots
  
  # Annotation attribution based on score positivity:
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


  # 2nd pass with lower threshold to annotate more cells:
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
  save(meta_annot_perso,file=paste0("./Grouped_objects/meta_annot_june2_",spl,".rd"))  
  
  
} 


#Figures
for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_june2_",spl,".rd"))
  
  sel_lin2<-c("iaDC","aDC","DC2_3","DC1","Macs","MNP_N","Mono","Neutro","pDC","Momacs_transition","Mono1","Mono3","Mono5")

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


###########################

# Continue the annotations in the other lineages to get an annotation sheet for each slide to load on the visualisation app

###########################

grp_dmrc<-function(sel_lin){
  spl="GIM23_InfROI1"
  load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_my_",spl,".rd"))
  load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
  bigmap<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
  big_dmrc<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%bigmap$cell),]
  rownames(big_dmrc)<-rownames(bigmap)
  
  
  for(spl in c("GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
    print(spl)
    load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))
    load(file=paste0("./Grouped_objects/meta_annot_my_",spl,".rd"))
    load(file=paste0("./Grouped_objects/meta_annot_",spl, ".rd"))
    set.seed(42)
    
    sub_ma_mnpn<-meta_annot[which(meta_annot_perso$my_annots%in%sel_lin),]
    sub_dmrc_mnpn<-dgcmtx_raw_counts[which(rownames(dgcmtx_raw_counts)%in%sub_ma_mnpn$cell),colnames(big_dmrc)] #colnames to ensure the gene order is the same before merging
    rownames(sub_dmrc_mnpn)<-rownames(sub_ma_mnpn)
    
    bigmap<-rbind(bigmap,sub_ma_mnpn)
    big_dmrc<-rbind(big_dmrc,sub_dmrc_mnpn)
  }  
  return( list(bigmap=bigmap, big_dmrc=big_dmrc))
}



##------------------------------##
##           T cells            ##
##------------------------------##


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

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  # First annotation round:
  bigmap$pos_n_CD8A<-ifelse(bigmap$sig_n_CD8A>=3,TRUE,FALSE)
  bigmap$pos_n_CD4<-ifelse(bigmap$sig_n_CD4>=3,TRUE,FALSE)
  bigmap$pos_n_CD3D<-ifelse(bigmap$sig_n_CD3D>=1,TRUE,FALSE)
  bigmap$pos_n_CXCL13<-ifelse(bigmap$sig_n_CXCL13>=2,TRUE,FALSE)
  bigmap$pos_CD8eff<-ifelse(bigmap$sig_CD8eff>=0.5,TRUE,FALSE)
  
  load(file=paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("G:/Mon Drive/UG_metacells/Xenium/SubsetsAppObj/meta_annot_june2_",spl,".rd"))
  load(file=paste0("/SubsetsAppObj/meta_annot_",spl, ".rd"))
  
  meta_annot_perso$currentscore<-NA
  
  # Annotation attribution based on score positivity:
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
  

  # 2nd pass with lower threshold to annotate more cells:
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
  meta_annot_perso$merge_annot[which(!is.na(meta_annot_perso$TILC_ann))]<-meta_annot_perso$TILC_ann[which(!is.na(meta_annot_perso$TILC_ann))]
  table(meta_annot_perso$merge_annot)
  
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  
  
  save(meta_annot_perso,file=paste0("./Grouped_objects/meta_annot_june3_",spl,".rd"))
  
  
} # End of TILC loop


###########################


##------------------------------##
##     B cells / Plasmacells    ##
##------------------------------##

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

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  # First annotation round:
  bigmap$pos_n_TNFRSF13B<-ifelse(bigmap$sig_n_TNFRSF13B>=7,TRUE,FALSE)
  bigmap$pos_PC<-ifelse(bigmap$sig_PC>=0.5,TRUE,FALSE)
  bigmap$pos_NaiveB<-ifelse(bigmap$sig_NaiveB>=0.5,TRUE,FALSE)
  bigmap$pos_GClike<-ifelse(bigmap$sig_GClike>=0.5,TRUE,FALSE)
  
  
  load(file=paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("G:/Mon Drive/UG_metacells/Xenium/SubsetsAppObj/meta_annot_june3_",spl,".rd"))
  load(file=paste0("/SubsetsAppObj/meta_annot_",spl, ".rd"))
  
  meta_annot_perso$currentscore<-NA
  
  # Annotation attribution based on score positivity:
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
  

  # 2nd pass with lower threshold to annotate more cells:
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
  
  meta_annot_perso$my_annots<-meta_annot_perso$merge_annot
  
  save(meta_annot_perso,file=paste0("./Grouped_objects/meta_annot_june4_",spl,".rd"))
  
  
} # End of B_PC loop

###########################



##------------------------------##
##         Fibroblasts          ##
##------------------------------##

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

for(spl in c("GIM23_InfROI1","GIM23_InfROI2","GIM23_NonInfl","GIM33_InfROI1","GIM33_InfROI2","GIM38_InfROI1","GIM38_InfROI2","GIM38_NonInfl")){
  
  # First annotation round:
  bigmap$pos_n_SFRP2<-ifelse(bigmap$sig_n_SFRP2>=2,TRUE,FALSE)
  
  bigmap$pos_SM<-ifelse(bigmap$sig_SM>=0.5,TRUE,FALSE)
  bigmap$pos_Infl_Fib<-ifelse(bigmap$sig_Infl_Fib>=0.5,TRUE,FALSE)
  bigmap$pos_SM_RERGL<-ifelse(bigmap$sig_SM_RERGL>=0.5,TRUE,FALSE)
  
  
  load(file=paste0("./SubsetsAppObj/dgcmtx_raw_",spl,".rd"))
  load(file=paste0("G:/Mon Drive/UG_metacells/Xenium/SubsetsAppObj/meta_annot_june4_",spl,".rd"))
  load(file=paste0("/SubsetsAppObj/meta_annot_",spl, ".rd"))
  
  meta_annot_perso$currentscore<-NA
  
  # Annotation attribution based on score positivity:
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
  
  save(meta_annot_perso,file=paste0("./Grouped_objects/meta_annot_june5_",spl,".rd"))
  
  
} # End of FIB loop



