library(shiny)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(Rfast)
library(dplyr)

load("./Grouped_objects/metadata_mac20_1_3sig.rd")
colnames(metadata)<-gsub("data_Siggrp","MacSig",colnames(metadata))
metadata$batch<-gsub("Mac(\\d)(\\d)[.].*","Mac\\1\\2",rownames(metadata))
load("./Grouped_objects/data_counts_mac20_1_3.rd")
load("./Grouped_objects/data_logcpm_mac20_1_3.rd")



pal=rev(rainbow(15))[4:15]
col_pal<-c("GM"="#008000","GM_TNF"="#008000","GM_IFNg"="#40B740","GM_TNF_IFNg"="#40B740","GM_(TNF)_IFNg"="#40B740",
           "GM_LPS"="#FF4A42","GM_LPS_TNF"="#FF4A42","GM_LPS_IFNg"="#8E4700","GM_LPS_TNF_IFNg"="#8E4700",
           "GM_IL4"="#06677A","GM_TNF_IL4"="#06677A","GM_(TNF)_IL4"="#06677A","GM_IL4_IFNg"="#11A9C4","GM_TNF_IL4_IFNg"="#11A9C4",
           "GM_LPS_IL4"="#2929A0","GM_LPS_TNF_IL4"="#2929A0","GM_LPS_IL4_IFNg"="#7070E0","GM_LPS_TNF_IL4_IFNg"="#7070E0" )
batchpalette<-c("Mac20"="#648FFF","Mac21"="#DC267F","Mac22"="#CE9413")

#allsum<-Rfast::colsums(log(1+data))    ##V1_1
allsum<-Rfast::colsums(l.data.cpm)  


metadata$condition<-factor(metadata$condition, levels=names(col_pal))


tl_data<-as.data.frame(as.matrix(t(l.data.cpm)))
# metadata$condition[match(rownames(minidf),rownames(metadata))] # Je sais plus l'ordre pour match, mais c'est le m?me donc autant faire :

merged_df <- merge(metadata[,c(1:2,8:ncol(metadata))], tl_data, by = "row.names")   # don't keep LPS/TNF/GM detail as it might interfer with genes
rownames(merged_df) <- merged_df$Row.names
merged_df <- merged_df[, -1] 
merged_df_complete<-merged_df
merged_df_complete




# PROGENy gene list:
model <- progeny::model_human_full
model_500 <- model %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 500)
model_500$pathway<-gsub("-","_",model_500$pathway)
# pthw<-names(table(model_500$pathway))
pthw<-c("NFkB","JAK_STAT","MAPK","Hypoxia")


out<-data.frame(row.names = rownames(merged_df), condition=merged_df$condition)

for(i in 1:length(pthw)){
  mods=pthw[i]
  i_y2<-model_500[which(model_500$pathway==mods&model_500$weight>0),]$gene
  merged_df$i_y2<-NA
  rm_2<-"Mac22.13,Mac22.6,Mac20.3.4,Mac20.5.6,Mac22.12,Mac21.1,Mac21.2,Mac21.3,Mac21.4,Mac21.5,Mac21.6,Mac21.7,Mac21.8,Mac21.9,Mac21.10,Mac21.11,Mac21.12,Mac21.13,Mac21.14,Mac21.15,Mac21.16"
  i_rm2<-unlist(strsplit(rm_2,split = ","))
  merged_df<-merged_df_complete[setdiff(rownames(merged_df_complete),i_rm2),]
  # X_order: revDGE_ord2_noIL4
  ord_x3=rev(c("GM","GM_TNF","GM_IFNg",
               "GM_LPS","GM_LPS_TNF","GM_LPS_IFNg","GM_LPS_TNF_IFNg") )
  colorpalette<-c("GM"="#2B5637","GM_TNF"="#008900","GM_IFNg"="#48C143",
                  "GM_LPS"="#FF6579","GM_LPS_TNF"="#CE2F2F","GM_LPS_IFNg"="#AA6747","GM_LPS_TNF_IFNg"="#682E0A")
  
  gnshow_y<-gsub(" ","",i_y2) 
  notfound_y<-setdiff(gnshow_y,rownames(l.data.cpm))
  gnshow_y2<-setdiff(gnshow_y,notfound_y)
  y_name=mods
  currentscore<-(Rfast::colsums(l.data.cpm[gnshow_y2,])/allsum )
  names(currentscore)<-colnames(l.data.cpm)
  currentscore<-currentscore[setdiff(names(currentscore),i_rm2)]
  crt_sc<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))
  names(crt_sc)<-names(currentscore)
  eval(parse(text=paste0("merged_df$i_y2<-crt_sc[match(rownames(merged_df),names(crt_sc))]")))
  
  
  brbx_2<-ggboxplot(merged_df, x = "condition", y = "i_y2", 
                    color = "condition",fill= "condition",alpha=0.2,palette = colorpalette, size = rel(1))+
    scale_x_discrete(limits = ord_x3)+
    theme(legend.position = "none", 
          axis.text.y = element_text(size = rel(3)),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
    ggtitle(paste0(y_name))+ylab(label = y_name)  +
    theme(plot.title = element_text(size = rel(2), face = "bold"))    
  
  if (min(range(merged_df$i_y2)) < 0 && max(range(merged_df$i_y2)) > 0) {
    brbx_2 <- brbx_2 + geom_hline(yintercept=0, linetype="dashed", color = "lightgrey")
  }
  
  print(brbx_2)
  eval(parse(text=paste0("out$score_",mods,"<-merged_df$i_y2")))
  
}

write.csv2(out,file="Gaelle_Fig_4D_S10A_S10D_invitro_datapoints.csv",quote=F,row.names=F)


########################################################
# Redo Fig clean
