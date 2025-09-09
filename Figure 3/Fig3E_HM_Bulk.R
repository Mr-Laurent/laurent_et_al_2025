library(ComplexHeatmap)
library(ggplot2)
library(scales)
library(reshape)
library(reshape2)
library(dplyr)
library(Rfast)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")


load("./Grouped_objects/metadata_mac20_1_2.rd")
load("./Grouped_objects/data_counts_mac20_1_2.rd")
load("./Grouped_objects/data_logcpm_mac20_1_3.rd")

#### If personnal lists to add :
load("./Grouped_objects/genelist_16nov22_Macs.rd")
automate_strsplit <- function(df, col, indices) {
  unlist(strsplit(paste0(df[[col]][indices], collapse = ","), split = ","))
}

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



# REMOVE BAD SAMPLES
rm_spl<-c("Mac22.13","Mac22.6","Mac20.3.4","Mac20.5.6","Mac22.12","Mac21.1","Mac21.2","Mac21.3","Mac21.4","Mac21.5","Mac21.6","Mac21.7",
          "Mac21.8","Mac21.9","Mac21.10","Mac21.11","Mac21.12","Mac21.13","Mac21.14","Mac21.15","Mac21.16")

# REMOVE IL4 COND TOO
rm_spl2<-c("Mac22.3","Mac22.4","Mac20.7","Mac20.8","Mac22.7","Mac22.8","Mac20.11","Mac20.12","Mac22.11","Mac22.12","Mac20.15","Mac20.16","Mac22.15","Mac22.16")




ldatacpm<-l.data.cpm[,-which(colnames(l.data.cpm)%in%rm_spl)] 
ldatacpm<-ldatacpm[,-which(colnames(ldatacpm)%in%rm_spl2)] 
lmdcpm<-metadata[-which(rownames(metadata)%in%rm_spl),]
lmdcpm<-lmdcpm[-which(rownames(lmdcpm)%in%rm_spl2),]


# Some genes are differently named : MARCH1 (scRNAseq) -> MARCHF1 (bulk)   # Might do it for the HM too ???

rownames(ldatacpm)[which(rownames(ldatacpm)=="MARCHF1")]<-"MARCH1"
rownames(ldatacpm)[which(rownames(ldatacpm)=="WARS1")]<-"WARS"
rownames(ldatacpm)[which(rownames(ldatacpm)=="TAMALIN")]<-"GRASP"
rownames(ldatacpm)[which(rownames(ldatacpm)=="POLR1F")]<-"TWISTNB"
rownames(ldatacpm)[which(rownames(ldatacpm)=="PELATON")]<-"SMIM25"
rownames(ldatacpm)[which(rownames(ldatacpm)=="PRECSIT")]<-"LINC00346"
rownames(ldatacpm)[which(rownames(ldatacpm)=="MAILR")]<-"AZIN1-AS1"
rownames(ldatacpm)[which(rownames(ldatacpm)=="PALM2AKAP2")]<-"PALM2-AKAP2"
rownames(ldatacpm)[which(rownames(ldatacpm)=="MIR9-1HG")]<-"C1orf61"
rownames(ldatacpm)[which(rownames(ldatacpm)=="TARP")]<-"TRGC1"


#### Do it with mean expression by condition ####
tldatacpm<-as.data.frame(t(ldatacpm))
tldatacpm$cond<- lmdcpm$condition[match(lmdcpm$name,rownames(tldatacpm))]

expression_data_long <- tidyr::gather(tldatacpm, Gene, Expression, -cond)  # also possible with pivot_loner or melt I guess ?
summary_data <- expression_data_long %>%
  group_by(cond, Gene) %>%
  summarise(Mean_Expression = mean(Expression))

tldatacpm2<-reshape2::dcast(summary_data, cond ~ Gene)
rownames(tldatacpm2) <- tldatacpm2$cond
tldatacpm2$cond <- NULL 
ldatacpm2<-as.matrix(t(tldatacpm2))



metad2<-data.frame(condition=colnames(ldatacpm2))

allsum2<-Rfast::colsums(ldatacpm2) 
for(Mod_score in Modsig$signatures){
  print(Mod_score)
  geneshow<-eval(parse(text=Mod_score))
  ilist<-intersect(geneshow,rownames(ldatacpm2))
  if(length(ilist)<2){currentscore<-(ldatacpm2[intersect(ilist,rownames(ldatacpm2)),]/allsum2 )
  eval(parse(text=paste0("metad2$lodacp_",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }else{currentscore<-(Rfast::colsums(ldatacpm2[ilist,])/allsum2 )
  eval(parse(text=paste0("metad2$lodacp_",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }}

rownames(metad2)<-metad2$condition
metad2$condition <- NULL

colorpalette<-c("GM"="#008000","GM_TNF"="#008000",
                "GM_IFNg"="#40B740","GM_TNF_IFNg"="#40B740","GM_(TNF)_IFNg"="#40B740",
                "GM_LPS"="#FF4A42","GM_LPS_TNF"="#FF4A42",
                "GM_LPS_IFNg"="#8E4700","GM_LPS_TNF_IFNg"="#8E4700" )
ord_sig_perso2_log<-paste0("lodacp_",c("Siggrp_X","Siggrp_XI","Siggrp_XII","Siggrp_I","Siggrp_III","Siggrp_IV","Siggrp_VI",
                                       "Siggrp_II","Siggrp_XVI","Siggrp_XVII","Siggrp_XV","Siggrp_XIV","Siggrp_XIII","Siggrp_IX","Siggrp_VIII","Siggrp_VII","Siggrp_V") )

metad3<-metad2[names(colorpalette)[names(colorpalette)%in%colnames(ldatacpm2)],]
metad3$cond<-rownames(metad3)
mltmeta2<-reshape::melt(metad3)

plot_reorderselectedrow_selectedcol2<-ggplot(mltmeta2, aes(x = cond, y = variable, fill = value)) +
  geom_tile() +  scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish, name = "Median Score") +
  scale_x_discrete(limits = rev(rownames(metad3))) +scale_y_discrete(limits = ord_sig_perso2_log )+
  theme_bw() +  theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position = "right")+ geom_hline(yintercept = 0.5 + 0:40, colour = "black", size = 0.05)


pdf("./Figure 3/Fig3E_HM_Bulk.pdf",height = 8,width = 8)

plot_reorderselectedrow_selectedcol2+scale_y_discrete(limits =rev(c("lodacp_Siggrp_V","lodacp_Siggrp_VII","lodacp_Siggrp_VIII","lodacp_Siggrp_XIII",
                                                                    "lodacp_Siggrp_XIV","lodacp_Siggrp_XV","lodacp_Siggrp_IX","lodacp_Siggrp_XII")),
                                                      labels=rev(c("V","VII","VIII","XIII","XIV","XV","IX","XII")),position = "right")+
  theme(axis.text.y =element_text(size=rel(1.5)),axis.text.x =element_text(hjust=0.5,angle=0, size=rel(1.5)),
        axis.title.x = element_blank(),axis.title.y = element_blank())+
  ggtitle("Relative enrichment of signatures (on mean log data cpm) by condition \n manual reordering (r+c)")+coord_flip()

dev.off()

