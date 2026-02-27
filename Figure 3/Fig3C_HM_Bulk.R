library(ComplexHeatmap)
library(ggplot2)
library(scales)
library(reshape)
library(reshape2)
library(dplyr)
library(Rfast)


load("./Grouped_objects/metadata_mac20_1_2.rd")
load("./Grouped_objects/data_counts_mac20_1_2.rd")
load("./Grouped_objects/data_logcpm_mac20_1_3.rd")

#### If personnal lists to add :
load("./Grouped_objects/Modsig_MoMacsEnr.rd")


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
for(Mod_score in names(Modsig)[1:17]){
  print(Mod_score)
  geneshow<-unlist(Modsig[Mod_score])
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
ord_sig_perso2_log<-paste0("lodacp_",c("Sig_X","Sig_XI","Sig_XII","Sig_I","Sig_III","Sig_IV","Sig_VI",
                                       "Sig_II","Sig_XVI","Sig_XVII","Sig_XV","Sig_XIV","Sig_XIII","Sig_IX","Sig_VIII","Sig_VII","Sig_V") )

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


pdf("./Figures_print/Fig3C_HM_Bulk.pdf",height = 8,width = 8)

# plot_reorderselectedrow_selectedcol2+scale_y_discrete(limits =rev(c("lodacp_Sig_V","lodacp_Sig_VII","lodacp_Sig_VIII","lodacp_Sig_XIII",
#                                                                     "lodacp_Sig_XIV","lodacp_Sig_XV","lodacp_Sig_IX","lodacp_Sig_XII")),
#                                                       labels=rev(c("V","VII","VIII","XIII","XIV","XV","IX","XII")),position = "right")+
#   theme(axis.text.y =element_text(size=rel(1.5)),axis.text.x =element_text(hjust=0.5,angle=0, size=rel(1.5)),
#         axis.title.x = element_blank(),axis.title.y = element_blank())+
#   ggtitle("Relative enrichment of signatures (on mean log data cpm) by condition \n manual reordering (r+c)")+coord_flip()
# 
# plot_reorderselectedrow_selectedcol2+scale_y_discrete(limits =rev(c("lodacp_Sig_VI","lodacp_Sig_V","lodacp_Sig_VII","lodacp_Sig_VIII","lodacp_Sig_XIII",
#                                                                     "lodacp_Sig_XIV","lodacp_Sig_XV","lodacp_Sig_IX","lodacp_Sig_XII")),
#                                                       labels=rev(c("VI 'Transition C'","V 'Macrophage'","VII 'MonoA'","VIII 'MonoB'","XIII 'Inflammatory A'","XIV 'Inflammatory B'","XV 'IFIM program'","IX 'Interferon'","XII 'Immunoregulation/Repair'")),position = "right")+
#   theme(axis.text.y =element_text(size=rel(1.5)),axis.text.x =element_text(hjust=0.5,angle=0, size=rel(1.5)),
#         axis.title.x = element_blank(),axis.title.y = element_blank())+scale_x_discrete(limits = (rownames(metad3)))+
#   ggtitle("Relative enrichment of signatures (on mean log data cpm) by condition \n manual reordering (r+c)")+coord_flip()

plot_reorderselectedrow_selectedcol2+
  scale_y_discrete(limits =rev(c("lodacp_Sig_VII","lodacp_Sig_VIII","lodacp_Sig_XIII","lodacp_Sig_XIV","lodacp_Sig_IX",
                             "lodacp_Sig_XV","lodacp_Sig_XII","lodacp_Sig_VI","lodacp_Sig_V")),
                   labels=rev(c("VII 'MonoA'","VIII 'MonoB'","XIII 'Inflammatory A'","XIV 'Inflammatory B'","IX 'Interferon'",
                            "XV 'IFIM program'","XII 'Immunoregulation/Repair'","VI 'Transition C'","V 'Macrophage'")),position = "right")+
  theme(axis.text.y =element_text(size=rel(1.5)),axis.text.x =element_text(hjust=1,vjust=0.5,angle=90, size=rel(1.5)),
        axis.title.x = element_blank(),axis.title.y = element_blank())+scale_x_discrete(limits = (rownames(metad3)))+
  ggtitle("Relative enrichment of signatures (on mean log data cpm) by condition \n manual reordering (r+c)")#+coord_flip()

dev.off()

