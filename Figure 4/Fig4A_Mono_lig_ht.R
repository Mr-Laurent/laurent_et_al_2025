library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(seriation)
library(dplyr)
library(Matrix)

reg=1e-6

load("./Grouped_objects/LigRec_MoMac_lineages_cytok_fev25.rd")
load("./Grouped_objects/all_cytokines_liste.rd")
load("./Grouped_objects/umitab_all_maj_lin_fev25.rd") # Sparse matrix of counts (row= individual cell, column= gene)
load("./Grouped_objects/dfall_IBD_all_df3CD.rd")      # Metadata by cell (name, metacell, disease, celltype ...)

#Adapt for Mono / Macs :
dfall_all_df3CD$Group_v2 <- gsub("\\[MoMac\\] (Mono\\d+)", "[Mono] \\1", dfall_all_df3CD$Group_v2)
dfall_all_df3CD$Group_v2 <- gsub("\\[MoMac\\] (Macro\\d+)", "[Macro] \\1", dfall_all_df3CD$Group_v2)  

# Compute mean expression by gene in each subtypes
ctyp_mapping <- dfall_all_df3CD$Group_v2[match(rownames(umitab_all_maj_lin), dfall_all_df3CD$names)]
ctyp_factors <- factor(ctyp_mapping)
ctyp_means <- do.call(rbind, by(seq_len(nrow(umitab_all_maj_lin)), ctyp_factors, function(idx) {
  Matrix::colMeans(as.matrix(umitab_all_maj_lin[idx, , drop = FALSE]))
}))



#List of ligands to plot
mono_lig_cytok<-c(as.vector( unique(interaction_stats$Ligand[interaction_stats$subtype1=="Mono"]) ),
                  "CSF1","CSF2","CCL7","CLU","CD274","CXCL12","EBI3","APOE","IGF1","PDGFB","EREG","ADM","VEGFA","PLAU","IL1RN",
                  "INHBA","PDGFC","VEGFB","NRG1","HBEGF","AREG","RNASET2","LRPAP1","TGFB1","SEMA4A","APLP2","THBS1","CXCL14","SEMA4D")

mo_ctyp=c("[Mono] Mono5","[Mono] Mono4","[Mono] Mono3","[Mono] Mono1")

m_all_ctyp<-as.matrix(as.data.frame(t(ctyp_means)))
m_lig_ctyp_mm_mo<-m_all_ctyp[mono_lig_cytok,mo_ctyp]
#save(m_lig_ctyp_mm_mo,file="m_lig_ctyp_mm_mo_fev25.rd")



###---###  USE  TGLKMEANS, ONLY ON LINUX, TO REORDER GENES BASED ON THEIR EXPRESSION IN MONO 1 AND 3 ONLY
# library(tglkmeans)
# load("m_lig_ctyp_mm_mo_fev25.rd")
# reg=1e-6
# M_LR_GRP<-m_lig_ctyp_mm_mo
# M_LR_GRP_normed<-(log2((reg+(M_LR_GRP))/(reg+rowMeans(M_LR_GRP))))
#
# mo1mo3<-M_LR_GRP_normed[,c("[Mono] Mono3","[Mono] Mono1")]
# 
# m_lig_ctyp_mm_k15 <- TGL_kmeans_tidy(M_LR_GRP_normed,k = 15,metric = "euclid",verbose = TRUE,seed = 42)
# mo1mo3_k15 <- TGL_kmeans_tidy(mo1mo3,k = 15,metric = "euclid",verbose = TRUE,seed = 42)
# save(m_lig_ctyp_mm_k15,file="m_lig_ctyp_mm_mo_k15_fev25.rd")
# save(mo1mo3_k15,file="mo1mo3_mo_k15_fev25.rd")
###---###  
load("./Grouped_objects/mo1mo3_mo_k15_fev25.rd")


makezsc<-function(M_LR_GRP){
  M_LR_GRP_normed<- t(apply(m_lig_ctyp_mm_mo,1,scale)) 
  colnames(M_LR_GRP_normed)<-colnames(M_LR_GRP)
  cormat_genes_M_LR_GRP=cor(t(M_LR_GRP_normed)) 
  M_LR_GRP_ordgw=rownames(M_LR_GRP_normed)[get_order(seriate(as.dist(1-cormat_genes_M_LR_GRP),method = "GW_complete")) ]  # Using the order of the leaf nodes in a dendrogram obtained by hierarchical clustering and reordered by the Gruvaeus and Wainer
  M_LR_GRP_melt<-reshape::melt(M_LR_GRP_normed) 
  
  setClass(Class="Myobj",
           representation(
             ordgw="character",
             melt="data.frame"
           ))
  return(new("Myobj",ordgw=M_LR_GRP_ordgw,melt=M_LR_GRP_melt))
}


# m_lig_ctyp_mm_mo has the mean expression of each gene for each monocyte subtype
# Then, makezsc() compute the enrichment of a gene in a subtype relative to the whole Mono subset 
#  It uses a z-score to scale and center the mean of the given subtype compared to the other means of the other mono subtypes
z_lig_ct_mm_mo<-makezsc(m_lig_ctyp_mm_mo)
colnames(z_lig_ct_mm_mo@melt)<-c("Ligand","subtype","value")

# First look at the heatmap reordered by tglkmeans :
ord_tglk_lig_ct_mm_mo<-mo1mo3_k15$cluster$id[order(mo1mo3_k15$cluster$clust)]
ggplot(z_lig_ct_mm_mo@melt, aes(x = Ligand, y = subtype, fill = value)) +
  geom_tile() + scale_fill_gradientn(colors = rev(brewer.pal(n = 8, name = "RdBu")), limits=c(-2, 2), oob=squish, name = "z-score of mean \nacross Mono means") +
  scale_y_discrete(limits=mo_ctyp )+ggtitle("Reorder tglk on mono1/mono3")+
  scale_x_discrete(limits=ord_tglk_lig_ct_mm_mo)+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3, size=rel(0.7)),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"),
                     legend.position = "right")+ geom_hline(yintercept = 0.5 + 0:length(unique(z_lig_ct_mm_mo@melt$subtype)), colour = "black", size = 0.05) 

# From this first clustering of cytokine expression based on Mono1 and Mono3 expression, I manually rearranged some genes for a smoother figure
genes_lst<-"CXCL10,SEMA4A,LTB,EBI3,CCL19,CD274,INHBA,SPP1,CCL2,CXCL5,CLU,CXCL9,CXCL1,CSF1,CXCL11,IL6,TNFSF14,CSF3,IL1A,TNFSF15,CSF2,IL1RN,CCL5,CCL22,CXCL8,TNF,CLCF1,CCL3,CCL4,CCL20,IL23A,IL1B,CCL4L2,CCL7,EREG,ADM,VEGFA,PDGFB,IL10,TNFSF12,HBEGF,APLP2,THBS1,SEMA4D,CXCL14,IL15,OSM,TNFSF13,TNFSF13B,AREG,CXCL16,PLAU,TGFB1,TNFSF10,NRG1,RNASET2,LRPAP1,IL16,IL18,PDGFC,VEGFB,CXCL12,IGF1,CCL18,APOE"


pdf("./Figures_print/Fig4A_Mono_lig_ht.pdf",width = 14, height = 5)
ggplot(z_lig_ct_mm_mo@melt, aes(x = Ligand, y = subtype, fill = value)) +
  geom_tile() + scale_fill_gradientn(colors = rev(brewer.pal(n = 8, name = "RdBu")), limits=c(-1.5, 1.5), oob=squish, name = "z-score of mean \nacross Mono means") +
  scale_y_discrete(limits=rev(mo_ctyp) )+ggtitle("tglk+manual reordering on all mono")+
  scale_x_discrete(limits=unlist(strsplit(genes_lst,split=",")))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3, size=rel(1.5)),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"),
                     legend.position = "right")+ geom_hline(yintercept = 0.5 + 0:length(unique(z_lig_ct_mm_mo@melt$subtype)), colour = "black", size = 0.05) 

dev.off()
