reticulate::use_miniconda('r-reticulate')
library(anndata)
library(heatmaply)
library(ComplexHeatmap)
library(plotly)

findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}
set.seed(42)
# Loading the metacell object will all MNP (MoMacs + DC + few conta)
mdata<-anndata::read_h5ad("./Grouped_objects/all__16nov22_delStro.h5ad") 
load("./Grouped_objects/mnp_HM_z_mymod_150.rd")  

remvlist<-which(colSums(mdata$X)==0) 
mdatarm<-mdata$X[,-remvlist] 
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))

#Normalise data before plotting
dgcmtxcounts<- as(t(mdatarm), "dgCMatrix")
data.cpm<-1e+06*sweep(dgcmtxcounts,2,colSums(as.matrix(dgcmtxcounts)),`/`)
preloco<-as.matrix(log2(1+(data.cpm/10)))
gnlist<-c("PTPRC,CD3D,CD7,MS4A1,PDGFRB,LUM,TNFRSF17,HLA-DRA,LYZ,C5AR1,S100A9,S100A12,FCGR3A,CD209,MAF,FOLR2,FCAR,MERTK,PLD3,DAB2,CD14,VSIG4,SELENOP,
          C1QA,C1QB,C1QC,TLR7,PHEX,GZMB,TCF4,IL3RA,LILRB4,FLT3,LY75,ZBTB46,LAMP3,CD200,CD274,CD80,CLEC9A,CADM1,XCR1,CCR6,BATF3,CD1C,FCER1A,CLEC10A")
gnshow<-unlist(strsplit(gnlist,split=","))                                                                          
gnshow<-gsub(" ","",gnshow)                                                                                                 
notfound<-setdiff(gnshow,rownames(dgcmtxcounts)) 


# Cluster cells with ComplexHeatmap hierarchical clustering then plot the genes using the same metacell order
set.seed(42)
# The column_km always makes different results, need to fix the seed + use "draw" to have a fixed heatmap
ht<-draw(Heatmap(t(hm_zmod), column_title = "Module clustering based on corrected z score", name = "mat",
                 row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod)))),show_row_dend = TRUE,
                 row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                 column_dend_reorder = F, cluster_columns = T,
                 column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km = 40,row_km = 0, use_raster = T))

metacells_kclust<- column_order(ht)
metacells_kclustrow<- row_order(ht)

htmxord2<-as.matrix(t(hm_zmod))[ as.vector(unlist(metacells_kclustrow)) , as.vector(unlist(metacells_kclust)) ]

loco<-preloco[setdiff(gnshow,notfound),colnames(htmxord2)]


plot<-heatmaply(loco, fontsize_row = rel(1),fontsize_col = 7,dendrogram = "none",Rowv=T,Colv=F,legendgroup="2nd",coloraxis = 'coloraxis2',
          scale_fill_gradient_fun =scale_fill_gradientn(colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis"))) 
plotly::save_image(plot, file="./Figures_print/FigS2A_Truthplot_MNP_enriched.pdf") 

