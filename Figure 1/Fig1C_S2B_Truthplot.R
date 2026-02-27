library(Matrix)
library(ComplexHeatmap)
library(dplyr)
library(tidyverse)
library(viridis)

# Load data
load("./Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v3.rd")
load("./Grouped_objects/Macs_dgCmc_16nov22_counts.rd")  
load("./Grouped_objects/genelist_16nov22_Macs.rd") 
load("./Grouped_objects/ht_Enrich_Macs.rd") 

pdf(file="./Figures_print/FigS2B_HM_modules_betterInApp.pdf",width = 8,height = 10)
print(ht)
dev.off()
#---------------------------------------#

colnames(genelist)<-c("Mods","Gens")
hm_zmod<-as.matrix(metadf_z[,c(154:303)])
df_zmod2<-as.data.frame(metadf_z[,c(2:3,154:303)])
colnames(df_zmod2)<-sub("([a-zA-Z_]+)(\\d+)([a-zA-Z]+)", "Module \\2", colnames(df_zmod2))  
preloco_nc<-as.matrix(log(1+dgcmtxcounts))# DU coup les matrices de genes sont pas corrigees, on peut pas
allsum<-Matrix::colSums(preloco_nc) 
row_mod_ord<- as.vector(unlist(row_order(ht)))
data.cpm<-1e+06*sweep(dgcmtxcounts,2,colSums(dgcmtxcounts),`/`)  # if you don't want another package
preloco<-as.matrix(log2(1+(data.cpm/10)))  # Divide by 10 the CPM for more readability

###----------------------------------------------------------### 
#### Function to transform list of genes to vector readable #### 
#### by the program, with exclusion of unknown genes        ####
###----------------------------------------------------------### 
gene_input2vect<-function(genelist,dgcmtx){
  gnshow<-unlist(strsplit(genelist,split=","))
  gnshow<-gsub(" ","",gnshow)
  notfound<-setdiff(gnshow,rownames(dgcmtx))
  print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  return(setdiff(gnshow,notfound) )
}
###----------------------------------------------------------### 
findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}      
###----------------------------------------------------------### 
gnlist<-"CXCL8,SERPINB2,MFSD2A,F3,IL1A,TNF,IL7R,IDO1,INHBA,CD274,IL3RA,P2RX7,SLAMF1,CRIM1,SLC39A8,SLC7A11,ADAM19,SPP1,DNAJB4,HSPA6,FCGR2B,CA2,SCD,BNIP3,RRAD,SLC2A1,IFI44L,CXCL10,CXCL9,IFIT2,OAS2,CXCL11,IL23A,CD40,ICAM4,SCN1B,SLAMF7,CCL20,CXCL3,CXCL2,S100A8,S100A9,VCAN,FCN1,RIPOR2,S100A12,APOBEC3A,SELL,FCAR,ALOX5AP,RETN,MCEMP1,OSM,ATP5PB,MRPS21,RNF7,NDUFB5,ATP5ME,C1QA,VSIG4,DNASE1L3,DAB2,PDK4,SLCO2B1,CMKLR1,CD209,FUCA1,MERTK,SELENOP,FOLR2,IGF1,AXL,LYVE1"

# Cut into highlighted gene by programs for easier manipulations:
gn_XIV<-c("CXCL8,SERPINB2,MFSD2A,F3,IL1A,TNF")
gn_XV<-c("IL7R,IDO1,INHBA,CD274,IL3RA,P2RX7,SLAMF1,CRIM1,SLC39A8,SLC7A11,ADAM19,SPP1")
gn_XVII<-c("DNAJB4,HSPA6,FCGR2B,CA2,SCD,BNIP3,RRAD,SLC2A1")
gn_IX<-c("IFI44L,CXCL10,CXCL9,IFIT2,OAS2,CXCL11")
gn_XIII<-c("IL23A,CD40,SLAMF7,CCL20,CXCL3,CXCL2")
gn_VIII<-c("S100A8,S100A9,VCAN,FCN1,RIPOR2,S100A12,APOBEC3A")
gn_VII<-c("SELL,FCAR,ALOX5AP,RETN,MCEMP1,OSM")
gn_XII<-c("IL10,CDK1,RAB7B,TPRA1,LY6K,BCL11A,PHLDB1,MMP19,EREG")
#gn_I<-c("CTSD,MS4A4A,MS4A6A,HRH2,ITGA4,CD36,ATOX1,CPVL,BST2")
gn_II<-c("ATP5PB,MRPS21,RNF7,NDUFB5,ATP5ME")
#gn_IV<-c("MNDA,VAMP5,LY96,C3AR1,CALM2")
gn_VI<-c("CSF1R,CD72,CD101,IL18")
gn_V<-c("C1QA,VSIG4,DNASE1L3,DAB2,PDK4,SLCO2B1,CMKLR1,CD209,FUCA1,MERTK,SELENOP,FOLR2,IGF1,AXL,LYVE1")

gnlist<-paste(gn_XIV,gn_XV,gn_XVII,gn_IX,gn_XIII,gn_VIII,gn_VII,gn_XII,gn_II,gn_VI,gn_V,sep=",")



ann_df<-data.frame(metacell=rownames(metadf_z),annotation="undefined")
ann_clu="6,40,2,33,39,38,37,34,35,36"
ann_clu_show<-unlist(strsplit(ann_clu,split=","))
ann_clu_show<-gsub(" ","",ann_clu_show)
mc_part1<-rownames(hm_zmod)[as.vector(unlist(column_order(ht)[as.numeric(ann_clu_show)])) ]
set.seed(42)
mc_part1_rdm<-sample(mc_part1)
ann_clu2="3,9,11,25,4,13,7,8,5,10,24,26,27,28,12,14,16,18,1,19,20,23,29,30,32,15,17,21,22,31"
ann_clu_show2<-unlist(strsplit(ann_clu2,split=","))
ann_clu_show2<-gsub(" ","",ann_clu_show2)
mc_part2<-rownames(hm_zmod)[as.vector(unlist(column_order(ht)[as.numeric(ann_clu_show2)])) ]
mc_ordhm<-append(mc_part1_rdm,mc_part2)


htmxord<-t(hm_zmod)[as.vector(unlist(row_order(ht))),mc_ordhm]
gnshow<-unlist(strsplit(gnlist,split=","))                                                                          
gnshow<-gsub(" ","",gnshow)                                                                                                 
notfound<-setdiff(gnshow,rownames(dgcmtxcounts)) 

# Normalized metacell count matrix with genes to show
loco<-preloco[setdiff(gnshow,notfound),colnames(htmxord)]

# Adapt the matrix to a format readable by ggplot
lim_ord=colnames(htmxord)
hmgg2<-loco%>% 
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  pivot_longer(-c(genes), names_to = "metacells", values_to = "counts")

reordplot<-ggplot(hmgg2,aes(x=metacells, y=genes, fill=counts)) + 
  geom_raster() + scale_x_discrete(limits=lim_ord)+ 
  scale_y_discrete(limits=rev(rownames(loco)))+
  scale_fill_gradientn(colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = rel(0.7)) )

pdf(file="./Figures_print/Fig1C_Truthplot_EnrichMomac.pdf",width = 8,height = 7.5)
print(reordplot)
dev.off()

#Save as picture aswell
ggsave(
  "./Figures_print/Fig1C_Truthplot_EnrichMomac.tiff",
  plot = print(reordplot),
  width = 8,
  height = 7.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)
