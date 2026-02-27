## -- ##

## Before running this code, use the CodeProp_XXXX_IBD_complet.R to create all cell distributions 

## ---------------------------------------------------------- ##


create_models<-function(dfall_lin,clean_lin,lin_name){ 
  
  cell_to_cluster=dfall_lin@df3CD$Group
  names(cell_to_cluster)=dfall_lin@df3CD$names
  cell_to_sample=dfall_lin@df3CD$sample
  names(cell_to_sample)=dfall_lin@df3CD$names
  
  u_clusters=unique(cell_to_cluster)
  u_samples=unique(cell_to_sample)
  dataset<-c()
  dataset$counts<-array(0,dim=c(length(u_samples),length(all_genes),length(u_clusters)),dimnames = list(u_samples,all_genes,u_clusters))
  
  # Compute "counts" : each matrix slide (a patient) is the sum of umi counts by gene for each MoMac cluster
  for (sampi in u_samples){
    # maski=(cell_to_sample==sampi)
    # tmp_counts=as.matrix(Matrix::t(fac2sparse(cell_to_cluster[maski]) %*% Matrix::t(t(clean_lin1$X)[,maski])))  
    ### My clean matrix has more cells than the df object (that was filtered) so use this way:
    mask_tho=dfall_lin@df3CD$names[which(dfall_lin@df3CD$sample==sampi)]
    tmp_counts=as.matrix(Matrix::t(fac2sparse(cell_to_cluster[mask_tho]) %*% Matrix::t(t(clean_lin1$X)[,mask_tho])))
    dataset$counts[sampi,rownames(tmp_counts),colnames(tmp_counts)]=tmp_counts
  }
  
  models=apply(dataset$counts,2:3,sum)
  models2=t(t(models)/colSums(models)) # also possible with models2 <- sweep(models, 2, colSums(models), "/")  but it's a bit less efficient from my benchmarks
  
  ncells_per_cluster<-setNames(as.vector(table(cell_to_cluster)), names(table(cell_to_cluster)))[u_clusters]

  cluster_to_subtype<-colnames(models2)
  names(cluster_to_subtype)<-colnames(models2)

  return(list(dataset=dataset,models=models2,ncells_per_cluster=ncells_per_cluster,cluster_to_subtype=cluster_to_subtype))
}

## ---------------------------------------------------------- ##


library(dplyr)
library(seriation)
library(dendsort)
library(matrixStats)
library(datasets)
library(circlize)
library(scales)
library(Matrix)
library(Rfast)
library(tidyverse)
library(abind)
library(reticulate)
reticulate::use_miniconda('r-reticulate') # it allows me to read h5ad file without errors from R anndata
use_python("C:/Users/XXXXX/AppData/Local/r-miniconda/envs/r-reticulate/python.exe", required = TRUE)
library(anndata)
library(anndataR)

load("./Grouped_objects/LigRec_obj/all_genes.rd")

load("./Grouped_objects/LigRec_obj/dfall_IBD_MoMac.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_Mono.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_Macro.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_Bcells.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_DC.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_Fibro.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_PC.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_Tcells.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_EndoC.rd")
load("./Grouped_objects/LigRec_obj/dfall_IBD_GliPerSM_detailed.rd")

setClass(Class="Myobj2", representation( df3CD="data.frame" ))
df3CD_o1<-c("190","196","193","158","209","69","n7","n1","n3","n17","n47","n6","n33","n37","128","138","187","181","n13") 

dfall_MoMac@df3CD$Group_v2<-paste0("[MoMac] ",dfall_MoMac@df3CD$Group)
dfall_Mono@df3CD$Group_v2<-paste0("[Mono] ",dfall_Mono@df3CD$Group)
dfall_Macro@df3CD$Group_v2<-paste0("[Macro] ",dfall_Macro@df3CD$Group)
dfall_Bcells@df3CD$Group_v2<-paste0("[Bcells] ",dfall_Bcells@df3CD$Group)
dfall_DC@df3CD$Group_v2<-paste0("[DC] ",dfall_DC@df3CD$Group)
dfall_Fibro@df3CD$Group_v2<-paste0("[Fibro] ",dfall_Fibro@df3CD$Group)
dfall_PC@df3CD$Group_v2<-paste0("[PC] ",dfall_PC@df3CD$Group)
dfall_Tcells@df3CD$Group_v2<-paste0("[Tcells] ",dfall_Tcells@df3CD$Group)
dfall_EndoC@df3CD$Group_v2<-paste0("[EndoC] ",dfall_EndoC@df3CD$Group)
dfall_GliPerSM@df3CD$Group_v2<-paste0("[GliPerSM] ",dfall_GliPerSM@df3CD$Group)

dfall_all_df3CD<-rbind(dfall_MoMac@df3CD,dfall_Bcells@df3CD,dfall_DC@df3CD,dfall_Fibro@df3CD,dfall_PC@df3CD,dfall_Tcells@df3CD,dfall_EndoC@df3CD,dfall_GliPerSM@df3CD)

save(dfall_all_df3CD,file="./Grouped_objects/LigRec_obj/dfall_IBD_all_df3CD.rd")

#MoMacs:
load("./Grouped_objects/LigRec_obj/df_mcall_Macs_select2.rd")
clean_lin1<-anndata::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totalldec_MacsfromMNP_select2.h5ad')
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_momac<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="MoMac")
MoMac_obj<-create_models(dfall_lin=dfall_MoMac,clean_lin=clean_lin1,lin_name="MoMac")
colnames(MoMac_obj$models)<-paste0("[MoMac] ",colnames(MoMac_obj$models))
save(MoMac_obj, file="./Grouped_objects/LigRec_obj/run_objects/MoMac_obj.rd")
#Monos:
Mono_obj<-create_models(dfall_lin=dfall_Mono,clean_lin=clean_lin1,lin_name="Mono")
colnames(Mono_obj$models)<-paste0("[Mono] ",colnames(Mono_obj$models))
save(Mono_obj, file="./Grouped_objects/LigRec_obj/run_objects/Mono_obj.rd")
#Macros:
Macro_obj<-create_models(dfall_lin=dfall_Macro,clean_lin=clean_lin1,lin_name="Macro")
colnames(Macro_obj$models)<-paste0("[Macro] ",colnames(Macro_obj$models))
save(Macro_obj, file="./Grouped_objects/LigRec_obj/run_objects/Macro_obj.rd")

#Bcells:
load("./Grouped_objects/LigRec_obj/df_mcall_Totdec_B.rd")
clean_lin1<-anndata::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totalldec_B.h5ad')
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_NNN<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="NNN")
Bcells_obj<-create_models(dfall_lin=dfall_Bcells,clean_lin=clean_lin1,lin_name="Bcells")
colnames(Bcells_obj$models)<-paste0("[Bcells] ",colnames(Bcells_obj$models))
save(Bcells_obj, file="./Grouped_objects/LigRec_obj/run_objects/Bcells_obj.rd")

#DC:
load("./Grouped_objects/LigRec_obj/df_mcall_DC_select2.rd")
clean_lin1<-anndata::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totalldec_DCfromMNP_select2.h5ad')
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_dc<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="DC")
DC_obj<-create_models(dfall_lin=dfall_DC,clean_lin=clean_lin1,lin_name="DC")
colnames(DC_obj$models)<-paste0("[DC] ",colnames(DC_obj$models))
save(DC_obj, file="./Grouped_objects/LigRec_obj/run_objects/DC_obj.rd")

#Fibro:
load("./Grouped_objects/LigRec_obj/df_mcall_Totdec_StrGli_Fibro.rd")
clean_lin1<-anndata::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totdec_StrGli_Fibro.h5ad')
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_fibro<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="Fibro")
Fibro_obj<-create_models(dfall_lin=dfall_Fibro,clean_lin=clean_lin1,lin_name="Fibro")
colnames(Fibro_obj$models)<-paste0("[Fibro] ",colnames(Fibro_obj$models))
save(Fibro_obj, file="./Grouped_objects/LigRec_obj/run_objects/Fibro_obj.rd")

#PC:
load("./Grouped_objects/LigRec_obj/df_Totdec_PlasmaCells.rd")
clean_lin1<-anndataR::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totdec_PlasmaCells.h5ad')
colnames(clean_lin1$X)<-clean_lin1$var_names
rownames(clean_lin1$X)<-clean_lin1$obs_names
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_NNN<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="NNN")
PC_obj<-create_models(dfall_lin=dfall_PC,clean_lin=clean_lin1,lin_name="PC")
colnames(PC_obj$models)<-paste0("[PC] ",colnames(PC_obj$models))
save(PC_obj, file="./Grouped_objects/LigRec_obj/run_objects/PC_obj.rd")

#Tcells:
load("./Grouped_objects/LigRec_obj/SubsetsAppObj/df_Totdec_Tcellsmini.rd")
clean_lin1<-anndataR::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totdec_mini_Tcells_oct24redo.h5ad')
colnames(clean_lin1$X)<-clean_lin1$var_names
rownames(clean_lin1$X)<-clean_lin1$obs_names
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_NNN<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="NNN")
Tcells_obj<-create_models(dfall_lin=dfall_Tcells,clean_lin=clean_lin1,lin_name="Tcells")
colnames(Tcells_obj$models)<-paste0("[Tcells] ",colnames(Tcells_obj$models))
save(Tcells_obj, file="./Grouped_objects/LigRec_obj/run_objects/Tcells_obj.rd")

#EndoC:
load("./Grouped_objects/LigRec_obj/df_Totdec_StrGli_EndoC.rd")
clean_lin1<-anndataR::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totdec_StrGli_EndoC.h5ad')
colnames(clean_lin1$X)<-clean_lin1$var_names
rownames(clean_lin1$X)<-clean_lin1$obs_names
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_NNN<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="NNN")
EndoC_obj<-create_models(dfall_lin=dfall_EndoC,clean_lin=clean_lin1,lin_name="EndoC")
colnames(EndoC_obj$models)<-paste0("[EndoC] ",colnames(EndoC_obj$models))
save(EndoC_obj, file="./Grouped_objects/LigRec_obj/run_objects/EndoC_obj.rd")

#GliPerSM:
load(paste0(googlepath,"SubsetsAppObj/df_Totdec_StrGli_GliPerSM.rd"))
clean_lin1<-anndataR::read_h5ad('./Grouped_objects/h5ad_files/all_clean_Totdec_StrGli_GliPerSM.h5ad')
colnames(clean_lin1$X)<-clean_lin1$var_names
rownames(clean_lin1$X)<-clean_lin1$obs_names
df_lin1_cd_inf<-df[which(df$sample%in%df3CD_o1),]
# ncm_NNN<-mean_pat_fct(cl_lin=clean_lin1,df_lin=df_lin1_cd_inf,lin_choice="NNN")
GliPerSM_obj<-create_models(dfall_lin=dfall_GliPerSM,clean_lin=clean_lin1,lin_name="GliPerSM")
colnames(GliPerSM_obj$models)<-paste0("[GliPerSM] ",colnames(GliPerSM_obj$models))
save(GliPerSM_obj, file="./Grouped_objects/LigRec_obj/run_objects/GliPerSM_detailed_obj.rd")


#################################

library(dplyr)
library(seriation)
library(dendsort)
library(matrixStats)
library(datasets)
library(circlize)
library(scales)
library(Matrix)
library(Rfast)
library(tidyverse)
library(anndata)
library(anndataR)
library(abind)



load("./Grouped_objects/LigRec_obj/run_objects_amasauce/all_genes.rd")
### Validated pairs from CellphoneDB
# https://github.com/ventolab/CellphoneDB/blob/dc8abd15b24e1d48c8862d54e8122603630a68ed/NatureProtocols2024_case_studies/v5.0.0/interaction_input.csv
tab_CpDB<-read.csv2(file="https://github.com/ventolab/CellphoneDB/raw/dc8abd15b24e1d48c8862d54e8122603630a68ed/NatureProtocols2024_case_studies/v5.0.0/interaction_input.csv",header = T,sep = ",")

tab_CpDB_LigRec<-tab_CpDB[which(tab_CpDB$directionality=="Ligand-Receptor"),]  
# 2508 of the 2911 interactions are Ligand-Receptor

## is_ppi : It's for protein-protein interactions
# ex. of a non ppi : Adrenaline_byPNMT_and_SLC18A2     (PNMT modifies Noradrenaline in adrenaline, so it's not direct interaction)
tab_CpDB_LigRec_ppi<-tab_CpDB_LigRec[which(tab_CpDB_LigRec$is_ppi=="True"),]  
# 1475 interactions left


## Adapt to R from https://github.com/saezlab/liana-py/issues/60
# replace + with _
tab_CpDB_int<-tab_CpDB_LigRec_ppi$interactors <- gsub('\\+', '_', tab_CpDB_LigRec_ppi$interactors)

#Replace the last "-" for a "&" as it's the limit between ligand name and receptor
morehyphen<-sapply(tab_CpDB_int, function(x) {
  sum(gregexpr("-", x)[[1]] > 0) > 1 # Count hyphens and check if more than 1
})
tab_CpDB_int[morehyphen]
# [1] "HLA-A-KIR3DL1"     "HLA-B-KIR3DL2"     "HLA-C-KIR2DL1"     "HLA-C-KIR2DL3"     "HLA-E-KLRC1"       "HLA-E-KLRC1_KLRD1" "HLA-E-KLRC2"      
# [8] "HLA-E-KLRC2_KLRD1" "HLA-E-KLRC3_KLRD1" "HLA-E-KLRK1"       "HLA-F-KIR3DL1"     "HLA-F-KIR3DL2"     "HLA-F-KIR3DS1"     "HLA-F-LILRB1"     
# [15] "HLA-F-LILRB2"      "HLA-G-LILRB1"      "HLA-G-LILRB2"      "MT-RNR2-FPR2"      "MT-RNR2-FPR3"      "ERVH48-1-SLC1A5"   "HLA-E-VSIR"       
# [22] "HLA-F-VSIR"  
tab_CpDB_int<-gsub("(.*?)-(.*?)$", "\\1&\\2", tab_CpDB_int)
split_interactors<-strsplit(tab_CpDB_int, '&')
sp_int<-do.call(rbind, split_interactors) # Make it into 2 columns instead of a list

expand_ligands_receptors <- function(data) {
  expanded_rows <- list()
  for (i in 1:nrow(data)) {
    # split ligands and receptors
    ligands <- unlist(strsplit(data[i, 1], "_"))  
    receptors <- unlist(strsplit(data[i, 2], "_")) 
    # Redo all ligand receptor combinations
    for (ligand in ligands) {
      for (receptor in receptors) {
        expanded_rows[[length(expanded_rows) + 1]] <- c(ligand, receptor)
      }
    }
  }
  return(do.call(rbind, expanded_rows))
}

# Expand the resource data
sp_int_xpd <- expand_ligands_receptors(sp_int)
df_sp_int<-as.data.frame(sp_int_xpd, stringsAsFactors = FALSE)
df_sp_int<-unique(df_sp_int)
colnames(df_sp_int)<-c("ligand","receptor")
# 1484 unique ligand/receptor pairs kept

# For information,some L/R had multiple L or R due to complex formign to bind one or the other
#subset on existing genes:
df_sp_int<-df_sp_int[df_sp_int$ligand%in%all_genes,]
df_sp_int<-df_sp_int[df_sp_int$receptor%in%all_genes,]
# 1372 L/R left

validated_pairs_CpDB<-df_sp_int

save(validated_pairs_CpDB, file="./Grouped_objects/LigRec_obj/run_objects/validated_pairs_CpDB.rd")



# Taken from LigRec_amasauce_v4_CpDB_MoMac_vsLin_seuil_0umi_25cell_pval.R

library(dplyr)
library(seriation)
library(dendsort)
library(matrixStats)
library(datasets)
library(circlize)
library(scales)
library(Matrix)
library(Rfast)
library(tidyverse)
library(anndata)
library(anndataR)
library(abind)


load("./Grouped_objects/LigRec_obj/run_objects/Mono_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/Macro_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/Bcells_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/DC_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/Fibro_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/PC_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/Tcells_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/EndoC_obj.rd")
load("./Grouped_objects/LigRec_obj/run_objects/GliPerSM_detailed_obj.rd")

df3CD_o1<-c("190","196","193","158","209","69","n7","n1","n3","n17","n47","n6","n33","n37","128","138","187","181","n13")  #"n49" is not in al the matrices, so we remove it
load("./Grouped_objects/LigRec_obj/all_genes.rd")

000000-------00-----000000----00000000--------0000----
00----00----0000----00----00-----00---------00----00--
00----00---00--00---00----00-----00-------------000---
000000-----000000---000000-------00-----00--------00--
00--------00----00--00----00-----00---------00----00--
00--------00----00--00----00-----00-----------0000----0

# Our Mono1 enriched samples are the following:
df3CD_mono1<-c("n17","n47","n33","n37","128","138","187","181","n13")
samps1=df3CD_mono1
samps2=setdiff(df3CD_o1,df3CD_mono1)
reg=1e-6
exprs_thresh=5e-6
samps=df3CD_o1
lineage_merge=T

# lineage_merge=F if I do it in 1 lineage only 
subtype_models=cbind(Mono_obj$models,Macro_obj$models,Bcells_obj$models,DC_obj$models,Fibro_obj$models,PC_obj$models,Tcells_obj$models,EndoC_obj$models,GliPerSM_obj$models)

cluster_to_subtype1<-colnames(subtype_models)
names(cluster_to_subtype1)<-colnames(subtype_models)

cluster_to_broader_subtype<-cluster_to_subtype1
#cluster_to_broader_subtype[grep("\\[MoMac\\]",cluster_to_broader_subtype)]<-"MoMacs"
cluster_to_broader_subtype[grep("\\[Mono\\]",cluster_to_broader_subtype)]<-"Mono"
cluster_to_broader_subtype[grep("\\[Macro\\]",cluster_to_broader_subtype)]<-"Macro"
cluster_to_broader_subtype[grep("\\[Bcells\\]",cluster_to_broader_subtype)]<-"Bcells"
cluster_to_broader_subtype[grep("\\[DC\\]",cluster_to_broader_subtype)]<-"DCs"
cluster_to_broader_subtype[grep("\\[Fibro\\]",cluster_to_broader_subtype)]<-"Fibros"
cluster_to_broader_subtype[grep("\\[PC\\]",cluster_to_broader_subtype)]<-"PCs"
cluster_to_broader_subtype[grep("\\[Tcells\\]",cluster_to_broader_subtype)]<-"Tcells"
cluster_to_broader_subtype[grep("\\[EndoC\\]",cluster_to_broader_subtype)]<-"EndoCs"
cluster_to_broader_subtype[grep("\\[GliPerSM\\]",cluster_to_broader_subtype)]<-"GliPerSMs"




load("./Grouped_objects/LigRec_obj/run_objects/validated_pairs_CpDB.rd")

validated_pairs<-validated_pairs_CpDB
## Keep only 1st Lig-Rec pairs from valitaded list  
cc2<-unique(validated_pairs[,1])
#cc2<-cc2[rowMaxs(subtype_models[cc2,])>exprs_thresh]
cc2<-cc2[matrixStats::rowMaxs(subtype_models[cc2,])>exprs_thresh]
rec<-unique(as.vector(validated_pairs[match(cc2,validated_pairs[,1]),2]))
# Pourquoi on suppr des lignes ? pas compris... Juste pour prendre la 1e paire qui vient I guess, et garder cells ou les Lig rec sont d?j? dans ce top ?
validated_pairs<-validated_pairs[validated_pairs[,1]%in%cc2&validated_pairs[,2]%in%rec,]
validated_pairs_matrix=matrix(NA,length(rec),length(cc2))
rownames(validated_pairs_matrix)=rec
colnames(validated_pairs_matrix)=cc2
for (i in 1:nrow(validated_pairs)){
  if (validated_pairs[i,2]%in%rec&(validated_pairs[i,1]%in%cc2)){
    validated_pairs_matrix[validated_pairs[i,2],validated_pairs[i,1]]=1
  }
}

rec_to_cc<-sapply(split(validated_pairs[,1],validated_pairs[,2]),paste,collapse=",")
cc_to_rec<-split(validated_pairs[,2],validated_pairs[,1])

combined_counts<-array(0,dim=c(length(df3CD_o1),length(all_genes),ncol(subtype_models)),
                       dimnames = list(df3CD_o1,all_genes, colnames(subtype_models) ))


## They all don't have all samples, so need to fill in blanks before combining them
for (i in df3CD_o1) {
  print(i)
  allcat=c("Mono_obj","Macro_obj","Bcells_obj","DC_obj","Fibro_obj","PC_obj","Tcells_obj","EndoC_obj","GliPerSM_obj")
  pb=txtProgressBar(width=length(allcat),min = 1,max=length(allcat),style=3)
  for(j in allcat){
    setTxtProgressBar(pb,value = which(allcat==j))
    eval(parse(text=paste0("j_o<-",j)))
    if (i %in% rownames(j_o$dataset$counts)) {
      j_counts <- j_o$dataset$counts[i, , ]
    } else {
      empty_mtx <- matrix(0, nrow = length(all_genes), ncol = length(colnames(j_o$models)))
      j_counts <- empty_mtx
    }
    colnames(j_counts)<-colnames(j_o$models)
    combined_counts[i, ,colnames(j_o$models)] <- j_counts
  }
} 


## Array of total counts from selected lineages by patients
corrected_counts<-pmax(combined_counts,0)
corrected_tot_mean=corrected_counts/apply(corrected_counts,1,sum)
n_expressed=(corrected_counts>0)[,c(validated_pairs[,1],validated_pairs[,2]),]

# Total_umi_frac doesn't have a good name, it's a mean expression of genes across subtypes by group (1 or 2)
total_umi_frac1=apply(corrected_tot_mean[samps1,,],2:3,mean)
total_umi_frac2=apply(corrected_tot_mean[samps2,,],2:3,mean)

if(lineage_merge==T){
  ############ Ignore if no subtypes
  pool_clusters=function(x){sapply(split(as.data.frame(t(x)),cluster_to_broader_subtype[colnames(x)],drop=F),colSums)}
  total_umi_frac_per_celltype1=pool_clusters(total_umi_frac1)
  total_umi_frac_per_celltype2=pool_clusters(total_umi_frac2)
}else{  ############ ignore if lineages are merged
  total_umi_frac_per_celltype1=(total_umi_frac1)
  total_umi_frac_per_celltype2=(total_umi_frac2)
}

n_expressed_pooled=array(0,dim = c(dim(n_expressed)[1:2],dim(total_umi_frac_per_celltype1)[2]),dimnames = list(dimnames(n_expressed)[[1]],dimnames(n_expressed)[[2]],dimnames(total_umi_frac_per_celltype1)[[2]]))
for (samp in dimnames(n_expressed)[[1]]){
  if(lineage_merge==T){ n_expressed_pooled[samp,,]=pool_clusters(n_expressed[samp,,])
  }else{ n_expressed_pooled[samp,,]=(n_expressed[samp,,]) }
}
n_patient_expressing1<-colSums(n_expressed_pooled[samps1,,]>0)
n_patient_expressing2<-colSums(n_expressed_pooled[samps2,,]>0)

m1=t(t(total_umi_frac_per_celltype1)/colSums(total_umi_frac_per_celltype1))
m2=t(t(total_umi_frac_per_celltype2)/colSums(total_umi_frac_per_celltype2))

total_umi_frac_per_celltype_ligands1<-total_umi_frac_per_celltype1[cc2,]
m1_rec<-m1[rec,]
total_umi_frac_per_celltype_ligands2<-total_umi_frac_per_celltype2[cc2,]
m2_rec<-m2[rec,]


000000-------00-----000000----00000000---------00-----
00----00----0000----00----00-----00-----------00------
00----00---00--00---00----00-----00----------00-------
000000-----000000---000000-------00-----00--00--00----
00--------00----00--00----00-----00---------00000000--
00--------00----00--00----00-----00-------------00----0

library(ggplot2)
library(reshape2)

load("./Grouped_objects/LigRec_obj/dfall_IBD_all_df3CD.rd")

#Adapt for Mono / Macs :

dfall_all_df3CD$Group_v2 <- gsub("\\[MoMac\\] (Mono\\d+)", "[Mono] \\1", dfall_all_df3CD$Group_v2)
dfall_all_df3CD$Group_v2 <- gsub("\\[MoMac\\] (Macro\\d+)", "[Macro] \\1", dfall_all_df3CD$Group_v2)  


load("./Grouped_objects/LigRec_obj/umitab_all_maj_lin.rd")
umitab_all_maj_lin<-Matrix::t(umitab_all_maj_lin)
tables_output_path = "./Figures_print/"

main_figures_path=tables_output_path 

get_two_subtypes_interactions_linorsub=function(subtype1,subtype2,rec_thresh=1e-5,lig_thresh=1e-6,n_expressed_thresh=2,tables_output_path,linorsub="sub"){   # Change number of atient to 2 as we have 9 vs 10 patients
  
  if(linorsub=="lin"){
    n_p_e1=n_patient_expressing1
    n_p_e2=n_patient_expressing2
    m1_rec_fct=m1_rec
    m2_rec_fct=m2_rec
    t_u_f_lig1=total_umi_frac_per_celltype_ligands1
    t_u_f_lig2=total_umi_frac_per_celltype_ligands2
  }else if(linorsub=="sub"){
    n_p_e1=n_patient_expressing_sub1
    n_p_e2=n_patient_expressing_sub2
    m1_rec_fct=m_sub1_rec
    m2_rec_fct=m_sub2_rec
    t_u_f_lig1=total_umi_frac_sub_ligands1
    t_u_f_lig2=total_umi_frac_sub_ligands2
  }else{print("linorsub not found")}
  reg=1e-14
  
  rec_crit1=m1_rec_fct[,subtype2]>rec_thresh&(n_p_e1[rownames(m1_rec_fct),subtype2]>n_expressed_thresh)&Matrix::rowSums(!is.na(validated_pairs_matrix))>0  # valid mtx > 0 is useless here as it's the same as checking if the gene is in rec
  rec_crit2=m2_rec_fct[,subtype2]>rec_thresh&(n_p_e2[rownames(m2_rec_fct),subtype2]>n_expressed_thresh)&Matrix::rowSums(!is.na(validated_pairs_matrix))>0
  
  ligmat1=matrix(t_u_f_lig1[,subtype1],nrow(m1_rec_fct),nrow(t_u_f_lig1),byrow = T)
  ligmat2=matrix(t_u_f_lig2[,subtype1],nrow(m2_rec_fct),nrow(t_u_f_lig2),byrow = T)
  
  recmat1=matrix(m1_rec_fct[,subtype2],nrow(m1_rec_fct),nrow(t_u_f_lig1))
  recmat2=matrix(m2_rec_fct[,subtype2],nrow(m2_rec_fct),nrow(t_u_f_lig2))
  
  intensity_score1=log10(reg+recmat1*ligmat1*validated_pairs_matrix)
  intensity_score2=log10(reg+recmat2*ligmat2*validated_pairs_matrix)
  lig_crit1=(ligmat1>lig_thresh)&matrix(n_p_e1[colnames(intensity_score1),subtype1]>n_expressed_thresh&Matrix::colSums(!is.na(validated_pairs_matrix[,colnames(intensity_score1)]))>0,nrow(intensity_score1),ncol(intensity_score1),byrow=T)
  lig_crit2=(ligmat2>lig_thresh)&
    matrix(n_p_e2[colnames(intensity_score2),subtype1]>n_expressed_thresh&
             Matrix::colSums(!is.na(validated_pairs_matrix[,colnames(intensity_score2)]))>0,
           nrow(intensity_score2),ncol(intensity_score2),byrow=T)
  
  
  pair_mask_mat=(lig_crit1|lig_crit2)&(rec_crit1|rec_crit2)&
    ifelse(is.na(validated_pairs_matrix),F,T)    
  
  z1=intensity_score1*ifelse(pair_mask_mat,pair_mask_mat,NA)    # For each TRUE, returns TRUE, and is FALSE : return NA 
  z2=intensity_score2*ifelse(pair_mask_mat,pair_mask_mat,NA)
  row_mask=pmax(Matrix::rowSums(pair_mask_mat,na.rm=T),0)>0   
  column_mask=pmax(Matrix::colSums(pair_mask_mat,na.rm=T),0)>0
  z1=z1[which(row_mask),which(column_mask),drop=F]  #only keep rows and columns with values
  z2=z2[which(row_mask),which(column_mask),drop=F]
  
  ligands=colnames(z1)
  receptors=rownames(z1)
  
  sample_to_patient<-df3CD_o1
  names(sample_to_patient)<-df3CD_o1
  
  
  
  inflamed_samples_v2_filtered<-df3CD_o1
  
  # Differential Expression fucntions:
  
  DE_ligand_receptor_interactions=function(mask_pat2_ligands,mask_pat1_ligands,clusters_ligands,mask_pat2_clusters_receptors,mask_pat1_clusters_receptors,ligands,receptors,nmin_umi_thresh=0,  
                                           nmin_cells_with_min_umi=20,reg=1e-14,nchunks=100,n_per_chunk=1000,noise_correction=F){
    
  
    mask_pat1_clusters_ligands=dfall_all_df3CD$names[dfall_all_df3CD$sample%in%samps1&dfall_all_df3CD$Group_v2%in%clusters_ligands]
    mask_pat2_clusters_ligands=dfall_all_df3CD$names[dfall_all_df3CD$sample%in%samps2&dfall_all_df3CD$Group_v2%in%clusters_ligands]
    
    ntot_umis_ligands=Matrix::colSums(umitab_all_maj_lin[,c(mask_pat2_ligands,mask_pat1_ligands)])
    u_ligands=umitab_all_maj_lin[,c(mask_pat2_ligands,mask_pat1_ligands)]
    u_ligands_clusters=umitab_all_maj_lin[,c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)]
    
    # Subset on ligand list + only ligands expressed >5 in at least 25 cells
    n_cells_with_numis_above_thresh_ligands=Matrix::rowSums(u_ligands_clusters>nmin_umi_thresh)
    names(n_cells_with_numis_above_thresh_ligands)=rownames(u_ligands_clusters)
    gene_mask_ligands=n_cells_with_numis_above_thresh_ligands>nmin_cells_with_min_umi&rownames(u_ligands_clusters)%in%c(ligands)
    if (sum(gene_mask_ligands)>1){
      u_ligands_clusters=u_ligands_clusters[gene_mask_ligands,,drop=F]
      u_ligands=u_ligands[gene_mask_ligands,,drop=F]
    
      u_receptors=umitab_all_maj_lin[,c(mask_pat2_clusters_receptors,mask_pat1_clusters_receptors)]
      u_receptors_clusters=umitab_all_maj_lin[,c(mask_pat2_clusters_receptors,mask_pat1_clusters_receptors)]
      n_cells_with_numis_above_thresh_receptors=Matrix::rowSums(u_receptors>nmin_umi_thresh)
      names(n_cells_with_numis_above_thresh_receptors)=rownames(u_receptors)
      gene_mask_receptors=n_cells_with_numis_above_thresh_receptors>nmin_cells_with_min_umi&rownames(u_receptors)%in%c(receptors)
      
      if (sum(gene_mask_receptors)>1){
        u_receptors=u_receptors[gene_mask_receptors,c(mask_pat2_clusters_receptors,mask_pat1_clusters_receptors),drop=F]
        message("Testing ",sum(gene_mask_ligands)," ligands vs. ",sum(gene_mask_receptors)," receptors")
        
        if(nrow(u_ligands_clusters)==1){
          obs_s_ligands=pmax(sum(u_ligands_clusters[,c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)]),0)
          obs_s_pat2_ligands=pmax(sum(u_ligands_clusters[,mask_pat2_clusters_ligands,drop=F]),0)
        }else{
          obs_s_ligands=pmax(Matrix::rowSums(u_ligands_clusters[,c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)]),0)
          obs_s_pat2_ligands=pmax(Matrix::rowSums(u_ligands_clusters[,mask_pat2_clusters_ligands,drop=F]),0)
        }
        obs_s_pat1_ligands=pmax(obs_s_ligands-obs_s_pat2_ligands,0)
        obs_m_pat2_ligands=obs_s_pat2_ligands/sum(ntot_umis_ligands[mask_pat2_ligands])
        obs_m_pat1_ligands=obs_s_pat1_ligands/sum(ntot_umis_ligands[mask_pat1_ligands])
        
        if(nrow(u_receptors_clusters)==1){
          obs_s_receptors=pmax(sum(u_receptors),0)
          obs_s_pat2_receptors=pmax(sum(u_receptors[,mask_pat2_clusters_receptors,drop=F]),0)
        }else{
          obs_s_receptors=pmax(Matrix::rowSums(u_receptors),0)
          obs_s_pat2_receptors=pmax(Matrix::rowSums(u_receptors[,mask_pat2_clusters_receptors,drop=F]),0)
        }
        obs_s_pat1_receptors=pmax(obs_s_receptors-obs_s_pat2_receptors,0)
        obs_m_pat2_receptors=obs_s_pat2_receptors/sum(obs_s_pat2_receptors)
        obs_m_pat1_receptors=obs_s_pat1_receptors/sum(obs_s_pat1_receptors)
        
        # Make log2FC between group 1 and 2
        obs_interaction_intensity_pat1=matrix(obs_m_pat1_ligands,length(obs_m_pat1_ligands),
                                              length(obs_m_pat1_receptors))*matrix(obs_m_pat1_receptors,length(obs_m_pat1_ligands),
                                                                                   length(obs_m_pat1_receptors),byrow = T,dimnames=list(names(obs_m_pat1_ligands),names(obs_m_pat1_receptors)))
        obs_interaction_intensity_pat2=matrix(obs_m_pat2_ligands,length(obs_m_pat2_ligands),
                                              length(obs_m_pat2_receptors))*matrix(obs_m_pat2_receptors,
                                                                                   length(obs_m_pat2_ligands),length(obs_m_pat2_receptors),byrow = T,dimnames=list(names(obs_m_pat2_ligands),names(obs_m_pat2_receptors)))
        
        obs_log2_fc=log2((reg+obs_interaction_intensity_pat1)/(reg+obs_interaction_intensity_pat2))
        obs_log2_fc_arr=array(obs_log2_fc,dim=c(dim(obs_log2_fc)[1],dim(obs_log2_fc)[2],n_per_chunk))
        
        ncounts_bigger=obs_interaction_intensity_pat1*0
        n1_ligands=length(mask_pat2_ligands)
        n1_receptors=length(mask_pat2_clusters_receptors)
        mat_ligands=matrix(c(rep(T,n1_ligands),rep(F,length(ntot_umis_ligands)-n1_ligands)),n_per_chunk,length(ntot_umis_ligands),byrow =T)
        mat_receptors=matrix(c(rep(T,n1_receptors),rep(F,ncol(u_receptors)-n1_receptors)),n_per_chunk,ncol(u_receptors),byrow =T)
        
        ntot_receptors=ncol(u_receptors)
        ntot_ligands=ncol(u_ligands)
        
    
        pb=txtProgressBar(min = 1,max=nchunks)
        # Sampling X times, this is the long step
        for (i in 1:nchunks){
          setTxtProgressBar(pb,value = i)
          
          mat_resampled_ligands=apply(mat_ligands,1,sample,ntot_ligands)
          
          s_bg_ligands=pmax(u_ligands%*%(mat_resampled_ligands*colnames(u_ligands)%in%c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)),0)
          s_bg_ligands_ntot=ntot_umis_ligands%*%mat_resampled_ligands
          s_bg_receptors=pmax(u_receptors%*%apply(mat_receptors,1,sample,ntot_receptors),0)
          
          s_fg_ligands=pmax(obs_s_ligands-s_bg_ligands,0)
          s_fg_ligands_ntot=sum(ntot_umis_ligands)-s_bg_ligands_ntot
          m_bg_ligands=t(t(s_bg_ligands)/t(s_bg_ligands_ntot))
          m_fg_ligands=t(t(s_fg_ligands)/t(s_fg_ligands_ntot))
          
          s_fg_receptors=pmax(obs_s_receptors-s_bg_receptors,0)
          m_bg_receptors=t(t(s_bg_receptors)/Matrix::colSums(s_bg_receptors))
          m_fg_receptors=t(t(s_fg_receptors)/Matrix::colSums(s_fg_receptors))
          
          arr_fg_ligands=array(m_fg_ligands,dim=c(dim(m_fg_ligands)[1],n_per_chunk,dim(m_fg_receptors)[1]))
          arr_fg_receptors=aperm(array(m_fg_receptors,dim=c(dim(m_fg_receptors)[1],n_per_chunk,dim(m_fg_ligands)[1])),c(3,2,1))
          arr_bg_ligands=array(m_bg_ligands,dim=c(dim(m_bg_ligands)[1],n_per_chunk,dim(m_bg_receptors)[1]))
          arr_bg_receptors=aperm(array(m_bg_receptors,dim=c(dim(m_bg_receptors)[1],n_per_chunk,dim(m_bg_ligands)[1])),c(3,2,1))
          
          log2_fc=aperm(log2((reg+arr_fg_ligands*arr_fg_receptors)/(reg+arr_bg_ligands*arr_bg_receptors)),c(1,3,2))
          ncounts_bigger=ncounts_bigger+apply(abs(log2_fc)>abs(obs_log2_fc_arr),1:2,sum)
        }
        
        p.value=ncounts_bigger/(i*n_per_chunk)
        adj.p.value=matrix(p.adjust(p.value,method = "BH"),nrow(p.value),ncol(p.value),dimnames = dimnames(p.value))
        
        if(nrow(u_ligands)==1){
          de_res=data.frame(ligand= rownames(u_ligands) ,receptor=rep(colnames(p.value),each=nrow(p.value)),counts_pat2_ligands=obs_s_pat2_ligands,counts_pat1_ligands=obs_s_pat1_ligands,freq_pat2_ligands=obs_m_pat2_ligands,freq_pat1_ligands=obs_m_pat1_ligands,n_cells_with_min_umis_above_thresh=n_cells_with_numis_above_thresh_ligands[gene_mask_ligands],log2_FC=as.vector(obs_log2_fc),p.value=as.vector(p.value),adj.p.value=as.vector(adj.p.value))
        }else{
          de_res=data.frame(ligand= rownames(p.value),receptor=rep(colnames(p.value),each=nrow(p.value)),counts_pat2_ligands=obs_s_pat2_ligands,counts_pat1_ligands=obs_s_pat1_ligands,freq_pat2_ligands=obs_m_pat2_ligands,freq_pat1_ligands=obs_m_pat1_ligands,n_cells_with_min_umis_above_thresh=n_cells_with_numis_above_thresh_ligands[gene_mask_ligands],log2_FC=as.vector(obs_log2_fc),p.value=as.vector(p.value),adj.p.value=as.vector(adj.p.value))
          return(de_res) 
        }
        
        
      }
      else{
        print("not enough receptors passing thrs")
        de_res=c()
      }
    }
    else{
      print("not enough ligands passing thrs")
      de_res=c()
    }
  }
  
  
  DE_ligand_receptor_interactions_patterns=function(clusters_ligands,clusters_receptors,samples,ligands,receptors,nmin_umi_thresh=0,
                                                    nmin_cells_with_min_umi=20,reg=1e-14,nchunks=100,n_per_chunk=1000,ncells_per_sample=1000,noise_correction=F,
                                                    min_n_cells=100){
    
    samp_by_sample=function(mask,cell_to_sample,ncells_per_sample){
      batches=unique(cell_to_sample[mask])
      mask2=c()
      for (b in batches){
        maskb=mask[cell_to_sample[mask]==b]
        if (length(maskb)>ncells_per_sample){
          maskb=sample(maskb,size = ncells_per_sample,replace = F)
        }
        mask2=c(mask2,maskb)
      }
      return(mask2)
    }
    
    mask_pat1_ligands=dfall_all_df3CD$names[dfall_all_df3CD$sample%in%samps1]
    mask_pat2_ligands=dfall_all_df3CD$names[dfall_all_df3CD$sample%in%samps2]
    
    ## Subsample to 2000 cells max by samples
    cell2spl<-dfall_all_df3CD$sample
    names(cell2spl)<-dfall_all_df3CD$names
    mask_pat1_ligands=samp_by_sample(mask_pat1_ligands,cell2spl,ncells_per_sample)
    mask_pat2_ligands=samp_by_sample(mask_pat2_ligands,cell2spl,ncells_per_sample)
    
    
    ## get cellnames of cells in samps1 / samps2  but only in the clusters selected
    mask_pat1_clusters_receptors=dfall_all_df3CD$names[dfall_all_df3CD$sample%in%samps1&dfall_all_df3CD$Group_v2%in%clusters_receptors]
    mask_pat2_clusters_receptors=dfall_all_df3CD$names[dfall_all_df3CD$sample%in%samps2&dfall_all_df3CD$Group_v2%in%clusters_receptors]
    mask_pat1_clusters_receptors=samp_by_sample(mask_pat1_clusters_receptors,cell2spl,ncells_per_sample)
    mask_pat2_clusters_receptors=samp_by_sample(mask_pat2_clusters_receptors,cell2spl,ncells_per_sample)
  
    if (min(c(length(mask_pat1_ligands),
              length(mask_pat2_ligands)))>=min_n_cells&min(c(length(mask_pat1_clusters_receptors),
                                                             length(mask_pat2_clusters_receptors)))>=min_n_cells){
      print("Ok for DE")
      
      return(DE_ligand_receptor_interactions(mask_pat1_ligands=mask_pat1_ligands,mask_pat2_ligands=mask_pat2_ligands,clusters_ligands=clusters_ligands,
                                             mask_pat1_clusters_receptors=mask_pat1_clusters_receptors,mask_pat2_clusters_receptors=mask_pat2_clusters_receptors,
                                             ligands = ligands,receptors = receptors,nmin_umi_thresh=nmin_umi_thresh,
                                             nmin_cells_with_min_umi=nmin_cells_with_min_umi,reg=reg,nchunks=nchunks,n_per_chunk=n_per_chunk,
                                             noise_correction=noise_correction))
    }else{print("if not passed")
      print(paste0("min_n_cells = ", min_n_cells))
      print(min(c(length(mask_pat1_ligands),length(mask_pat2_ligands))))
      print(min(c(length(mask_pat1_clusters_receptors),length(mask_pat2_clusters_receptors))))
      return(NULL)}
  }

  if(linorsub=="lin"){
    c2s<-cluster_to_broader_subtype
  }else if(linorsub=="sub"){
    c2s<-cluster_to_subtype1
  }
  
  de_res=DE_ligand_receptor_interactions_patterns(clusters_ligands=names(c2s)[c2s==subtype1],clusters_receptors=names(c2s)[c2s==subtype2],samples=inflamed_samples_v2_filtered,
                                                  ligands = ligands,receptors = receptors,nchunks = 100,nmin_cells_with_min_umi=25,ncells_per_sample=2000,nmin_umi_thresh = 0,n_per_chunk = 1e3, # changed nmin_umi_thresh to 25, nchunks =10 to be faster
                                                  noise_correction=F,min_n_cells=50)    
  
  
  
  if (!is.null(de_res)){
    rownames(de_res)=paste(de_res$ligand,de_res$receptor,sep="_")
    write.csv(de_res[order(de_res$log2_FC),],file=paste(tables_output_path,"/DE_inf_pat1_vs_pat2_",gsub("/","_",subtype1),"_",gsub("/","_",subtype2),".csv",sep=""))
  }
  
  
  mat=z1-z2
  ord1=order(Matrix::rowSums(mat,na.rm=T),decreasing = T)
  ord2=order(Matrix::colSums(mat,na.rm=T))
  
  
  mat_conv=as.data.frame(as.table(mat))
  mat_conv=mat_conv[,c(2:1,3)]
  colnames(mat_conv)=c("Ligand","Receptor","Log10_ratio_pat1_pat2")
  
  stats_mat=cbind(subtype1,subtype2,mat_conv[,1:2],log2_ratio_score1_score2=log2(10)*mat_conv[,3],
                  log10_exprs_ligand_pat1=log10(reg+as.vector(t_u_f_lig1[mat_conv[,1],subtype1])),
                  log10_exprs_ligand_pat2=log10(reg+as.vector(t_u_f_lig2[mat_conv[,1],subtype1])),
                  n_expressed_ligand_pat1=n_p_e1[as.character(mat_conv[,1]),subtype1],
                  n_expressed_ligand_pat2=n_p_e2[as.character(mat_conv[,1]),subtype1],
                  log10_exprs_receptor_pat1=log10(reg+m1_rec_fct[as.character(mat_conv[,2]),subtype2]),
                  log10_exprs_receptor_pat2=log10(reg+m2_rec_fct[as.character(mat_conv[,2]),subtype2]),
                  n_expressed_receptor_pat1=n_p_e1[as.character(mat_conv[,2]),subtype2],
                  n_expressed_receptor_pat2=n_p_e2[as.character(mat_conv[,2]),subtype2],
                  log10_score1=as.vector(z1),log10_score2=as.vector(z2))
  
  stats_mat=stats_mat[!is.na(stats_mat[,5]),]
  if (!is.null(de_res)){
    stats_mat$p.value=de_res[paste(stats_mat$Ligand,stats_mat$Receptor,sep="_"),"p.value"]
    stats_mat$adj.p.value=p.adjust(stats_mat$p.value)
  }else{
    stats_mat$p.value=NA
    stats_mat$adj.p.value=NA
  }
  rownames(stats_mat)=paste(stats_mat[,1],stats_mat[,2],stats_mat[,3],stats_mat[,4],sep="_")
  return(list(stats_mat=stats_mat,intensity_score1=intensity_score1,intensity_score2=intensity_score2,
              pair_mask=pair_mask_mat,plot_list=list(mat=mat,ord1=ord1,ord2=ord2)))
}



lineages<-c("Bcells","DCs","EndoCs","Fibros","GliPerSMs","Mono","Macro","PCs","Tcells" )
st=expand.grid(c("Mono","Macro"),lineages)
intensity_score1=c()
intensity_score2=c()
pair_mask=c()
interaction_stats=c()

# For each pair Lineage ligands vs Lineage Receptors, do the DE test 
res_l=list()
for (i in 1:nrow(st)){ 
  subtype1=as.character(st[i,1])
  subtype2=as.character(st[i,2])
  key_s=paste(subtype1,subtype2,sep="_")
  print(paste0(i,"/",nrow(st),"  ",key_s))
  
  res_l[[key_s]]=get_two_subtypes_interactions_linorsub(subtype1=subtype1,subtype2=subtype2,tables_output_path = tables_output_path,linorsub = "lin")
  interaction_stats=rbind(interaction_stats,res_l[[key_s]][["stats_mat"]])
  intensity_score1=c(intensity_score1,res_l[[key_s]][["intensity_score1"]])
  intensity_score2=c(intensity_score2,res_l[[key_s]][["intensity_score2"]])
  pair_mask=c(pair_mask,res_l[[key_s]][["pair_mask"]])
}

save(list = c("res_l","interaction_stats","intensity_score1","intensity_score2","pair_mask"),file=paste(tables_output_path,"LigRec_MoMac_lineages_v2.rd",sep=""))  # V2 has Mono vs Mono and Mac vs Mac

#?# save(list = c("res_l","interaction_stats","intensity_score1","intensity_score2","pair_mask"),file=paste(tables_output_path,"res_l_testOct.rd",sep=""))  # V2 has Mono vs Mono and Mac vs Mac


# Plot the Ligand-receptor results
bin_plot2_ggplot <- function(subtype1, subtype2, mat, ord1, ord2, stats, fdr_thresh=1e-2, figure_path) {
  mat_ord <- mat[ord1, ord2, drop = FALSE]
  # convert mat for ggplot2
  mat_df <- melt(mat_ord)
  colnames(mat_df) <- c("Receptor","Ligand","Value")
  mat_df$Ligand <- factor(mat_df$Ligand, levels = colnames(mat)[ord2])  # if mat has ligands in colnames
  mat_df$Receptor <- factor(mat_df$Receptor, levels = rownames(mat)[ord1])
  
  stats2 <- stats[stats$subtype1 == subtype1 & stats$subtype2 == subtype2, ]
  
  # If we have p values, create Asterisk dataframe + plot it, else ignore
  if(!is.na(stats2$adj.p.value )){
    pmask <- stats2$adj.p.value < fdr_thresh
    
    if(count(pmask==TRUE)>0){ # If it pass the threshold, put asterisks on the graph, if not put all as NA
      stats_asterisks <- data.frame(
        Receptor = factor(stats2[pmask, ]$Receptor, levels = rownames(mat)[ord1]),
        Ligand = factor(stats2[pmask, ]$Ligand, levels = colnames(mat)[ord2]),
        label = "*"
      )
      
      mat_df$star <- ifelse(
        with(mat_df, paste(Receptor, Ligand) %in% paste(stats_asterisks$Receptor, stats_asterisks$Ligand)),
        "*", NA
      )
    }else{
      stats_asterisks <- data.frame(
        Receptor = NA, 
        Ligand = NA, 
        label = NA
      )
      mat_df$star <- ifelse(
      with(mat_df, paste(Receptor, Ligand) %in% paste(stats_asterisks$Receptor, stats_asterisks$Ligand)),
      "*", NA
    )}
    
    
    p <- ggplot(mat_df, aes(x = Receptor, y = Ligand, fill = Value)) +
      geom_tile() +
      scale_fill_gradient2(low = "black", mid = "gray", high = "red", midpoint = 0, na.value = "white",limits=c(-2,2), oob=squish) +
      geom_text(aes(label = star), size = 5.5,fontface = 'bold', color = "black") +
      geom_text(aes(label = star), size = 5, color = "white") +
      # geom_text(data = stats_asterisks, aes(x = Receptor, y = Ligand, label = label), size = 5, color = "white") +
      labs(x = paste(subtype2, "Receptors"), y = paste(subtype1, "Ligands")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = rel(1.5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1,size = rel(1.5)),
        axis.title.x = element_text(margin = margin(t = 15),size = rel(1)),
        axis.title.y = element_text(margin = margin(r = 15),size = rel(1)),
        panel.border = element_rect(colour = "black", fill=NA, size=rel(2))
      )+scale_x_discrete( expand = c(0, 0))+scale_y_discrete( expand = c(0, 0)) +  # expand = no space btw plot and axis
      coord_fixed()  
    
  }else{
    
    p <- ggplot(mat_df, aes(x = Receptor, y = Ligand, fill = Value)) +
      geom_tile() +
      scale_fill_gradient2(low = "black", mid = "gray", high = "red", midpoint = 0, na.value = "white",limits=c(-2,2), oob=squish) +
      labs(x = paste(subtype2, "Receptors"), y = paste(subtype1, "Ligands")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = rel(1.5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1,size = rel(1.5)),
        axis.title.x = element_text(margin = margin(t = 15),size = rel(1)),
        axis.title.y = element_text(margin = margin(r = 15),size = rel(1)),
        panel.border = element_rect(colour = "black", fill=NA, size=rel(2))
      )+scale_x_discrete( expand = c(0, 0))+scale_y_discrete(expand = c(0, 0))   +   # limits=rev,
      coord_fixed() 
  }
  
  # Save plots as PDF
  ggsave(filename = file.path(figure_path, paste("ggbin_LIN_", gsub("/","_",subtype1), gsub("/","_",subtype2), ".pdf", sep = "_")),
         plot = p, width = pmax(2 + 0.25 * nrow(mat),3.5), height = pmax(1.75 + 0.25 * ncol(mat),3.25) )

  
  return(p)
}


for (i in 1:nrow(st)){
  subtype1=as.character(st[i,1])
  subtype2=as.character(st[i,2])
  key_s=paste(subtype1,subtype2,sep="_")
  print(key_s)
  if(!is.null(res_l[[key_s]])){
    bin_plot2_ggplot(subtype1,subtype2,mat=res_l[[key_s]]$plot_list[["mat"]],ord1=res_l[[key_s]]$plot_list[["ord1"]],ord2=res_l[[key_s]]$plot_list[["ord2"]],stats=interaction_stats,figure_path = main_figures_path)
  }
}


