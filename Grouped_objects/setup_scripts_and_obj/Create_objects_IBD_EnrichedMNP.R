#### ~ Functions ~ ####
#' Get metacells names from the heatmap of the metacells object
#'
#' @param data Metacell z matrix
#' @param metadata metadata with cells corresponding to metacells
#' @param HMkc Number of k groups of columns
#' @param HMkr Number of k groups of rows
#' @param rmclu list of clusters to remove, ex: "1,2,4"
#' @param name Name of the file
#' @return A matrix of the infile
#' @export
#' 
#' 
cells_from_mc_hm<-function(data=NULL,metadata=NULL,HMkc=40,HMkr=0,rmclu=NULL,name=NULL){
  hm_zmod_rm<-data
  print(paste0("computing heatmap of ",HMkc,"cols and ",HMkr,"rows"))
  set.seed(42)
  ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                    row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                    row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                    column_dend_reorder = F, cluster_columns = T,
                    column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = HMkc,row_km = HMkr, use_raster = T) )
  
  metacells_kclust<- column_order(ht) 
  metacells_kclustrow<- row_order(ht) 
  for (i in 1:length(metacells_kclust)){
    if (i == 1) {
      clu <- t(t(rownames(hm_zmod_rm)[metacells_kclust[[i]]]))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("Metacell", "Cluster")
    } else {
      clu <- t(t(rownames(hm_zmod_rm)[metacells_kclust[[i]]]))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }
  rm(clu)
  out<-as.data.frame(out)
  
  rm_clu<-paste0("cluster",unlist(strsplit(rmclu,split=",")))
  Mc_to_subset<- out$Metacell[which(out$Cluster%in%rm_clu)]
  print("metacells identified")
  metadata$mc2<-paste0("mc",metadata$mc)
  write.csv(x = metadata$names[which(metadata$mc2%in%Mc_to_subset)], file=paste0("cells_to_keep_",name,".csv") ,quote = F,row.names = T)
  
  print("csv written")
}

#' Get the gene to use in module analysis and plot them
#' They are the most variable+expressed genes, following loess curve
#' Adapted from get_highly_variable_genes in https://github.com/effiken/martin_et_al_cell_2019
#' 
#' @param varmean_df dataframe of mean and variance of the genes
#' @param vm_MeanThres threshold for the log(mean) of expression
#' @param vm_VarmeanThres threshold for the curve start
#' @param x1min min value for the plot
#' @param x1max max value for the plot
#' @return TRUE/FALSE genes kept
#' @export
#' 
#' 
vm_genes_in_mod<-function(varmean_df=NULL,
                          vm_MeanThres=-1,
                          vm_VarmeanThres=1,
                          x1min=-1,
                          x1max=3.5){
  #Compute loess curve to get the variable genes
  x=log10(varmean_df$m)
  breaks=seq(min(x),max(x),.2)
  lv=log2(varmean_df$v/varmean_df$m)
  z=sapply(split(lv,cut(x,breaks)),min,na.rm=T)
  maskinf=is.infinite(z)
  z=z[!maskinf]
  b=breaks[-length(breaks)]
  b=b[!maskinf]
  lo=loess(z~b)
  # Plot of the loess curve for variable gene selection 
  plot(log10(varmean_df$m),log2(varmean_df$v/varmean_df$m),xlab="Log10(mean)",ylab="log2(var/mean)",panel.first=grid())
  x1=seq(x1min,x1max,l=100)
  lline=predict(lo,newdata =x1)
  lines(x1,lline+as.numeric(vm_VarmeanThres),col=2)
  abline(v=vm_MeanThres,col=2)
  lline2=predict(lo,newdata =log10(varmean_df$m))
  
  ## Selection of the genes 
  # Use the same threshold but this time to subset the genes to keep from the matrix
  geneModuleMask<-log10(varmean_df$m)>as.numeric(vm_MeanThres)&log2(varmean_df$v/varmean_df$m)>lline2+as.numeric(vm_VarmeanThres)
  return(geneModuleMask)
  
}



#' Compute an enrichment score of a single gene or a list of genes
#' The score is the sum of log(1+genes counts) of the signature
#' Divided by the sum of log(1+gene counts) of all the genes of a cell 
#' 
#' @param program A list of programs/signatures containing the genes that make them
#' @param loop_select number or name selected to loop in the program/signature list names
#' @param concat TRUE/FALSE for when the gene list is a concatenated sentence: "XXX,XXX,XXX"
#' @param count_mtx count matrix with gene expression by cell/metacell (columns are cell/metacell names)
#' @param center TRUE/FALSE variable to scale and center the result or not
#' @param metadf A dataframe for metadata/scores where rows have the cell/metacell names
#' @param pseudocount 1 by default, but can be changed
#' @return metadf, modified with the program scores in additional columns
#' @export
#' 
#' 
enrich_score_fast<-function(program=NULL,loop_select=1,concat=FALSE,count_mtx=NULL,center=TRUE,metadf=NULL,pseudocount=1){
  allsum<-Rfast::colsums(log(pseudocount+count_mtx))
  for(prgm in loop_select){
    print(prgm)
    if(concat==FALSE){
      gnshow<-program[prgm]
    }else{
      gnshow<-unlist(strsplit(as.character(program[prgm]),split=","))
      gnshow<-gsub(" ","",gnshow)
    }
    
    if(length(gnshow)<2){
      values<-log(pseudocount+count_mtx[gnshow,])/allsum
      if(center==TRUE){
        eval(parse(text=paste0("metadf$My_Mod_",prgm,"scaled<-as.vector(scale(values, center = TRUE, scale = TRUE))")))
      }else{eval(parse(text=paste0("metadf$My_Mod_",prgm,"score<-values"))) }
    }else{
      values<-Rfast::colsums(log(pseudocount+count_mtx[gnshow,]))/allsum
      if(center==TRUE){
        eval(parse(text=paste0("metadf$My_Mod_",prgm,"scaled<-as.vector(scale(values, center = TRUE, scale = TRUE))")))
      }else{eval(parse(text=paste0("metadf$My_Mod_",prgm,"score<-values"))) }
    }
    
    
  }
  return(metadf)
}



##---##---##---##---##---##---##---##---##---##---##---##---##---##---##---##


#### Step 1: Create MNP enriched metacell object ####
## Load the matrix (.rd or .h5ad) to merge the 7 samples enriched in MNP together before creating a metacell object 


# My python environment: (anndata 0.8.0, matplotlib 3.7.3, metacells 0.5, numpy 1.20.3, pandas 1.3.1, scipy 1.7.0, seaborn 0.11.2)
``` 01_Merge_enrichedMNP_h5.py ```
# /// Script python 
# /// 01_Merge_enrichedMNP_h5.py
#
# Input:  umitabcr5_n2.h5ad, umitabcr5_n4.h5ad, umitabcr5_n8.h5ad, umitabcr5_n18.h5ad, umitabcr5_n34.h5ad, umitabcr5_n48.h5ad, umitabcr5_n50.h5ad
#
# Output: ad16nov22_reAug24.h5ad 
#         (merged AnnData object with n_obs × n_vars = 119472 × 36591)
# ///


``` 02_Metacell_enrichedMNP.py ```
# /// Script python 
# /// 02_Metacell_enrichedMNP.py
#
# Input:  ad16nov22_reAug24.h5ad   (merged AnnData object with n_obs × n_vars = 119472 × 36591)
#
# Output: all_clean_all_16nov22.h5ad  (AnnData object of Enriched MNP cells) 
#         all__16nov22.h5ad           (AnnData object of Enriched MNP metacells)
#         
# ///

##----------------------------------------------------------------## 


######################### At some point, I made a lot of analysis of a "supergut" model to discover general signatures in MNP: Mc12mayPBMC with PBMC+Tissue cells from Nantes, Sinai, Powrie 

load("./genelist_150mod_myelo_v2_12mayPBMC.rd")

#### Create df metadata for Enriched MNP ####
# It assigns cell name to their metacell and sample

library(anndata)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj")

mdata<- anndata::read_h5ad('all__16nov22.h5ad')
clean<- anndata::read_h5ad('all_clean_all_16nov22.h5ad')

df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$sample<-sub("[^-].*[-]","",df$names)
df$sample<-sub("t$","",sub("^ad","",df$sample) )
df$sample<-sub("^cr5_s","",df$sample)
# remove non-expressed genes
remvlist<-which(colSums(mdata$X)==0)
mdatarm<-mdata$X[,-remvlist] 
# 1545 29198
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))

dgcmtxcounts<-as(t(mdatarm),"dgCMatrix")
gnlist <- as.list(setNames(genelist$V2, genelist$V1))
metadf<-data.frame(nCount_RNA=as.numeric(colSums(dgcmtxcounts)),nFeature_RNA=as.numeric(colSums(dgcmtxcounts>0)),row.names=colnames(dgcmtxcounts))
# Compute the scores by module:
metadf_z<-enrich_score_fast(program=gnlist,loop_select=c(1:length(gnlist)),concat=TRUE,
                            count_mtx=as.matrix(dgcmtxcounts),metadf=metadf)

hm_zmod_rm<-metadf_z[,-c(1:2)]

set.seed(42)
ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                  row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns = T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = 40,row_km = 0, use_raster = T) )

htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  

metacells_kclust<- column_order(ht) 
# Problem if only one metacell is in the cluster : it doesn't read as a matrix  
# Use rownames(hm_zmod_rm)[metacells_kclust[[i]]]
# instead of 
#     rownames(hm_zmod_rm[metacells_kclust[[i]],])
for (i in 1:length(metacells_kclust)){
  if (i == 1) {
    clu <- t(t(rownames(hm_zmod_rm)[metacells_kclust[[i]]]))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("Metacell", "Cluster")
  } else {
    clu <- t(t(rownames(hm_zmod_rm)[metacells_kclust[[i]]]))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
out<-as.data.frame(out)

## From the expressed genes, we identified clusters 2,3,4,5,6 and 37 as stromal cells, that we need to remove
Mc_to_exclude<- out$Metacell[which(out$Cluster%in%c("cluster2","cluster3","cluster4","cluster5","cluster6","cluster37"))]
df$mc<-paste0("mc",df$mc)
write.csv(x = df$names[which(df$mc%in%Mc_to_exclude)], file="strom_to_excl.csv",quote = F,row.names = T)


##----------------------------------------------------------------## 

``` 03_enrichedMNP_excl_stromal.py ```
# /// Script python 
# /// 03_enrichedMNP_excl_stromal.py
#
# Input:  strom_to_excl.csv           (cell names of identified contaminants to exclude)
#         all_clean_all_16nov22.h5ad  (AnnData object of Enriched MNP cells)
#
# Output: all_clean_all_16nov22_delStro.h5ad  (AnnData object of Enriched MNP cells) 
#         all__16nov22_delStro.h5ad           (AnnData object of Enriched MNP metacells)
#         
# ///

##----------------------------------------------------------------## 



setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj')
library(anndata)
library(chameleon)
library(pheatmap)
library(pracma)
library(stats)
library(dplyr)

mdata<- anndata::read_h5ad('./all__16nov22_delStro.h5ad')  # 1328 metacells × 30043 g
clean<-read_h5ad('./all_clean_all_16nov22_delStro.h5ad')   # 26922 cells × 30043 g

umis <- as.matrix(mdata$X)

fractions <- umis / rowSums(umis)
log_fractions <- log2(1e-5 + fractions)
feature_log_fractions <- log_fractions[,mdata$var$top_feature_gene]
dim(feature_log_fractions)
# 1328  762  gene expression level each   (V4)
set.seed(42)

# sort metacells in 50 clusters
k_means <- stats::kmeans(feature_log_fractions, centers=50)
cluster_of_metacells <- as.integer(k_means$cluster)
mdata$obs$cluster <- cluster_of_metacells

df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$sample<-sub("[^-].*[-]","",df$names)
df$sample<-sub("t$","",sub("^ad","",df$sample) )
df$sample<-sub("^cr5_s","",df$sample)

df2<-cbind(mc=mdata$obs_names,clust=mdata$obs$cluster)
df2<-as.data.frame(df2)
df$clust<-NA
df["clust"] <- lapply("clust", function(x) df2[[x]][match(df$mc, df2$mc)])
df$ori<-NA
df[df$sample %>% grep(pattern = "^n"),]$ori<-"nantes"

save(df, file="df_mcall_50clust_16nov22_delStro.rd")
# 26922 cells with associated metacell, sample and cluster

remvlist<-which(colSums(mdata$X)==0)
mdatarm<-mdata$X[,-remvlist] 
# 1328 27733
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))

load("./genelist_150mod_myelo_v2_12mayPBMC.rd")
dgcmtxcounts<-as(t(mdatarm),"dgCMatrix")

gnlist <- as.list(setNames(genelist$V2, genelist$V1))
metadf<-data.frame(nCount_RNA=as.numeric(colSums(dgcmtxcounts)),nFeature_RNA=as.numeric(colSums(dgcmtxcounts>0)),row.names=colnames(dgcmtxcounts))
# Compute the scores by module:
metadf_z<-enrich_score_fast(program=gnlist,loop_select=c(1:length(gnlist)),concat=TRUE,
                            count_mtx=as.matrix(dgcmtxcounts),metadf=metadf)

save(dgcmtxcounts,file="./mnp_mc_16nov22_counts.rd")
save(metadf_z,file="./mnp_z_mymod_150_metadata.rd")

nmod<-seq(1:150)
hmtxz<-matrix()
for(i in nmod){
  hmtxz<-( (metadf_z[,paste0("My_Mod_",nmod,"scaled")]%>%(function(x) { x[is.na(x)] <- 0; return(x) })) )
  colnames(hmtxz)[nmod]<-paste0("corr_z_Mod_",nmod)
}

matcorm<-as.matrix(hmtxz)
matcorm <- apply(matcorm, 2 ,as.numeric) # Need to be sure that it's in numeric or else I'll see the   " Error: `col` should have names to map to values in `mat` " message 
rownames(matcorm)<-rownames(hmtxz)
hm_zmod<-matcorm
save(hm_zmod,file="./mnp_HM_z_mymod_150.rd")

## From the object MNP enriched, select clusters of interest to write a CSV of cells to keep
cells_from_mc_hm(data=hm_zmod,metadata=df,HMkc=40,HMkr=0,rmclu="33,34,35,36,37,38,39",name="DC_inMNPenrich")
# It writes cells_to_keep_DC_inMNPenrich.csv
cells_from_mc_hm(data=hm_zmod,metadata=df,HMkc=40,HMkr=0,rmclu="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,32",name="Macs_inMNPenrich")
# It writes cells_to_keep_Macs_inMNPenrich.csv




##----------------------------------------------------------------## 

``` 04_Momacs_mc_enrichedMNP.py ```
# /// Script python 
# /// 04_Momacs_mc_enrichedMNP.py
#
# Input:  cells_to_keep_Macs_inMNPenrich.csv  (cell names of identified momacs to include)
#         all_clean_all_16nov22_delStro.h5ad  (AnnData object of Enriched MNP cells)
#
# Output: all_clean_all_16nov22_Macs.h5ad     (AnnData object of Enriched Momacs cells) 
#         all__16nov22_Macs.h5ad              (AnnData object of Enriched Momacs metacells: 802mc)
#         
# ///


``` 05_DCs_mc_enrichedMNP.py ```
# /// Script python 
# /// 05_DCs_mc_enrichedMNP.py
#
# Input:  cells_to_keep_DC_inMNPenrich.csv  (cell names of identified DCs to include)
#         all_clean_all_16nov22_delStro.h5ad  (AnnData object of Enriched MNP cells)
#
# Output: all_clean_all_16nov22_DC.h5ad     (AnnData object of Enriched Momacs cells) 
#         all__16nov22_DC.h5ad              (AnnData object of Enriched Momacs metacells: 802mc)
#         
# ///

##----------------------------------------------------------------## 



setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj')
library(anndata)
library(chameleon)
library(pheatmap)
library(pracma)
library(stats)
library(dplyr)
library(parallelDist)


Createobj_mnpsub<-function(object,vm_VarmeanThres){
  mdata<-anndata::read_h5ad(paste0('all__16nov22_',object,'.h5ad'))
  clean<-anndata::read_h5ad(paste0('all_clean_all_16nov22_',object,'.h5ad'))
  
  # count quick normalization
  umis <- as.matrix(mdata$X)
  fractions <- umis / rowSums(umis)
  log_fractions <- log2(1e-5 + fractions)
  feature_log_fractions <- log_fractions[,mdata$var$top_feature_gene]
  
  
  set.seed(42)
  # Cluster metacells based on normalized counts
  k_means <- stats::kmeans(feature_log_fractions, centers=50)
  cluster_of_metacells <- as.integer(k_means$cluster)
  mdata$obs$cluster <- cluster_of_metacells
  
  mdata$write_h5ad(paste0('all_16nov22_',object,'_50clust_metacells.h5ad'))
  
  # save metadata 
  df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
  df$sample<-sub("[^-].*[-]","",df$names)
  df$sample<-sub("t$","",sub("^ad","",df$sample) )
  df$sample<-sub("^cr5_s","",df$sample)
  
  df2<-cbind(mc=mdata$obs_names,clust=mdata$obs$cluster)
  df2<-as.data.frame(df2)
  df$clust<-NA
  df["clust"] <- lapply("clust", function(x) df2[[x]][match(df$mc, df2$mc)])
  
  df$ori<-NA
  df[df$sample %>% grep(pattern = "^n"),]$ori<-"nantes"
  
  remvlist<-which(colSums(mdata$X)==0)
  mdatarm<-mdata$X[,-remvlist] 
  rownames(mdatarm)<-paste0("mc",rownames(mdatarm))
  dgcmtxcounts<-as(t(mdatarm),"dgCMatrix")
  
  
  # Subset on metacells with > 20 counts and genes with > 1 expression
  ds<-t(mdata$X[,which(colSums(mdata$X)>20)])
  ds<-ds[which(Matrix::rowSums(ds)>1),]
  ds<-ds[-which(rownames(ds)=="MALAT1"),]  # MALAT1 is too much of an outlier and sometimes break the gene selection, so we remove it
  
  ds_mean<-Matrix::rowMeans(ds)
  ds_var<-apply(ds, 1, var)
  
  varmean_df=data.frame(m=ds_mean,v=ds_var,gene=rownames(ds))
  rownames(varmean_df)=rownames(ds)
  # Compute loess curve to get the variable genes
  # Selection of the threshold (Adapt with the plot, or the number of genes you want)
  # Plot of the loess curve for variable gene selection 
  # Selection of the genes 
  geneModuleMask<-vm_genes_in_mod(varmean_df=varmean_df,
                                  vm_MeanThres=-1,
                                  vm_VarmeanThres=vm_VarmeanThres,
                                  x1min=-1,x1max=3.5)
 
  effivar<-rownames(ds)[which(geneModuleMask==T)]
  ## Computing gene-gene correlations on a reduced matrix with only the selected genes,
  ## and k-mean clustering to select the wanted number of clusters
  nclust<-"150"
  ds_select<-ds[effivar,]
  res2<-cor(as.matrix(t(ds_select)), method = c("pearson"))
  res2_dist=parDist(res2,method = 'euclidean')
  res2_clust=hclust(res2_dist,method = 'complete')
  res2_cut<-cutree(res2_clust, k = nclust)
  
  # Dataframe of the associated module for each gene
  dfcut<-as.data.frame(names(res2_cut))
  dfcut$related_genes_module<-res2_cut
  colnames(dfcut)<-c("X","related_genes_module")
  rownames(dfcut)<-rownames(dfcut)
  
  ## Writting a csv file with genelist by module:
  newmod<-c()
  for(i in 1:nclust){ newmod<-append(newmod, paste0("My_Mod_",i))}
  for( i in levels(as.factor(dfcut$related_genes_module))){
    assign(paste0("My_Mod_",i),dfcut$X[which(dfcut$related_genes_module==i)])
  }
  allmod<-c()
  alllist<-list()
  for(i in 1:nclust){
    submod<-data.frame(Mod=paste0(newmod[i]),Gen=paste0(eval(as.name(newmod[i]))))
    allmod<-rbind(allmod,submod)
    alllist[i]<-list(paste0( paste0(eval(as.name(newmod[i])),collapse=",") ) )
  }
  write.table(unlist(alllist),file=paste0("Modulelist_",object,"_",nclust,"mod_mTh",
                                          gsub("-","m",as.character(inVarMean_MeanThresh)),
                                          "_vmTh",gsub("-","m",as.character(vm_VarmeanThres)),
                                          ".csv"),row.names = newmod,col.names = F)
  
  
  genelist<-c()
  for(i in 1:nclust){
    genelist<-rbind(genelist,data.frame(Mod=paste0(newmod[i]),Gens=alllist[[i]]  )  )
  }
  save(genelist,file=paste0("./genelist_16nov22_",object,".rd"))
  save(df, file=paste0("./df_mcall_50clust_16nov22_",object,".rd"))
  
  
  gnlist <- as.list(setNames(genelist$Gens, genelist$Mod))
  metadf<-data.frame(nCount_RNA=as.numeric(colSums(dgcmtxcounts)),nFeature_RNA=as.numeric(colSums(dgcmtxcounts>0)),row.names=colnames(dgcmtxcounts))
  # Compute the scores by module:
  metadf_z<-enrich_score_fast(program=gnlist,loop_select=c(1:length(gnlist)),concat=TRUE,
                              count_mtx=as.matrix(dgcmtxcounts),metadf=metadf)
  save(dgcmtxcounts,file=paste0("./../",object,"_dgCmc_16nov22_counts.rd"))
  save(metadf_z,file=paste0("./../",object,"_z_mymod_150_metadata_goodlog.rd"))
  
  #### Create Macs ht objects (Enrich) ####
  library(ggplot2)
  library(ComplexHeatmap)
  # Make the "ht" object in Enriched MNP cohort:
  hm_zmod<-metadf_z[,3:152]
  # Create the ComplexHeatmap object to order clusters
  set.seed(42)
  ht = draw(Heatmap(t(hm_zmod), column_title = "Module clustering based on corrected z score", name = "mat",
                    row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod)))),show_row_dend = TRUE,
                    row_dend_reorder = F, cluster_rows = T, row_order=NULL,show_column_names = FALSE,
                    column_dend_reorder = F, cluster_columns = T,
                    column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km = 40,row_km =0, use_raster = T) )
  save(ht,file=paste0("./../ht_Enrich_",object,".rd"))
  
}


Createobj_mnpsub(object = "Macs",vm_VarmeanThres = 1.25)
# dim(feature_log_fractions): 802x731 
# 15885 Macs cells in df
# geneModuleMask is 1899 genes

Createobj_mnpsub(object = "DC",vm_VarmeanThres = 1.5)
# dim(feature_log_fractions): 202x458 
# 5190 DC cells in df
# geneModuleMask is 3213 genes

##----------------------------------------------------------------## 

