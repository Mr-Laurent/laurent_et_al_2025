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

#' Compute an enrichment score of a single gene or a list of genes
#' The score is the sum of log(1+genes counts) of the signature
#' Divided by the sum of log(1+gene counts) of all the genes of a cell
#' Adapted for heavy matrices (stay in sparse matrix) 
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
enrich_score_heavy<-function(program=NULL,loop_select=1,concat=FALSE,count_mtx=NULL,center=TRUE,metadf=NULL,pseudocount=1){
  allsum <- colSums({
    log_sparse <- dgcmtxcounts       # copy sparse structure
    log_sparse@x <- log(log_sparse@x + pseudocount)
    log_sparse
  }) + log(pseudocount) * (nrow(dgcmtxcounts) - colSums(dgcmtxcounts != 0))
  
  for(prgm in loop_select){
    print(prgm)
    if(concat==FALSE){
      gnshow<-program[prgm]
    }else{
      gnshow<-unlist(strsplit(as.character(program[prgm]),split=","))
      gnshow<-gsub(" ","",gnshow)
    }
    
    if(length(gnshow)<2){
      # values<-log(pseudocount+count_mtx[gnshow,])/allsum   ADAPT FOR LARGE SPARSE:
      row_sel <- count_mtx[gnshow, , drop = FALSE]   # dgCMatrix subset
      log_sparse <- row_sel
      log_sparse@x <- log(log_sparse@x + pseudocount)
      # add zero contribution
      out <- as.numeric(log_sparse) # sparse values expanded into dense vector
      out[out == 0] <- log(pseudocount)
      values<-out/allsum
      if(center==TRUE){
        eval(parse(text=paste0("metadf$My_Mod_",prgm,"scaled<-as.vector(scale(values, center = TRUE, scale = TRUE))")))
      }else{eval(parse(text=paste0("metadf$My_Mod_",prgm,"score<-values"))) }
    }else{
      #values<-Rfast::colsums(log(pseudocount+count_mtx[gnshow,]))/allsum    ADAPT FOR LARGE SPARSE:
      row_sel <- count_mtx[gnshow[which(gnshow%in%rownames(count_mtx))], , drop = FALSE] # dgCMatrix subset, keep only genes found
      log_sparse <- row_sel
      log_sparse@x <- log(log_sparse@x + pseudocount)
      # add zero contribution
      zero_contrib <- log(pseudocount) * (nrow(row_sel) - colSums(row_sel != 0))
      out <- colSums(log_sparse) + zero_contrib
      values<-out/allsum
      if(center==TRUE){
        eval(parse(text=paste0("metadf$My_Mod_",prgm,"scaled<-as.vector(scale(values, center = TRUE, scale = TRUE))")))
      }else{eval(parse(text=paste0("metadf$My_Mod_",prgm,"score<-values"))) }
    }
    
    
  }
  return(metadf)
}

#' Quick tool to simplify the R object making following the metacell creation
#' 
#' @param object Name of the subset used
#' @return dgcmtxcounts and metadf_z, the count matrix and the metadata with program enrichment scores
#' @export
#' 
#' 
make_objs<-function(object=""){
  colnames(genelist)<-c("Mod","Gens")
  mdata<- anndata::read_h5ad(paste0('all__',object,'.h5ad'))
  clean<-read_h5ad(paste0('all_clean_',object,'.h5ad'))
  
  remvlist<-which(colSums(mdata$X)==0)
  mdatarm<-mdata$X[,-remvlist] 
  rownames(mdatarm)<-paste0("mc",rownames(mdatarm))
  
  dgcmtxcounts<-as(t(mdatarm),"dgCMatrix")
  gnlist <- as.list(setNames(genelist$Gens, genelist$Mod))
  metadf<-data.frame(nCount_RNA=as.numeric(Matrix::colSums(dgcmtxcounts)),nFeature_RNA=as.numeric(Matrix::colSums(dgcmtxcounts>0)),row.names=colnames(dgcmtxcounts))
  # Compute the scores by module:
  metadf_z<-enrich_score_heavy(program=gnlist,loop_select=c(1:length(gnlist)),concat=TRUE,  # Use enrich_score_fast if the matrices are not too big
                               count_mtx=dgcmtxcounts,metadf=metadf)
  return(list(dgcmtxcounts=dgcmtxcounts, metadf_z=metadf_z))
}


##---##---##---##---##---##---##---##---##---##---##---##---##---##---##---##



#### Step 1: Create a merged object from all total lamina propria samples (n=51) #### 

``` 01_Merge_totalLP_h5.py ```
python /SCRATCH-BIRD/users/tlaurent/metacells_cr5/anndatatot_merge.py alltot_v2_s_n_dec22
# /// Script python
# /// 01_Merge_totalLP_h5.py 
#
# Input: All the .h5ad files from the Nantes + Sinai patients with total LP
#
# Output: adalltot_v2_s_n_dec22.h5ad (merged AnnData object with n_obs x n_vars = 489621 x 21932)
#
# ///

mdata<- anndata::read_h5ad('E:/BiRD_user_tlaurent/metacells_cr5/A_mcs/adalltot_v2_s_n_dec22.h5ad')

``` 02_Metacell_totalLP.py ```
# /// Script python
# /// 02_Metacell_totalLP.py 
#
# Input: adalltot_v2_s_n_dec22.h5ad (merged AnnData object with n_obs x n_vars = 489621 x 21932)
#
# Output: all_clean_all_alltot_v2_s_n_dec22.h5ad  (AnnData object of Total LP cells)
#         all__alltot_v2_s_n_dec22.h5ad           (AnnData object of Total LP metacells)
#
# ///

#-----------------------------------------------------------#



setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')
library(anndata)
library(chameleon)
library(pheatmap)
library(pracma)
library(stats)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plyr)
library(ComplexHeatmap)
library(viridis)
library(gplots)
library(circlize)
library(pals)
library(heatmaply)
library(scales)
library(Matrix)

#If anndata doesn't work, try anndataR
mdata<- anndata::read_h5ad('./all__alltot_v2_s_n_dec22.h5ad')
clean<- anndata::read_h5ad('./all_clean_all_alltot_v2_s_n_dec22.h5ad')
mcnam<-"alltot_v2_s_n_dec22"

umis <- as.matrix(mdata$X)

fractions <- umis / rowSums(umis)
log_fractions <- log2(1e-5 + fractions)
feature_log_fractions <- log_fractions[,mdata$var$top_feature_gene]
dim(feature_log_fractions)
# 9091 788 
set.seed(42)

# sort metacells in 100 clusters
k_means <- stats::kmeans(feature_log_fractions, centers=100)
cluster_of_metacells <- as.integer(k_means$cluster)
mdata$obs$cluster <- cluster_of_metacells
# save(mdata, file="all_alltot_s_n_dec22_100clust_metacells.h5ad")



# Make a dataframe of metadata by cell (cell names, metacell, sample...)
df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$sample<-sub("[^-].*[-]","",df$names)
df$sample<-sub("[^_].*[_]","",df$names)
df$sample<-sub("t$","",sub("^ad","",df$sample) )

df2<-cbind(mc=mdata$obs_names,clust=mdata$obs$cluster)
df2<-as.data.frame(df2)
df$clust<-NA
df["clust"] <- lapply("clust", function(x) df2[[x]][match(df$mc, df2$mc)])
df$ori<-NA
df[df$sample %>% grep(pattern = "^n"),]$ori<-"nantes"
df[is.na(df$ori),]$ori<-"sinai"

df$chem<-NA
df[df$sample %>% grep(pattern = "^n"),]$chem<-"V3 chem."
df[is.na(df$chem),]$chem<-"V2 chem."
df[df$sample %>% grep(pattern = "237"),]$chem<-"V3 chem."
df[df$sample %>% grep(pattern = "238"),]$chem<-"V3 chem."
df[df$sample%in%c("67","68","69"),]$chem<-"V1 chem."

save(df, file=paste0("df_mcall_100clust_",mcnam,".rd"))

load("G:/Mon Drive/UG_metacells/McExplorer_cell_lineages/df_mcall_100clust_alltot_v2_s_n_dec22.rd")

load("./genelist_150mod_myelo_v2_12mayPBMC.rd")
# remove non-expressed genes
remvlist<-which(colSums(mdata$X)==0) 
mdatarm<-mdata$X[,-remvlist] 
# 9091 19875
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))
#gene names are in: rownames(mdata$var)
#metacell names are in: rownames(mdata$obs)

dgcmtxcounts<-as(t(mdatarm),"dgCMatrix")
gnlist <- as.list(setNames(genelist$V2, genelist$V1))
metadf<-data.frame(nCount_RNA=as.numeric(Matrix::colSums(dgcmtxcounts)),nFeature_RNA=as.numeric(Matrix::colSums(dgcmtxcounts>0)),row.names=colnames(dgcmtxcounts))
# Compute the scores by module:
metadf_z<-enrich_score_heavy(program=gnlist,loop_select=c(1:length(gnlist)),concat=TRUE,
                            count_mtx=dgcmtxcounts,metadf=metadf)


save(dgcmtxcounts,file=paste0("./total_mc_",mcnam,"_counts.rd"))
save(metadf_z,file=paste0("./total_z_",mcnam,"_150_metadata.rd"))


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
save(hm_zmod,file=paste0("./total_HMz_",mcnam,"_150_metadata.rd"))


#### Step 2: from the gene modules, we identify lineages to sub divide: ####

# Write the csv files of cells to keep by lineage
# cells_from_mc_hm(data=hm_zmod,metadata=metadf_z,HMkc=40,HMkr=0,rmclu="30,31,32,33,34",name="Stromal_inTotall_dec")  # Need more precision
cells_from_mc_hm(data=hm_zmod,metadata=metadf_z,HMkc=40,HMkr=0,rmclu="8,9,10",name="B_inTotall_dec")
# cells_from_mc_hm(data=hm_zmod,metadata=metadf_z,HMkc=40,HMkr=0,rmclu="6,7,11,12,13,14,15,16,17,18,19",name="T_ILC_inTotall_dec")  # Need more precision
# cells_from_mc_hm(data=hm_zmod,metadata=metadf_z,HMkc=40,HMkr=0,rmclu="20,21,22,23,24,25,26,35,36,37,38,39,40",name="PlasmaC_inTotall_dec") # Need more precision
# cells_from_mc_hm(data=hm_zmod,metadata=metadf_z,HMkc=40,HMkr=0,rmclu="1,2,3,4,5",name="MNP_inTotall_dec")  # Need more precision
## Update: take alspo pDC and Act DC, remove Mast cells in the MNP group:

set.seed(42)

ht = draw(Heatmap(t(hm_zmod), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod)))),show_row_dend = TRUE,
                  row_dend_reorder =F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns =T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km = 40,row_km =0, use_raster = T) )

htmxord<-as.matrix(t(hm_zmod))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  
mc2kp<-append(colnames(htmxord)[c(1:which(colnames(htmxord)=="mc1126"))],c(
  colnames(htmxord)[c(which(colnames(htmxord)=="mc8898"):which(colnames(htmxord)=="mc1135"))],
  colnames(htmxord)[c(which(colnames(htmxord)=="mc2955"):which(colnames(htmxord)=="mc3008"))], 
  colnames(htmxord)[c(which(colnames(htmxord)=="mc2876"):which(colnames(htmxord)=="mc2861"))] ) )

mc2kp<-setdiff(mc2kp,c("mc1109","mc1130","mc1131","mc1137"))

load("df_mcall_100clust_alltot_v2_s_n_dec22.rd")
df$mc2<-paste0("mc",df$mc)
c2kp<-df$names[which(df$mc2%in%mc2kp)]

write.csv(x = c2kp, file=paste0("cells_to_keep_MNP_tot_select2.csv") ,quote = F,row.names = T)



#### Step 3: Build an MNP metacell object from the total LP ####
#-----------------------------------------------------------#

``` 03_Metacell_Subset_in_totalLP.py ```
# /// Script python
# /// 03_Metacell_Subset_in_totalLP.py 
#
# Input: all_clean_all_alltot_v2_s_n_dec22.h5ad  (AnnData object of Total LP cells)
#        cells_to_keep_MNP_tot_select2.csv
#
# Output: all_clean_Totalldec_MNP_select2.h5ad  (AnnData object of MNP in Total LP cells)
#         all__Totalldec_MNP_select2.h5ad       (AnnData object of MNP in Total LP metacells)
#
# ///
python 03_Metacell_Subset_in_totalLP.py \
  --clean all_clean_all_alltot_v2_s_n_dec22.h5ad \
  --excl cells_to_keep_MNP_tot_select2.csv \
  --target_size 5000 \
  --clean_out all_clean_Totalldec_MNP_select2.h5ad \
  --metacells_out all__Totalldec_MNP_select2.h5ad

#-----------------------------------------------------------#


setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

library(anndata)
library(chameleon)
library(pheatmap)
library(pracma)
library(stats)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plyr)
library(ComplexHeatmap)
library(viridis)
library(gplots)
library(circlize)
library(pals)
library(heatmaply)
library(scales)
library(Matrix)


load("./df_mcall_100clust_alltot_s_n_dec22.rd")
load("./genelist_150mod_myelo_v2_12mayPBMC.rd")
result<-make_objs(object="Totalldec_MNP_select2")
dgcmtxcounts<-result$dgcmtxcounts
metadf_z<-result$metadf_z
save(dgcmtxcounts,file="MNP_select2_dgCmc_counts.rd")  
save(metadf_z, file="MNP_select2_z_mymod_150_metadata_goodlog.rd")

# Make df object to recall cells to subset for Macs then DC :
df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$sample<-sub("[^-].*[-]","",df$names)
df$sample<-sub("[^_].*[_]","",df$names)
df$sample<-sub("t$","",sub("^ad","",df$sample) )

df$ori<-NA
df[df$sample %>% grep(pattern = "^n"),]$ori<-"nantes"
df[is.na(df$ori),]$ori<-"sinai"
table(df$ori)
df$chem<-NA
df[df$sample %>% grep(pattern = "^n"),]$chem<-"V3 chem."
df[is.na(df$chem),]$chem<-"V2 chem."
df[df$sample %>% grep(pattern = "237"),]$chem<-"V3 chem."
df[df$sample %>% grep(pattern = "238"),]$chem<-"V3 chem."
df[df$sample%in%c("67","68","69"),]$chem<-"V1 chem."

save(df, file=paste0("df_mcall_MNP_select2.rd"))

hm_zmod<-metadf_z[,3:152]
# Select clusters to keep:
set.seed(42)
ht = draw(Heatmap(t(hm_zmod), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod)))),show_row_dend = TRUE,
                  row_dend_reorder =F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns =T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km = 40,row_km =0, use_raster = T) )

htmxord<-as.matrix(t(hm_zmod))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  

mc2kpmac<-c( 
  colnames(htmxord)[c(which(colnames(htmxord)=="mc395"):which(colnames(htmxord)=="mc282"))],
  colnames(htmxord)[c(which(colnames(htmxord)=="mc304"):which(colnames(htmxord)=="mc376"))], 
  colnames(htmxord)[c(which(colnames(htmxord)=="mc289"):which(colnames(htmxord)=="mc451"))] )
mc2kpdc<-c( colnames(htmxord)[c(which(colnames(htmxord)=="mc280"):which(colnames(htmxord)=="mc402"))] )

df$mc2<-paste0("mc",df$mc)

c2kpmac<-df$names[which(df$mc2%in%mc2kpmac)]
c2kpdc<-df$names[which(df$mc2%in%mc2kpdc)]

write.csv(x = c2kpmac, file=paste0("cells_to_keep_MacsfromMNP_tot_select2.csv") ,quote = F,row.names = T)
write.csv(x = c2kpdc, file=paste0("cells_to_keep_DCfromMNP_tot_select2.csv") ,quote = F,row.names = T)



#### Step 4: Build Macs and DC metacell objects from the total LP ####
#-----------------------------------------------------------#
# /// Script python
# /// 03_Metacell_Subset_in_totalLP.py 
# For Momacs: 
python 03_Metacell_Subset_in_totalLP.py \
  --clean all_clean_all_alltot_v2_s_n_dec22.h5ad \
  --excl cells_to_keep_MacsfromMNP_tot_select2.csv \
  --target_size 5000 \
  --clean_out all_clean_Totalldec_MacsfromMNP_select2.h5ad \
  --metacells_out all__Totalldec_MacsfromMNP_select2.h5ad


#-----------------------------------------------------------#
# /// Script python
# /// 03_Metacell_Subset_in_totalLP.py 
# For DC: 
python 03_Metacell_Subset_in_totalLP.py \
  --clean all_clean_all_alltot_v2_s_n_dec22.h5ad \
  --excl cells_to_keep_DCfromMNP_tot_select2.csv \
  --target_size 5000 \
  --clean_out all_clean_Totalldec_DCfromMNP_select2.h5ad \
  --metacells_out all__Totalldec_DCfromMNP_select2.h5ad

#-----------------------------------------------------------#

# Dimensions:   # Momacs: 246 x 20129        # DC: 146 x 20129


setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

library(anndata)
library(chameleon)
library(pheatmap)
library(pracma)
library(stats)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plyr)
library(ComplexHeatmap)
library(viridis)
library(gplots)
library(circlize)
library(pals)
library(heatmaply)
library(scales)
library(Matrix)

# Load the genelist from Enriched MNP analysis
load("./../genelist_16nov22_Macs.rd")
result<-make_objs(object="Totalldec_MacsfromMNP_select2")
dgcmtxcounts<-result$dgcmtxcounts
metadf_z<-result$metadf_z
save(dgcmtxcounts,file="MacsfromMNP_select2_dgCmc_counts.rd")  
save(metadf_z, file="MacsfromMNP_select2_z_mymod_150_metadata_goodlog.rd")
  
# Same for DC:
load("./../genelist_16nov22_DC.rd")
result<-make_objs(object="Totalldec_DCfromMNP_select2")
dgcmtxcounts<-result$dgcmtxcounts
metadf_z<-result$metadf_z
save(dgcmtxcounts,file="DCfromMNP_select2_dgCmc_counts.rd")  
save(metadf_z, file="DCfromMNP_select2_z_mymod_150_metadata_goodlog.rd")





#### Step 5: Build Fibro metacell object from the total LP ####
setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

load("total_z_alltot_v2_s_n_dec22_150_metadata.rd")
load("total_mc_alltot_v2_s_n_dec22_counts.rd") 

library(ComplexHeatmap)
library(scales)
library(ggplot2)
library(heatmaply)

hm_zmod_rm<-metadf_z[,3:152]
set.seed(42)
ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                  row_dend_reorder =F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns =T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = 40,row_km =0, use_raster = T) )
htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  
# Fibroblast selection:
mc2kp<- colnames(htmxord)[c(which(colnames(htmxord)=="mc3277"):which(colnames(htmxord)=="mc1985"))]

load("df_mcall_100clust_alltot_v2_s_n_dec22.rd")
df$mc2<-paste0("mc",df$mc)
c2kp<-df$names[which(df$mc2%in%mc2kp)]

write.csv(x = c2kp, file=paste0("cells_to_keep_StroGli_tot_select2.csv") ,quote = F,row.names = T)

#-----------------------------------------------------------#
# /// Script python
# /// 03_Metacell_Subset_in_totalLP.py 
# For Stromal/Glial: 
python 03_Metacell_Subset_in_totalLP.py \
--clean all_clean_all_alltot_v2_s_n_dec22.h5ad \
--excl cells_to_keep_StroGli_tot_select2.csv \
--target_size 5000 \
--clean_out all_clean_Totalldec_StroGli.h5ad \
--metacells_out all__Totalldec_StroGli.h5ad

# 757 floats
#-----------------------------------------------------------#



setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

nameobj="Totalldec_StroGli"
nclust=150

library(anndata)
library(Matrix)
library(parallelDist)

## Read the Cell object ad the Metacell object :
clean<-anndata::read_h5ad(paste0('all_clean_',nameobj,'.h5ad'))
mdata<-read_h5ad(paste0('all__',nameobj,'.h5ad'))
remvlist<-which(Matrix::colSums(mdata$X)<3)
mdatarm<-mdata$X[,-remvlist] #3860 19016   to reduce only to usefull genes
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))

## Compute the variable genes to keep, play on the 2 thresholds to gets around 2500 genes
ds=t(mdatarm)
ds_mean<-Matrix::rowMeans(ds)
ds_var<-apply(ds, 1, var)
varmean_df=data.frame(m=ds_mean,v=ds_var,gene=rownames(ds))
rownames(varmean_df)=rownames(ds)

geneModuleMask<-vm_genes_in_mod(varmean_df=varmean_df,
                          vm_MeanThres=-0.25,
                          vm_VarmeanThres=1,
                          x1min=-1,
                          x1max=3.5)
  
table(geneModuleMask) 
# FALSE  TRUE
# 14207  2775  with Mth=-0.25, vmTh=1
effivar<-rownames(varmean_df)[which(geneModuleMask==T)]
min(rowSums(ds[which(rownames(ds)%in%effivar),])) #
ds_select<-ds[effivar,]
## computing gene-gene correlations, and k clusters from hierarch clust to select the wanted number of clusters
res2<-cor(as.matrix(t(ds_select)), method = c("pearson"))
res2_dist=parallelDist::parDist(res2,method = 'euclidean')
res2_clust=hclust(res2_dist,method = 'complete')
res2_cut<-cutree(res2_clust, k = nclust)  #
# Dataframe of the associated module for each gene
dfcut<-as.data.frame(names(res2_cut))
dfcut$related_genes_module<-res2_cut
colnames(dfcut)<-c("X","related_genes_module")
rownames(dfcut)<-rownames(dfcut)

## Compute module scores 
newmod<-c()
for(i in 1:nclust){ newmod<-append(newmod, paste0("My_Mod_",i))}
for( i in levels(as.factor(dfcut$related_genes_module))){ assign(paste0("My_Mod_",i),dfcut$X[which(dfcut$related_genes_module==i)])}
allmod<-c()
alllist<-list()
for(i in 1:nclust){ submod<-data.frame(Mod=paste0(newmod[i]),Gen=paste0(eval(as.name(newmod[i]))))
allmod<-rbind(allmod,submod)
alllist[i]<-list(paste0( paste0(eval(as.name(newmod[i])),collapse=",") ) ) }
genelist<-data.frame(Mods=newmod,Gens=unlist(alllist))
save(genelist,file=paste0("genelist_",nclust,"m_",nameobj,".rd"))

## Save matrix and metadata
result<-make_objs(object=nameobj)

dgcmtxcounts<-result$dgcmtxcounts
metadf_z<-result$metadf_z
save(dgcmtxcounts,file=paste0(nameobj,"_dgCmc_counts.rd")) 
save(metadf_z,file=paste0(nameobj,"_z_",nclust,"m_metadata.rd"))


#-----------------------------------------------------------#
## Make a html truthplot (expression heatmap) to select Fibroblasts ##
library(heatmaply)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

nameobj="Totalldec_StroGli"
nclust=150
load(paste0(nameobj,"_dgCmc_counts.rd"))
load(paste0(nameobj,"_z_",nclust,"m_metadata.rd"))
load(paste0("genelist_",nclust,"m_",nameobj,".rd"))

hm_zmod<-metadf_z[,3:152]
findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}

thmz<-t(hm_zmod)
set.seed(42)
ht = draw(Heatmap(thmz, column_title = "Module clustering based on corrected z score", name = "mat",
                  show_row_dend = TRUE,show_row_names = F,
                  row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns = T,show_column_names = F,
                  column_km = 40,row_km =0, use_raster = T) )

htmxord<-thmz[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  
save(htmxord,file=paste0(nameobj,"_m",nclust,".rd"))

geneshow="CD7,CD2,CD3D,STMN1,KIAA0101,PCNA,TNFRSF4,TNFRSF18,CTLA4,IL2RA,FOXP3,CCR7,SELL,LEF1,CST7,GZMA,GZMK,CD8A,CD8B,CD69,ID2,ANKRD28,GPR171,TNFAIP3,ANXA1,IL7R,TRDC,TYROBP,FCER1G,PTGDR,KRT81,KRT86,KIT,PCDH9,ALDOC,LINC00299,CMC1, XCL2, IL2RB, KLRD1,KLRF1, CLIC3, MCTP2,MZB1,TNFRSF17,SEC11C,DERL3,XBP1,IGHA2,IGHM,IGHG1,IGHG3,BANK1, VPREB3, CD24, ARHGAP24, FCRLA, RALGPS2, TNFRSF13C, SPIB,MS4A1, CD79B, CD19,IGHD, FCER2,TCL1A, CD72,CD27,CLECL1,TNFRSF13B,FCGR3A,CD14,CD68, SLC40A1, STAB1, SEPP1, CSF1R, MS4A4A, SLCO2B1, MAFB,MS4A7,C1QA,C1QB,C1QC,IL1RN,S100A8,S100A9,TREM1,KYNU,C1orf54,CLEC9A,CADM1,CPNE3,XCR1,IDO1,CLEC10A,FCER1A,CD1C,CD1D,GPR157,LAMP3,DAPP1,FSCN1,CCL19,CCL22,EBI3,GZMB,TCF4,LILRA4,IRF7,CLEC4C,IL1RL1,CPA3,MS4A2,TPSB2,TPSAB1,ADCYAP1,GATA2,HPGDS,HPGD,VWF, RAMP3, NPDC1, JAM2, PLVAP, NOTCH4, HSPG2, ESAM, CYYR1, CD93, ICAM2, S1PR1, RAMP2,CD36, CA4, TMEM88, FLT1,ACKR1, SELP,DUSP23,PROX1,MMRN1,CCL21,IGF2,LYVE1,TFF3,RGS5,NDUFA4L2,C11orf96,ACTG2,MYH11,CXCL14,ADH1B,CTSK,MMP2,LUM,PTGDS,TCF21,ADAMDEC1,CCL13,CCL8,ADAM28,HAPLN1, ABCA8,CFD,THY1,CHI3L1,BGN,PDPN,TNFRSF12A,MPZ,NRXN1,SCN7A,SEMA3B,SOX2,MYOT,GFRA3,TUBB2B,GPM6B,PLP1,XKR4,ALDH1A1,CRYAB,NTM,ANK3,LGI4,S100B"
gnshow<-unlist(strsplit(geneshow,split=","))
gnshow<-gsub(" ","",gnshow)
notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
loco<-as.matrix(log(1+dgcmtxcounts[setdiff(gnshow,notfound),colnames(htmxord)]))

heatmaply(loco, fontsize_row = 7,fontsize_col = 7,dendrogram = "none",Rowv=T,Colv=F,legendgroup="2nd",coloraxis = 'coloraxis2',scale_fill_gradient_fun =scale_fill_gradientn(
  colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis")),file = paste0(nameobj,"_m",nclust,".html") ) 

pdf(file=paste0(nameobj,"_m",nclust,".pdf"))
set.seed(42)
ht = draw(Heatmap(thmz, column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod)))),show_row_dend = TRUE,
                  row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns = T,show_column_names = F,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km = 40,row_km =0, use_raster = T) )
dev.off()

# Stromal/glial markers : https://github.com/effiken/martin_et_al_cell_2019/blob/master/input/gene_lists/gene_list_figure_2g.txt
geneshow="VWF, RAMP3, NPDC1, JAM2, PLVAP, NOTCH4, HSPG2, ESAM, CYYR1, CD93, ICAM2, S1PR1, RAMP2,CD36,CA4,TMEM88,FLT1,SELP,SELE,ACKR1,DUSP23,CXorf36,KDR,PROX1,MMRN1,RELN,ANGPT2,CCL21,LYVE1,TFF3,RGS5,NDUFA4L2,ACTG2,MYH11,CXCL14,ADH1B,CTSK,MMP2,LUM,PTGDS,CCL2,TCF21,ADAMDEC1,CCL13,CCL8,ADAM28,HAPLN1,ABCA8,CFD,THY1,CHI3L1,BGN,PDPN,CXCL2,CXCL8,TNFRSF12A,SOX2,CADM1, SORCS1, LGI4, ALDH1A1, CAB39L, NTM, ANK3, PRIMA1, HAND2, GPM6B, CRYAB, SLC22A17, ERBB3, CDH19, TMEM71, NRXN1, S100B, MPZ, MYOT, SCN7A, PLP1, CNP, SPP1,AP1S2,CLU"
gnshow<-unlist(strsplit(geneshow,split=","))
gnshow<-gsub(" ","",gnshow)
notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
loco<-as.matrix(log(1+dgcmtxcounts[setdiff(gnshow,notfound),colnames(htmxord)]))

heatmaply(loco, fontsize_row = 7,fontsize_col = 7,dendrogram = "none",Rowv=T,Colv=F,legendgroup="2nd",coloraxis = 'coloraxis2',scale_fill_gradient_fun =scale_fill_gradientn(
  colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis")),file = paste0(nameobj,"_m",nclust,"_Liste2.html") ) 

# Select the Fibroblast cells in the whole Stromal/glial subset:
mc2kp<- colnames(htmxord)[append(c(which(colnames(htmxord)=="mc402"):which(colnames(htmxord)=="mc247")),
                                 c(which(colnames(htmxord)=="mc411"):which(colnames(htmxord)=="mc728")) )]

clean<-anndata::read_h5ad(paste0('all_clean_',nameobj,'.h5ad'))
df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$mc2<-paste0("mc",df$mc)
c2kp<-df$names[which(df$mc2%in%mc2kp)]  #7106 cells
write.csv(x = c2kp, file=paste0("cells_to_keep_Fibro_fromStrGli_fromTot.csv") ,quote = F,row.names = T)

# Select the Endothelial cells in the whole Stromal/glial subset:
mc2kp <- colnames(htmxord)[c(which(colnames(htmxord)=="mc729"):which(colnames(htmxord)=="mc533")) ]
c2kp<-df$names[which(df$mc2%in%mc2kp)]  #4480 cells
write.csv(x = c2kp, file=paste0("cells_to_keep_EndoC_fromStrGli_fromTot.csv") ,quote = F,row.names = T)

# Select the Glial + Pericytes + Smooth Muscle-like cells in the whole Stromal/glial subset:
mc2kp <- append(colnames(htmxord)[c(which(colnames(htmxord)=="mc133"):which(colnames(htmxord)=="mc289")) ],
                colnames(htmxord)[c(which(colnames(htmxord)=="mc307"):which(colnames(htmxord)=="mc347")) ] )
c2kp<-df$names[which(df$mc2%in%mc2kp)]  #2471 cells
write.csv(x = c2kp, file=paste0("cells_to_keep_GliPerSM_fromStrGli_fromTot.csv") ,quote = F,row.names = T)

#-----------------------------------------------------------#
# /// Script python
# /// 03_Metacell_Subset_in_totalLP.py 
# For Fibro: 
python 03_Metacell_Subset_in_totalLP.py \
--clean all_clean_Totalldec_StroGli.h5ad \
--excl cells_to_keep_Fibro_fromStrGli_fromTot.csv \
--target_size 5000 \
--clean_out all_clean_Totdec_StrGli_Fibro.h5ad \
--metacells_out all_mc_Totdec_StrGli_Fibro.h5ad

# 284 floats
#-----------------------------------------------------------#


setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')
nameobj="Totdec_StrGli_Fibro"
nclust=150

library(anndata)
library(Matrix)
library(parallelDist)

## Read the Cell object ad the Metacell object :
clean<-anndata::read_h5ad(paste0('all_clean_',nameobj,'.h5ad'))
mdata<-read_h5ad(paste0('all__',nameobj,'.h5ad'))
remvlist<-which(Matrix::colSums(mdata$X)<3)
mdatarm<-mdata$X[,-remvlist]
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))

## Compute the variable genes to keep, play on the 2 thresholds to gets around 2500 genes
ds=t(mdatarm)
ds_mean<-Matrix::rowMeans(ds)
ds_var<-apply(ds, 1, var)
varmean_df=data.frame(m=ds_mean,v=ds_var,gene=rownames(ds))
rownames(varmean_df)=rownames(ds)

geneModuleMask<-vm_genes_in_mod(varmean_df=varmean_df,
                                vm_MeanThres=-0.25,
                                vm_VarmeanThres=1,
                                x1min=-1,
                                x1max=3.5)
table(geneModuleMask) 
# FALSE  TRUE
# 13491  2130  with Mth=-0.25, vmTh=1
effivar<-rownames(varmean_df)[which(geneModuleMask==T)]
min(rowSums(ds[which(rownames(ds)%in%effivar),])) #
ds_select<-ds[effivar,]
## computing gene-gene correlations, and k clusters from hierarch clust to select the wanted number of clusters
res2<-cor(as.matrix(t(ds_select)), method = c("pearson"))
res2_dist=parallelDist::parDist(res2,method = 'euclidean')
res2_clust=hclust(res2_dist,method = 'complete')
res2_cut<-cutree(res2_clust, k = nclust)  #
# Dataframe of the associated module for each gene
dfcut<-as.data.frame(names(res2_cut))
dfcut$related_genes_module<-res2_cut
colnames(dfcut)<-c("X","related_genes_module")
rownames(dfcut)<-rownames(dfcut)

## Compute module scores 
newmod<-c()
for(i in 1:nclust){ newmod<-append(newmod, paste0("My_Mod_",i))}
for( i in levels(as.factor(dfcut$related_genes_module))){ assign(paste0("My_Mod_",i),dfcut$X[which(dfcut$related_genes_module==i)])}
allmod<-c()
alllist<-list()
for(i in 1:nclust){ submod<-data.frame(Mod=paste0(newmod[i]),Gen=paste0(eval(as.name(newmod[i]))))
allmod<-rbind(allmod,submod)
alllist[i]<-list(paste0( paste0(eval(as.name(newmod[i])),collapse=",") ) ) }
genelist<-data.frame(Mods=newmod,Gens=unlist(alllist))
save(genelist,file=paste0("genelist_",nclust,"m_",nameobj,".rd"))

## Save matrix and metadata
result<-make_objs(object=nameobj)   

dgcmtxcounts<-result$dgcmtxcounts
metadf_z<-result$metadf_z
save(dgcmtxcounts,file=paste0(nameobj,"_dgCmc_counts.rd")) 
save(metadf_z,file=paste0(nameobj,"_z_",nclust,"m_metadata.rd"))



#-----------------------------------------------------------#
## Make the final Fibro heatmap + save ht object for the app
library(heatmaply)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

nameobj="Totdec_StrGli_Fibro"
nclust=150
load(paste0(nameobj,"_dgCmc_counts.rd"))
load(paste0(nameobj,"_z_",nclust,"m_metadata.rd"))
load(paste0("genelist_",nclust,"m_",nameobj,".rd"))

hm_zmod<-metadf_z[,3:152]
findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}

thmz<-t(hm_zmod)
set.seed(42)
ht = draw(Heatmap(thmz, column_title = "Module clustering based on corrected z score", name = "mat",
                  show_row_dend = TRUE,show_row_names = F,
                  row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns = T,show_column_names = F,
                  column_km = 40,row_km =0, use_raster = T) )

save(ht, file="./ht_Total_Fibro.rd")




#### Step 6: Build Tcells (+ILC+NK) metacell object from the total LP ####

setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

load("total_z_alltot_v2_s_n_dec22_150_metadata.rd")
load("total_mc_alltot_v2_s_n_dec22_counts.rd") 

library(ComplexHeatmap)
library(scales)
library(ggplot2)
library(heatmaply)
hm_zmod<-metadf_z[,3:152]
hm_zmod_rm<-hm_zmod
set.seed(42)

ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                  row_dend_reorder =F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns =T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = 40,row_km =0, use_raster = T) )
#-# hm_env$ht_pos = ht_pos_on_device
htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  

load("df_mcall_100clust_alltot_v2_s_n_dec22.rd")
df$mc2<-paste0("mc",df$mc)

# Tcells - ILC - NK :
mc2kp<-append(colnames(htmxord)[c(which(colnames(htmxord)=="mc425"):which(colnames(htmxord)=="mc366"))],
              c(colnames(htmxord)[c(which(colnames(htmxord)=="mc764"):which(colnames(htmxord)=="mc8355"))],
                colnames(htmxord)[c(which(colnames(htmxord)=="mc966"):which(colnames(htmxord)=="mc1668"))] ) )  
mc2kp<-setdiff(mc2kp,colnames(htmxord)[c(which(colnames(htmxord)=="mc1139"):which(colnames(htmxord)=="mc4971"))])
mc2kp<-setdiff(mc2kp,colnames(htmxord)[which(colnames(htmxord)=="mc3343")])
c2kp<-df$names[which(df$mc2%in%mc2kp)]

write.csv(x = c2kp, file=paste0("cells_to_keep_Tcells_tot_select2.csv") ,quote = F,row.names = T)


#-----------------------------------------------------------#
# /// Script python
# /// 03_Metacell_Subset_in_totalLP.py 
# For Fibro: 
python 03_Metacell_Subset_in_totalLP.py \
--clean all_clean_all_alltot_v2_s_n_dec22.h5ad \
--excl cells_to_keep_Tcells_tot_select2.csv \
--target_size 80000 \   # trget size is higher due to a large number of cells, we need bigger metacells
--clean_out all_clean_mini_Totdec_Tcells.h5ad \
--metacells_out all_mc_mini_Totdec_Tcells.h5ad

#-----------------------------------------------------------#


setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')
library(anndata)
library(Matrix)
library(parallelDist)
library(Seurat)

cell_type="Tcells"
clean<-anndata::read_h5ad(paste0('all_clean_mini_Totdec_',cell_type,'.h5ad'))
mdata<-read_h5ad(paste0('all_mc_mini_Totdec_',cell_type,'.h5ad'))
remvlist<-which(Rfast::colsums(mdata$X)<2)  
mdatarm<-mdata$X[,-remvlist] 
rownames(mdatarm)<-paste0("mc",rownames(mdatarm))
total<-CreateSeuratObject(t(mdatarm), min.cells = 0, min.features = 0)

ds=t(mdatarm)
ds_mean<-Matrix::rowMeans(ds)
ds_var<-apply(ds, 1, var)

varmean_df=data.frame(m=ds_mean,v=ds_var,gene=rownames(ds))
rownames(varmean_df)=rownames(ds)

# Var genes computation : aim to 2500 - 3000 var genes for 150 modules 
# (low thresholds : inVarMean_MeanThresh=-0.5
#                   inVarMean_varmeanThresh=2 )
x=log10(varmean_df$m)
breaks=seq(min(x),max(x),.2)
lv=log2(varmean_df$v/varmean_df$m)
z=sapply(split(lv,cut(x,breaks)),min,na.rm=T)
maskinf=is.infinite(z)
z=z[!maskinf]
b=breaks[-length(breaks)]
b=b[!maskinf]
lo=loess(z~b)

# Automatization of the Threshold selection for varmean : (not lower than 0.5)
inVarMean_MeanThresh=-1
inVarMean_varmeanThresh=3
desired_genes = 2250
while (TRUE) {
  # Calculate geneModuleMask using current thresholds
  lline2 = predict(lo, newdata = log10(varmean_df$m))
  geneModuleMask = log10(varmean_df$m)>as.numeric(inVarMean_MeanThresh) &
    log2(varmean_df$v/varmean_df$m)>lline2+as.numeric(inVarMean_varmeanThresh)
  # Check number of "variable" genes
  num_genes = length(which(geneModuleMask==TRUE))
  print(num_genes)
  # Check if number of genes is within desired range
  if (num_genes >= 2000 & num_genes <= 2500) {
    print(paste0("VarMean Thrs=",inVarMean_varmeanThresh ))
    break
  }
  # Update thresholds based on number of genes, limit at 0.5
  if (num_genes < desired_genes & inVarMean_varmeanThresh> 0.4) {
    inVarMean_varmeanThresh = inVarMean_varmeanThresh - 0.05
  } else {
    
    print(paste0("VarMean Thrs=",inVarMean_varmeanThresh ))
    break
  }
}

# Selecting the number of modules :
if(num_genes<1000){nclust=50}else if(num_genes >= 1000 & num_genes <= 2000)
{nclust=100}else if(num_genes > 2000 & num_genes <= 3000)
{nclust=150}else
{nclust=200}


effivar<-rownames(varmean_df)[which(geneModuleMask==T)]
min(rowSums(ds[which(rownames(ds)%in%effivar),])) # 
ds_select<-ds[effivar,]

# computing gene-gene correlations, and k clusters from hierarch clust to select the wanted number of clusters


res2<-cor(as.matrix(t(ds_select)), method = c("pearson"))
res2_dist=parDist(res2,method = 'euclidean')
res2_clust=hclust(res2_dist,method = 'complete')
res2_cut<-cutree(res2_clust, k = nclust)

# Dataframe of the associated module for each gene
dfcut<-as.data.frame(names(res2_cut))
dfcut$related_genes_module<-res2_cut
colnames(dfcut)<-c("X","related_genes_module")
rownames(dfcut)<-rownames(dfcut)


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

write.table(unlist(alllist),file=paste0("Modlist_",cell_type,"_",nclust,"mod_",num_genes,"g_Totdec.csv"),row.names = newmod,col.names = F)


genelist<-data.frame(Mods=newmod,Gens=unlist(alllist))
save(genelist,file=paste0("genelist_",nclust,"m_",cell_type,"_Totdec.rd"))
ds<-as.matrix(total@assays$RNA@counts)


for(i in newmod){
  print(i)
  ilist<-unlist(strsplit(genelist$Gens[which(genelist$Mods==i)],split = ","))
  if(length(ilist)<2){values<-(log(1+ds[intersect(ilist,rownames(ds)),])/Rfast::colsums(log(1+ds)) )
  eval(parse(text=paste0("total@meta.data$",i,"score<-values")))
  }else{values<-(Rfast::colsums(log(1+ds[intersect(ilist,rownames(ds)),]))/Rfast::colsums(log(1+ds)) )
  eval(parse(text=paste0("total@meta.data$",i,"score<-values")))
  }}   # 10 modules : 30sec pour Matrix::colSums(total@..) , 10sec pour Rfast::colsums(ds)

for(i in newmod){
  eval(parse(text=paste0("total@meta.data$",i,"scaled<-as.vector(scale(total@meta.data$",i,"score, center = TRUE, scale = TRUE))")))
}

dgcmtxcounts<- as(total@assays$RNA@counts, "dgCMatrix")
save(dgcmtxcounts,file=paste0(cell_type,"mini_dgCmc_counts_",nclust,"m_Totdec.rd")) 
metadf_z<-as.data.frame(total@meta.data)
save(metadf_z,file=paste0(cell_type,"mini_z_",nclust,"m_metadata_Totdec.rd")) 

df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$sample<-sub("[^-].*[-]","",df$names)
df$sample<-sub("[^_].*[_]","",df$names)
df$sample<-sub("t$","",sub("^ad","",df$sample) )

df$mc2<-paste0("mc",df$mc)
save(df, file=paste0("df_",cell_type,"_Totdec.rd"))

############################################

library(heatmaply)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
setwd("G:/Mon Drive/UG_metacells/McExplorer_cell_lineages/")

cell_type="Tcells"
cell_type="Totdec_Tcellsmini"
datapath<-"G:/Mon Drive/UG_metacells/McExplorer_cell_lineages"

mtx_file=paste(datapath,"SubsetsAppObj",list.files(path=paste0(datapath,"/SubsetsAppObj/"),pattern = paste0(cell_type,"_dgCmc_counts_\\d+m.rd")),sep="/")
metadf_z_file=paste(datapath,"SubsetsAppObj",list.files(path=paste0(datapath,"/SubsetsAppObj/"),pattern = paste0(cell_type,"_z_\\d+m_metadata.rd")),sep="/")
genelist_file=paste(datapath,"SubsetsAppObj",list.files(path=paste0(datapath,"/SubsetsAppObj/"),pattern = paste0("genelist_\\d+m_",cell_type,".rd")),sep="/")
load(file=mtx_file)
load(file=metadf_z_file)
load(file=genelist_file)
nbrmod <- gsub(".*_(\\d+)m_.*", "\\1", genelist_file)

hm_zmod<-as.matrix(metadf_z[,(eval(as.numeric(nbrmod))+4):((eval(as.numeric(nbrmod))*2)+3)])

findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}

thmz<-t(hm_zmod)
set.seed(42)
ht = draw(Heatmap(thmz, column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod)))),show_row_dend = TRUE,
                  row_dend_reorder = F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns = T,show_column_names = F,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km = 40,row_km =0, use_raster = T) )
#save(ht,file="G:/Mon Drive/UG_metacells/Figures/ht_Total_Tcells.rd")
save(ht,file=paste0("./ht_",cell_type,".rd"))
#### Clean Heatmap
pdf(file = paste0("./HM_40clust_",cell_type,".pdf"),width=11,height = 9)
rownames(thmz)<-gsub("scaled","",rownames(thmz))
set.seed(42)
draw(Heatmap(thmz, column_title_gp =gpar(fontsize=rel(800/length(colnames(hm_zmod)))) , name = "mat",column_title = c(1:40),
             row_names_gp =gpar(fontsize=rel(700/length(colnames(hm_zmod)))),show_row_dend = F,show_column_dend = F,
             row_dend_reorder = F, cluster_rows = F, row_order=row_order(ht), # cluster_rows become F and we give the old HM order 
             column_dend_reorder = F, cluster_columns = T,row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
             column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod)))),column_km =40,row_km = 0, use_raster = T) )
dev.off()

htmxord<-thmz[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ]  
save(htmxord,file=paste0("SubsetsAppObj/htmxord_",cell_type,"_",nbrmod,"m.rd"))




#### Step 7: Build other metacell object from the total LP ####


setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

load("total_z_alltot_v2_s_n_dec22_150_metadata.rd")
load("total_mc_alltot_v2_s_n_dec22_counts.rd") 

library(ComplexHeatmap)
library(scales)
library(ggplot2)
library(heatmaply)
hm_zmod<-metadf_z[,3:152]
hm_zmod_rm<-hm_zmod
set.seed(42)

ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                  row_dend_reorder =F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns =T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = 40,row_km =0, use_raster = T) )

htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ] 
save(htmxord, file="htmxord_mc_alltot_v2_s_n_dec22.rd")

load("df_mcall_100clust_alltot_v2_s_n_dec22.rd")
# Plasma cells metacells:
df$mc2<-paste0("mc",df$mc)

mc2kp<-append(colnames(htmxord)[c(which(colnames(htmxord)=="mc1686"):which(colnames(htmxord)=="mc8629"))],
              c(colnames(htmxord)[c(which(colnames(htmxord)=="mc8540"):which(colnames(htmxord)=="mc1664"))],
                colnames(htmxord)[c(which(colnames(htmxord)=="mc3700"):which(colnames(htmxord)=="mc8421"))] ) )
c2kp<-df$names[which(df$mc2%in%mc2kp)]
write.csv(x = c2kp, file=paste0("cells_to_keep_PlasmaCells_tot_select2.csv") ,quote = F,row.names = T)

# Mast cells metacells:
mc2kp<-colnames(htmxord)[c(which(colnames(htmxord)=="mc1069"):which(colnames(htmxord)=="mc8809"))] 
c2kp<-df$names[which(df$mc2%in%mc2kp)]  #1377 cells
write.csv(x = c2kp, file=paste0("cells_to_keep_Mast_tot_select2.csv") ,quote = F,row.names = T)


##
###
##


setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')

load("total_z_alltot_v2_s_n_dec22_150_metadata.rd")
load("total_mc_alltot_v2_s_n_dec22_counts.rd") 
nameobj="Totalldec_StroGli"
nclust=150

load("Totalldec_StroGli_dgCmc_counts.rd")) 
load("Totalldec_StroGli_z_150m_metadata.rd"))

library(anndata)
library(Matrix)
library(parallelDist)

## Read the Cell object ad the Metacell object :
clean<-anndata::read_h5ad(paste0('all_clean_Totalldec_StroGli.h5ad'))
mdata<-read_h5ad(paste0('all__Totalldec_StroGli.h5ad'))


library(ComplexHeatmap)
library(scales)
library(ggplot2)
library(heatmaply)
hm_zmod<-metadf_z[,3:152]
hm_zmod_rm<-hm_zmod
set.seed(42)

ht = draw(Heatmap(t(hm_zmod_rm), column_title = "Module clustering based on corrected z score", name = "mat",
                  row_names_gp =gpar(fontsize=rel(400/length(colnames(hm_zmod_rm)))),show_row_dend = TRUE,
                  row_dend_reorder =F, cluster_rows = T, row_order=NULL,
                  column_dend_reorder = F, cluster_columns =T,
                  column_names_gp = gpar(fontsize=rel(400/length(rownames(hm_zmod_rm)))),column_km = 40,row_km =0, use_raster = T) )

htmxord<-as.matrix(t(hm_zmod_rm))[ as.vector(unlist(row_order(ht))) , as.vector(unlist(column_order(ht))) ] 
clean<-anndata::read_h5ad(paste0('all_clean_',nameobj,'.h5ad'))
df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$mc2<-paste0("mc",df$mc)
mc2kp <- colnames(htmxord)[c(which(colnames(htmxord)=="mc729"):which(colnames(htmxord)=="mc533")) ]
c2kp<-df$names[which(df$mc2%in%mc2kp)]  #4480 cells
write.csv(x = c2kp, file=paste0("cells_to_keep_EndoC_fromStrGli_fromTot.csv") ,quote = F,row.names = T)


##
###
##

# Once the cells from lineages are selected, use
``` 03_Metacell_Subset_in_totalLP.py ``` 
# To make the metacells and
``` Rscript Make_appfiles_subMC_Totdec.R ```
# To make the objects





