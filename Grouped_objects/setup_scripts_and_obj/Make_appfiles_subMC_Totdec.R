args = commandArgs(trailingOnly=TRUE)
cell_type=args[1]
setwd('G:/Mon Drive/UG_metacells/Figures paper Aout/Grouped_objects/setup_scripts_and_obj/')
library(anndata)
library(Seurat)
library(Rfast)
library(reticulate)
use_python("~/miniforge3/bin/python")   

clean<-anndata::read_h5ad(paste0('all_clean_Totdec_',cell_type,'.h5ad'))
mdata<-read_h5ad(paste0('all_mc_Totdec_',cell_type,'.h5ad'))


remvlist<-which(Rfast::colsums(mdata$X)<2)   # est ce que ça m'empêche de voir les gènes en HM après ? bon en même temps si c'est pour voir 0 expression osef

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

library(parallelDist)
res2<-cor(as.matrix(t(ds_select)), method = c("pearson"))
res2_dist=parDist(res2,method = 'euclidean')
res2_clust=hclust(res2_dist,method = 'complete')
res2_cut<-cutree(res2_clust, k = nclust)

# Dataframe of the associated module for each gene
dfcut<-as.data.frame(names(res2_cut))
dfcut$related_genes_module<-res2_cut
colnames(dfcut)<-c("X","related_genes_module")
rownames(dfcut)<-rownames(dfcut)

library(RColorBrewer)
pal=rev(rainbow(12))

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
  }}   # 10 modules: 30sec with Matrix::colSums(total@..) , 10sec with Rfast::colsums(ds)

for(i in newmod){
  eval(parse(text=paste0("total@meta.data$",i,"scaled<-as.vector(scale(total@meta.data$",i,"score, center = TRUE, scale = TRUE))")))
}


dgcmtxcounts<- as(total@assays$RNA@counts, "dgCMatrix")

save(dgcmtxcounts,file=paste0(cell_type,"_dgCmc_counts_",nclust,"m_Totdec.rd")) 



metadf<-as.data.frame(total@meta.data)
for(i in newmod){
  eval(parse(text=paste0("metadf$",i,"score<-NA")))
  eval(parse(text=paste0("metadf$",i,"score[match(colnames(total),rownames(metadf))]<-total@meta.data$",i,"score")))
}

metadf_z<-metadf
save(metadf_z,file=paste0(cell_type,"_z_",nclust,"m_metadata_Totdec.rd")) 

df<-as.data.frame(cbind(names=clean$obs_names,mc=clean$obs$metacell))
df$sample<-sub("[^-].*[-]","",df$names)
df$sample<-sub("[^_].*[_]","",df$names)
df$sample<-sub("t$","",sub("^ad","",df$sample) )

df$mc2<-paste0("mc",df$mc)
save(df, file=paste0("df_",cell_type,"_Totdec.rd"))

