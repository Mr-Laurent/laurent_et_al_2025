library(RcisTarget)
library(DT)
library(htmlwidgets)
library(ggplot2)
library(Rfast)


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
sigall<-list()
for(i in ls(pattern = "^grp_")){
  eval(parse(text=paste0("sigall$Sig",gsub("_ids$","",i),"<- automate_strsplit(genelist, 'Gens', ",i,")" )))
}

## RcisTarget rankings of top motifs in 500bp up - 100bp down around TSS of genes in the given program

data(motifAnnotations_hgnc)
motifRankings <- importRankings("./Grouped_objects/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")

unqmod<-c()
for(i in names(sigall) ){
  geneLists <- as.vector( unlist(sigall[i]) )
  t<-try( motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,motifAnnot=motifAnnotations_hgnc))
  if("try-error" %in% class(t)){}else{
    # Motif enrichment analysis:
    motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                             motifAnnot=motifAnnotations_hgnc)
    motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
    resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:200,]  #Save the 200 top motifs as HTML  
    y<-DT::datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
                     escape = FALSE, # To show the logo
                     filter="top", options=list(pageLength=5))
    DT::saveWidget(y, paste0("./Figure 3/Rcis_results/",i,"_TF_MacsSig_500bpu100bpd.html"))
    write.csv2(motifEnrichmentTable_wGenes[1:200,c(3,4,5,7,9)],file=paste("./Figure 3/Rcis_results/MacsSig_500bpup100bpdown_",i,".csv",sep=""), row.names = TRUE )
  }}




## RcisTarget rankings of top motifs in 10kbp up - 10kbp down
motifRankings <- importRankings("./Grouped_objects/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
unqmod<-c()
for(i in names(sigall) ){
  geneLists <- as.vector( unlist(sigall[i]) )
  t<-try( motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,motifAnnot=motifAnnotations_hgnc))
  if("try-error" %in% class(t)){}else{
    # Motif enrichment analysis:
    motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                             motifAnnot=motifAnnotations_hgnc)
    motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
    resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:200,]
    y<-DT::datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
                     escape = FALSE, 
                     filter="top", options=list(pageLength=5))
    DT::saveWidget(y, paste0("./Figure 3/Rcis_results/",i,"_TF_MacsSig_10kbu10kbd.html"))
    write.csv2(motifEnrichmentTable_wGenes[1:200,c(3,4,5,7,9)],file=paste("./Figure 3/Rcis_results/MacsSig_10kbup10kbdown_",i,".csv",sep=""), row.names = TRUE )
    
  }}

