# Post-op recurrence of Crohn (REMIND cohort)
# Ngollo, Allez paper (J Crohn Colitis 2022)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186582
# [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

setwd("~/Grouped_objects/RNAseq_REMIND/")
library(affyio)
library(affy)
library(Biobase)
library(GEOquery)
library(annotate)
library(org.Hs.eg.db)
library(AnnotationForge)
library(human.db0)
library(hthgu133pluspm.db)
library(hgu133plus2.db)
library(dplyr)

gse<- getGEO('GSE186582')

# Access the expression data:
matex<-exprs(gse[[1]])
GSE1<-gse[[1]]
# dim(matex)
# 20186   489

## Translate Affymetrix probes to gene names:
ae.annots <- AnnotationDbi::select(
  x       = hgu133plus2.db,
  keys    = rownames(GSE1),
  columns = c("PROBEID", "SYMBOL"),
  keytype = "PROBEID"
)


## If there are duplicates, use:

# dup.ids <- ae.annots$PROBEID[duplicated(ae.annots$PROBEID)] %>% 
#   unique %>%
#   sort
# An example probe that has multiple gene-level mappings
# ae.annots[ae.annots$PROBEID == dup.ids[1], ]
# 
# collapser <- function(x){
#   x %>% unique %>% sort %>% paste(collapse = "|")
# }
# 
# ae.annots <- AnnotationDbi::select(
#   x       = hgu133plus2.db ,
#   keys    = rownames(GSE1),
#   columns = c("PROBEID", "SYMBOL"),
#   keytype = "PROBEID"
# ) %>%
#   group_by(PROBEID) %>%
#   summarise_each(funs(collapser)) %>%
#   ungroup
# 
# all(ae.annots$PROBEID == rownames(GSE1))


rownames(matex)<-ae.annots$SYMBOL
# there are 25 NA, nned to be removed:
matex<-matex[!is.na(ae.annots$SYMBOL),]
# No duplicates: table(duplicated(rownames(matex)) )

colnames(pData(phenoData(gse[[1]])))
# data_processing: The Affymetrix raw data were processed to obtain a log2 expression value for each gene probe set using the robust multichip average (RMA) method implemented in R

# Make a metadata objects with only useful info:
minipheno<-pData(phenoData(gse[[1]]))[,c(1,2,8,38:43)]
colnames(matex)<-minipheno$title[match(colnames(matex),minipheno$geo_accession)]
save(matex,file="./matex_AllezNgollo_postop.rd")


##---------------------## 
# 'table(minipheno$`location:ch1`,minipheno$source_name_ch1)
# Ctrl : 25 Ileum          (IC, CTRL)
# M0I:  196 Infl Ileum     (M0)
# M0M:  147 Ileum Margin   (MI)
# M6:   121 Neoterm Ileum  (M6)'

#remove last part to get common samples:
mtp<-table(minipheno$`location:ch1`,sub("_[0-9,A-Z]*$", "", minipheno$title))
# 221 total patients for the 498 samples
tmtp<-t(mtp)
pat_pos_M0I_M6 <- sum(rowSums(tmtp[, c("M0I", "M6")] == 1) == 2)  # 116 patients with M0I and M6
pat_pos_M0I_M0M_M6 <- sum(rowSums(tmtp[, c("M0I", "M0M", "M6")] == 1) == 3)  # 81 patients with M0I, M0M and M6
# /!\ Patients with and without underscore are not the same, it's not a typo: P17 is a Male and P_17 a Female
# Don't: minipheno$patient<-sub("^([^_]*_[^_]*).*", "\\1", minipheno$title), Do:
minipheno$patient<-sub("_[0-9,A-Z]*$", "", minipheno$title)

load("./../Modgene_Macs_enr.rd")

minipheno$clin_rec<-minipheno$`rutgeertrec:ch1`
minipheno$clin_rec[minipheno$`rutgeerts:ch1`%in%c("0","1")]<-"Rem"


## Compute the enrichment score of each program in a similar fashion as in scRNAseq:
allsum<-Rfast::colsums(matex) # need same order than for score computation
for(i in 1:length(names(Modgene))){
  print(i)
  geneshow<-as.vector(unlist(Modgene[i]))
  ilist<-intersect(geneshow,rownames(matex))
  if(length(ilist)<2){currentscore<-(matex[intersect(ilist,rownames(matex)),minipheno$title]/allsum )
  eval(parse(text=paste0("minipheno$data_Sig",names(Modgene)[i],"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }else{currentscore<-(Rfast::colsums(matex[ilist,minipheno$title])/allsum )
  eval(parse(text=paste0("minipheno$data_Sig",names(Modgene)[i],"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }}

save(minipheno,file="./minipheno_AllezNgollo_postop_goodRemRec_sig.rd")

