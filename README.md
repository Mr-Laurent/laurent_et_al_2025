Future repository for figure codes

Packages on R version 4.0.3

affy_1.68.0
affyio_1.60.0
anndata_0.7.5.6 
AnnotationForge_1.32.0
Biobase_2.50.0
ComplexHeatmap_2.6.2
dplyr_1.0.7
DT_0.19
factoextra_1.0.7 
GEOquery_2.58.0  
ggplot2_3.3.5
ggprism_1.0.3
ggpubr_0.4.0
ggrastr_0.2.3
heatmaply_1.2.1
hgu133plus2.db_3.2.3
hthgu133pluspm.db_12.1.0
htmlwidgets_1.6.4
human.db0_3.12.0
Matrix_1.6-5 
AnnotationForge_1.32.0 
parallelDist_0.2.6 
patchwork_1.3.0
plotly_4.10.4
progeny_1.12.0
purrr_0.3.4
RcisTarget_1.10.0
reshape_0.8.8
reshape2_1.4.4
reticulate_1.24 
Rfast_2.0.7
rstatix_0.7.2
scales_1.3.0 
tglkmeans ???
tidyverse_1.3.1
viridis_0.6.1






Objects

Fig 1 --------------------------

load("Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v3.rd") # metadf_z     dataframe with module scores and scaled scores by metacells
load("Grouped_objects/df_mcall_50clust_16nov22_Macs.rd")  # df contains the metadatas: cell ID, metacell ID, sample of origin, cluster and lab of origin
load("Grouped_objects/Macs_dgCmc_16nov22_counts.rd")      # dgcmtxcounts is the dgCMatrix (condensed) with counts by metacells
load("Grouped_objects/genelist_16nov22_Macs.rd")          # genelist has the gene lists associated to each of the 150 modules 
load("Grouped_objects/ht_Enrich_Macs.rd")                 # ht has the cluster + metacells order defined with ComplexHeatmap  

dfinfo<-read.csv2(file="./Grouped_objects/sample_annots.csv",sep=",",header=T)   # metadata by patient (name, tissue, status...)
load("./Grouped_objects/MacsTot_z_mymod_150_metadata_goodlog_v3.rd")             # module scores by metacells + nCount, nFeature
load("./Grouped_objects/ht_Total_Macs.rd")                                       # ht has the cluster + metacells order defined with ComplexHeatmap  (in Momacs from the total cohort)
load("./Grouped_objects/df_mcall_Macs_select2.rd")                               # metadata by cell (sample, metacell asociated, disease, tissue...)

load("Grouped_objects/Modsig_MoMacsEnr.rd")               # Contains the genelists associated to the 17 mo-macs programs (and which modules were merged to obtain them)

load("Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v5.rd")               # module scores by metacells + nCount, nFeature + program scores by regulon

Fig 2 ----------------------------

df = read.csv("./Grouped_objects/momacs_cell_metadata.csv")                      # df contains the metadata of TAURUS cohort: cell ID, metacell ID, sample of origin, disease, annotation...
annots = read.csv("./Grouped_objects/Buckley_momacs_annot_241217.csv")           # annotation associated to each metacell of the TAURUS cohort

load("Grouped_objects/minipheno_AllezNgollo_postop_goodRemRec_sig.rd")           # minipheno contains the metadata of REMIND cohort: sample ID, timepoint, sample of origin, treatment, program score...
load("Grouped_objects/matex_AllezNgollo_postop.rd")                              # Normalized count matrix of genes (rows) by sample (columns)

load("Grouped_objects/zonegradient_v4_33roi1_500cells.rd")                       # cell count distribution in each selected zone
load("Grouped_objects/zonegradient_v4_33roi2_500cells.rd")
load("Grouped_objects/zonegradient_v4_38roi1_500cells.rd")
load("Grouped_objects/zonegradient_v4_38roi2_500cells.rd")

load("Grouped_objects/freqs3_no69_noadj_nolog.rd")                               # Cell subtype frequencies for each atient of the total cohort

load(file=paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))              # Xenium slide count matrix
load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))        # Xenium slide cell annotation of subtypes

Fig 3 ----------------------------

load(file="./Grouped_objects/df_modauc3_234regfrom490.rd")                       # Comes from #SCENIC_Heatmap_regAllGenes.R, regulon program score by cell
regulons <- readRDS("./Grouped_objects/2.6_regulons_asGeneSet.Rds")              # List of genes by regulon

load("./Grouped_objects/metadata_mac20_1_2.rd")
load("./Grouped_objects/data_counts_mac20_1_2.rd")
load("./Grouped_objects/data_logcpm_mac20_1_3.rd")

load("./Grouped_objects/norm_ut_a_ml.rd")
load("./Grouped_objects/dfall_IBD_all_df3CD.rd")

motifRankings <- importRankings("./Grouped_objects/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
motifRankings <- importRankings("./Grouped_objects/hg38__refseq-r80__10kp_up_and_down_tss.mc9nr.feather")

Fig 4 ----------------------------
load("./Grouped_objects/LigRec_MoMac_lineages_cytok_fev25.rd")
load("./Grouped_objects/all_cytokines_liste.rd")
load("./Grouped_objects/umitab_all_maj_lin_fev25.rd")
load("./Grouped_objects/mo1mo3_mo_k15_fev25.rd")

load("./Grouped_objects/data_counts_mac43to46.rd") 
load("./Grouped_objects/data_logcpm_mac43to46.rd")
obj<-readRDS("./Grouped_objects/MM250_DGEseq_DDS.rds")

Fig S1
Fig S2 ---------------------------                                     ###################### A AJOUTER AU UNCLOUD
load("./Grouped_objects/MacsfromMNP_select2_dgCmc_counts.rd")

Fig S4 ---------------------------
df = read.csv("./Grouped_objects/momacs_cell_metadata.csv")
annots = read.csv("./Grouped_objects/Buckley_momacs_annot_241217.csv")

Fig S5 ---------------------------

load(file=paste0("./Grouped_objects/dgcmtx_raw_",spl,".rd"))  ########################### Not uploaded yet
load(file=paste0("./Grouped_objects/meta_annot_june1_",spl,".rd"))
load(paste0(googlepath,"./Grouped_objects/meta_annot_",spl, ".rd"))  

