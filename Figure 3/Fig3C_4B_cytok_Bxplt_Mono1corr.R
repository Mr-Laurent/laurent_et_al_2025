library(Matrix)
library(dplyr)
library(ggplot2)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")

load("./Grouped_objects/norm_ut_a_ml.rd")
load("./Grouped_objects/dfall_IBD_all_df3CD.rd")
load("./Grouped_objects/ord_pval_pos_mono1.rd") # subtypes sorted by higher to lower correlation to Mono1 (top 29) created in Fig2H_corr.R
load("./Grouped_objects/ord_pval_pos_mono3.rd") # created in Fig2H_corr.R

#Adapt for Mono / Macs :
dfall_all_df3CD$Group_v2 <- gsub("\\[MoMac\\] (Mono\\d+)", "[Mono] \\1", dfall_all_df3CD$Group_v2)
dfall_all_df3CD$Group_v2 <- gsub("\\[MoMac\\] (Macro\\d+)", "[Macro] \\1", dfall_all_df3CD$Group_v2)  

genes_2see<-c("IFNG","CSF2","TNF","CXCL10","CXCR3","IL23A","CCL2","CXCL1")
norm_ut_df <- as.data.frame(as.matrix(norm_ut_a_ml[,genes_2see]) )
norm_ut_df$names <- rownames(norm_ut_a_ml[,genes_2see])
merged_data <- inner_join(norm_ut_df, dfall_all_df3CD, by = "names")

mean_expression <- merged_data %>%
  group_by(sample, Group_v2) %>%
  dplyr::summarise(across(genes_2see, mean, na.rm = TRUE))

ctyp_sample_norm_cytok<-as.data.frame(mean_expression)


df3CD_mono1hi<-c("n17","n47","n33","n37","128","138","187","181","n13")
ctyp_sample_norm_cytok <- ctyp_sample_norm_cytok %>%
  mutate(Sample_Type = ifelse(sample %in% df3CD_mono1hi, "Mono1hi", "Mono1lo"))

New_names_m1 <- setNames(ord_pval_pos_mono1, ord_pval_pos_mono1)
New_names_m1["[DC] Act. DC"] <- "[DC] aDC"
New_names_m1["[DC] Act. DC Infl."] <- "[DC] iaDCb"
New_names_m1["[DC] Infl. DC3"] <- "[DC] iaDCa"
New_names_m1["[Bcells] Atypical_Mem"] <- "[Bcells] Atypical_Mem_Bcells"
New_names_m1["[Tcells] CM-like"] <- "[Tcells] CM-like_Tcells"
New_names_m1["[Tcells] Activated"] <- "[Tcells] Activated_Tcells"
New_names_m1["[Tcells] Naive/CM"] <- "[Tcells] Naive/CM_Tcells"  
New_names_m1["[Bcells] GC_like" ] <- "[Bcells] GC_like_Bcells" 
New_names_m1["[Tcells] Effector_CD8s"] <- "[Tcells] Effector_CD8s_Tcells" 
New_names_m1["[PC] IgG" ] <- "[PC] IgG_Plasma_cells" 
New_names_m1["[Bcells] Naive"] <- "[Bcells] Naive_Bcells"
New_names_m1["[Bcells] GC_like" ] <- "[Bcells] GC_like_Bcells" 
New_names_m1_nl<-gsub("^\\[.*\\] ","",New_names_m1) # names without lineage category
New_names_m1_nl[which(New_names_m1_nl=="NK-like"):which(New_names_m1_nl=="iaDCa")]<-paste0("**",New_names_m1_nl[which(New_names_m1_nl=="NK-like"):which(New_names_m1_nl=="iaDCa")],"**") # Make the enriched one as bold text


pdf("./Figure 3/Fig3C_4B_cytok_Bxplt_Mono1corr.pdf",width = 4, height = 7)
the_gene=genes_2see[1]
# First figure to get the names only, next figures are the boxplots only by gene
eval(parse(text=paste0("print(ggplot(ctyp_sample_norm_cytok[which(ctyp_sample_norm_cytok$Group_v2%in%c(ord_pval_pos_mono1,'[Mono] Mono1') ),], aes(x =Group_v2, y =",the_gene,"))+
  geom_boxplot(lwd=1, color='black', fill='grey', alpha=0.5)+scale_x_discrete(limits = c(ord_pval_pos_mono1,'[Mono] Mono1'), labels=c(as.vector(New_names_m1_nl),'Mono1') ) +xlab('')+
  geom_jitter(stat='identity',inherit.aes = TRUE,shape=16, position=position_jitter(0.2),size=1)+theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 01),axis.text.y = ggtext::element_markdown())+ ylab('Mean ",the_gene," normalized\nexpression by sample')+
  coord_flip()+geom_vline(xintercept = 19.5,linetype='dashed')+ggtitle('",the_gene," in Subtypes positively correlated to Mono1') )"
)))
for(the_gene in genes_2see){
  eval(parse(text=paste0("print(ggplot(ctyp_sample_norm_cytok[which(ctyp_sample_norm_cytok$Group_v2%in%c(ord_pval_pos_mono1,'[Mono] Mono1') ),], aes(x =Group_v2, y =",the_gene,"))+
  geom_boxplot(lwd=1, color='black', fill='grey', alpha=0.5)+scale_x_discrete(limits = c(ord_pval_pos_mono1,'[Mono] Mono1'), labels=c(as.vector(New_names_m1_nl),'Mono1') ) +xlab('')+
  geom_jitter(stat='identity',inherit.aes = TRUE,shape=16, position=position_jitter(0.2),size=1)+theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 01))+ ylab('Mean ",the_gene," normalized\nexpression by sample')+
  coord_flip()+guides(y = 'none')+geom_vline(xintercept = 19.5,linetype='dashed')+ggtitle('",the_gene," in Subtypes positively correlated to Mono1')+
                         ylim(0,max(ctyp_sample_norm_cytok$",the_gene,"[which(ctyp_sample_norm_cytok$Group_v2%in%c(ord_pval_pos_mono1,'[Mono] Mono1',ord_pval_pos_mono3) ) ] ))+scale_y_continuous(trans='log1p')  )"  # Used if we plot in Mono1 correlated subtypes or Mono3 
  )))
  # eval(parse(text=paste0("print(ggplot(ctyp_sample_norm_cytok[which(ctyp_sample_norm_cytok$Group_v2%in%c(ord_pval_pos_mono3,'[Mono] Mono3') ),], aes(x =Group_v2, y =",the_gene,"))+
  # geom_boxplot(lwd=1, color='black', fill='grey', alpha=0.5)+scale_x_discrete(limits = c(ord_pval_pos_mono3,'[Mono] Mono3') ) +xlab('')+
  # geom_jitter(stat='identity',inherit.aes = TRUE,shape=16, position=position_jitter(0.2),size=1)+theme_minimal()+ 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 01))+ ylab('Mean ",the_gene," normalized\nexpression by sample')+
  # coord_flip()+guides(y = 'none')+geom_vline(xintercept = 21.5,linetype='dashed')+ggtitle('",the_gene," in Subtypes positively correlated to Mono3')+
  #                        ylim(0,max(ctyp_sample_norm_cytok$",the_gene,"[which(ctyp_sample_norm_cytok$Group_v2%in%c(ord_pval_pos_mono1,'[Mono] Mono1',ord_pval_pos_mono3) ) ] ))+scale_y_continuous(trans='log1p')  )"
  # )))
}
dev.off()
