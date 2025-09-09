
# Parts from "Subtypes_cor_mono1_3_cytok.R

library(dplyr)
library(reshape2)
library(Matrix)
library(ggplot2)
library(plyr)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
load("Grouped_objects/freqs3_no69_noadj_nolog.rd")  # subtypes frequency by patients


###------------------------------------------------------------------### 
#### Function for spearman correlation between subtypes proportions #### 
###------------------------------------------------------------------### 
spearman_pvalues <- function(data) {
  p_values <- matrix(NA, ncol(data), ncol(data))
  colnames(p_values) <- rownames(p_values) <- colnames(data)
  for (i in 1:(ncol(data)-1)) {
    for (j in (i+1):ncol(data)) {
      test_result <- cor.test(data[, i], data[, j], method = "spearman")
      p_values[i, j] <- test_result$p.value
      p_values[j, i] <- test_result$p.value
    }
  } 
  return(p_values)
}
###-------------------------------------------------------------### 



cormat_sp=cor((freqs3),use="pairwise.complete.obs", method="spearman") 
p_values3 <- spearman_pvalues(freqs3)

# Correlation of subtypes to Mono1 distribution
mlt_spe_cor<-as.data.frame(reshape2::melt(cormat_sp["[MoMac] Mono1",-which(colnames(cormat_sp)=="[MoMac] Mono1")])  )
mlt_spe_pval<-as.data.frame(reshape2::melt(p_values3["[MoMac] Mono1",-which(colnames(cormat_sp)=="[MoMac] Mono1")]) )
mlt_spe_cor$pval <- as.numeric(mlt_spe_pval$value)
mlt_spe_cor$mlogpval <- -log10(as.numeric(mlt_spe_pval$value))
colnames(mlt_spe_cor)<-c("spear_cor","pval","mlogpval")
mlt_spe_cor$group<-rownames(mlt_spe_cor)
#Add info of significance + correlation direction
mlt_spe_cor <- mlt_spe_cor %>%
  mutate(pval_sig = if_else(pval < 0.05, "Yes", "No"))
mlt_spe_cor <- mlt_spe_cor %>%
  mutate(pos_cor = if_else(spear_cor < 0, "Neg", "Pos"))

# Correlation of subtypes to Mono3 distribution
mlt3_spe_cor<-as.data.frame(reshape2::melt(cormat_sp["[MoMac] Mono3",-which(colnames(cormat_sp)=="[MoMac] Mono3")])  )
mlt3_spe_pval<-as.data.frame(reshape2::melt(p_values3["[MoMac] Mono3",-which(colnames(cormat_sp)=="[MoMac] Mono3")]) )
mlt3_spe_cor$pval <- as.numeric(mlt3_spe_pval$value)
mlt3_spe_cor$mlogpval <- -log10(as.numeric(mlt3_spe_pval$value))
colnames(mlt3_spe_cor)<-c("spear_cor","pval","mlogpval")
mlt3_spe_cor$group<-rownames(mlt3_spe_cor)
mlt3_spe_cor <- mlt3_spe_cor %>%
  mutate(pval_sig = if_else(pval < 0.05, "Yes", "No"))
mlt3_spe_cor <- mlt3_spe_cor %>%
  mutate(pos_cor = if_else(spear_cor < 0, "Neg", "Pos"))


ord_pval<-rownames(mlt_spe_cor[order(mlt_spe_cor$pval,decreasing = T),] )

mono1_pos_sig<-mlt_spe_cor$group[mlt_spe_cor$pval_sig=="Yes"&mlt_spe_cor$pos_cor=="Pos"]
mono3_pos_sig<-mlt3_spe_cor$group[mlt3_spe_cor$pval_sig=="Yes"&mlt3_spe_cor$pos_cor=="Pos"]

mlt_spe_cor$mo1mo3<-"Nothing"
mlt_spe_cor$mo1mo3[which(mlt_spe_cor$group%in%intersect(mono1_pos_sig,mono3_pos_sig))]<-"Sig. in both"
mlt_spe_cor$mo1mo3[which(mlt_spe_cor$group%in%setdiff(mono1_pos_sig,mono3_pos_sig))]<-"Sig. in Mono1"
mlt_spe_cor$mo1mo3[which(mlt_spe_cor$group%in%setdiff(mono3_pos_sig,mono1_pos_sig))]<-"Sig. in Mono3"


mlt_spe_cor$group[which(mlt_spe_cor$group=="[Mast]")] <- "[Mast] Mast"
mlt_spe_cor$group<-gsub("^\\[.*\\] ","",mlt_spe_cor$group)
old.cluster.ids<-c("Act. DC","Infl. DC3","Act. DC Infl.",
                   "Naive","Atypical_Mem","GC_like","Mem_1_NFkB","Mem_2","Mem_3",
                   "IgG","IgA","IgM",
                   "CM-like","Activated","Naive/CM","Effector_CD8s","gd-like_GZMB+","gd-like_GZMB-")     

new.cluster.ids <- c("aDC","iaDCa","iaDCb",
                     "Naive_Bcells","Atypical_Mem_Bcells","GC_like_Bcells","Mem_1_NFkB_Bcells","Mem_2_Bcells","Mem_3_Bcells",
                     "IgG_Plasma_cells","IgA_Plasma_cells","IgM_Plasma_cells",
                     "CM-like_Tcells","Activated_Tcells","Naive/CM_Tcells","Effector_CD8s_Tcells","T_gd-like_GZMB+","T_gd-like_GZMB-")

mlt_spe_cor$group2<- plyr::mapvalues(x = mlt_spe_cor$group, from = old.cluster.ids, to = new.cluster.ids)

mlt_spe_cor_pos<-mlt_spe_cor[which(mlt_spe_cor$pos_cor=="Pos"),]
ord_pval_pos<-mlt_spe_cor_pos[order(mlt_spe_cor_pos$pval,decreasing = T),"group2"] 



pdf("./Figure 2/Fig2H_mo1mo3_corr.pdf",width = 8, height = 7)

# Plot the Top positively-cortrelated subtypes to Mono1 distribution
print(ggplot(mlt_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord_pval_pos ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 19.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto Mono1 value",color = "Spearman\ncorrelation\nto Mono1 value")
)


mlt_spe_cor_neg<-mlt_spe_cor[which(mlt_spe_cor$pos_cor=="Neg"),]
ord_pval_neg<-mlt_spe_cor_neg[order(mlt_spe_cor_neg$pval,decreasing = T),"group2"] 

# Plot the Top negatively-cortrelated subtypes to Mono1 distribution
print(ggplot(mlt_spe_cor_neg, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord_pval_neg ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 24.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto Mono1 value",color = "Spearman\ncorrelation\nto Mono1 value")
)

mlt3_spe_cor$mo1mo3<-"Nothing"
mlt3_spe_cor$mo1mo3[which(mlt3_spe_cor$group%in%intersect(mono1_pos_sig,mono3_pos_sig))]<-"Sig. in both"
mlt3_spe_cor$mo1mo3[which(mlt3_spe_cor$group%in%setdiff(mono1_pos_sig,mono3_pos_sig))]<-"Sig. in Mono1"
mlt3_spe_cor$mo1mo3[which(mlt3_spe_cor$group%in%setdiff(mono3_pos_sig,mono1_pos_sig))]<-"Sig. in Mono3"

mlt3_spe_cor$group[which(mlt3_spe_cor$group=="[Mast]")] <- "[Mast] Mast"
mlt3_spe_cor$group<-gsub("^\\[.*\\] ","",mlt3_spe_cor$group)

old.cluster.ids<-c("Act. DC","Infl. DC3","Act. DC Infl.",
                   "Naive","Atypical_Mem","GC_like","Mem_1_NFkB","Mem_2","Mem_3",
                   "IgG","IgA","IgM",
                   "CM-like","Activated","Naive/CM","Effector_CD8s","gd-like_GZMB+","gd-like_GZMB-")     

new.cluster.ids <- c("aDC","iaDCa","iaDCb",
                     "Naive_Bcells","Atypical_Mem_Bcells","GC_like_Bcells","Mem_1_NFkB_Bcells","Mem_2_Bcells","Mem_3_Bcells",
                     "IgG_Plasma_cells","IgA_Plasma_cells","IgM_Plasma_cells",
                     "CM-like_Tcells","Activated_Tcells","Naive/CM_Tcells","Effector_CD8s_Tcells","T_gd-like_GZMB+","T_gd-like_GZMB-")

mlt3_spe_cor$group2<- plyr::mapvalues(x = mlt3_spe_cor$group, from = old.cluster.ids, to = new.cluster.ids)

mlt3_spe_cor_pos<-mlt3_spe_cor[which(mlt3_spe_cor$pos_cor=="Pos"),]
ord_pval_pos<-mlt3_spe_cor_pos[order(mlt3_spe_cor_pos$pval,decreasing = T),"group2"] 

# Plot the Top positively-cortrelated subtypes to Mono3 distribution
print(ggplot(mlt3_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord_pval_pos ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 21.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto Mono3 value",color = "Spearman\ncorrelation\nto Mono3 value")
)

mlt3_spe_cor_neg<-mlt3_spe_cor[which(mlt3_spe_cor$pos_cor=="Neg"),]
ord_pval_neg<-mlt3_spe_cor_neg[order(mlt3_spe_cor_neg$pval,decreasing = T),"group2"]
# Plot the Top negatively-cortrelated subtypes to Mono3 distribution
print(ggplot(mlt3_spe_cor_neg, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord_pval_neg ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 34.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto Mono3 value",color = "Spearman\ncorrelation\nto Mono3 value")
)

dev.off()



