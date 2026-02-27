# Parts from "Subtypes_cor_mono1_3_cytok.R
library(dplyr)
library(reshape2)
library(Matrix)
library(ggplot2)
library(plyr)
library(openxlsx)

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

mlt_spe_cor_pos<-mlt_spe_cor[which(mlt_spe_cor$pos_cor=="Pos"),]
ord_pval_pos<-rownames(mlt_spe_cor_pos[order(mlt_spe_cor_pos$pval,decreasing = T),] )
ord_pval_pos_mono1<-ord_pval_pos
save(ord_pval_pos_mono1,file="./Grouped_objects/ord_pval_pos_mono1.rd")

# Correlation of subtypes to IRM/Mono3 distribution
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

mlt3_spe_cor_pos<-mlt3_spe_cor[which(mlt3_spe_cor$pos_cor=="Pos"),]
ord3_pval_pos<-rownames(mlt3_spe_cor_pos[order(mlt3_spe_cor_pos$pval,decreasing = T),] )
ord_pval_pos_mono3<-ord3_pval_pos
save(ord_pval_pos_mono3,file="./Grouped_objects/ord_pval_pos_mono3.rd")

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
                   "DC2_KIT","DC2_CD207",
                   "Naive","Atypical_Mem","GC_like","Mem_1_NFkB","Mem_2","Mem_3",
                   "IgG","IgA","IgM",
                   "CM-like","Activated","Naive/CM","Effector_CD8s","gd-like_GZMB+","gd-like_GZMB-",
                   "TRM-like_NFkB-induced","CXCR6_TRM-like","CD8s_MAIT_like","Effector_CD8s_GZMK_low",
                   "Mono1","Mono2","Mono3","Mono4","Mono5",
                   "Macro6","Macro7","Macro8","Macro9",
                   "Fib1_Infl. FRC-like","Fib2_Infl.","Fib3_WNT5B+ SOX6+","Fib4_WNT2B+","Fib5_WNT2B+ HGFhi","Fib6_WNT2B+ MME+",
                   "Effector_Tregs","Tregs","Ig-skewed veinous ACKR1+","Prolif Vein. ACKR1+","Arterial CD36-CA4-",
                   "Mature Vein. ACKR1+","Arterial CD36+CA4-","Capillary 1","Capillary 2","Lymphatics",
                   "Vein. ACKR1+",
                   "Plasmablast IgG","Plasmablast IgA")    

# last name update:
new.cluster.ids <- c("aDCb","iaDCa","iaDCb",
                     "DC2 KIT","DC2 CD207",
                     "B Naive","B Atyp Mem","B GC-like","B Mem 1 NFkB","B Mem 2","B Mem 3",
                     "PC IgG","PC IgA","PC IgM",
                     "T CM-like","T Act","T Naive/CM","T CD8 Eff","Tgd-like GZMB+","Tgd-like GZMB-",
                     "TRM-like NFkB","TRM-like CXCR6","T CD8 MAIT-like","T CD8 Eff GZMKlow",
                     "IFIM","IFN Mono","IRM","Early Mono","Late Mono",
                     "IFN Macro","Infl Macro","FOLR2- Macro","FOLR2+ Macro",
                     "Fib Infl FRC-like","Fib Infl","Fib WNT5B SOX6","Fib WNT2B","Fib WNT2B HGFhigh","Fib WNT2B MME",
                     "Treg Eff","Treg","EC Vein Ig-sk","EC vein cycl","EC Art CD36-",
                     "EC Vein Mature","EC Art CD36+","EC Capil 1","EC Capil 2","EC Lymph",
                     "EC Vein",
                     "PB IgG","PB IgA")

mlt_spe_cor$group2<- plyr::mapvalues(x = mlt_spe_cor$group, from = old.cluster.ids, to = new.cluster.ids)

mlt_spe_cor_pos<-mlt_spe_cor[which(mlt_spe_cor$pos_cor=="Pos"),]
ord_pval_pos<-mlt_spe_cor_pos[order(mlt_spe_cor_pos$pval,decreasing = T),"group2"] 

pdf("./Figures_print/Fig2F_mo1mo3_corr.pdf",width = 12, height = 6.5)
print(ggplot(mlt_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
        geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
        geom_bar(lwd=0.7, stat="identity",width = 0.8)+scale_x_discrete(limits = rev(ord_pval_pos) ) + 
        # scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
        # scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
        scale_fill_manual(values="#AAAAAA")+
        scale_color_manual(values="black")+
        theme_minimal()+geom_vline(xintercept = 10.5,linetype="dashed")+ 
        ylab("-log10(p-value)")+xlab("")+ 
        theme(legend.position = "none",axis.title.x =element_blank(),panel.grid.major = element_blank(),
              axis.text.y = element_text(size=rel(3.25),margin = margin(0,8,0,0, unit = "points")),
              axis.text.x = element_text(size=rel(2.25),angle = 90, hjust = 1, vjust = 0.5),
              axis.ticks.length.x = unit(5,"points"),axis.ticks.length.y = unit(5,"points"),
              axis.line.y= element_line(colour="black",size=rel(1.75)),axis.ticks.y = element_line(colour="black",size=rel(1.75)),
              plot.margin = unit(c(40, 4, 4, 4),"points"))+
        labs(fill = "Spearman\ncorrelation\nto Mono1 value",color = "Spearman\ncorrelation\nto Mono1 value")+ 
        scale_y_continuous(expand = c(0, 0)) # removes space between x axis and annotations so the y line starts at 0 not below
)
dev.off()


pdf("./Figures_print/Fig2F_S9A_mo1mo3_corr.pdf",width = 8, height = 7)

# Plot the Top positively-cortrelated subtypes to Mono1 distribution
print(ggplot(mlt_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord_pval_pos ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 19.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto Mono1 value",color = "Spearman\ncorrelation\nto Mono1 value")
)


print(ggplot(mlt_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
        geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
        geom_bar(lwd=0.7, stat="identity",width = 0.8)+scale_x_discrete(limits = rev(ord_pval_pos) ) + 
        # scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
        # scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
        scale_fill_manual(values="#AAAAAA")+
        scale_color_manual(values="black")+
        theme_minimal()+geom_vline(xintercept = 10.5,linetype="dashed")+ 
        ylab("-log10(p-value)")+xlab("")+ 
        theme(legend.position = "none",axis.title.x =element_blank(),panel.grid.major = element_blank(),
              axis.text.y = element_text(size=rel(3.25),angle = 90, hjust = 0, vjust = 1,margin = margin(0,8,0,0, unit = "points")),
              axis.text.y.right =element_text(hjust=0.5),
              axis.title.y.right = element_text(angle = 90),
              axis.text.x = element_text(size=rel(2.25),angle = 90, hjust = 1, vjust = 0.5),
              axis.ticks.length.x = unit(5,"points"),axis.ticks.length.y = unit(5,"points"),
              axis.line.y= element_line(colour="black",size=rel(1.75)),axis.ticks.y = element_line(colour="black",size=rel(1.75)),
              plot.margin = unit(c(40, 4, 4, 4),"points"))+
        labs(fill = "Spearman\ncorrelation\nto IFIM value",color = "Spearman\ncorrelation\nto IFIM value")+ 
        scale_y_continuous(expand = c(0, 0)) # ,position = "right" # removes space between x axis and annotations so the y line starts at 0 not below
)



mlt_spe_cor_pos[rev(match(ord_pval_pos,mlt_spe_cor_pos$group2)),c(8,1:3,5)]

mlt_spe_cor_neg<-mlt_spe_cor[which(mlt_spe_cor$pos_cor=="Neg"),]
ord_pval_neg<-mlt_spe_cor_neg[order(mlt_spe_cor_neg$pval,decreasing = T),"group2"] 

# Plot the Top negatively-cortrelated subtypes to Mono1 distribution
print(ggplot(mlt_spe_cor_neg, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord_pval_neg ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 24.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto IFIM value",color = "Spearman\ncorrelation\nto IFIM value")
)

print(ggplot(mlt_spe_cor_neg, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
        geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
        geom_bar(lwd=0.7, stat="identity",width = 0.8)+scale_x_discrete(limits = rev(ord_pval_neg) ) + 
        # scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
        # scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
        scale_fill_manual(values="#555555")+
        scale_color_manual(values="black")+
        theme_minimal()+geom_vline(xintercept = 12.5,linetype="dashed")+ 
        ylab("-log10(p-value)")+xlab("")+ 
        theme(legend.position = "none",axis.title.x =element_blank(),panel.grid.major = element_blank(),
              axis.text.y = element_text(size=rel(3.25),angle = 90, hjust = 0, vjust = 1,margin = margin(0,8,0,0, unit = "points")),
              axis.text.y.right =element_text(hjust=0.5),
              axis.title.y.right = element_text(size=rel(1.5),angle = 90),
              axis.text.x = element_text(size=rel(1.25),angle = 90, hjust = 1, vjust = 0.5),
              axis.ticks.length.x = unit(5,"points"),axis.ticks.length.y = unit(5,"points"),
              axis.line.y= element_line(colour="black",size=rel(1.75)),axis.ticks.y = element_line(colour="black",size=rel(1.75)),
              plot.margin = unit(c(40, 4, 4, 4),"points"))+
        labs(fill = "Spearman\ncorrelation\nto IFIM value",color = "Spearman\ncorrelation\nto IFIM value")+ 
        scale_y_continuous(expand = c(0, 0),position = "right") # removes space between x axis and annotations so the y line starts at 0 not below
)


mlt3_spe_cor$mo1mo3<-"Nothing"
mlt3_spe_cor$mo1mo3[which(mlt3_spe_cor$group%in%intersect(mono1_pos_sig,mono3_pos_sig))]<-"Sig. in both"
mlt3_spe_cor$mo1mo3[which(mlt3_spe_cor$group%in%setdiff(mono1_pos_sig,mono3_pos_sig))]<-"Sig. in Mono1"
mlt3_spe_cor$mo1mo3[which(mlt3_spe_cor$group%in%setdiff(mono3_pos_sig,mono1_pos_sig))]<-"Sig. in Mono3"

mlt3_spe_cor$group[which(mlt3_spe_cor$group=="[Mast]")] <- "[Mast] Mast"
mlt3_spe_cor$group<-gsub("^\\[.*\\] ","",mlt3_spe_cor$group)

mlt3_spe_cor$group2<- plyr::mapvalues(x = mlt3_spe_cor$group, from = old.cluster.ids, to = new.cluster.ids)

mlt3_spe_cor_pos<-mlt3_spe_cor[which(mlt3_spe_cor$pos_cor=="Pos"),]
ord3_pval_pos<-mlt3_spe_cor_pos[order(mlt3_spe_cor_pos$pval,decreasing = T),"group2"] 

# Plot the Top positively-cortrelated subtypes to Mono3 distribution
print(ggplot(mlt3_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord3_pval_pos ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 21.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto IRM value",color = "Spearman\ncorrelation\nto IRM value")
)

print(ggplot(mlt3_spe_cor_pos, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
        geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
        geom_bar(lwd=0.7, stat="identity",width = 0.8)+scale_x_discrete(limits = rev(ord3_pval_pos) ) + 
        # scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
        # scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
        scale_fill_manual(values="#AAAAAA")+
        scale_color_manual(values="black")+
        theme_minimal()+geom_vline(xintercept = 7.5,linetype="dashed")+ 
        ylab("-log10(p-value)")+xlab("")+ 
        theme(legend.position = "none",axis.title.x =element_blank(),panel.grid.major = element_blank(),
              axis.text.y = element_text(size=rel(3.25),angle = 90, hjust = 0, vjust = 1,margin = margin(0,8,0,0, unit = "points")),
              axis.text.y.right =element_text(hjust=0.5),
              axis.title.y.right = element_text(size=rel(1.5),angle = 90),
              axis.text.x = element_text(size=rel(1.25),angle = 90, hjust = 1, vjust = 0.5),
              axis.ticks.length.x = unit(5,"points"),axis.ticks.length.y = unit(5,"points"),
              axis.line.y= element_line(colour="black",size=rel(1.75)),axis.ticks.y = element_line(colour="black",size=rel(1.75)),
              plot.margin = unit(c(40, 4, 4, 4),"points"))+
        labs(fill = "Spearman\ncorrelation\nto IRM value",color = "Spearman\ncorrelation\nto IRM value")+ 
        scale_y_continuous(expand = c(0, 0),position = "right") # removes space between x axis and annotations so the y line starts at 0 not below
)
mlt3_spe_cor_neg<-mlt3_spe_cor[which(mlt3_spe_cor$pos_cor=="Neg"),]
ord3_pval_neg<-mlt3_spe_cor_neg[order(mlt3_spe_cor_neg$pval,decreasing = T),"group2"]
# Plot the Top negatively-cortrelated subtypes to Mono3 distribution
print(ggplot(mlt3_spe_cor_neg, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
  geom_bar(lwd=0.5, stat="identity",  alpha=0.5)+scale_x_discrete(limits = ord3_pval_neg ) + 
  scale_fill_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive")  )+ 
  scale_color_manual(values=c("Neg"="#7899EE","Pos"="#EF5513"),labels=c("Negative","Positive") )+
  theme_minimal()+geom_vline(xintercept = 34.5,linetype="dashed")+ geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
  coord_flip()+ ylab("-log10(p-value)")+xlab("")+ 
  labs(fill = "Spearman\ncorrelation\nto IRM value",color = "Spearman\ncorrelation\nto IRM value")
)

print(ggplot(mlt3_spe_cor_neg, aes(x =group2, y =mlogpval, color=pos_cor, fill=pos_cor))+
        geom_hline(yintercept = -log10(0.05),linetype="dashed",alpha=0.1)+
        geom_bar(lwd=0.7, stat="identity",width = 0.8)+scale_x_discrete(limits = rev(ord3_pval_neg) ) + 
        scale_fill_manual(values="#555555")+
        scale_color_manual(values="black")+
        theme_minimal()+geom_vline(xintercept = 3.5,linetype="dashed")+ 
        ylab("-log10(p-value)")+xlab("")+ 
        theme(legend.position = "none",axis.title.x =element_blank(),panel.grid.major = element_blank(),
              axis.text.y = element_text(size=rel(3.25),angle = 90, hjust = 0, vjust = 1,margin = margin(0,8,0,0, unit = "points")),
              axis.text.y.right =element_text(hjust=0.5),
              axis.title.y.right = element_text(size=rel(1.5),angle = 90),
              axis.text.x = element_text(size=rel(1.25),angle = 90, hjust = 1, vjust = 0.5),
              axis.ticks.length.x = unit(5,"points"),axis.ticks.length.y = unit(5,"points"),
              axis.line.y= element_line(colour="black",size=rel(1.75)),axis.ticks.y = element_line(colour="black",size=rel(1.75)),
              plot.margin = unit(c(40, 4, 4, 4),"points"))+
        labs(fill = "Spearman\ncorrelation\nto IRM value",color = "Spearman\ncorrelation\nto IRM value")+ 
        scale_y_continuous(expand = c(0, 0),position = "right") # removes space between x axis and annotations so the y line starts at 0 not below
)

dev.off()

# Save the values in an excel file:
wb <- createWorkbook()
addWorksheet(wb, "Spearman_corr")

writeData(wb, "Spearman_corr", mlt_spe_cor_pos[rev(match(ord_pval_pos,mlt_spe_cor_pos$group2)),c(8,1:3,5)],
          startCol = 1,startRow = 1)
writeData(wb, "Spearman_corr", mlt_spe_cor_neg[rev(match(ord_pval_neg,mlt_spe_cor_neg$group2)),c(8,1:3,5)],
          startCol = 7,startRow = 1)
writeData(wb, "Spearman_corr", mlt3_spe_cor_pos[rev(match(ord3_pval_pos,mlt3_spe_cor_pos$group2)),c(8,1:3,5)],
          startCol = 13,startRow = 1)
writeData(wb, "Spearman_corr", mlt3_spe_cor_neg[rev(match(ord3_pval_neg,mlt3_spe_cor_neg$group2)),c(8,1:3,5)],
          startCol = 19,startRow = 1)

saveWorkbook(wb, "./Figure 2/Fig2F_S9A_mo1mo3_corr.xlsx", overwrite = TRUE)
