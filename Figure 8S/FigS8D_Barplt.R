library(ggplot2)
library(ggpubr)
library(patchwork)
library(Rfast)
library(xlsx)
library(edgeR)
library(ggrepel)
library(factoextra)
library(dplyr)
library(scales)
library(openxlsx)
library(rstatix)



Make_prism_barplot<-function(merged_df,title_name){
  df_sum <- merged_df %>%
    group_by(condition) %>%
    summarise(mean  = mean(value, na.rm = TRUE),
              se    = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
              upper = mean + se,  lower = mean - se, .groups = "drop"  )
  
  ggplot() + geom_col( data = df_sum,
                       aes(x = condition, y = mean, fill = condition, color = condition),
                       alpha = 1,width = 0.8,size=1,
                       position = position_dodge(width = 0.8)) +scale_color_manual(values = rep("black",8))+
    
    geom_jitter(  data = merged_df,size = 2,
                  aes(x = condition, y = value, color = condition ), # shape = batch, 
                  position = position_jitterdodge(
                    jitter.width = 0.20,jitter.height = 0,dodge.width = 0.8 )) +
    geom_errorbar(data = df_sum, aes(x = condition, ymin = lower, ymax = upper), width = 0.4, size = rel(1))+
    scale_fill_manual(values = palette2) +
    scale_x_discrete(limits = ord_x3, labels= labels) +
    ggtitle(title_name) + 
    ggprism::theme_prism(base_size = 16)+ labs(y= title_name)+
    theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),
          axis.text.y = element_text(size=rel(1.5),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(1),angle=90,hjust=1,vjust=0.5),
          axis.ticks.length = unit(15,"points"),
          plot.margin = unit(c(20, 4, 4, 4),"points"))
  
}


ord_x3=c("GM_LPS","GM_LPS_TNF","GM_IFN_TNF","GM_LPS_TNF_IFN") 
palette2<-c("GM_LPS"="white","GM_LPS_TNF"="white","GM_IFN_TNF"="white","GM_LPS_TNF_IFN"="white")
labels=c("GM_LPS","GM_LPS_TNF","GM_IFNg_TNF","GM_LPS_TNF_IFNg")


# Table S14: Data points for Figure S8D
pct_pc<-xlsx::read.xlsx("./Figure 8S/FigS8D_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 1:5)
rownames(pct_pc)<-pct_pc$Experiment.ID
pct_pc$Experiment.ID<-NULL
colnames(pct_pc)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(pct_pc)) )

pct_p<-xlsx::read.xlsx("./Figure 8S/FigS8D_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 7:11)
rownames(pct_p)<-pct_p$Experiment.ID
pct_p$Experiment.ID<-NULL
colnames(pct_p)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(pct_p)) )

pct_c<-xlsx::read.xlsx("./Figure 8S/FigS8D_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 13:17)
rownames(pct_c)<-pct_c$Experiment.ID
pct_c$Experiment.ID<-NULL
colnames(pct_c)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(pct_c)) )


data_set<-reshape2::melt(as.matrix(pct_pc))
colnames(data_set)<-c("sample","variable","value")

## Here all conditios use the same samples, so even if non parametric, it's paired, so use Friedman's not KW :
Fried_rstatix<-data_set %>%rstatix::friedman_test(value ~ variable |sample)
pwc_label <- bquote(paste("pwc: ", bold("Friedman's test")) )  #pwc= PairWise Comparison
test_label <- get_test_label(Fried_rstatix, detailed = TRUE)
combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))

# Fried is stat, perform paired Wilcoxon signed-rank test ?  (Dunn is not for paired)
# It's best test but we only have 3 samples so we will never get anything lower than p=0.25, it's useless to do a stat test with n=3

wb_f <- createWorkbook()
addWorksheet(wb_f, "Fried_Prop")
writeData(wb_f, "Fried_Prop","Friedman test on % of PD-L1+ CD123+ cells" , startRow = 1)
writeData(wb_f, "Fried_Prop", Fried_rstatix[,c(2:5)], startRow = 2)


pdf("./Figures_print/FigS8D_Barplt_pct.pdf",width =7,height = 6)

merged_df<-reshape2::melt(pct_pc)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value)+5 ))  #+5 here to have nicer plot
Make_prism_barplot(merged_df,title_name="% of PD-L1+ CD123+ cells")+
  scale_y_continuous(
    breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)))+ # to make a clean y axis cut 
       scale_fill_manual(values = alpha(c("white", "white", "white", "white"), 0) )
#+
  # stat_pvalue_manual(dunn_results, label = "padj_round", tip.length = 0,hide.ns = TRUE,coord.flip = F,size=rel(5.5),bracket.size = rel(1))  # NOT DUNN
  # 


data_set<-reshape2::melt(as.matrix(pct_p))
colnames(data_set)<-c("sample","variable","value")
Fried_rstatix<-data_set %>%rstatix::friedman_test(value ~ variable |sample)
pwc_label <- bquote(paste("pwc: ", bold("Friedman's test")) )  #pwc= PairWise Comparison
test_label <- get_test_label(Fried_rstatix, detailed = TRUE)
combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))

writeData(wb_f, "Fried_Prop","Friedman test on % of PD-L1+ cells" , startRow = 4)
writeData(wb_f, "Fried_Prop", Fried_rstatix[,c(2:5)], startRow = 5)

merged_df<-reshape2::melt(pct_p)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="% of PD-L1+ cells")+
  scale_y_continuous(
    breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )+ scale_fill_manual(values = alpha(c("white", "white", "white", "white"), 0) )


data_set<-reshape2::melt(as.matrix(pct_c))
colnames(data_set)<-c("sample","variable","value")
Fried_rstatix<-data_set %>%rstatix::friedman_test(value ~ variable |sample)
pwc_label <- bquote(paste("pwc: ", bold("Friedman's test")) )  #pwc= PairWise Comparison
test_label <- get_test_label(Fried_rstatix, detailed = TRUE)
combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))

writeData(wb_f, "Fried_Prop","Friedman test on % of CD123+ cells" , startRow = 7)
writeData(wb_f, "Fried_Prop", Fried_rstatix[,c(2:5)], startRow = 8)


merged_df<-reshape2::melt(pct_c)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="% of CD123+ cells")+
  scale_y_continuous(
    breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )+ scale_fill_manual(values = alpha(c("white", "white", "white", "white"), 0) )

dev.off()


saveWorkbook(wb_f, "./Figure 8S/FigS8D_PDL1CD123.xlsx", overwrite = TRUE)

