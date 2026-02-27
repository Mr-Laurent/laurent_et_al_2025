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




Make_prism_barplot<-function(merged_df,title_name){
  df_sum <- merged_df %>%
    group_by(condition) %>%
    summarise(mean  = mean(value, na.rm = TRUE),
              se    = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
              upper = mean + se, .groups = "drop"  )
  
  ggplot() + geom_col( data = df_sum,
                       aes(x = condition, y = mean, fill = condition, color = condition),
                       alpha = 1,width = 0.8,size=1,
                       position = position_dodge(width = 0.8)) +scale_color_manual(values = rep("black",8))+
    
    geom_jitter(  data = merged_df,size = 2,
                  aes(x = condition, y = value, color = condition ), # shape = batch, 
                  position = position_jitterdodge(
                    jitter.width = 0.20,jitter.height = 0,dodge.width = 0.8 )) +
    geom_errorbar(data = df_sum, aes(x = condition, ymin = mean, ymax = upper), width = 0.4, size = rel(1))+
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




ord_x3=c("GM","GM_TNF","GM_IFNg","GM_TNF_IFNg","nothing",
         "GM_LPS","GM_LPS_TNF","GM_LPS_IFNg","GM_LPS_TNF_IFNg") 
palette2<-c("GM"="#E2F0DA","GM_TNF"="#C6DFB4","GM_IFNg"="#A9D08E","GM_TNF_IFNg"="#548335","nothing"="#black",
            "GM_LPS"="#FCE4D6","GM_LPS_TNF"="#F8CBAD","GM_LPS_IFNg"="#F4B085","GM_LPS_TNF_IFNg"="#C65911")
labels=c("GM","GM_TNF","GM_IFNg","GM_TNF_IFNg","",
         "GM_LPS","GM_LPS_TNF","GM_LPS_IFNg","GM_LPS_TNF_IFNg")

# Table S13: Data points for Figure 3E
cxcl9_data<-xlsx::read.xlsx("./Figure 3/Fig3E_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 1:9)
rownames(cxcl9_data)<-cxcl9_data$Experiment.ID
cxcl9_data$Experiment.ID<-NULL
colnames(cxcl9_data)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(cxcl9_data)) )

il23_data<-xlsx::read.xlsx("./Figure 3/Fig3E_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 11:19)
rownames(il23_data)<-il23_data$Experiment.ID
il23_data$Experiment.ID<-NULL
colnames(il23_data)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(il23_data)) )

ccl7_data<-xlsx::read.xlsx("./Figure 3/Fig3E_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 21:29)
rownames(ccl7_data)<-ccl7_data$Experiment.ID
ccl7_data$Experiment.ID<-NULL
colnames(ccl7_data)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(ccl7_data)) )

pdf("./Figures_print/Fig3E_Barplt_genes.pdf",width =4,height = 7)

merged_df<-reshape2::melt(cxcl9_data)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="CXCL9")+
  scale_y_continuous( breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )

merged_df<-reshape2::melt(il23_data)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="IL-23")+
  scale_y_continuous( breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )

merged_df<-reshape2::melt(ccl7_data)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="CCL7")+
  scale_y_continuous( breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )

dev.off()




# Table S15: Data points for Figure S8E
cxcl10_data<-xlsx::read.xlsx("./Figure 8S/FigS8E_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 1:9)
rownames(cxcl10_data)<-cxcl10_data$Experiment.ID
cxcl10_data$Experiment.ID<-NULL
colnames(cxcl10_data)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(cxcl10_data)) )

il1b_data<-xlsx::read.xlsx("./Figure 8S/FigS8E_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 11:19)
rownames(il1b_data)<-il1b_data$Experiment.ID
il1b_data$Experiment.ID<-NULL
colnames(il1b_data)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(il1b_data)) )

il1a_data<-xlsx::read.xlsx("./Figure 8S/FigS8E_datapoints.xlsx","Feuil1", rowIndex = 2:5, colIndex = 21:29)
rownames(il1a_data)<-il1a_data$Experiment.ID
il1a_data$Experiment.ID<-NULL
colnames(il1a_data)<-gsub("\\.","_", gsub("GM.CSF","GM",colnames(il1a_data)) )


pdf("./Figures_print/FigS8E_Barplt_genes.pdf",width =4,height = 7)

merged_df<-reshape2::melt(cxcl10_data)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="CXCL10")+
  scale_y_continuous(
    breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )+geom_hline(yintercept = max(merged_df$value) ,lty="11", color = "black",size=rel(1) )    # Add limit of detection

merged_df<-reshape2::melt(il1b_data)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="IL-1b")+
  scale_y_continuous(
    breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )

merged_df<-reshape2::melt(il1a_data)
colnames(merged_df)<-c("condition","value")
breaks <- pretty(c(0, max(merged_df$value) ))
Make_prism_barplot(merged_df,title_name="IL-1a")+
  scale_y_continuous(
    breaks = breaks,
    limits = c(0, max(breaks)),expand = expansion(mult = c(0, 0)) # to make a clean y axis cut
  )

dev.off()




