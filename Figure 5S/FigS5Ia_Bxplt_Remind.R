library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(openxlsx)


# Load the metadata and count matrix from Ngollo cohort
load("./Grouped_objects/RNAseq_REMIND/minipheno_AllezNgollo_postop_goodRemRec_sig.rd")
load("./Grouped_objects/genelist_16nov22_Macs.rd")

miniphenoM0I<-minipheno[which(minipheno$`location:ch1`%in%c("M0I")),]
signatures <- c("data_Siggrp_I","data_Siggrp_II","data_Siggrp_III","data_Siggrp_IV","data_Siggrp_V",
                "data_Siggrp_VI","data_Siggrp_VII","data_Siggrp_VIII","data_Siggrp_IX","data_Siggrp_X",   
                "data_Siggrp_XI","data_Siggrp_XII","data_Siggrp_XIII","data_Siggrp_XIV","data_Siggrp_XV",  
                "data_Siggrp_XVI","data_Siggrp_XVII")
levels(signatures)<-signatures



##---------------------------------------##
## Boxplots w/ Rutgeerts scores (CTM0I)  ##
##---------------------------------------##

bx_rutgeerts<-function(data_set,val2,col_rut,lgd2,plot_title,comp_rut,stat){
  eval(parse(text=paste0("min_y_value<-min(data_set$",val2,")" )))
  eval(parse(text=paste0("max_y_value<-max(data_set$",val2,")" )))
  plot<-ggplot(data_set, aes(x =new_rut, y =!!sym(val2), color =new_rut, fill=new_rut)) +
    geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
    scale_x_discrete(limits = names(col_rut) )+
    geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
    ggprism::theme_prism(base_size = 20)+
    labs(y= sig_rec[lgd2]) +
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.line = element_line(size = rel(1.6)),
          axis.text = element_text(face="bold",size=rel(1.5), color="black"),
          axis.ticks = element_line(size = rel(1.6)),
          axis.ticks.length =unit(16,"pt"),
          axis.title = element_text(face="bold",size=rel(1.8), color="black"),
          axis.text.x = element_text(margin = margin(t = 8)),
          axis.text.y = element_text(margin = margin(r = 8)),
          panel.grid.minor = element_blank(),panel.grid.major = element_blank() )+
    scale_color_manual(values=col_rut, name="Rutgeerts\nendocopic\nrecurrence\nscore")+
    scale_fill_manual(values=col_rut, name="Rutgeerts\nendocopic\nrecurrence\nscore")+
    ggtitle(plot_title)+
    coord_cartesian(ylim = c(min_y_value, max_y_value+(0.1*(max_y_value-min_y_value) )))#Get 10% more space above the plot so the stat text is not cut
  if(stat==T){
    #It does a wilcoxon test to compare means between 2 groups
    return(plot+stat_compare_means(comparisons = comp_rut,size = rel(10),vjust=-0.2,
                                   bracket.size = rel(2)) )
  }else{return(plot)}
}  


##---------------------------------------##


col_rut<-c("i2"="#ba7000",
           "3_4"="#801c01")
comp_rut=list(c("i2","3_4"))

all_sig_comp_BR<-function(sub_set){
  eval(parse(text=paste0("data_set=minipheno",sub_set)))
  data_set$new_rut<-data_set$`rutgeerts:ch1`
  data_set<-data_set[-which(data_set$new_rut%in%c("0","1")),]
  data_set$new_rut[which(data_set$new_rut%in%c("i2a","i2b"))]<-"i2"
  data_set$new_rut[which(data_set$new_rut%in%c("i3","i4"))]<-"3_4"
  n <- length(signatures)
  for (j in 1:n) {
    y <- signatures[j]
    lgd2=paste0("Signature ",gsub("data_Siggrp_","",y))
    print(paste0(lgd2))
    eval(parse(text=paste0("min_y_value<-min(data_set$",y,")" )))
    eval(parse(text=paste0("max_y_value<-max(data_set$",y,")" )))
    print(bx_rutgeerts(data_set =data_set ,val2=y,lgd2=lgd2,plot_title=sub_set,
                       col_rut=col_rut,comp_rut=comp_rut,stat=T))
    
    eval(parse(text=paste0("data_set$currentscore<-as.numeric(data_set$",y,")" )))
    results_wilcoxon<-wilcox.test(data_set$currentscore ~ data_set$new_rut)
    # Save results for spreadsheet:
    writeData(wb, "Wilcox_Signatures", paste0(lgd2," - Wilcoxon test:"), startRow = xl_row, colNames = FALSE)
    xl_row <- xl_row + 1
    writeData(wb, "Wilcox_Signatures", names(results_wilcoxon[1]),startRow = xl_row)
    writeData(wb, "Wilcox_Signatures", names(results_wilcoxon[3]),startCol = 2,startRow = xl_row)
    xl_row <- xl_row + 1
    writeData(wb, "Wilcox_Signatures", results_wilcoxon[c(1,3)],startRow = xl_row)
    xl_row <- xl_row + 1  # 2 lines as separator row
    
  }
}

sig_rec<-c("Signature I"="I_Transition A","Signature II"="II_OXPHOS","Signature III"="Program III","Signature IV"="IV_Transition B","Signature V"="V_Macrophage","Signature VI"="VI_Transition C",
           "Signature VII"="VII_Mono A","Signature VIII"="VIII_Mono B","Signature IX"="IX_Interferon","Signature X"="Program X","Signature XI"="Program XI","Signature XII"="XII_Immunoreg/repair",
           "Signature XIII"="XIII_Inflammatory A","Signature XIV"="XIV_Inflammatory B","Signature XV"="XV_Inflammatory Mono1","Signature XVI"="XVI_Stress","Signature XVII"="XVII_Hypoxia")

wb <- createWorkbook()
addWorksheet(wb, "Wilcox_Signatures")
xl_row <- 1

pdf(paste0("./Figures_print/FigS5Ia_Bxplt_Remind.pdf"),width =6,height = 7.5)
all_sig_comp_BR("M0I")
dev.off()

saveWorkbook(wb, "./Figure 5S/FigS5Ia_Bxplt_prgm_stats.xlsx", overwrite = TRUE)



