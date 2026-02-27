library(ggplot2)
library(dplyr)
library(xlsx)
library(rstatix)
library(ggpubr)
library(openxlsx)


data_set=xlsx::read.xlsx2("./Figure 3/Fig3D_datapoints.xlsx",sheetName = "Feuil1",startRow = 2)
data_set$value<-as.numeric(data_set$value)
data_set$condition<-factor(data_set$condition,levels=c("GM TNF IFNg","GM LPS", "GM LPS TNF","GM LPS TNF IFNg"))

## It's non parametric, it's not always the same patients so KW then Dunn post hoc not Friedman
KW_rstatix<-data_set %>%rstatix::kruskal_test(value ~ condition)
pwc_label <- bquote(paste("pwc: ", bold("Wilcoxon test")) )  #pwc= PairWise Comparison
test_label <- get_test_label(KW_rstatix, detailed = TRUE)
combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))

# KW positive, perform Dunn
dunn_results<-data_set %>%rstatix::dunn_test(value ~ condition, p.adjust.method = "holm")

dunn_results <- dunn_results %>% add_y_position()
dunn_results$y.position<-c(0,0,0,0,90,97)  # better if I manually tell them where to be
dunn_results$padj_round<-round(dunn_results$p.adj,4)

wb <- createWorkbook()
addWorksheet(wb, "Dunn_Prop")
writeData(wb, "Dunn_Prop", KW_rstatix[,c(2:5)], startRow = 2)
writeData(wb, "Dunn_Prop", dunn_results[,c(2:9)], startRow = nrow(KW_rstatix)+4)
saveWorkbook(wb, "./Figure 3/Fig3D_PDL1CD123_culture.xlsx", overwrite = TRUE)



pdf(paste0("./Figures_print/Fig3D_PDL1CD123_culture_illu.pdf"),width = 3,height = 5.5)

  df_sum <- data_set %>%
    group_by(condition) %>%
    summarise(mean  = mean(value, na.rm = TRUE),
              se    = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
              upper = mean + se,
              lower = mean - se, .groups = "drop"  )
  
  ggplot() + geom_col( data = df_sum, aes(x = condition, y = mean, fill = condition, color = condition),
                       alpha = 1,width = 0.8,size=1, position = position_dodge(width = 0.8)) +scale_color_manual(values = rep("black",7))+
    geom_jitter(data = data_set,aes(x = condition, y = value),shape=21, position=position_jitter(width = 0.3),size=2.5, fill="black",stroke=1.25) +

    geom_errorbar(data = df_sum, aes(x = condition, ymin =eval(parse(text=paste0("mean"))), ymax = upper), width = 0.4, size = rel(1))+
    scale_fill_manual(values = c("GM TNF IFNg"="#447228","GM LPS"="#FBDECD", "GM LPS TNF"="#F5C09D","GM LPS TNF IFNg"="#B84410")) +
    scale_x_discrete(limits = levels(data_set$condition)) +
    scale_y_continuous(expand = c(0, 0),n.breaks = 6)+
    ggtitle("PD-L1+ CD123+") + 
    ggprism::theme_prism(base_size = 16)+ labs(y= "% cells")+
    theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),
          axis.text.y = element_text(size=rel(1.75),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(0.5),angle=90,hjust=0.5,vjust=0.5),
          axis.ticks.length.y = unit(15,"points"),axis.ticks.length.x = unit(7.5,"points"),
          plot.margin = unit(c(20, 4, 4, 4),"points"))+
    coord_cartesian(ylim = c(0, 100),clip = 'off')+
    # stat_compare_means(comparisons = my_comparisons, size = 4,bracket.size = rel(0.5), vjust = 0.5,tip.length = 0)+
    stat_pvalue_manual(dunn_results, label = "padj_round", tip.length = 0,hide.ns = TRUE,coord.flip = F,size=rel(5.5),bracket.size = rel(1)) 
    
dev.off()
