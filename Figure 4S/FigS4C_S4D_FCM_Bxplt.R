# Fig S4C/D  02/12/2025
#From Table S6:

library(ggplot2)
library(dplyr)
library(xlsx)
library(rstatix)
library(ggpubr)


palette<-c("#999999","#000000")
palette2<-c("black","white")


prismplot<-function(y_title,y_val,p_v_lab,labels,plt_title,psiz=25){
  
  data_set=xlsx::read.xlsx2("./Figure 4S/FigS4C_S4D_Datapoints.xlsx",sheetName = plt_title)
  
  if(length(which(data_set$Sample==""))>0){ data_set<-data_set[-which(data_set$Sample==""),]}
  
  data_set$value<-as.numeric(data_set$value)
  data_set$category<-factor(data_set$condition, levels=c("CD uninflamed","CD inflamed"))
  
  #Parform Wilcoxon rank sum test (a.k.a. Mann-Whitney U test)
  wilc_rstatix<-data_set %>%rstatix::wilcox_test(value ~ category)
  df_pvals <- wilc_rstatix %>% add_xy_position() %>%  add_significance("p") %>%mutate(category="i2") 
  df_pvals$y.position<-max(data_set$value)*1.1  # Automatic position is varying too much, use max + 10% space
  df_pvals$p_round<-round(df_pvals$p,4)
  df_pvals$xmin="CD uninflamed"
  df_pvals$xmax="CD inflamed"
  
  long_df2=data_set
  min_y=min(data_set$value)
  max_y=max(data_set$value)*1.1
  
  print(ggplot(long_df2, aes(x =category, y =!!sym(y_val), fill=category))+
          ggtitle(label=paste0(plt_title,"\n"))+
          stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
          geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
          scale_x_discrete(limits = levels(long_df2$category) , labels= labels)+
          scale_color_manual(values = palette)+ scale_y_continuous(n.breaks = 5) +#breaks = c(0.25,0.5,0.75,1), expand = c(0, 0) ) +
          ggprism::theme_prism(base_size = 16)+ labs(y= y_title)+
          stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = category), fatten = 4)+ # Visualize the median on top, allows to adapt color
          scale_color_manual(values = palette2 ) +
          geom_jitter(shape=21, position=position_jitter(width = 0.20),size=4, fill="white",stroke=1.5) +
          stat_pvalue_manual(df_pvals, label = p_v_lab, tip.length = 0,hide.ns = FALSE, bracket.size=rel(1.4), label.size = rel(psiz))+
          theme(legend.position = "none",axis.title.x =element_blank(),
                axis.text.y = element_text(size=rel(2.5),margin = margin(0,8,0,0, unit = "points")),
                # axis.text.x = element_text(size=rel(1)),
                axis.text.x = element_blank(),
                axis.ticks.length.y = unit(15,"points"), axis.ticks.length.x = unit(0,"points"),
                plot.margin = unit(c(20, 4, 4, 4),"points"))+   #add margin above plot so asterisks are not cropped
          coord_cartesian(ylim = c(0, max_y),clip = 'off') )    #cut the plot at 1
}

pdf("./Figures_print/FigS4C_FCM_Bxplt.pdf",width =4,height =7)
prismplot(y_title="% of mo-macs",y_val="value",plt_title="Sinai",
          p_v_lab="p.signif", labels=c("CD\nuninflamed","CD\ninflamed"))
prismplot(y_title="% of mo-macs",y_val="value",plt_title="Sinai",
          p_v_lab="p_round", psiz=10, labels=c("CD\nuninflamed","CD\ninflamed"))
prismplot(y_title="% of mo-macs",y_val="value",plt_title="Nantes",
          p_v_lab="p.signif", labels=c("CD\nuninflamed","CD\ninflamed"))
prismplot(y_title="% of mo-macs",y_val="value",plt_title="Nantes",
          p_v_lab="p_round", psiz=10, labels=c("CD\nuninflamed","CD\ninflamed"))
dev.off()






#-----------------------------------------------#

# Fig S4D



library(ggplot2)
library(dplyr)
library(xlsx)
library(rstatix)
library(ggpubr)

palette<-c("#999999","#000000")
palette2<-c("black","white")


prismplot<-function(y_title,y_val,p_v_lab,labels,plt_title,psiz=25){
  
  data_set=xlsx::read.xlsx2("./Figure 4S/FigS4D_Datapoints.xlsx",sheetName = plt_title)
  if(length(which(data_set$Sample==""))>0){ data_set<-data_set[-which(data_set$Sample==""),]}
  data_set$value<-as.numeric(data_set$value)
  data_set$category<-factor(data_set$condition, levels=c("CD uninflamed","CD inflamed"))
  wilc_rstatix<-data_set %>%rstatix::wilcox_test(value ~ category)
  df_pvals <- wilc_rstatix %>% add_xy_position() %>%  add_significance("p") %>%mutate(category="i2")
  df_pvals$y.position<-max(data_set$value)*1.1  # Automatic position is varying too much, use max + 10% space
  df_pvals$p_round<-round(df_pvals$p,4)
  
  long_df2=data_set
  min_y=min(data_set$value)
  max_y=max(data_set$value)*1.1
  
  print(ggplot(long_df2, aes(x =category, y =!!sym(y_val), fill=category))+
          ggtitle(label=paste0(plt_title,"\n"))+
          stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
          geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
          scale_x_discrete(limits = levels(long_df2$category) , labels= labels)+
          scale_color_manual(values = palette)+ scale_y_continuous(n.breaks = 5) +#breaks = c(0.25,0.5,0.75,1), expand = c(0, 0) ) +
          ggprism::theme_prism(base_size = 16)+ labs(y= y_title)+
          stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = category), fatten = 4)+ # Visualize the median on top, allows to adapt color
          scale_color_manual(values = palette2 ) +
          geom_jitter(shape=21, position=position_jitter(width = 0.20),size=4, fill="white",stroke=1.5) +
          stat_pvalue_manual(df_pvals, label = p_v_lab, tip.length = 0,hide.ns = FALSE, bracket.size=rel(1.4), label.size = rel(psiz))+
          theme(legend.position = "none",axis.title.x =element_blank(),
                axis.text.y = element_text(size=rel(2.5),margin = margin(0,8,0,0, unit = "points")),
                # axis.text.x = element_text(size=rel(1)),
                axis.text.x = element_blank(),
                axis.ticks.length.y = unit(15,"points"), axis.ticks.length.x = unit(0,"points"),
                plot.margin = unit(c(20, 4, 4, 4),"points"))+   #add margin above plot so asterisks are not cropped
          coord_cartesian(ylim = c(0, max_y),clip = 'off') )    #cut the plot at 1
}

pdf("./Figures_print/FigS4D_FCM_Bxplt.pdf",width =4,height = 7)
prismplot(y_title="% of mo-macs",y_val="value",plt_title="Nantes",
          p_v_lab="p.signif", labels=c("CD\nuninflamed","CD\ninflamed"))
prismplot(y_title="% of mo-macs",y_val="value",plt_title="Nantes",
          p_v_lab="p_round", psiz=10, labels=c("CD\nuninflamed","CD\ninflamed"))
dev.off()

