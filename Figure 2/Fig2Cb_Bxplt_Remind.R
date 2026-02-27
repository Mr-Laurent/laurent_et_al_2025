library(Rfast)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(rstatix)


# Load the metadata and count matrix from Ngollo cohort
load("./Grouped_objects/RNAseq_REMIND/minipheno_AllezNgollo_postop_goodRemRec_sig.rd")
load("./Grouped_objects/genelist_16nov22_Macs.rd")

miniphenoM0I<-minipheno[which(minipheno$`location:ch1`%in%c("M0I")),]
signatures <- c("data_Siggrp_XV")
levels(signatures)<-signatures



##---------------------------------------##
## Boxplots w/ Rutgeerts scores (CTM0I)  ##
##---------------------------------------##

bx_rutgeerts<-function(data_set,val2,col_rut,lgd2,plot_title,comp_rut,stat){
  plot<-ggplot(data_set, aes(x =new_rut, y =!!sym(val2), color =new_rut, fill=new_rut)) +
    geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
    scale_x_discrete(limits = names(col_rut) )+
    geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
    ggprism::theme_prism(base_size = 25)+
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30))+
    labs(y= lgd2) +
    scale_color_manual(values=col_rut, name="Rutgeerts\nendocopic\nrecurrence\nscore")+
    scale_fill_manual(values=col_rut, name="Rutgeerts\nendocopic\nrecurrence\nscore")+
    ggtitle(plot_title)
  if(stat==T){
    #It does a wilcoxon test to compare means between 2 groups
    return(plot+stat_compare_means(comparisons = comp_rut,size = 6,
                                   bracket.size = rel(2)) )
  }else{return(plot)}
}  


##---------------------------------------##


col_rut<-c("i2"="#ba7000",
           "i3_4"="#801c01")
comp_rut=list(c("i2","i3_4"))


all_sig_comp_BR<-function(sub_set){
  eval(parse(text=paste0("data_set=minipheno",sub_set)))
  data_set$new_rut<-data_set$`rutgeerts:ch1`
  data_set<-data_set[-which(data_set$new_rut%in%c("0","1")),]
  data_set$new_rut[which(data_set$new_rut%in%c("i2a","i2b"))]<-"i2"
  data_set$new_rut[which(data_set$new_rut%in%c("i3","i4"))]<-"i3_4"
  n <- length(signatures)
  for (j in 1:n) {
    y <- signatures[j]
    # lgd2=paste0("Signature ",gsub("data_Siggrp_","",y))
    lgd2="XV_Infl._III"
    print(paste0(lgd2))
    print(bx_rutgeerts(data_set=data_set,val2=y,lgd2=lgd2,plot_title=sub_set,
                       col_rut=col_rut,comp_rut=comp_rut,stat=T))
  }
}




# pdf("./Figure 2/Fig2Da_XV_Bxplt.pdf",width =6,height = 7.5)
all_sig_comp_BR("M0I")
# dev.off()


data_set=miniphenoM0I
data_set$category<-data_set$`rutgeerts:ch1`
data_set<-data_set[-which(data_set$category%in%c("0","1")),]
data_set$category[which(data_set$category%in%c("i2a","i2b"))]<-"i2"
data_set$category[which(data_set$category%in%c("i3","i4"))]<-"i3_4"
# mpM0I<-data_set[,12:29]
# write.csv2(mpM0I,file="Fig2Dremind_datapoints.csv",quote=F,row.names=T)




palette<-c("#999999","#000000")
palette2<-c("black","white")
data_set$category<-factor(data_set$category, levels=c("i2","i3_4"))

wilc_rstatix<-data_set %>%rstatix::wilcox_test(data_Siggrp_XV ~ category)
df_pvals <- wilc_rstatix %>% add_xy_position() %>%  add_significance("p") %>%mutate(category="i2") 
df_pvals$y.position<-df_pvals$y.position+0.05 
df_pvals$p_round<-round(df_pvals$p,4)

prismplot3<-function(long_df2,y_title,y_val,p_v_lab,min_y,max_y, psiz=25){
  ggplot(long_df2, aes(x =category, y =!!sym(y_val), fill=category))+
    stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
    geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
    scale_x_discrete(limits = levels(long_df2$category) )+  # , labels=c("U","I","R")
    scale_color_manual(values = palette)+ scale_y_continuous(n.breaks = 4) +#breaks = c(0.25,0.5,0.75,1)) +
    ggprism::theme_prism(base_size = 16)+ labs(y= y_title)+
    stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = category), fatten = 4)+ # Visualize the median on top, allows to adapt color
    scale_color_manual(values = palette2 ) +
    geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
    stat_pvalue_manual(df_pvals, label = p_v_lab, tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(psiz))+
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.text.y = element_text(size=rel(2.75),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(1)),
          axis.ticks.length = unit(15,"points"),
          plot.margin = unit(c(40, 4, 4, 4),"points"))+   #add margin above plot so asterisks are not cropped
    coord_cartesian(ylim = c(min_y, max_y),clip = 'off')#+  #cut the plot at 1
}

pdf("./Figures_print/Fig2Cb_XV_Bxplt.pdf",width =4,height = 8)
prismplot3(long_df2=data_set,y_title="XV_IFIM prog",y_val="data_Siggrp_XV",
           p_v_lab="p.signif",min_y=min(data_set$data_Siggrp_XV),max_y=max(data_set$data_Siggrp_XV)+0.05)
prismplot3(long_df2=data_set,y_title="XV_IFIM prog",y_val="data_Siggrp_XV",
           p_v_lab="p_round",psiz=10,min_y=min(data_set$data_Siggrp_XV),max_y=max(data_set$data_Siggrp_XV)+0.05)
dev.off()

