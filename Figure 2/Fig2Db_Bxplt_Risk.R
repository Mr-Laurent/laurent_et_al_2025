library(Rfast)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(rstatix)


# Load the metadata and count matrix from RISK cohort
load("Grouped_objects/r_dsgn_InfUninf_cd68norm_Momacsscores.rd")
r_dsgn_Inf<-r_dsgn[which(r_dsgn$status%in%c("Inflamed")),]
r_dsgn_Uninf<-r_dsgn[which(r_dsgn$status%in%c("Uninflamed")),]

r_dsgn_r<-r_dsgn[r_dsgn$response%in%c("NR","R"),]
r_dsgn_r_Inf<-r_dsgn_r[which(r_dsgn_r$status%in%c("Inflamed")),]
r_dsgn_r_Uninf<-r_dsgn_r[which(r_dsgn_r$status%in%c("Uninflamed")),]
signatures<-c("Siggrp_XV")

##------------------------------##
## Boxplots w/ Signature scores ##
##------------------------------##

bx_rep<-function(data_set,sub_set,val2,col_rep,lgd2,plot_title,comp_rep,stat){
  #val_2<-paste0("data",sub_set,"_",val2)
  val_2<-val2
  plot<-ggplot(data_set, aes(x =response, y =!!sym(val_2), color =response, fill=response)) +
    geom_boxplot(lwd=2, alpha=0.5,outlier.alpha = 0)+
    scale_x_discrete(limits = names(col_rep) )+
    geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
    ggprism::theme_prism(base_size = 25)+
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30))+
    labs(y= lgd2) +
    scale_color_manual(values=col_rep, name="Response")+
    scale_fill_manual(values=col_rep, name="Response")+
    ggtitle(plot_title)
  if(stat==T){
    #It does a wilcoxon test to compare means between 2 groups
    return(plot+stat_compare_means(comparisons = comp_rep,size = 6,
                                   bracket.size = rel(2)) )
  }else{return(plot)}
}


##---------------------------------------##

#### Signature Boxplot R vs NR ####

col_rep<-c("R"="#20A1DA","NR"="#BB2820" )
comp_rep=list(c("R","NR"))

all_sig_comp_BR<-function(sub_set){
  eval(parse(text=paste0("data_set=r_dsgn_r",sub_set)))
  # signatures<-colnames(r_dsgn_r)[grep(paste0("dataCD68n_Siggrp"),colnames(r_dsgn_r))]
  signatures<-"dataCD68n_Siggrp_XV"
  n <- length(signatures)
  for (j in 1:n) {
    y <- signatures[j]
    # lgd2=paste0("Signature ",gsub(paste0("dataCD68n_Siggrp_"),"",y))
    lgd2="XV_Infl._III"
    print(lgd2)
    print(bx_rep(data_set =data_set,sub_set=sub_set,val2=y,lgd2=lgd2,plot_title=sub_set,
                 col_rep=col_rep,comp_rep=comp_rep,stat=T))
  }
}




# pdf(paste0("./Figure 2/Fig2Db_XV_Bxplt.pdf"),width =6,height = 7.5)
# all_sig_comp_BR("_Inf")
# dev.off()

data_set=r_dsgn_r_Inf
# r_inf<-data_set[,10:27]
# write.csv2(r_inf,file="Gaelle_Fig2Drisk_datapoints.csv",quote=F,row.names=T)
 



palette<-c("#999999","#000000")
palette2<-c("black","white")
data_set$category<-factor(data_set$response, levels=c("R","NR"))

wilc_rstatix<-data_set %>%rstatix::wilcox_test(dataCD68n_Siggrp_XV ~ category)
df_pvals <- wilc_rstatix %>% add_xy_position() %>%  add_significance("p") %>%mutate(category="R") 
df_pvals$y.position<-df_pvals$y.position+0.05 
df_pvals$p_round<-round(df_pvals$p,4)
  
  
prismplot3<-function(long_df2,y_title,y_val,p_v_lab,min_y,max_y, psiz=25){
  ggplot(long_df2, aes(x =category, y =!!sym(y_val), fill=category))+
    stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
    geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
    scale_x_discrete(limits = levels(long_df2$category) )+  # , labels=c("U","I","R")
    scale_color_manual(values = palette)+ scale_y_continuous(n.breaks = 5) +#breaks = c(0.25,0.5,0.75,1)) +
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

pdf("./Figures_print/Fig2Db_XV_Bxplt.pdf",width =4,height = 8)
prismplot3(long_df2=data_set,y_title="XV_IFIM prog",y_val="dataCD68n_Siggrp_XV",
           p_v_lab="p.signif",min_y=min(data_set$dataCD68n_Siggrp_XV),max_y=max(data_set$dataCD68n_Siggrp_XV)+0.05)
prismplot3(long_df2=data_set,y_title="XV_IFIM prog",y_val="dataCD68n_Siggrp_XV",
           p_v_lab="p_round",psiz=10,min_y=min(data_set$dataCD68n_Siggrp_XV),max_y=max(data_set$dataCD68n_Siggrp_XV)+0.05)
dev.off()

