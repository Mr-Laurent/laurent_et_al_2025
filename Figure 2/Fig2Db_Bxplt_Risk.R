library(Rfast)
library(ggplot2)
library(ggprism)
library(ggpubr)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
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
  data_set=r_dsgn_r
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

pdf(paste0("./Figure 2/Fig2Db_XV_Bxplt.pdf"),width =6,height = 7.5)
all_sig_comp_BR("_Inf")
dev.off()


