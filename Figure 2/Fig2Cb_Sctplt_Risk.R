library(Rfast)
library(ggplot2)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
# Load the metadata and count matrix from RISK cohort
load("Grouped_objects/r_dsgn_InfUninf_cd68norm_Momacsscores.rd")
r_dsgn_Inf<-r_dsgn[which(r_dsgn$status%in%c("Inflamed")),]
r_dsgn_Uninf<-r_dsgn[which(r_dsgn$status%in%c("Uninflamed")),]


##---------------------------------------##
##--##  Scatter with All spl, black  ##--##
##---------------------------------------##
sctr_black_all<-function(data_set,subset,val1,val2,lgd1,lgd2,plot_title){
  val_1<-paste0("dataCD68n_",val1)
  val_2<-paste0("dataCD68n_",val2)
  vec_val1<-eval(parse(text=paste0("r_dsgn$",val_1)))
  vec_val2<-eval(parse(text=paste0("r_dsgn$",val_2)))
  df2<-data_set[c(val_1,val_2)]
  colnames(df2)<-c("x","y")
  # Spearman as it's non-parametric
  spr_tst<-cor.test(df2$x,df2$y,method = "spearman") 
  if(spr_tst$p.value==0){pval<-"< 2.2e-16"}else if(spr_tst$p.value<=0.001){
    pval<-format(spr_tst$p.value, scientific = T,digits = 3) }else{pval<-round(spr_tst$p.value,5)}
  # Clinical outcome:
  return(ggplot(data_set, aes(x=!!sym(val_1), y=!!sym(val_2))) +geom_point(size=rel(3.5), stroke=0)+
           theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
           theme(plot.title= element_text(size=rel(1.8)),plot.subtitle= element_text(size=rel(1.8),face="bold"),
                 axis.title.x = element_text(size=rel(1.8), hjust = 0.5,),
                 axis.title.y = element_text(size=rel(1.8), vjust=2,hjust = 0.5),
                 axis.text.x = element_text(size=rel(2)),
                 axis.text.y = element_text(size=rel(2)),
                 legend.position="bottom", legend.text = element_text(size=rel(1.5)),
                 legend.title = element_text(size=rel(1.5)))+
           labs(y= lgd2, x = lgd1)+xlim(min(vec_val1),max(vec_val1))+ylim(min(vec_val2),max(vec_val2))+
           ggtitle(plot_title,label =bquote("Spearman: "~rho~"= "~.(round(spr_tst$estimate,3))~", "~italic(p[val])~"= "~.(pval)))
  )
}  

##---------------------------------------##

#### Figure in black, all RISK samples #### 

signatures<-c("Siggrp_VIII","Siggrp_XV")
all_sig_comp_SBA<-function(sub_set){
  eval(parse(text=paste0("data_set=r_dsgn",sub_set)))
  n <- length(signatures)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      x <- signatures[i]
      y <- signatures[j]
      lgd1=paste0("Signature ",gsub("Siggrp_","",x))
      lgd2=paste0("Signature ",gsub("Siggrp_","",y))
      print(paste0(lgd1," ",lgd2))
      
      print(sctr_black_all(data_set =data_set, subset=sub_set,val1 = x,val2=y,lgd1=lgd1,lgd2=lgd2,plot_title=paste0("CD68n",sub_set)))
    }} 
}

pdf(paste0("./Figure 2/Fig2Cb_Sct_Risk.pdf"),width =6,height = 6)
todolist<-c("_Uninf","_Inf")
for (sub_set in todolist) {
  all_sig_comp_SBA(sub_set)
}
dev.off()






