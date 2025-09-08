library(Rfast)
library(ggplot2)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
# Load the metadata and count matrix from Ngollo cohort
load("Grouped_objects/minipheno_AllezNgollo_postop_goodRemRec_sig.rd")
load("Grouped_objects/matex_AllezNgollo_postop.rd")   
load("Grouped_objects/genelist_16nov22_Macs.rd")    
minipheno$title<-gsub("^(\\w)(\\w)(_.*)$", "\\1_\\2\\3", minipheno$title)
setdiff(minipheno$title,colnames(matex))
allsum<-Rfast::colsums(matex) # need same order as the score computation


val1<-"data_Siggrp_VIII"
val2<-"data_Siggrp_XV"
lgd1=paste0("Signature ",gsub("data_Siggrp_","",val1))
lgd2=paste0("Signature ",gsub("data_Siggrp_","",val2))
pal=rev(rainbow(15))[4:15]

# Subset categories
miniphenoM0I<-minipheno[which(minipheno$`location:ch1`%in%c("M0I")),]
miniphenoCTRLM0IM0M<-minipheno[which(minipheno$`location:ch1`%in%c("M0I","Ctrl","M0M")),]
miniphenoCTRL<-minipheno[which(minipheno$`location:ch1`%in%c("Ctrl")),]


todolist<-c("CTRL","M0I") 
pdf(paste0("./Figure 2/Fig2Ca_Scatter.pdf"),width =6,height = 6)
for (sub_set in todolist) {
  eval(parse(text=paste0("data_set=minipheno",sub_set)))
  plot_title=sub_set
  vec_val1<-eval(parse(text=paste0("miniphenoCTRLM0IM0M$",val1)))
  vec_val2<-eval(parse(text=paste0("miniphenoCTRLM0IM0M$",val2)))
  df2<-data_set[c(val1,val2)]
  colnames(df2)<-c("x","y")
  spr_tst<-cor.test(df2$x,df2$y,method = "spearman") 
  if(spr_tst$p.value==0){pval<-"< 2.2e-16"}else if(spr_tst$p.value<=0.001){
    pval<-format(spr_tst$p.value, scientific = T,digits = 3) }else{pval<-round(spr_tst$p.value,5)}
  print(ggplot(data_set, aes(x=!!sym(val1), y=!!sym(val2))) +geom_point(size=rel(3.5), stroke=0)+
          theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
          theme(plot.title= element_text(size=rel(1.8)),plot.subtitle= element_text(size=rel(1.8),face="bold"),
                axis.title.x = element_text(size=rel(1.8), hjust = 0.5,),
                axis.title.y = element_text(size=rel(1.8), vjust=2,hjust = 0.5),
                axis.text.x = element_text(size=rel(2)),
                axis.text.y = element_text(size=rel(2)),
                legend.position="bottom", legend.text = element_text(size=rel(1.5)),
                legend.title = element_text(size=rel(1.5)))+
          labs(y= lgd2, x = lgd1)+xlim(min(vec_val1),max(vec_val1))+ylim(min(vec_val2),max(vec_val2))+
          ggtitle(plot_title,label =bquote("Spearman: "~rho~"= "~.(round(spr_tst$estimate,3))~", "~italic(p[val])~"= "~.(pval))) )
}
  dev.off()


# Get the stats :
  
sub_set="CTRL"
eval(parse(text=paste0("data_set=minipheno",sub_set)))
plot_title=sub_set
vec_val1<-eval(parse(text=paste0("miniphenoCTRLM0IM0M$",val1)))
vec_val2<-eval(parse(text=paste0("miniphenoCTRLM0IM0M$",val2)))
df2<-data_set[c(val1,val2)]
colnames(df2)<-c("x","y")
spr_tst_u<-cor.test(df2$x,df2$y,method = "spearman") 

sub_set="M0I"
eval(parse(text=paste0("data_set=minipheno",sub_set)))
plot_title=sub_set
vec_val1<-eval(parse(text=paste0("miniphenoCTRLM0IM0M$",val1)))
vec_val2<-eval(parse(text=paste0("miniphenoCTRLM0IM0M$",val2)))
df2<-data_set[c(val1,val2)]
colnames(df2)<-c("x","y")
spr_tst_i<-cor.test(df2$x,df2$y,method = "spearman") 

  
format(p.adjust(c(format(spr_tst_u$p.value, scientific = T,digits = 3),
                    format(spr_tst_i$p.value, scientific = T,digits = 3)), 
                  method = "bonferroni"), scientific = T,digits = 3)  
  
# "8.22e-03" "0.00e+00"
  
