library(Matrix)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggplot2)
library(progeny)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
load("Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v3.rd")
load("Grouped_objects/Macs_dgCmc_16nov22_counts.rd") 
df<-read.csv("Grouped_objects/EnrMacs_Mono1merged_annot_june24.csv") # Load metacells annotation




#---------------------------------------#


hm_zmod<-as.matrix(metadf_z[,c(154:303)])
df_zmod<-as.data.frame(metadf_z[,c(2:3,154:303)])
colnames(df_zmod)<-sub("([a-zA-Z_]+)(\\d+)([a-zA-Z]+)", "Module_\\2", colnames(df_zmod))  
preloco_nc<-as.matrix(log(1+dgcmtxcounts))
allsum<-Matrix::colSums(preloco_nc) 



###----------------------------------------------------------### 
####    Function to have multiple pvalue levels as stars     ### 
get_significance_label <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )}
###----------------------------------------------------------### 
#### Function to get all the genes as vector for each list  #### 
findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}  
###----------------------------------------------------------### 

palette_1=c('Mono-like 1'='#874037','Mono-like 2'='#C3A34B','Mono-like 3'='#f77774','Mono-like 4'='#ffa3a1','Mono-like 5'='#77c799',
            'Macrophage-like 6'='#c967eb','Macrophage-like 7'='#5fc2ed','Macrophage-like 8'='#4F88B9','Macrophage-like 9'='#5C538B')
bx_order=c('Mono-like 1','Mono-like 2','Mono-like 3','Mono-like 4','Mono-like 5','Macrophage-like 6','Macrophage-like 7','Macrophage-like 8','Macrophage-like 9')

# "Symbol stats - comparison vs mean"
# Figres were first done in the app, but to have a similar look, I'm modifying the code here:
bx_hei=600
bx_wid=550
bx_siz=2
vjust_stat=0.8
y_txt_size=3.75
bx_stat_siz=1.6
clr_ln_score="#444444"

ann_df<-data.frame(metacell=rownames(metadf_z),annotation="undefined")

# Get the genes enriched in the selected pathways from the PROGENy data:
model <- progeny::model_human_full

model_500 <- model %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 500)

model_500$pathway<-gsub("-","_",model_500$pathway)
pthw<-names(table(model_500$pathway))

pthw<-c("NFkB","JAK_STAT","MAPK","Hypoxia")


  
df_zmod2<-df_zmod
df_zmod2$Annotation<-df$annotation[match(rownames(df_zmod),df$metacell )]
df_zmod2$currentscore<-NA
rn_slct_dfz2<-rownames(df_zmod2[which(df_zmod2$Annotation%in%bx_order),])



#Get min and max score values across all score to homogenize the plots:
minval=0
maxval=0
for(mods in pthw){
  ilist<-intersect(model_500[which(model_500$pathway==mods&model_500$weight>0),]$gene,rownames(ds))
  values<- scale( Matrix::colSums(preloco_nc[ilist,rn_slct_dfz2,drop=F])/allsum[rn_slct_dfz2] )[,]
  if(min(values)<minval){minval<-min(values)}
  if(max(values)>minval){maxval<-max(values)}
}  



pdf(file="./Figure 4/Fig4D_S9A_Bxplt_prgm_progeny.pdf",width = 8, height = 5.5)
for(mods in pthw){
  ilist<-intersect(model_500[which(model_500$pathway==mods&model_500$weight>0),]$gene,rownames(ds))
  values<- scale( Matrix::colSums(preloco_nc[ilist,rn_slct_dfz2,drop=F])/allsum[rn_slct_dfz2] )[,]  # Need [,] to drop the attributes of scale(x) so it works for after
  df_zmod3<-df_zmod2[rn_slct_dfz2,]  
  eval(parse(text=paste0("df_zmod3$currentscore[match(rownames(df_zmod3),names(values))]<-values")))
  brplt_title<-paste0("Program ",gsub("Sig_","",mods))
  stat_results <- compare_means(currentscore ~ Annotation, data=df_zmod3, method = "wilcox.test", p.adjust.method = "bonferroni")
  
  # Stats : check K-W result, then do wilcoxon + correction
  KW_rstatix<-df_zmod3 %>%rstatix::kruskal_test(currentscore ~ Annotation)
  pwc_label <- bquote(paste("pwc: ", bold("Wilcoxon test"),"; p.adjust: ", bold("Bonferroni")) )  #pwc= PairWise Comparison
  test_label <- get_test_label(KW_rstatix, detailed = TRUE)
  combined_label <- bquote(paste(.(test_label), " | ", .(pwc_label)))
  
  results_wilcoxon<- df_zmod3 %>% group_by(Annotation) %>%
    summarise(p_value = wilcox.test(currentscore, mu = mean(df_zmod3$currentscore))$p.value, .groups = "drop") %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
  
  # preplot containing stats
  plot_w4<-ggplot(data = df_zmod3, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation ))+
    ggtitle(paste0(brplt_title))+
    geom_text(data = results_wilcoxon, aes(x = Annotation, y = maxval + 0.2*bx_stat_siz,
                                           label = get_significance_label(p_adj)),
              size = 5*bx_stat_siz, color="black", vjust = 0)+scale_x_discrete(limits=bx_order)
  
  
  palette_en_cours2<-palette_1[bx_order]
  
  # Add a space between mono and macs with an empty column
  bx_order_space=c('Mono-like 1','Mono-like 2','Mono-like 3','Mono-like 4','Mono-like 5','nothing','Macrophage-like 6','Macrophage-like 7','Macrophage-like 8','Macrophage-like 9')
  
  # print the final plot with boxplots+stats
  print( plot_w4+geom_hline(yintercept = mean(df_zmod3$currentscore), linetype = "dashed",color=clr_ln_score, size=rel(0.7)*bx_siz) +
           geom_boxplot(width=0.8, size=0.7*bx_siz,outlier.alpha = 0.5)+ #fill="white",
           xlab("Cluster")+theme_bw()+
           scale_fill_manual(values = alpha(palette_en_cours2, .4)) +
           scale_color_manual(values = palette_en_cours2) +coord_cartesian(ylim = c(minval, (maxval+ 0.2*bx_stat_siz)))+
           theme(axis.text.x = element_blank(),
                 # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = rel(1.8),margin = margin(t = 6)),
                 axis.title.x = element_blank(),
                 axis.text.y = element_text(size = rel(y_txt_size),margin = margin(r = 4)),
                 axis.title.y=element_blank(),
                 plot.title = element_text(size = rel(1.5), face = "bold"),
                 plot.subtitle = element_text(size = rel(1)),
                 panel.border = element_blank(),
                 legend.key.size = unit(1*bx_siz,'cm'),
                 legend.position = "none",
                 axis.ticks = element_line(size=rel(2.5)),axis.ticks.length = unit(0.2, "cm"),
                 axis.ticks.x = element_blank(),
                 legend.text = element_text(size=rel(1)*bx_siz),
                 axis.line = element_line(colour = "black",size = rel(1.25)*bx_siz)) +
           geom_jitter(height = 0,size=0.5*bx_siz, width = 0.2,alpha=0.5) +
           geom_vline(xintercept=6, linetype = "21",color='#CCC', size=rel(0.5)*bx_siz)+
           
           scale_x_discrete(expand=c(0.1,0), limits=bx_order_space,breaks=bx_order, drop=FALSE) # Add a space between mono and macs with an empty column
  )
  
}

dev.off()



