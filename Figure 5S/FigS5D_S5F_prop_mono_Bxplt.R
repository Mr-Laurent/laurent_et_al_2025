# S5D: "Proportion of Monocytes among MoMacs"
library(dplyr)
library(purrr)
library(tidyverse)
library(ggpubr)


df = read.csv("./Grouped_objects/momacs_cell_metadata.csv")
annots = read.csv("./Grouped_objects/Buckley_momacs_annot_241217.csv")


df = df[df$metacell_name != "Outliers",]
df$annot = unlist(lapply(X = df$metacell_name, FUN = function(x) {return(annots$annotation[annots$metacell == x])} ))
df = df[df$annot != "undefined" ,]   # Remove the "undefined" annotation

# df objects by Disease condition
df_UC = df[df$Disease == "UC",]
df_CD = df[df$Disease == "CD",]
df_Hy = df[df$Disease == "Healthy",]
df_list_Dis = list(df_UC = df_UC, df_CD = df_CD, df_Hy = df_Hy)

# Subdivide them by Remission status
df_list_DisRem = list()
for (i in 1:length(df_list_Dis)) {
  my_df = df_list_Dis[i]
  df_i = data.frame(my_df[[1]])
  
  df_i_Resp  = df_i[df_i$Remission == "Responder",]
  df_i_NResp = df_i[df_i$Remission == "Non_responder",]
  df_i_Empty = df_i[df_i$Remission == "",]
  df_i_NotAv = df_i[df_i$Remission == "Not_avail",]
  
  df_list_i = list(df_i_Resp = df_i_Resp, df_i_NResp = df_i_NResp, df_i_Empty = df_i_Empty, df_i_NotAv = df_i_NotAv)
  names(df_list_i) = c(paste0(names(my_df),"_Resp"),
                       paste0(names(my_df),"_NResp"),
                       paste0(names(my_df),"_Empty"),
                       paste0(names(my_df),"_NotAv"))
  df_list_DisRem = append(df_list_DisRem, df_list_i)
}

# subdivide them by Treatment condition
df_list_DisRemTtt = list()
for (i in 1:length(df_list_DisRem)) {
  my_df = df_list_DisRem[i]
  df_i = data.frame(my_df[[1]])
  
  df_i_Pre  = df_i[df_i$Treatment == "Pre",]
  df_i_Post = df_i[df_i$Treatment == "Post",]
  df_i_None = df_i[df_i$Treatment == "None",]
  
  df_list_i = list(df_i_Pre = df_i_Pre, df_i_Post = df_i_Post, df_i_None = df_i_None)
  names(df_list_i) = c(paste0(names(my_df),"_Pre"),
                       paste0(names(my_df),"_Post"),
                       paste0(names(my_df),"_None"))
  df_list_DisRemTtt = append(df_list_DisRemTtt, df_list_i)
}




### Now only pre-ttt samples, separated according to condition (H, NInfl, Infl) and site (Ileum, Ascending_Colon, Rectum...) ----
df_Healthy = df[df$Disease == "Healthy",]
df_CD = df[df$Disease == "CD",]
df_PreCD = df_CD[df_CD$Treatment == "Pre",]

df_PreCD_NInfl = df_PreCD[df_PreCD$Inflammation == "Non-inflamed",]
df_PreCD_Infl  = df_PreCD[df_PreCD$Inflammation == "Inflamed",]

df_Healthy_Ileum = df_Healthy[df_Healthy$Site == "Terminal_Ileum",]
df_Healthy_ColRc = df_Healthy[df_Healthy$Site != "Terminal_Ileum",]
df_PreCD_NInfl_Ileum = df_PreCD_NInfl[df_PreCD_NInfl$Site == "Terminal_Ileum",]
df_PreCD_NInfl_ColRc = df_PreCD_NInfl[df_PreCD_NInfl$Site != "Terminal_Ileum",]
df_PreCD_Infl_Ileum = df_PreCD_Infl[df_PreCD_Infl$Site == "Terminal_Ileum",]
df_PreCD_Infl_ColRc = df_PreCD_Infl[df_PreCD_Infl$Site != "Terminal_Ileum",]

df_list_HPreCD = list(df_Healthy_Ileum = df_Healthy_Ileum, df_Healthy_ColRc = df_Healthy_ColRc,
                      df_PreCD_NInfl_Ileum = df_PreCD_NInfl_Ileum, df_PreCD_NInfl_ColRc = df_PreCD_NInfl_ColRc,
                      df_PreCD_Infl_Ileum = df_PreCD_Infl_Ileum, df_PreCD_Infl_ColRc = df_PreCD_Infl_ColRc)

# Create a summary (with cell counts for each cell annotation within each combination of conditions)
summ_list_HPreCD = list()
for (i in 1:length(df_list_HPreCD)) {   # For each of the condition dfs :
  my_df = df_list_HPreCD[i]
  df_i = data.frame(my_df[[1]])
  df_i = group_by(df_i, sample_id, annot)
  
  annot_vect = c("Mono1", "Mono2", "Mono3", "Mono4", "Mono5", "Macro6", "Macro7", "Macro8", "Macro9")
  
  summ_df = ungroup(summarise(df_i, cell_count = n()))
  summ_df$annot = factor(summ_df$annot, levels = annot_vect)
  summ_df = list(summ_df)
  names(summ_df) = gsub("df_", replacement = "", names(my_df))
  summ_list_HPreCD = append(summ_list_HPreCD, summ_df)
}

#### Code to add Patient and HBI tiles   ----
HBI_mapping = c(CD1 = 13, CD2 = 14, CD3 = 9, CD4 = 6, CD5 = 11, CD6 = 12, CD7 = 6, CD8 = 8, CD9 = 10, CD10 = 8,
                CD11 = 10, CD12 = 8, CD13 = 14, CD14 = 15, CD15 = 6, CD16 = 9, CD17 = NA)
for (i in 1:length(summ_list_HPreCD)) {
  my_df = summ_list_HPreCD[[i]]
  my_df$Patient = sapply(X = my_df$sample_id, FUN = function(x) {unique(df$Patient[df$sample_id == x])})
  if (!is_empty(grep(pattern = "H", x = my_df$Patient))) {
    my_df$HBI = 0
  } else {
    my_df = my_df %>% mutate(HBI = HBI_mapping[Patient])
  }
  summ_list_HPreCD[[i]] = my_df
}


### Make df of Monocyte proportion by patient ----
mono_mac_prop_df = data.frame()
for (i in 1:length(summ_list_HPreCD)) {
  my_df = summ_list_HPreCD[[i]]
  
  if (nrow(my_df) != 0) {
    wider_df = my_df %>% pivot_wider(names_from = annot, values_from = cell_count) %>% data.frame()
    wider_df = wider_df[, c(-2, -3)]
    wider_df[is.na(wider_df)] = 0
    wider_df = wider_df[rowSums(wider_df[,-1]) > 15,]   # Remove all samples with less than 15 momacs
    
    Mono_cols  = grep(pattern = "Mono",  x = colnames(wider_df))
    Macro_cols = grep(pattern = "Macro", x = colnames(wider_df))
    wider_df$Mono = rowSums(wider_df[,Mono_cols])
    wider_df$Macro = rowSums(wider_df[,Macro_cols])
    wider_df$Prop_mono = wider_df$Mono / (wider_df$Mono + wider_df$Macro)
    
    my_df = wider_df[, colnames(wider_df) %in% c("sample_id", "Prop_mono")]
    colnames(my_df)[1] = "condition"
    my_df$condition = names(summ_list_HPreCD)[i]
    
    mono_mac_prop_df = rbind(mono_mac_prop_df, my_df)
  }
}

mono_mac_prop_df$condition = factor(mono_mac_prop_df$condition, levels = c("Healthy_Ileum", "Healthy_ColRc", "PreCD_NInfl_Ileum", "PreCD_NInfl_ColRc", "PreCD_Infl_Ileum", "PreCD_Infl_ColRc"))

# PLOT only with ileal samples
select_rows_bool = grepl(x = mono_mac_prop_df$condition, pattern = "Ileum")
mono_mac_prop_df = mono_mac_prop_df[select_rows_bool ,]
mono_mac_prop_df$condition = factor(mono_mac_prop_df$condition, levels = c("Healthy_Ileum", "PreCD_NInfl_Ileum", "PreCD_Infl_Ileum"))



# For IFIM proportionsin momacs:
df_Ileum = df[df$Site == "Terminal_Ileum",]
df_Healthy = df_Ileum[df_Ileum$Disease == "Healthy",]
df_CD = df_Ileum[df_Ileum$Disease == "CD",]
df_PreCD = df_CD[df_CD$Treatment == "Pre",]

df_PreCD_NInfl = df_PreCD[df_PreCD$Inflammation == "Non-inflamed",]
df_PreCD_Infl  = df_PreCD[df_PreCD$Inflammation == "Inflamed",]
df_PreCD_Infl_NResp = df_PreCD_Infl[df_PreCD_Infl$Remission == "Non_responder",]
df_PreCD_Infl_Resp  = df_PreCD_Infl[df_PreCD_Infl$Remission == "Responder",]

df_list_HPreCD = list(df_Healthy = df_Healthy, df_PreCD_NInfl = df_PreCD_NInfl,
                      df_PreCD_Infl_NResp = df_PreCD_Infl_NResp, df_PreCD_Infl_Resp = df_PreCD_Infl_Resp)

summ_list_HPreCD = list()
for (i in 1:length(df_list_HPreCD)) {   # For each of the condition dfs :
  my_df = df_list_HPreCD[i]
  df_i = data.frame(my_df[[1]])
  df_i = group_by(df_i, sample_id, annot)

  annot_vect = c("Mono1", "Mono2", "Mono3", "Mono4", "Mono5", "Macro6", "Macro7", "Macro8", "Macro9")

  summ_df = ungroup(summarise(df_i, cell_count = n()))
  summ_df$annot = factor(summ_df$annot, levels = annot_vect)
  summ_df = list(summ_df)
  names(summ_df) = gsub("df_", replacement = "", names(my_df))
  summ_list_HPreCD = append(summ_list_HPreCD, summ_df)
}

### Now separate monos and macs
mono_summ_list = list()
macs_summ_list = list()
for (i in 1:length(summ_list_HPreCD)) {
  my_df = summ_list_HPreCD[[i]]
  mono_summ_list[[i]] = my_df[my_df$annot %in% c("Mono1", "Mono2", "Mono3", "Mono4", "Mono5"),]
  names(mono_summ_list)[i] = names(summ_list_HPreCD)[i]
  macs_summ_list[[i]] = my_df[my_df$annot %in% c("Macro6", "Macro7", "Macro8", "Macro9"),]
  names(macs_summ_list)[i] = names(summ_list_HPreCD)[i]
}



long_df = data.frame()
for (i in 1:length(summ_list_HPreCD)) {
  my_df = summ_list_HPreCD[[i]]

  if (nrow(my_df) != 0) {
    wider_df = my_df %>% pivot_wider(names_from = annot, values_from = cell_count) %>% data.frame()
    wider_df[is.na(wider_df)] = 0
    wider_df = wider_df[rowSums(apply(wider_df[,-c(1,2,3)],2,as.numeric) ) > 15,]   # Remove all samples with less than 15 momacs

    if (!("Mono1" %in% colnames(wider_df))) {
      wider_df$Mono1_in_Monos  = 0
      wider_df$Mono1_in_Momacs = 0
    } else {
      Mono_cols  = grep(pattern = "Mono",  x = colnames(wider_df))
      wider_df$Mono1_in_Monos  = wider_df$Mono1 / rowSums(wider_df[,Mono_cols])
      wider_df$Mono1_in_Momacs = wider_df$Mono1 / rowSums(wider_df[,-c(1,2,3, length(colnames(wider_df)))])
    }
    wider_df = wider_df[, colnames(wider_df) %in% c("sample_id", "Mono1_in_Monos", "Mono1_in_Momacs")]
    my_df = wider_df %>% pivot_longer(cols = colnames(wider_df)[-1], names_to = "prop_type", values_to = "proportion")

    colnames(my_df)[1] = "condition"
    my_df$condition = names(summ_list_HPreCD)[i]

    long_df = rbind(long_df, my_df)
  }
}



long_df$condition = factor(long_df$condition, levels = c("Healthy", "PreCD_NInfl", "PreCD_Infl_NResp", "PreCD_Infl_Resp"))
in_momacs = long_df[long_df$prop_type == "Mono1_in_Momacs",]


data_set<-mono_mac_prop_df
palette<-c("#EEEEEE","#AAAAAA","#000000")
palette2<-c("black","black","white")
data_set$category<-factor(data_set$condition, levels=c("Healthy_Ileum","PreCD_NInfl_Ileum","PreCD_Infl_Ileum"))


data_set%>%rstatix::kruskal_test(Prop_mono ~ category)
# No dunn post-hoc test to do, KW not signif


prismplot3_nostat<-function(long_df2,y_title,y_val,min_y,max_y,labels=NULL){
  ggplot(long_df2, aes(x =category, y =!!sym(y_val), fill=category))+
    stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
    geom_boxplot(lwd=1.15, color="black", fill=palette, fatten = NULL, outlier.colour = "white",width=0.6)+
    scale_x_discrete(limits = levels(long_df2$category) , labels= labels)+
    scale_color_manual(values = palette)+ scale_y_continuous(n.breaks = 5) +#breaks = c(0.25,0.5,0.75,1)) +
    ggprism::theme_prism(base_size = 16)+ labs(y= y_title)+
    stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = category), fatten = 4)+ # Visualize the median on top, allows to adapt color
    scale_color_manual(values = palette2 ) +
    geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
    # NO STAT # stat_pvalue_manual(df_pvals, label = p_v_lab, tip.length = 0,hide.ns = TRUE, bracket.size=rel(1.4), label.size = rel(25))+
    theme(legend.position = "none",axis.title.x =element_blank(),
          axis.text.y = element_text(size=rel(1.5),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(1)),
          axis.ticks.length = unit(15,"points"),
          plot.margin = unit(c(40, 4, 4, 4),"points"))+   #add margin above plot so asterisks are not cropped
    coord_cartesian(ylim = c(min_y, max_y),clip = 'off')#+  #cut the plot at 1
}

pdf("./Figures_print/FigS5D_Bxplt_Monoprop.pdf",width =4,height = 6)
prismplot3_nostat(long_df2=data_set,y_title="% of Mono in Momacs",y_val="Prop_mono",
                  min_y=0,max_y=max(data_set$Prop_mono))
dev.off()



data_set<-in_momacs
palette<-c("#EEEEEE","#AAAAAA","#444444","#000000")
palette2<-c("black","black","white","white")
data_set$category<-factor(data_set$condition, levels=c("Healthy","PreCD_NInfl","PreCD_Infl_Resp","PreCD_Infl_NResp"))


data_set%>%rstatix::kruskal_test(proportion ~ category)
# No dunn post-hoc test to do, KW not signif


pdf("./Figures_print/FigS5F_Bxplt_IFIMprop.pdf",width =4,height = 6)
prismplot3_nostat(long_df2=data_set,y_title="% of IFIM in Momacs",y_val="proportion",
                  min_y=0,max_y=max(data_set$proportion))
dev.off()


