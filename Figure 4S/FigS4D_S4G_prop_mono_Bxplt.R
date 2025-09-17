# S4D: "Proportion of Monocytes among MoMacs"
library(dplyr)
library(purrr)
library(tidyverse)
library(ggpubr)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
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

# Rstatix list of comparisons (Wilcoxon):
my_comparisons <- list( c("Healthy_Ileum", "PreCD_NInfl_Ileum"),
                        c("PreCD_NInfl_Ileum", "PreCD_Infl_Ileum"),
                        c("Healthy_Ileum", "PreCD_Infl_Ileum") )

color_vect = c("Healthy_Ileum" = "#ABCBD6", "PreCD_NInfl_Ileum" = "#C97185", "PreCD_Infl_Ileum" = "#C22D4F")
pdf("./Figure 4S/FigS4D_prop_mono_Bxplt.pdf", width = 7, height = 6)

my_comparisons <- list( c("Healthy_Ileum","PreCD_NInfl_Ileum"), c("PreCD_NInfl_Ileum","PreCD_Infl_Ileum"), c("Healthy_Ileum","PreCD_Infl_Ileum") )

ggplot(mono_mac_prop_df, aes(x =condition, y =Prop_mono, color =condition, fill=condition))+
  geom_boxplot(lwd=1.6, color=color_vect, fill=color_vect, alpha=0.5)+
  scale_x_discrete(limits = c("Healthy_Ileum","PreCD_NInfl_Ileum","PreCD_Infl_Ileum"),
                   labels = c("Healthy\nIleum","PreCD_NInfl\nIleum","PreCD_Infl\nIleum") )+
  geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
  scale_color_manual(values = color_vect)+
  ggprism::theme_prism(base_size = 16)+
  theme(legend.position = "none",axis.title.x =element_blank())+
  labs(y= "Monocyte frequency\nin Macrophage subset")+
  stat_compare_means(comparisons = my_comparisons,size = 6,bracket.size = rel(1.5), method = "wilcox.test", p.adjust.method = "bonferroni") +
  stat_compare_means(label.y = max(mono_mac_prop_df$Prop_mono)+0.35,size = 8) # size is for K-W pvalue title 

# It's not stat, so make the same without bars for the paper :
ggplot(mono_mac_prop_df, aes(x =condition, y =Prop_mono, color =condition, fill=condition))+
  geom_boxplot(lwd=1.6, color=color_vect, fill=color_vect, alpha=0.5)+
  scale_x_discrete(limits = c("Healthy_Ileum","PreCD_NInfl_Ileum","PreCD_Infl_Ileum"),
                   labels = c("Healthy\nIleum","Pre-ttt\nNon-infl.\nIleum","Pre-ttt\nInflamed\nIleum") )+
  geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
  scale_color_manual(values = color_vect)+
  ggprism::theme_prism(base_size = 16)+
  theme(legend.position = "none",axis.title.x =element_blank())+
  labs(y= "Proportion of monocytes\nin mo-macs")
dev.off()



##############################
## Fig 4G:



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

my_comparisons <- list( c("Healthy", "PreCD_NInfl"),
                        c("Healthy", "PreCD_Infl_NResp"),
                        c("Healthy", "PreCD_Infl_Resp"),
                        c("PreCD_NInfl", "PreCD_Infl_NResp"),
                        c("PreCD_NInfl", "PreCD_Infl_Resp"),
                        c("PreCD_Infl_NResp", "PreCD_Infl_Resp") )

pdf("./Figure 4S/FigS4G_prop_mono_Bxplt.pdf", width = 7, height = 6)

ggplot(in_momacs, aes(x = condition, y = proportion, fill = condition, color = condition)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, size = rel(0.5)) + # I remove outlier dots as I have jitter, it will avoid misunderstanding
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", fatten = 2) + # Add black median line
  theme_minimal() +
  theme(strip.text = element_blank(), axis.text.y = element_text(size=rel(0.7)), axis.title.y = element_blank()  ) +    #,panel.border = element_rect(colour = "black", fill=NA, size=1) )+  # To remove subplot titles, needs to be put after theme_minimal
  labs(title = "Proportion of Mono1 among all MoMacs in Ileum samples",    # In TL's plot : "(only patients with >15 MoMacs)"
       x = "Category",
       y = "Mono / Momacs proportion") +
  scale_fill_manual(values = color_vect) +
  scale_x_discrete(limits = levels(in_momacs$condition)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = color_vect) +
  geom_jitter(height = 0, size = rel(1), width = 0.25, alpha = 0.5, aes(color = prop_type)) +
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 6, angle= 90)) +
  stat_compare_means(comparisons = my_comparisons, size = rel(3), bracket.size = rel(0.5)) + # Add pairwise comparisons p-value, size is for p-values
  stat_compare_means(label.y = 0.75, size = rel(4)) # size is for K-W pvalue title 

# K-W not stat: show only the boxplot without stats
color_vect = c("Healthy" = "#ABCBD6", "PreCD_NInfl" = "#C97185", "PreCD_Infl_NResp" = "#C22D4F", "PreCD_Infl_Resp" = "#E27D4F" )
ggplot(in_momacs, aes(x =condition, y =proportion, color =condition, fill=condition))+
  geom_boxplot(lwd=1.6, color=color_vect, fill=color_vect, alpha=0.5)+
  scale_x_discrete(limits = c("Healthy", "PreCD_NInfl", "PreCD_Infl_NResp", "PreCD_Infl_Resp"),
                   labels = c("Healthy\nIleum", "Pre-ttt\nNon-infl.\nIleum", "Pre-ttt\nInfl. Ileum\nNon-resp.", "Pre-ttt\nInfl. Ileum\nResp.") )+
  geom_jitter(shape=16, position=position_jitter(0.2),size=2.5)+
  scale_color_manual(values = color_vect)+
  ggprism::theme_prism(base_size = 16)+
  theme(legend.position = "none",axis.title.x =element_blank())+
  labs(y= "Proportion of Mono1\nin mo-macs")

dev.off()
