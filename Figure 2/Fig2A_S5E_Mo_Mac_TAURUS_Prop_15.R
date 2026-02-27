library(ComplexHeatmap)
library(parallelDist)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(tidyverse)
library(ggpubr)
library(rstatix)


df = read.csv("./Grouped_objects/momacs_cell_metadata.csv")
annots = read.csv("./Grouped_objects/Buckley_momacs_annot_241217.csv")


df = df[df$metacell_name != "Outliers",]
df$annot = unlist(lapply(X = df$metacell_name, FUN = function(x) {return(annots$annotation[annots$metacell == x])} ))
df = df[df$annot != "undefined" ,]   # Remove the "undefined" annotation



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

# Renaming samples
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



mono_summ_list = list()
macs_summ_list = list()
for (i in 1:length(summ_list_HPreCD)) {
  my_df = summ_list_HPreCD[[i]]
  mono_summ_list[[i]] = my_df[my_df$annot %in% c("Mono1", "Mono2", "Mono3", "Mono4", "Mono5"),]
  names(mono_summ_list)[i] = names(summ_list_HPreCD)[i]
  macs_summ_list[[i]] = my_df[my_df$annot %in% c("Macro6", "Macro7", "Macro8", "Macro9"),]
  names(macs_summ_list)[i] = names(summ_list_HPreCD)[i]
}


# Reordering
mono_summ_list[[1]] = mono_summ_list[[1]][order(factor(mono_summ_list[[1]]$sample_id, levels = c("CID005775-1", "CID005779-1", "CID005709-1"))),]
mono_summ_list[[2]] = mono_summ_list[[2]][order(factor(mono_summ_list[[2]]$sample_id, levels = c("CID005711-1", "CID005778-1", "CID005777-1", "CID005780-1",
                                                                                                 "CID005710-1", "CID005781-1", "CID005712-1", "CID005782-1", "CID005776-1"))),]
mono_summ_list[[3]] = mono_summ_list[[3]][order(factor(mono_summ_list[[3]]$sample_id, levels = c("CID004751-1", "CID004708-1"))),]
mono_summ_list[[4]] = mono_summ_list[[4]][order(factor(mono_summ_list[[4]]$sample_id, levels = c("CID003378-1", "CID003673-1", "CID004720-2", "CID004732-1",
                                                                                                 "CID004733-1", "CID004734-1", "CID003384-1", "CID004719-1", "CID005760-1", "CID004718-1", "CID003380-1", "CID004741-1", "CID003383-1",
                                                                                                 "CID003385-1", "CID004754-1", "CID004740-1", "CID004727-1", "CID005750-1", "CID003678-1", "CID003379-1", "CID004711-1", "CID006561-1", "CID005749-1"))),]
mono_summ_list[[5]] = mono_summ_list[[5]][order(factor(mono_summ_list[[5]]$sample_id, levels = c("CID006558-1", "CID003677-1", "CID006575-1", "CID003381-1",
                                                                                                 "CID005759-1", "CID005747-1", "CID005763-1", "CID005771-1", "CID004731-1", "CID003386-1", "CID005755-1", "CID003671-1"))),]
mono_summ_list[[6]] = mono_summ_list[[6]][order(factor(mono_summ_list[[6]]$sample_id, levels = c("CID003672-1", "CID004739-1", "CID005756-1", "CID003679-1",
                                                                                                 "CID005765-1", "CID004709-1", "CID005748-1", "CID004753-1", "CID004725-1", "CID004752-1", "CID005772-1", "CID005766-1", "CID004710-1", "CID004726-1", "CID005764-1"))),]

macs_summ_list[[1]] = macs_summ_list[[1]][order(factor(macs_summ_list[[1]]$sample_id, levels = c("CID005775-1", "CID005779-1", "CID005709-1"))),]
macs_summ_list[[2]] = macs_summ_list[[2]][order(factor(macs_summ_list[[2]]$sample_id, levels = c("CID005711-1", "CID005778-1", "CID005777-1", "CID005780-1",
                                                                                                 "CID005710-1", "CID005781-1", "CID005712-1", "CID005782-1", "CID005776-1"))),]
macs_summ_list[[3]] = macs_summ_list[[3]][order(factor(macs_summ_list[[3]]$sample_id, levels = c("CID004751-1", "CID004708-1"))),]
macs_summ_list[[4]] = macs_summ_list[[4]][order(factor(macs_summ_list[[4]]$sample_id, levels = c("CID003378-1", "CID003673-1", "CID004720-2", "CID004732-1",
                                                                                                 "CID004733-1", "CID004734-1", "CID003384-1", "CID004719-1", "CID005760-1", "CID004718-1", "CID003380-1", "CID004741-1", "CID003383-1",
                                                                                                 "CID003385-1", "CID004754-1", "CID004740-1", "CID004727-1", "CID005750-1", "CID003678-1", "CID003379-1", "CID004711-1", "CID006561-1", "CID005749-1"))),]
macs_summ_list[[5]] = macs_summ_list[[5]][order(factor(macs_summ_list[[5]]$sample_id, levels = c("CID006558-1", "CID003677-1", "CID006575-1", "CID003381-1",
                                                                                                 "CID005759-1", "CID005747-1", "CID005763-1", "CID005771-1", "CID004731-1", "CID003386-1", "CID005755-1", "CID003671-1"))),]
macs_summ_list[[6]] = macs_summ_list[[6]][order(factor(macs_summ_list[[6]]$sample_id, levels = c("CID003672-1", "CID004739-1", "CID005756-1", "CID003679-1",
                                                                                                 "CID005765-1", "CID004709-1", "CID005748-1", "CID004753-1", "CID004725-1", "CID004752-1", "CID005772-1", "CID005766-1", "CID004710-1", "CID004726-1", "CID005764-1"))),]


# color_vect = c("Mono1" = "#874037", "Mono2" = "#C3A34B", "Mono3" = "#F77774", "Mono4" = "#FFA3A1", "Mono5" = "#B4DEC6",
#                "Macro6" = "#C967EB", "Macro7" = "#A4DCEB", "Macro8" = "#4F88B9", "Macro9" = "#5C538B")
color_vect = c("Mono1" = "#874037", "Mono2" = "#C3A34B", "Mono3" = "#F77774", "Mono4" = "#FFA3A1", "Mono5" = "#B4DEC6")

name_i="PreCD_Infl_Ileum"

plot_df = mono_summ_list[[name_i]]


pdf(file="./Figures_print/Fig2A_TAURUS_mono_prop.pdf",width=6,height = 6)
print(ggplot(plot_df,aes(x= sample_id, y=cell_count,fill= annot)) +  geom_col(position = 'fill')+
        scale_fill_manual(values=color_vect)+  
        theme_classic()+labs( x = ' ', y = ' ',fill = 'annot')+ 
        scale_x_discrete(limits =unique(plot_df$sample_id),labels=as.character(c(6:17)))+ 
        theme(axis.text.y = element_text(size=rel(1.8)),axis.text.x = element_text(size=rel(1.5),angle = 90, hjust = 1, vjust = 0.5))
)
dev.off()


# Fig S5E
name_i="Healthy_Ileum"
plot_df = mono_summ_list[[name_i]]
pdf(file="./Figures_print/FigS5E_TAURUS_mono_prop.pdf",width=6,height = 6)
plt1<-ggplot(plot_df,aes(x= sample_id, y=cell_count,fill= annot)) +  geom_col(position = 'fill')+
        scale_fill_manual(values=color_vect)+  
        theme_classic()+labs( x = ' ', y = ' ',fill = 'annot')+ 
        scale_x_discrete(limits =unique(plot_df$sample_id),labels=as.character(c(1:3)))+ 
        theme(axis.text.y = element_text(size=rel(1.8)),axis.text.x = element_text(size=rel(1.5),angle = 90, hjust = 1, vjust = 0.5)
)
name_i="PreCD_NInfl_Ileum"
plot_df = mono_summ_list[[name_i]]
plt2<-ggplot(plot_df,aes(x= sample_id, y=cell_count,fill= annot)) +  geom_col(position = 'fill')+
        scale_fill_manual(values=color_vect)+  
        theme_classic()+labs( x = ' ', y = ' ',fill = 'annot')+ 
        scale_x_discrete(limits =unique(plot_df$sample_id),labels=as.character(c(4:5)))+ 
        theme(axis.text.y = element_text(size=rel(1.8)),axis.text.x = element_text(size=rel(1.5),angle = 90, hjust = 1, vjust = 0.5)
)
print(plt1+plt2+plot_layout(design = "AAABB"))
dev.off()
