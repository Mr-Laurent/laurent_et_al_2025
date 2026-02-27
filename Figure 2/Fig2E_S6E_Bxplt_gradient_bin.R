library(Matrix)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(rstatix) # for stat_compare
library(tidyr)  # For pivot_wider


`%ni%` <- Negate(`%in%`)

load("./Grouped_objects/Xenium/zonegradient_v4_33roi1_500cells.rd")
load("./Grouped_objects/Xenium/zonegradient_v4_33roi2_500cells.rd")
load("./Grouped_objects/Xenium/zonegradient_v4_38roi1_500cells.rd")
load("./Grouped_objects/Xenium/zonegradient_v4_38roi2_500cells.rd")

all_zones<-dplyr::bind_cols(zonecounts_33roi1[rownames(zonecounts_33roi1),],zonecounts_33roi2[rownames(zonecounts_33roi1),],
                            zonecounts_38roi1[rownames(zonecounts_33roi1),],zonecounts_38roi2[rownames(zonecounts_33roi1),] )
all_zones_prop<-100*sweep(all_zones,2,colSums(all_zones),`/`) 

# colSums(all_zones)[which(colSums(all_zones)<500)]
# sub_grad_33_4a sub_grad_33_5a sub_grad_33_5c sub_grad_33_5d sub_grad_38_4a sub_grad_38_5a sub_grad_38_3b sub_grad_38_4b sub_grad_38_5c sub_grad_38_4d sub_grad_38_5d 
#            456            459            419            440            469            389            463            441            262            397            355 
## /!\ 11 of the 40 zones are under 500 cells







Make_prism_boxplot<-function(merged_df, title_name){
  ggplot(merged_df, aes(x =condition, y =value, fill=condition))+ 
    stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
    geom_boxplot(lwd=1.15, color="black", fatten = NULL, outlier.colour = "white",width=0.6)+
    scale_x_discrete(limits = levels(merged_df$condition) , labels= labels)+
    scale_fill_manual(values = palette2)+ scale_y_continuous(n.breaks = 5) +
    ggprism::theme_prism(base_size = 16)+ ggtitle(title_name) + 
    stat_summary(fun = median,geom = "crossbar",width = 0.57, fatten = 2)+ 
    geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
    theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),
          axis.text.y = element_text(size=rel(2),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=0.5),
          axis.ticks.length = unit(15,"points"),
          plot.margin = unit(c(20, 4, 4, 4),"points"))
  
}

Make_prism_barplot<-function(merged_df,title_name,yminerrbar){
  df_sum <- merged_df %>%
    group_by(condition) %>%
    summarise(mean  = mean(value, na.rm = TRUE),
              se    = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
              upper = mean + se,
              lower = mean - se, .groups = "drop"  )
  
  ggplot() + geom_col( data = df_sum,
                       aes(x = condition, y = mean, fill = condition, color = condition),
                       alpha = 1,width = 0.8,size=1,
                       position = position_dodge(width = 0.8)) +scale_color_manual(values = rep("black",7))+
    geom_jitter(data = merged_df,aes(x = condition, y = value),shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
    geom_errorbar(data = df_sum, aes(x = condition, ymin =eval(parse(text=paste0(yminerrbar))), ymax = upper), width = 0.4, size = rel(1))+
    scale_fill_manual(values = palette2) +
    scale_x_discrete(limits = levels(merged_df$condition), labels= labels) +
    ggtitle(title_name) + 
    ggprism::theme_prism(base_size = 16)+ labs(y= title_name)+
    theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),
          axis.text.y = element_text(size=rel(2),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=0.5),
          axis.ticks.length.y = unit(15,"points"),axis.ticks.length.x = unit(7.5,"points"),
          plot.margin = unit(c(20, 4, 4, 4),"points"))
  
}

pdf(paste0("./Figures_print/Fig2E_subtype_distrib_prism.pdf"),width = 3.5,height = 4.5)
palette2<-rep("black",5)
for(i in rownames(all_zones_prop) ){
  df<-all_zones_prop[i,which(gsub("^.*[_](.*)[_].*[_].*","\\1",colnames(all_zones_prop))=="grad")]
  df_long <- df %>%
    mutate(sample_id = rownames(df)) %>% #1:n()) %>%
    pivot_longer( cols = -sample_id,names_to = "sample_name", values_to = "value" ) %>%
    mutate( condition = (sub(".*_(\\d+)[a-z]$", "\\1", sample_name))     )
  
  df_long$condition<-factor(df_long$condition, levels=c("1","2","3","4","5"))
  # print(Make_prism_boxplot(merged_df=df_long, title_name=i) )
  # print(Make_prism_barplot(merged_df=df_long, title_name=i,yminerrbar="mean")+
  #   coord_cartesian(ylim = c(0, max(df_long$value)*1.1),clip = 'on')+
  #   scale_y_continuous(expand = c(0, 0)) )
  print(Make_prism_barplot(merged_df=df_long, title_name=i,yminerrbar="mean")+  # Keep this version in paper
    coord_cartesian(ylim = c(0, max(df_long$value)*1.1),clip = 'off') )
   
 }
dev.off()


### Ratio Monos/Macs
pdf(paste0("./Figures_print/Fig2E_S6E_grdt_ratio_prism.pdf"),width = 9,height = 7)

all_zones2<-all_zones[c("Macs","Momacs_transition","Mono","Mono1","Mono3","Mono5"),]
all_zones2<-rbind(all_zones2,"Monos"=colSums(all_zones2[c("Momacs_transition","Mono","Mono1","Mono3","Mono5"),]) )
all_zones2<-all_zones2[-which(rownames(all_zones2)%in%c("Momacs_transition","Mono","Mono1","Mono3","Mono5")),]
Mono_ratio_p1=(all_zones2["Monos",])/(all_zones2["Macs",])
df_l_mora <- Mono_ratio_p1 %>%
  mutate(sample_id = rownames(Mono_ratio_p1)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  
plot1<-ggplot(df_l_mora, aes(x = factor(group), y = value,  color=sample_id,
                             fill = factor(sample_id,level=c("Macs","Monos") ))) +scale_y_continuous(trans = "log1p")+
  geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
  labs(title = "ratio Mono/Macs in 'gradients'", x="", y = "ratio MoMacs") + #x = "1=Ulceration --> 5=Deeper in the tissue",
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none")+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values = c("Monos"="#CF6910"))+scale_color_manual(values = c("Monos"="#453020"))

plot1_p<-ggplot(df_l_mora, aes(x = factor(group), y = value,  fill=factor(group) )) +
  scale_y_continuous(trans = "log1p",n.breaks = 6)+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=rep("black",5), fatten = NULL, outlier.colour = "white",width=0.6)+
  ggprism::theme_prism(base_size = 16)+
  scale_color_manual(values = rep("black",5) ) +
  stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = sample_id), fatten = 4)+
  geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
  scale_color_manual(values = rep("white",5) ) + 
  labs(title = "ratio Mono/Macs\nin 'gradients'", x="", y = "ratio Mono/Macs") + #x = "1=Ulceration --> 5=Deeper in the tissue",
  theme(legend.position = "none",axis.title.x =element_blank(),
        axis.text.y = element_text(size=rel(1.5),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(1.5)),
        axis.ticks.length = unit(15,"points"),
        plot.margin = unit(c(20, 4, 4, 4),"points"))

# Transformed ratio (to have same distance from 0 if ratio is under 1:1)
log_ratio <- sapply(1:ncol(all_zones2), function(i) {
  monos <- all_zones2["Monos",i]
  macs <- all_zones2["Macs",i]
  if (monos > macs) {
    log_ratio <- log10(monos/macs)
  } else {
    log_ratio <- -log10(macs/monos)
  }
  return(log_ratio)
})
log_ratio<-as.data.frame(t(log_ratio))
colnames(log_ratio)<-colnames(Mono_ratio_p1)
rownames(log_ratio)<-"log ratio"

df_l_mora <- log_ratio %>%
  mutate(sample_id = rownames(log_ratio)) %>% 
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  

plot2<-ggplot(df_l_mora, aes(x = factor(group), y = value,  color=sample_id,
                             fill = factor(sample_id,level=c("log ratio") ))) + 
  geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
  labs(title = "log10 ratio Mono/Macs in 'gradients'", x = "", y = "log10 ratio Mono/Macs") +
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none",
                     text = element_text(size = 28),axis.title.x = element_text(size=18),
                     title = element_text(size=18),
                     axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)) )+
  scale_fill_manual(values = c("log ratio"="#CF6910"))+scale_color_manual(values = c("log ratio"="#453020"))+
  geom_hline(yintercept = 0,linetype="dashed",color="#6D6D6D")

plot2_p<-ggplot(df_l_mora, aes(x = factor(group), y = value,  fill=factor(group) )) +
  scale_y_continuous(n.breaks = 5)+
  geom_hline(size=rel(1),yintercept = 0,linetype="dashed",color="#6D6D6D")+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=rep("black",5), fatten = NULL, outlier.colour = "white",width=0.6)+
  ggprism::theme_prism(base_size = 16)+
  scale_color_manual(values = rep("black",5) ) +
  stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = sample_id), fatten = 4)+
  geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
  scale_color_manual(values = rep("white",5) ) + 
  labs(title = "log10 ratio Mono/Macs\nin 'gradients'", x="", y = "log10 ratio Mono/Macs") + #x = "1=Ulceration --> 5=Deeper in the tissue",
  theme(legend.position = "none",axis.title.x =element_blank(),
        axis.text.y = element_text(size=rel(1.5),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(1.5)),
        axis.ticks.length = unit(15,"points"),
        plot.margin = unit(c(20, 4, 4, 4),"points"))

print(plot1+plot2)
print(plot1_p+plot2_p)

## Ratio Mono1/Mono3

all_zones2<-all_zones[c("Mono1","Mono3"),]
Mono_ratio_p1=(1+all_zones2["Mono1",])/(1+all_zones2["Mono3",])
df_l_mora <- Mono_ratio_p1 %>%
  mutate(sample_id = rownames(Mono_ratio_p1)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  
plot1<-ggplot(df_l_mora, aes(x = factor(group), y = value, color=sample_id, 
                             fill = factor(sample_id ))) +scale_y_continuous(trans = "log1p")+
  geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
  labs(title = "ratio Mono1/Mono3 in 'gradients'", x = "1=Ulceration --> 5=Deeper in the tissue", y = "(+1) ratio MoMacs") +
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none",
                     text = element_text(size = 28),axis.title.x = element_text(size=18),
                     title = element_text(size=18),
                     axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)))+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values = alpha(c("Mono1"="#8E8EA3"),4))+scale_color_manual(values = c("Mono1"="#555568"))




# Transformed ratio (to have same distance from 0 if ratio is under 1:1)
log_ratio <- sapply(1:ncol(all_zones2), function(i) {
  val1 <- all_zones2["Mono1",i]+1   # add +1 to avoid infinite ratios ?
  val2 <- all_zones2["Mono3",i]+1   # add +1 to avoid infinite ratios ?
  if (val1 > val2) {log_ratio <- log10(val1/val2)
  } else {  log_ratio <- -log10(val2/val1)
  }
  return(log_ratio)
})
log_ratio<-as.data.frame(t(log_ratio))
colnames(log_ratio)<-colnames(Mono_ratio_p1)
rownames(log_ratio)<-"log ratio"

df_l_mora <- log_ratio %>%
  mutate(sample_id = rownames(log_ratio)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  
plot2<-ggplot(df_l_mora, aes(x = factor(group), y = value, color=sample_id,
                             fill = factor(sample_id,level=c("log ratio") ))) + #scale_y_continuous(trans = "log1p")+
  geom_boxplot(width=0.8, size=1.5, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.5), alpha=0.5 )+
  labs(title = "log10 ratio Mono1/Mono3 in 'gradients'", x = "1=Ulceration --> 5=Deeper in the tissue", y = "log10 (+1) ratio Mono1/Mono3") +
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1.5)),legend.title = element_blank(),legend.position = "none",
                     text = element_text(size = 28),axis.title.x = element_text(size=18),
                     title = element_text(size=18),
                     axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)) )+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values = alpha(c("log ratio"="#8E8EA3"),4))+scale_color_manual(values = c("log ratio"="#555568"))+
  geom_hline(yintercept = 0,size=rel(1),linetype="dashed",color="#ADADAD")
print(plot1+plot2)


plot3_p<-ggplot(df_l_mora, aes(x = factor(group), y = value,  fill=factor(group) )) +
  scale_y_continuous(n.breaks = 5)+
  geom_hline(size=rel(1),yintercept = 0,linetype="dashed",color="#6D6D6D")+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=rep("black",5), fatten = NULL, outlier.colour = "white",width=0.6)+
  ggprism::theme_prism(base_size = 16)+
  scale_color_manual(values = rep("black",5) ) +
  stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = sample_id), fatten = 4)+
  geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
  scale_color_manual(values = rep("white",5) ) + 
  labs(title = "log10 ratio IFIM/IRM\nin 'gradients'", x="", y = "log10 ratio IFIM/IRM") + #x = "1=Ulceration --> 5=Deeper in the tissue",
  theme(legend.position = "none",axis.title.x =element_blank(),
        axis.text.y = element_text(size=rel(1.5),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(1.5)),
        axis.ticks.length = unit(15,"points"),
        plot.margin = unit(c(20, 4, 4, 4),"points"))

# write.csv2(df_l_mora,file="FigS5E_log10IFIM_IRM_ratio_datapoints.csv",quote=F,row.names=F)



## Ratio Mono1/Mono5

all_zones2<-all_zones[c("Mono1","Mono5"),]
Mono_ratio_p1=(1+all_zones2["Mono1",])/(1+all_zones2["Mono5",])
df_l_mora <- Mono_ratio_p1 %>%
  mutate(sample_id = rownames(Mono_ratio_p1)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  
plot1<-ggplot(df_l_mora, aes(x = factor(group), y = value,  color=sample_id,
                             fill = factor(sample_id ))) +scale_y_continuous(trans = "log1p")+
  geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
  labs(title = "ratio Mono1/Mono5 in 'gradients'", x = "1=Ulceration --> 5=Deeper in the tissue", y = "(+1) ratio MoMacs") +
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none",
                     text = element_text(size = 28),axis.title.x = element_text(size=18),
                     title = element_text(size=18),
                     axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)))+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values = alpha(c("Mono1"="#8E8EA3"),4))+scale_color_manual(values = c("Mono1"="#555568"))


# Transformed ratio (to have same distance from 0 if ratio is under 1:1)
log_ratio <- sapply(1:ncol(all_zones2), function(i) {
  val1 <- all_zones2["Mono1",i]+1
  val2 <- all_zones2["Mono5",i]+1
  if (val1 > val2) {log_ratio <- log10(val1/val2)
  } else {  log_ratio <- -log10(val2/val1)
  }
  return(log_ratio)
})
log_ratio<-as.data.frame(t(log_ratio))
colnames(log_ratio)<-colnames(Mono_ratio_p1)
rownames(log_ratio)<-"log ratio"

df_l_mora <- log_ratio %>%
  mutate(sample_id = rownames(log_ratio)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  
plot2<-ggplot(df_l_mora, aes(x = factor(group), y = value,  color=sample_id,
                             fill = factor(sample_id,level=c("log ratio") ))) + #scale_y_continuous(trans = "log1p")+
  geom_boxplot(width=0.8, size=1.5, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.5), alpha=0.5 )+
  labs(title = "log10 ratio Mono1/Mono5 in 'gradients'", x = "1=Ulceration --> 5=Deeper in the tissue", y = "log10 (+1) ratio Mono1/Mono3") +
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1.5)),legend.title = element_blank(),legend.position = "none",
                     text = element_text(size = 28),axis.title.x = element_text(size=18),
                     title = element_text(size=18),
                     axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)))+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values = alpha(c("log ratio"="#8E8EA3"),4))+scale_color_manual(values = c("log ratio"="#555568"))+
  geom_hline(yintercept = 0,size=rel(1),linetype="dashed",color="#ADADAD")

print(plot1+plot2)

plot4_p<-ggplot(df_l_mora, aes(x = factor(group), y = value,  fill=factor(group) )) +
  scale_y_continuous(n.breaks = 5)+
  geom_hline(size=rel(1),yintercept = 0,linetype="dashed",color="#6D6D6D")+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=rep("black",5), fatten = NULL, outlier.colour = "white",width=0.6)+
  ggprism::theme_prism(base_size = 16)+
  scale_color_manual(values = rep("black",5) ) +
  stat_summary(fun = median,geom = "crossbar",width = 0.57,aes(color = sample_id), fatten = 4)+
  geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
  scale_color_manual(values = rep("white",5) ) + 
  labs(title = "log10 ratio IFIM/Late-Mono\nin 'gradients'", x="", y = "log10 ratio IFIM/Late-Mono") + #x = "1=Ulceration --> 5=Deeper in the tissue",
  theme(legend.position = "none",axis.title.x =element_blank(),
        axis.text.y = element_text(size=rel(1.5),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(1.5)),
        axis.ticks.length = unit(15,"points"),
        plot.margin = unit(c(20, 4, 4, 4),"points"))

print(plot3_p+plot4_p)

# write.csv2(df_l_mora,file="FigS5E_log10IFIM_Late_ratio_datapoints.csv",quote=F,row.names=F)



dev.off()






#Redo in small format for 22E Mono/Mac ration:

pdf(paste0("./Figures_print/Fig2E_momac_ratio_prism.pdf"),width = 3.5,height = 4.5)

all_zones2<-all_zones[c("Macs","Momacs_transition","Mono","Mono1","Mono3","Mono5"),]
all_zones2<-rbind(all_zones2,"Monos"=colSums(all_zones2[c("Momacs_transition","Mono","Mono1","Mono3","Mono5"),]) )
all_zones2<-all_zones2[-which(rownames(all_zones2)%in%c("Momacs_transition","Mono","Mono1","Mono3","Mono5")),]
Mono_ratio_p1=(all_zones2["Monos",])/(all_zones2["Macs",])
df_l_mora <- Mono_ratio_p1 %>%
  mutate(sample_id = rownames(Mono_ratio_p1)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id, names_to = "sample_name", values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  


# Transformed ratio (to have same distance from 0 if ratio is under 1:1)
log_ratio <- sapply(1:ncol(all_zones2), function(i) {
  monos <- all_zones2["Monos",i]
  macs <- all_zones2["Macs",i]
  if (monos > macs) { log_ratio <- log10(monos/macs) } 
  else { log_ratio <- -log10(macs/monos) }
  return(log_ratio)
})
log_ratio<-as.data.frame(t(log_ratio))
colnames(log_ratio)<-colnames(Mono_ratio_p1)
rownames(log_ratio)<-"log ratio"

df_l_mora <- log_ratio %>%
  mutate(sample_id = rownames(log_ratio)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  

ggplot(df_l_mora, aes(x = factor(group), y = value,  fill=factor(group) )) +
  scale_y_continuous(n.breaks = 5)+
  geom_hline(size=rel(1),yintercept = 0,linetype="dashed",color="#6D6D6D")+
  stat_boxplot(geom = "errorbar", width = 0.3,size=1) +
  geom_boxplot(lwd=1.15, color="black", fill=rep("black",5), fatten = NULL, outlier.colour = "white",width=0.8)+
  ggprism::theme_prism(base_size = 16)+
  scale_color_manual(values = rep("black",5) ) +
  stat_summary(fun = median,geom = "crossbar",width = 0.75,aes(color = sample_id), fatten = 4)+
  geom_jitter(shape=21, position=position_jitter(width = 0.20),size=3, fill="white",stroke=1.5) +
  scale_color_manual(values = rep("white",5) ) + 
  labs(title = "log10 ratio Mono/Macs\nin 'gradients'", x="", y = "log10 ratio Mono/Macs") + #x = "1=Ulceration --> 5=Deeper in the tissue",
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),
        axis.text.y = element_text(size=rel(2),margin = margin(0,8,0,0, unit = "points")),
        axis.text.x = element_text(size=rel(2)),
        axis.ticks.length.y = unit(15,"points"),axis.ticks.length.x = unit(7.5,"points"),
        plot.margin = unit(c(20, 4, 4, 4),"points"))



dev.off()
