library(Matrix)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(rstatix) # for stat_compare
library(tidyr)  # For pivot_wider


`%ni%` <- Negate(`%in%`)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
load("Grouped_objects/Xenium/zonegradient_v4_33roi1_500cells.rd")
load("Grouped_objects/Xenium/zonegradient_v4_33roi2_500cells.rd")
load("Grouped_objects/Xenium/zonegradient_v4_38roi1_500cells.rd")
load("Grouped_objects/Xenium/zonegradient_v4_38roi2_500cells.rd")

all_zones<-dplyr::bind_cols(zonecounts_33roi1[rownames(zonecounts_33roi1),],zonecounts_33roi2[rownames(zonecounts_33roi1),],
                            zonecounts_38roi1[rownames(zonecounts_33roi1),],zonecounts_38roi2[rownames(zonecounts_33roi1),] )
all_zones_prop<-100*sweep(all_zones,2,colSums(all_zones),`/`) 

# colSums(all_zones)[which(colSums(all_zones)<500)]
# sub_grad_33_4a sub_grad_33_5a sub_grad_33_5c sub_grad_33_5d sub_grad_38_4a sub_grad_38_5a sub_grad_38_3b sub_grad_38_4b sub_grad_38_5c sub_grad_38_4d sub_grad_38_5d 
#            456            459            419            440            469            389            463            441            262            397            355 
## /!\ 11 of the 40 zones are under 500 cells



my_palette <- setNames(rep("#666666", length(rownames(all_zones_prop)) ), rownames(all_zones_prop ))
my_palette["Neutro"] <- "#4D4DCA"
my_palette["Infl Fib"] <- "#C967EA"
my_palette["Mono1"] <- "#ad2718"
my_palette["CD8 eff."] <- "#1B458E"
my_palette["iaDC"] <- "#632523"
my_palette["aDC"] <- "#D45F03"
my_palette["Submucosa Myofibro"] <- "#BE8ED4"
my_palette["Naive B cell"] <- "#FFCC56"
my_palette["addNKlike"] <- "#2BB542"
my_palette["GClike B cell"] <- "#DDA493"
my_palette["Plasma Cell"] <- "#FF8061"
levels(my_palette)<-names(my_palette)

pdf(paste0("./Figure 2/Fig2FGI_subtype_distrib.pdf"),width = 4,height = 6)

for(i in rownames(all_zones_prop) ){
  df<-all_zones_prop[i,which(gsub("^.*[_](.*)[_].*[_].*","\\1",colnames(all_zones_prop))=="grad")]
  df_long <- df %>%
    mutate(sample_id = rownames(df)) %>% #1:n()) %>%
    pivot_longer( cols = -sample_id,names_to = "sample_name", values_to = "value" ) %>%
    mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name))     ) # Get final number to group 1 together, 2 tgthr etc
  
  print( ggplot(df_long, aes(x = factor(group), y = value, fill = sample_id)) +
           scale_fill_manual(values =my_palette[i])+
           geom_boxplot(width=0.8, size=rel(1.2), alpha=0.8,outlier.alpha = 0) + 
           geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.5), alpha=0.5 )+
           labs(title = paste0("% of ",i,"\nin 'gradients'"),x="",  y = "% of captured cells") +  #x = "1=Ulceration --> 5=Deeper in the tissue",
           theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.position = "none",
                              text = element_text(size = 28),axis.title.x = element_text(size=18),
                              title = element_text(size=18),
                              axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)) )#+scale_fill_manual(values = c("Mono1"="#874037","Mono3"="#F77774","Mono5"="#50AE86"))
  ) 
}
dev.off()

### Ratio Monos/Macs

pdf(paste0("./Figure 2/Fig2E_S5E_gradient_ratio.pdf"),width = 9,height = 7)

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
  mutate(sample_id = rownames(log_ratio)) %>% #1:n()) %>%
  pivot_longer(cols = -sample_id,
               names_to = "sample_name",
               values_to = "value" ) %>%
  mutate( group = as.numeric(sub(".*_(\\d+)[a-z]$", "\\1", sample_name)) )  
plot2<-ggplot(df_l_mora, aes(x = factor(group), y = value,  color=sample_id,
                             fill = factor(sample_id,level=c("log ratio") ))) + #scale_y_continuous(trans = "log1p")+
  geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
  labs(title = "log10 ratio Mono/Macs in 'gradients'", x = "", y = "log10 ratio Mono/Macs") +
  theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none",
                     text = element_text(size = 28),axis.title.x = element_text(size=18),
                     title = element_text(size=18),
                     axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)) )+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values = c("log ratio"="#CF6910"))+scale_color_manual(values = c("log ratio"="#453020"))+
  geom_hline(yintercept = 0,linetype="dashed",color="#6D6D6D")

print(plot1+plot2)





## Ratio Mono1/Mono3

all_zones2<-all_zones[c("Mono1","Mono3"),]
# Mono_ratio_p1=(all_zones2["Mono1",])/(all_zones2["Mono3",])
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




## Ratio Mono1/Mono5

all_zones2<-all_zones[c("Mono1","Mono5"),]
# Mono_ratio_p1=(all_zones2["Mono1",])/(all_zones2["Mono5",])
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


dev.off()







