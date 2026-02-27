## First, for each slide, we compute 
library(ggplot2)
library(ggpubr)

spl="GIM33_InfROI1"

load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)

meta_annot_perso$CXCL9<-dgcmtx_raw_counts[,"CXCL9"]
meta_annot_perso$CXCL10<-dgcmtx_raw_counts[,"CXCL10"]

# V4 gradients further from ulceration   (Boxplot_prop_subtypes_Zones_select.R)
rect_1a <- data.frame (xmin=2800, xmax=3050, ymin=13500, ymax=13750)
rect_2a <- data.frame (xmin=2750, xmax=3000, ymin=13750, ymax=14000)
rect_3a <- data.frame (xmin=2700, xmax=2950, ymin=14000, ymax=14250)
rect_4a <- data.frame (xmin=2650, xmax=2900, ymin=14250, ymax=14500)
rect_5a <- data.frame (xmin=2600, xmax=2850, ymin=14500, ymax=14750)

rect_1b <- data.frame (xmin=1050, xmax=1300, ymin=13350, ymax=13600)
rect_2b <- data.frame (xmin=1100, xmax=1350, ymin=13600, ymax=13850)
rect_3b <- data.frame (xmin=1150, xmax=1400, ymin=13850, ymax=14100)
rect_4b <- data.frame (xmin=1200, xmax=1450, ymin=14100, ymax=14350)
rect_5b <- data.frame (xmin=1250, xmax=1500, ymin=14350, ymax=14600)

rect_data <- bind_rows(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,  
                       rect_1b,rect_2b,rect_3b,rect_4b,rect_5b, .id = "rect_id")

rect_list <- list(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,rect_1b,rect_2b,rect_3b,rect_4b,rect_5b)
subsets <- list()
nb_select=500  # Subset for 500 cells 
for (i in seq_along(rect_list)) { 
  rect <- rect_list[[i]]
  subset <- subset(meta_annot_perso, x>=rect$xmin & x<=rect$xmax & y>=rect$ymin & y<=rect$ymax)
  if (nrow(subset) > nb_select) { # Subset nb_select rows (if possible)
    subset <- subset[sample(nrow(subset), nb_select), ]
  }
  subsets[[i]] <- subset }

subset_names <- c("sub_grad_33_1a", "sub_grad_33_2a", "sub_grad_33_3a",
                  "sub_grad_33_4a", "sub_grad_33_5a",
                  "sub_grad_33_1b", "sub_grad_33_2b", "sub_grad_33_3b",
                  "sub_grad_33_4b", "sub_grad_33_5b")
names(subsets) <- subset_names

lapply(subsets,function(x) dim(x))
tib_meta<-bind_rows(lapply(names(subsets), function(x) tibble(grp = x, num = subsets[[x]])))  #paste all subsets metadata in one tibble with a column info of bin
tib_meta$grp2<-paste0("bin_",gsub("(.*)_(.*)[a-z]$","\\2",tib_meta$grp))

save(tib_meta,file=paste0("./Grouped_objects/Xenium/tib_meta_v4_33roi1.rd"))





spl="GIM33_InfROI2"

load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)

meta_annot_perso$CXCL9<-dgcmtx_raw_counts[,"CXCL9"]
meta_annot_perso$CXCL10<-dgcmtx_raw_counts[,"CXCL10"]

# V4 gradients further from ulceration   (Boxplot_prop_subtypes_Zones_select.R)
rect_1a <- data.frame (xmin=12000, xmax=12250, ymin=18300, ymax=18550)
rect_2a <- data.frame (xmin=11900, xmax=12150, ymin=18050, ymax=18300)
rect_3a <- data.frame (xmin=11800, xmax=12050, ymin=17800, ymax=18050)
rect_4a <- data.frame (xmin=11700, xmax=11950, ymin=17550, ymax=17800)
rect_5a <- data.frame (xmin=11600, xmax=11850, ymin=17300, ymax=17550)

rect_1b <- data.frame (xmin=12350, xmax=12600, ymin=18200, ymax=18450)
rect_2b <- data.frame (xmin=12350, xmax=12600, ymin=17950, ymax=18200)
rect_3b <- data.frame (xmin=12350, xmax=12600, ymin=17700, ymax=17950)
rect_4b <- data.frame (xmin=12350, xmax=12600, ymin=17450, ymax=17700)
rect_5b <- data.frame (xmin=12350, xmax=12600, ymin=17200, ymax=17450)

rect_data <- bind_rows(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,  
                       rect_1b,rect_2b,rect_3b,rect_4b,rect_5b, .id = "rect_id")

rect_list <- list(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,rect_1b,rect_2b,rect_3b,rect_4b,rect_5b)
subsets <- list()
for (i in seq_along(rect_list)) { 
  rect <- rect_list[[i]]
  subset <- subset(meta_annot_perso, x>=rect$xmin & x<=rect$xmax & y>=rect$ymin & y<=rect$ymax)
  if (nrow(subset) > nb_select) { # Subset nb_select rows (if possible)
    subset <- subset[sample(nrow(subset), nb_select), ]
  }
  subsets[[i]] <- subset }

subset_names <- c("sub_grad_33_1c", "sub_grad_33_2c", "sub_grad_33_3c",
                  "sub_grad_33_4c", "sub_grad_33_5c",
                  "sub_grad_33_1d", "sub_grad_33_2d", "sub_grad_33_3d",
                  "sub_grad_33_4d", "sub_grad_33_5d")
names(subsets) <- subset_names

lapply(subsets,function(x) dim(x))
tib_meta<-bind_rows(lapply(names(subsets), function(x) tibble(grp = x, num = subsets[[x]])))  #paste all subsets metadata in one tibble with a column info of bin
tib_meta$grp2<-paste0("bin_",gsub("(.*)_(.*)[a-z]$","\\2",tib_meta$grp))

save(tib_meta,file=paste0("./Grouped_objects/Xenium/tib_meta_v4_33roi2.rd"))


spl="GIM38_InfROI1"

load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)

meta_annot_perso$CXCL9<-dgcmtx_raw_counts[,"CXCL9"]
meta_annot_perso$CXCL10<-dgcmtx_raw_counts[,"CXCL10"]

# V4 gradients further from ulceration   (Boxplot_prop_subtypes_Zones_select.R)
rect_1a <- data.frame (xmin=16200, xmax=16450, ymin=3900, ymax=4150)
rect_2a <- data.frame (xmin=16200, xmax=16450, ymin=3650, ymax=3900)
rect_3a <- data.frame (xmin=16200, xmax=16450, ymin=3400, ymax=3650)
rect_4a <- data.frame (xmin=16200, xmax=16450, ymin=3150, ymax=3400)
rect_5a <- data.frame (xmin=16200, xmax=16450, ymin=2900, ymax=3150)

rect_1b <- data.frame (xmin=15500, xmax=15750, ymin=4000, ymax=4250)
rect_2b <- data.frame (xmin=15400, xmax=15650, ymin=3750, ymax=4000)
rect_3b <- data.frame (xmin=15300, xmax=15550, ymin=3500, ymax=3750)
rect_4b <- data.frame (xmin=15200, xmax=15450, ymin=3250, ymax=3500)
rect_5b <- data.frame (xmin=15100, xmax=15350, ymin=3000, ymax=3250)

rect_data <- bind_rows(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,  
                       rect_1b,rect_2b,rect_3b,rect_4b,rect_5b, .id = "rect_id")

rect_list <- list(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,rect_1b,rect_2b,rect_3b,rect_4b,rect_5b)
subsets <- list()
for (i in seq_along(rect_list)) { 
  rect <- rect_list[[i]]
  subset <- subset(meta_annot_perso, x>=rect$xmin & x<=rect$xmax & y>=rect$ymin & y<=rect$ymax)
  if (nrow(subset) > nb_select) { # Subset nb_select rows (if possible)
    subset <- subset[sample(nrow(subset), nb_select), ]
  }
  subsets[[i]] <- subset }

subset_names <- c("sub_grad_38_1a", "sub_grad_38_2a", "sub_grad_38_3a",
                  "sub_grad_38_4a", "sub_grad_38_5a",
                  "sub_grad_38_1b", "sub_grad_38_2b", "sub_grad_38_3b",
                  "sub_grad_38_4b", "sub_grad_38_5b")
names(subsets) <- subset_names

lapply(subsets,function(x) dim(x))
tib_meta<-bind_rows(lapply(names(subsets), function(x) tibble(grp = x, num = subsets[[x]])))  #paste all subsets metadata in one tibble with a column info of bin
tib_meta$grp2<-paste0("bin_",gsub("(.*)_(.*)[a-z]$","\\2",tib_meta$grp))

save(tib_meta,file=paste0("./Grouped_objects/Xenium/tib_meta_v4_38roi1.rd"))





spl="GIM38_InfROI2"

load(paste0("./Grouped_objects/Xenium/dgcmtx_raw_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)

meta_annot_perso$CXCL9<-dgcmtx_raw_counts[,"CXCL9"]
meta_annot_perso$CXCL10<-dgcmtx_raw_counts[,"CXCL10"]

# V4 gradients further from ulceration   (Boxplot_prop_subtypes_Zones_select.R)
rect_1a <- data.frame (xmin=24400, xmax=24650, ymin=4100, ymax=4350)
rect_2a <- data.frame (xmin=24650, xmax=24900, ymin=4300, ymax=4550)
rect_3a <- data.frame (xmin=24900, xmax=25150, ymin=4500, ymax=4750)
rect_4a <- data.frame (xmin=25150, xmax=25400, ymin=4700, ymax=4950)
rect_5a <- data.frame (xmin=25400, xmax=25650, ymin=4900, ymax=5150)

rect_1b <- data.frame (xmin=24950, xmax=25200, ymin=3500, ymax=3750)
rect_2b <- data.frame (xmin=25200, xmax=25450, ymin=3350, ymax=3600)
rect_3b <- data.frame (xmin=25450, xmax=25700, ymin=3200, ymax=3450)
rect_4b <- data.frame (xmin=25700, xmax=25950, ymin=3050, ymax=3300)
rect_5b <- data.frame (xmin=25950, xmax=26200, ymin=2900, ymax=3150)

rect_data <- bind_rows(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,  
                       rect_1b,rect_2b,rect_3b,rect_4b,rect_5b, .id = "rect_id")

rect_list <- list(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,rect_1b,rect_2b,rect_3b,rect_4b,rect_5b)
subsets <- list()
for (i in seq_along(rect_list)) { 
  rect <- rect_list[[i]]
  subset <- subset(meta_annot_perso, x>=rect$xmin & x<=rect$xmax & y>=rect$ymin & y<=rect$ymax)
  if (nrow(subset) > nb_select) { # Subset nb_select rows (if possible)
    subset <- subset[sample(nrow(subset), nb_select), ]
  }
  subsets[[i]] <- subset }

subset_names <- c("sub_grad_38_1c", "sub_grad_38_2c", "sub_grad_38_3c",
                  "sub_grad_38_4c", "sub_grad_38_5c",
                  "sub_grad_38_1d", "sub_grad_38_2d", "sub_grad_38_3d",
                  "sub_grad_38_4d", "sub_grad_38_5d")
names(subsets) <- subset_names

lapply(subsets,function(x) dim(x))
tib_meta<-bind_rows(lapply(names(subsets), function(x) tibble(grp = x, num = subsets[[x]])))  #paste all subsets metadata in one tibble with a column info of bin
tib_meta$grp2<-paste0("bin_",gsub("(.*)_(.*)[a-z]$","\\2",tib_meta$grp))

save(tib_meta,file=paste0("./Grouped_objects/Xenium/tib_meta_v4_38roi2.rd"))

#---------------------------------------------------------------------------------------------------------#

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
    # geom_jitter(  data = merged_df,size = 2,
    #               aes(x = condition, y = value, color = condition ), # shape = batch, 
    #               position = position_jitterdodge(
    #                 jitter.width = 0.20,jitter.height = 0,dodge.width = 0.8 )) +
    geom_errorbar(data = df_sum, aes(x = condition, ymin =eval(parse(text=paste0(yminerrbar))), ymax = upper), width = 0.4, size = rel(1))+
    scale_fill_manual(values = rep("black",7)) +
    scale_x_discrete(limits = levels(merged_df$condition), labels= labels) +
    ggtitle(title_name) + 
    ggprism::theme_prism(base_size = 16)+ labs(y= title_name)+
    theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),
          axis.text.y = element_text(size=rel(2),margin = margin(0,8,0,0, unit = "points")),
          axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=0.5),
          axis.ticks.length.y = unit(15,"points"),axis.ticks.length.x = unit(7.5,"points"),
          plot.margin = unit(c(20, 4, 4, 4),"points"))
  
}



load("./Grouped_objects/Xenium/tib_meta_v4_33roi1.rd")
tib_m33r1<-aggregate(tib_meta$num$CXCL9, by=list(tib_meta$grp), FUN=mean)
load("./Grouped_objects/Xenium/tib_meta_v4_33roi2.rd")
tib_m33r2<-aggregate(tib_meta$num$CXCL9, by=list(tib_meta$grp), FUN=mean)
load("./Grouped_objects/Xenium/tib_meta_v4_38roi1.rd")
tib_m38r1<-aggregate(tib_meta$num$CXCL9, by=list(tib_meta$grp), FUN=mean)
load("./Grouped_objects/Xenium/tib_meta_v4_38roi2.rd")
tib_m38r2<-aggregate(tib_meta$num$CXCL9, by=list(tib_meta$grp), FUN=mean)

infl_gn<-do.call("rbind", list(tib_m33r1, tib_m33r2, tib_m38r1,tib_m38r2))
colnames(infl_gn)<-c("sample_id","value")
infl_gn$bin<-paste0("",gsub("(.*)_(.*)[a-z]$","\\2",infl_gn$sample_id))

# write.csv2(infl_gn,file="Gaelle_Fig4C_CXCL9_datapoints.csv",quote=F,row.names=F)

pdf(paste0("./Figures_print/Fig4C_bin_cxcl9_10.pdf"),width = 3.5,height = 4.5)
# ggplot(infl_gn, aes(x = factor(bin), y = value,  color="all",
#                              fill = "all")) +scale_y_continuous(trans = "log1p")+
#   geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
#   labs(title = "CXCL9 expression in 'gradients'", x = "", y = "expression") +
#   theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none",
#                      text = element_text(size = 28),axis.title.x = element_text(size=18),
#                      title = element_text(size=18),
#                      axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)))+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
#   scale_fill_manual(values = alpha("#8E8EA3",4))+scale_color_manual(values = "#555568")

infl_gn$condition<-infl_gn$bin
print(Make_prism_barplot(merged_df=infl_gn, title_name="CXCL9",yminerrbar="mean")+
        coord_cartesian(ylim = c(0, max(infl_gn$value)*1.1),clip = 'off') )

load("./Grouped_objects/Xenium/tib_meta_v4_33roi1.rd")
tib_m33r1<-aggregate(tib_meta$num$CXCL10, by=list(tib_meta$grp), FUN=mean)
load("./Grouped_objects/Xenium/tib_meta_v4_33roi2.rd")
tib_m33r2<-aggregate(tib_meta$num$CXCL10, by=list(tib_meta$grp), FUN=mean)
load("./Grouped_objects/Xenium/tib_meta_v4_38roi1.rd")
tib_m38r1<-aggregate(tib_meta$num$CXCL10, by=list(tib_meta$grp), FUN=mean)
load("./Grouped_objects/Xenium/tib_meta_v4_38roi2.rd")
tib_m38r2<-aggregate(tib_meta$num$CXCL10, by=list(tib_meta$grp), FUN=mean)

infl_gn<-do.call("rbind", list(tib_m33r1, tib_m33r2, tib_m38r1,tib_m38r2))
colnames(infl_gn)<-c("sample_id","value")
infl_gn$bin<-paste0("",gsub("(.*)_(.*)[a-z]$","\\2",infl_gn$sample_id))

# write.csv2(infl_gn,file="Gaelle_Fig4C_CXCL10_datapoints.csv",quote=F,row.names=F)

# ggplot(infl_gn, aes(x = factor(bin), y = value,  color="all",
#                     fill = "all")) +scale_y_continuous(trans = "log1p")+
#   geom_boxplot(width=0.8, size=0.95, alpha=0.8,outlier.alpha = 0) + geom_point(position = position_jitterdodge(jitter.width =0.2 ),size=rel(1.2), alpha=0.5 )+
#   labs(title = "CXCL10 expression in 'gradients'", x = "", y = "expression") +
#   theme_pubr()+theme(panel.grid.major = element_line(size=rel(1)),legend.title = element_blank(),legend.position = "none",
#                      text = element_text(size = 28),axis.title.x = element_text(size=18),
#                      title = element_text(size=18),
#                      axis.line = element_line(size = rel(1.5)), axis.ticks = element_line(size = rel(1.5)))+#guides(fill=guide_legend(nrow=2,byrow=TRUE))+
#   scale_fill_manual(values = alpha("#8E8EA3",4))+scale_color_manual(values = "#555568")
infl_gn$condition<-infl_gn$bin
print(Make_prism_barplot(merged_df=infl_gn, title_name="CXCL10",yminerrbar="mean")+
        coord_cartesian(ylim = c(0, max(infl_gn$value)*1.1),clip = 'off') )
dev.off()





