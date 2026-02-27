library(Matrix)
library(ggplot2)
library(reticulate)
library(dplyr)
library(scales)
library(viridis)
library(shiny)
library(patchwork)
library(ggpubr)
library(ggrastr)
`%ni%` <- Negate(`%in%`)


for(spl in c("GIM23_NonInfl","GIM23_InfROI1","GIM23_InfROI2","GIM38_NonInfl")){
  
  load(file=paste0("./Grouped_objects/Xenium/dgcmtx_",spl,".rd")) # Take the one already normalized like in the App
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
  load(file=paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
  set.seed(42)
  
  
  meta_data<-meta_annot
  meta_data$cell_id<-rownames(meta_data)
  meta_data$currentscore<-NA
  input_geneask="=EPCAM+AGR3+GPX2+MMP7"
  legd<-paste0("Additive expression")
  gn2sum<-unlist(strsplit(gsub("^=","",input_geneask),split="\\+"))
  meta_data$currentscore<-Matrix::rowSums(dgcmtxcounts[,which(colnames(dgcmtxcounts)%in%gn2sum),drop=F]) 
  
  max_cur_sco<-max(meta_data$currentscore )
  max_cur_sco<-6.5
  brplt_title<-paste0("projection of current score in ",spl," slide")
  
  plot_sc3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
    geom_point(size = 0.2) +
    theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(title = brplt_title, color = legd )+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  plot_sc3 <- plot_sc3+coord_fixed()
  pal_hiro_min<-c("#000005","#230058","#5F0C6E","#B12362","#E56157","#ED894D","#FECF75","#DEDEDE")
  plot_score3 <- plot_sc3+scale_color_gradientn(colors = rev(pal_hiro_min), na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
  
  pdf(paste0("./Figures_print/Fig2E_Epith_zones_",spl,".pdf"),width=13,height=12)
  print(plot_score3) 
  dev.off()
}





## Now for the inflamed slides where we selected zones:

do_zone_epith<-function(){
  rect_data <- bind_rows(rect_1a,rect_2a,rect_3a,rect_4a,rect_5a,  
                         rect_1b,rect_2b,rect_3b,rect_4b,rect_5b, .id = "rect_id")
  
  
  meta_data<-meta_annot
  meta_data$cell_id<-rownames(meta_data)
  meta_data$currentscore<-NA
  input_geneask="=EPCAM+AGR3+GPX2+MMP7"
  legd<-paste0("Additive expression")
  gn2sum<-unlist(strsplit(gsub("^=","",input_geneask),split="\\+"))
  meta_data$currentscore<-Matrix::rowSums(dgcmtxcounts[,which(colnames(dgcmtxcounts)%in%gn2sum),drop=F]) 
  
  max_cur_sco<-max(meta_data$currentscore )
  max_cur_sco<-6.5
  brplt_title<-paste0("projection of current score in ",spl," slide")
  
  plot_sc3<-ggplot(meta_data, aes(x = x, y = y, color = currentscore)) +
    geom_point(size = 0.2) +
    theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(title = brplt_title, color = legd )+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  plot_sc3 <- plot_sc3+coord_fixed()
  pal_hiro_min<-c("#000005","#230058","#5F0C6E","#B12362","#E56157","#ED894D","#FECF75","#DEDEDE")
  plot_score3 <- plot_sc3+scale_color_gradientn(colors = rev(pal_hiro_min), na.value ="#EDEDED", limits = c(0, max_cur_sco ) )
  
  pdf(paste0("./Figures_print/Fig2E_Epith_zones_",spl,".pdf"),width=13,height=12)
  print(plot_score3 + 
    geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              color = "black", alpha = 0.1, inherit.aes = FALSE,size=rel(0.7)) )
  dev.off()
  
}








spl="GIM33_InfROI1"
load(file=paste0("./Grouped_objects/Xenium/dgcmtx_",spl,".rd")) # Take the one already normalized like in the App
load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(file=paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)
# V4 gradients further from ulceration
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

do_zone_epith()


spl="GIM33_InfROI2"
load(file=paste0("./Grouped_objects/Xenium/dgcmtx_",spl,".rd")) # Take the one already normalized like in the App
load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(file=paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)
# V4 gradients further from ulceration
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

do_zone_epith()


spl="GIM38_InfROI1"
load(file=paste0("./Grouped_objects/Xenium/dgcmtx_",spl,".rd")) # Take the one already normalized like in the App
load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(file=paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)
# V4 gradients further from ulceration
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

do_zone_epith()


spl="GIM38_InfROI2"
load(file=paste0("./Grouped_objects/Xenium/dgcmtx_",spl,".rd")) # Take the one already normalized like in the App
load(file=paste0("./Grouped_objects/Xenium/meta_annot_june5_",spl,".rd"))
load(file=paste0("./Grouped_objects/Xenium/meta_annot_",spl, ".rd"))
set.seed(42)
# V4 gradients further from ulceration
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

do_zone_epith()

