library(ggplot2)
library(ggpubr)

setwd("G:/Mon Drive/UG_metacells/Figures paper Aout/")
load("Grouped_objects/MacsEnr_z_mymod_150_metadata_goodlog_v5.rd") # Already has signature scores compared to v3
load("Grouped_objects/Macs_dgCmc_16nov22_counts.rd")
load("Grouped_objects/genelist_16nov22_Macs.rd") 



# Simple scatter plot comparing the distribution of 2 programs for each metacell:
versus_sig_plot<-function(x_ax,y_ax){
  print(ggplot(metadf_z, aes(x=!!sym(x_ax) , y=!!sym(y_ax) )) +
          geom_point(size=rel(2.5),colour="#2C8587" )+theme_minimal()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=rel(1.5))) )
}
## If a 3rd program needs to be tested, can be color coded:
# versus_3Dsig_plot<-function(x_ax,y_ax,z_ax){
#   print(ggplot(metadf_z, aes(x=!!sym(x_ax) , y=!!sym(y_ax), color=!!sym(z_ax) )) +
#           geom_point(size=rel(2.5),)+theme_minimal()+
#           theme(panel.border = element_rect(colour = "black", fill=NA, size=rel(1.5)))+
#           scale_color_gradientn(colours = rev(rainbow(15))[4:15])+theme_dark()+ labs(color = "score")+
#           ggtitle(paste("projection of ",z_ax," enrichment in Macs metacells",sep="")) )
# }


# Spearman is best as it's not linear :
plot1<-(versus_sig_plot("Signature_XV","Signature_XII") )
# plot2<-(versus_sig_plot("Signature_XV","Signature_XII")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),method="pearson") )
plot3<-(versus_sig_plot("Signature_XV","Signature_XII")+stat_cor(method="spearman",cor.coef.name = "rho",label.x.npc = "center",label.y.npc = "top" ) )

plot4<-(versus_sig_plot("Signature_XV","Signature_VI") )
# plot5<-(versus_sig_plot("Signature_XV","Signature_VI")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),method="pearson") )
plot6<-(versus_sig_plot("Signature_XV","Signature_VI")+stat_cor(method="spearman",cor.coef.name = "rho",label.x.npc = "center",label.y.npc = "top") )

pdf(file="./Figure 1/Fig1I_Sctplt_sig_green.pdf",width = 6,height = 6)
print(plot1)
# print(plot2)
print(plot3)
print(plot4)
# print(plot5)
print(plot6)
dev.off()

