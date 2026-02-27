# Fig S5G
library(ComplexHeatmap)
library(ggplot2)
library(Rfast)
library(ggrepel)
library(ggpubr)
library(scales)


dfinfo<-read.csv2(file="./Grouped_objects/Buckley_momacs_annot_241217.csv",sep=",",header=T)

object="100mods_momacs_ThomasBuckley"
load(paste0("./Grouped_objects/Taurus/genelist_",object,".rd"))
load(paste0("./Grouped_objects/Taurus/metadata_",object,".rd"))
load(paste0("./Grouped_objects/Taurus/counts_",object,".rd"))
load(paste0("./Grouped_objects/Taurus/ht_",object,".rd"))

load("./Grouped_objects/Modsig_MoMacsEnr.rd")


ds<-as.matrix(dgcmtxcounts)
allsum<-Rfast::colsums(log(1+ds))

for(Mod_score in names(Modsig)[1:17] ){
  print(Mod_score)
  geneshow<-unlist(Modsig[Mod_score])
  print(setdiff(geneshow,rownames(ds)) ) # Only keeps genes detected in the count matrix
  ilist<-intersect(geneshow,rownames(ds))
  if(length(ilist)<2){currentscore<-log(1+ds[ilist,])/allsum 
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }else{currentscore<-(Rfast::colsums(log(1+ds[ilist,,drop=F]))/allsum )
  eval(parse(text=paste0("metadf_z$",Mod_score,"<-as.vector(scale(currentscore, center = TRUE, scale = TRUE))")))
  }
}

# Simple scatter plot comparing the distribution of 2 programs for each metacell:
versus_sig_plot<-function(x_ax,y_ax){
  print(ggplot(metadf_z, aes(x=!!sym(x_ax) , y=!!sym(y_ax) )) +
          geom_point(size=rel(2.5),colour="#2C8587" )+theme_minimal()+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=rel(1.5))) )
}

current_df <- metadf_z[dfinfo$metacell[which(dfinfo$annotation=="Mono1")], ]
other_df <- metadf_z[dfinfo$metacell[which(dfinfo$annotation!="Mono1")], ]

# Spearman is best as it's not linear :
plot3<-(versus_sig_plot("Sig_XV","Sig_XII")+stat_cor(method="spearman",cor.coef.name = "rho",label.x.npc = "center",label.y.npc = "top" ) )


pdf(file="./Figures_print/FigS5G_Sctplt_sig.pdf",width = 6,height = 6)
print( ggplot(other_df, aes(x =Sig_XV, y =Sig_XII )) + geom_point(size=rel(2.5), color = "gray")+
         geom_point(size=rel(2.5),data = current_df, aes(x =Sig_XV, y =Sig_XII), color = "red") +theme_minimal()+
         theme(panel.border = element_rect(colour = "black", fill=NA, size=rel(1.5))) ) 
print(plot3) # prints the rho and p value
dev.off()

