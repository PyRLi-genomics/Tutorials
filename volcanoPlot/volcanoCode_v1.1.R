library("dplyr")
library("tidyverse")
library("genefilter")
library("ggplot2")
library("grDevices")
library("ggrepel")

#Read in data frame with DE genes


volcanoPlot<-function(datainput, 
                      upColor="#bb0c00", 
                      downColor="blue2", 
                      unchangedColor="grey",
                      labelGenes=FALSE, # TRUE : Label genes
                      genes2label=NULL, # Deafult label top 10 genes, or custom number of genes, or gene names c("gene1", "gene2"...)
                      geneName.col="gene_name",
                      x.col="log2FoldChange", 
                      y.col="padj", 
                      plot.title="VolcanoPlot",
                      x.label="log2FoldChange",
                      y.label="-log10(padj)",
                      point.size=2,
                      transparency=1,
                      x.lower.limit=NULL,
                      x.Upper.limit=NULL,
                      x.break=NULL,
                      y.lower.limit=NULL,
                      y.Upper.limit=NULL,
                      y.break=NULL,
                      axis_text_size=5,
                      legend_title_size=12,
                      legend_text_size=8,
                      plot_title_size=12,
                      plot_title_face="bold",
                      label_size=5,
                      xy_axis_labels=14,
                      lfc_cut_off=1){
  
  # Data table formating
  datainput<-datainput[,c(geneName.col,x.col,y.col)]
  colnames(datainput)<-c("geneName","log2FoldChange","padj")
  datainput<-datainput[complete.cases(datainput),]
  DEstatus=ifelse(datainput["log2FoldChange"]>lfc_cut_off & datainput["padj"]<0.05,"Upregulated",
                  ifelse(datainput["log2FoldChange"] < (-lfc_cut_off) & datainput["padj"]<0.05,"Downregulated", "Not Significant"))
  DEstatus<-as.data.frame(DEstatus)
  colnames(DEstatus)<-"DEstatus"
  datainput<-cbind(datainput,DEstatus)
  x.col="log2FoldChange"
  y.col="padj"
  geneName.col="geneName"
  
  #datainput<-transform(datainput,DEstatus=ifelse(datainput[[x.col]]>1 & datainput[[y.col]]<0.05,"Upregulated",
  #                                    ifelse(datainput[[x.col]] < (-1) & datainput[[y.col]]<0.05,"Downregulated", "Not Significant")))
  #print(head(datainput,10))
  if (labelGenes==TRUE){
    #print(genes2label)
    if(is.null(genes2label)){
      message ("Top 10 genes with highest foldchange are labelled ")
      genes2label_list<- datainput %>% filter(.,abs(!!as.symbol(x.col)) > lfc_cut_off & (!!as.symbol(y.col)) < 0.05) %>% 
        arrange(desc(abs(!!as.symbol(x.col)))) %>% 
        slice_head(n=10) %>% 
        pull(geneName)
      print("Genes labeled in plot are:")
      print(paste(genes2label_list,sep=" "))
    }else if (is.vector(genes2label) & is.character(genes2label)){
      message ("Available genes from the genes2labels will be included in label")
      genes2label_list<-genes2label
      print("Genes labeled in plot are:")
      print(paste(genes2label_list,sep=" "))
  #}else if (!is.vector(genes2label & is.numeric(genes2label))){
    }else if (is.numeric(genes2label)){
      message (paste0("Top ", genes2label ," genes with highest foldchange are labelled"))
      genes2label_list<- datainput %>% filter(.,abs(!!as.symbol(x.col)) > lfc_cut_off & (!!as.symbol(y.col)) < 0.05) %>% 
        arrange(desc(abs(!!as.symbol(x.col)))) %>% 
        slice_head(n=genes2label) %>% 
        pull(geneName)
      print("Genes labeled in plot are:")
      print(paste(genes2label_list,sep=" "))
    }
  }else{
    genes2label_list=c()
  }
  
    
  CompletegeneList<-datainput[[geneName.col]]
  datainput$geneLabel <- ifelse(CompletegeneList %in%genes2label_list,CompletegeneList, NA )

  datainput<-transform(datainput,log10padj=-log10(datainput[[y.col]]))
  
  if(!is.null(x.lower.limit) & !is.null(x.Upper.limit) & !is.null(x.break)){
    x.axis.scale=sort(c(lfc_cut_off,-lfc_cut_off,seq(x.lower.limit,x.Upper.limit,x.break)))
    x_low_lim=x.lower.limit
    x_Up_lim=x.Upper.limit
  }else{
    x.axis.scale=sort(c(lfc_cut_off,-lfc_cut_off,seq(floor(min(datainput[x.col])*1.5), ceiling(max(datainput[x.col])*1.5), (((ceiling(max(datainput[x.col]))*1.5)-(floor((min(datainput[x.col]))*1.5)))/10))))
    x_low_lim=floor(min(x.axis.scale))
    x_Up_lim=ceiling(max(x.axis.scale))
  }
 # print(paste(floor(min(datainput[x.col])*1.5), ceiling(max(datainput[x.col])*1.5), (((ceiling(max(datainput[x.col]))*1.5)-(floor((min(datainput[x.col]))*1.5)))/10)))
#  print(paste(min(datainput["log10padj"])*1.15, max(datainput["log10padj"])*1.15, ((max(datainput["log10padj"])*1.15)-min(datainput["log10padj"])*1.15)/20))
  
  if(!is.null(y.lower.limit) & !is.null(y.Upper.limit) & !is.null(y.break)){
    y.axis.scale=seq(y.lower.limit,y.Upper.limit,y.break)
  }else{
    y.axis.scale=seq(min(datainput["log10padj"])*1.15, max(datainput["log10padj"])*1.15, ((max(datainput["log10padj"])*1.15)-min(datainput["log10padj"])*1.15)/20)
  }
  
  colourVector<-c(downColor, unchangedColor, upColor)
  names(colourVector)<-c("Downregulated", "Not Significant", "Upregulated")
  DEobsrvd<-sort(unique(datainput$DEstatus))
  colorIN<-colourVector[DEobsrvd]
  names(colorIN)<-NULL
  
 # print(y.axis.scale)
 
  ####Plot Code
  ggplot(data = datainput, aes_string(x = x.col, y = "log10padj", col = "DEstatus", label = "geneLabel")) +
    geom_vline(xintercept = c(-lfc_cut_off, lfc_cut_off), col = "blue", linetype = 'dashed', alpha=0.5) +
    geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', alpha=0.5) + 
    geom_point(size = point.size, alpha=transparency) + 
    scale_color_manual(values = colorIN,
                       labels = DEobsrvd) +
    #scale_color_manual(values = c(downColor, unchangedColor, upColor), # to set the colours of our variable  
    #                   labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = c(0, ceiling(max(y.axis.scale))), xlim = c(x_low_lim, x_Up_lim)) + # since some genes can have minuslog10padj of inf, we set these limits
    #coord_cartesian(ylim = c(0, 80), xlim = c(-10,10)) +
    labs(x = x.label, y = y.label) + 
    scale_x_continuous(breaks = x.axis.scale) + # to customise the breaks in the x axis
    ggtitle(plot.title) + # Plot title 
    geom_text_repel(max.overlaps = Inf,size=label_size) + # To show all labels 
    theme_classic() +
    theme(text = element_text(size=rel(axis_text_size)), #change fontsize of axis ticks and label
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=legend_title_size), #change legend title font size
          legend.text = element_text(size=legend_text_size),
          title =element_text(size=plot_title_size, face=plot_title_face), #change legend text font size)
          axis.title = element_text(family = "Helvetica", size = xy_axis_labels, colour = "black"))
}
#)

#df<-read.csv("../Sanja/DE_results_01/CLuster_2_nicdFrac_vs_cntFrac.csv")

# #volcanoPlot(df,
#             labelGenes=TRUE,
#             geneName.col="geneID",
#             x.col="avg_log2FC", 
#             y.col="p_val_adj", 
#             plot.title="VolcanoPlot : MSC2",
#             x.label="log2FoldChange",
#             y.label="-log10(padj)",
#             x.Upper.limit=2.5,
#             x.lower.limit=-1.5,
#             x.break=2,
#             y.lower.limit=0,
#             y.Upper.limit=100,
#             y.break=50,
#             point.size=2,
#             transparency=0.5,
#             legend_text_size=12,
#             label_size=7,
#             xy_axis_labels=18,
#             genes2label=5)







