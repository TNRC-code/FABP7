  # Install seurat version 4.4.0
  # remotes::install_version('Seurat',version='4.4.0',lib='C:/Users/Devin/Desktop/FABP7/lib')
  # NOTE: DO NOT USE A DIFFERENT VERSION OF SEURAT (e.g. Seurat v5)
    
  library(Seurat,lib.loc = '../lib/')
  library(ggplot2)
  library(RColorBrewer)
  library(pals)
  library(patchwork)
  library(ggpubr)
  library(clustree)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  

  ### ---------------------------------------------------------------------- ###
  ### FIGURE 6                                                               ###
  ### ---------------------------------------------------------------------- ###
  
  ### Astrocytes
  ast <- readRDS('ast_final2_040323.RDS')
  plot_order <- c('Control','Periplaque','Lesion-core','Rim')
  
  ### Astrocyte overview (Fig. 6d)
  clus_id <- levels(ast)
  set.seed(3)
  cols <- setNames(sample(stepped(n=length(clus_id))),clus_id)
  Idents(ast) <- 'seurat_clusters'
  DefaultAssay(ast) <- 'SCT'
  DimPlot(ast,label=TRUE) + NoLegend() + coord_fixed() +
    scale_color_manual(values=cols) +
    annotate('text',x=Inf,y=Inf,label=paste(ncol(ast),'cells'),vjust=1,hjust=1) +
    annotate('text',x=Inf,y=-Inf,label='5 MS Patients',vjust=-2.5,hjust=1) +
    annotate('text',x=Inf,y=-Inf,label='3 Controls',vjust=-1,hjust=1) +
    ggtitle(label='Astrocytes') +
    theme(text = element_text(size = 12)) +
    theme(plot.title=element_text(size=14,face='plain',hjust=0.5)) 
  
  ### FABP7 expression (Fig. 6e)
  # Statistics
  Idents(ast) <- 'Lesion_Compartment'
  ast <- PrepSCTFindMarkers(ast,assay='SCT') # ensure correct data 
  stats1 <- sapply(plot_order,function(i) {
    mrk <- FindMarkers(ast,features='FABP7',ident.1=i,ident.2='Control',min.pct=0,logfc.threshold=0)
    p <- mrk$p_val_adj
    p <- signif(p,digits=2)
    if(p > 0.05) {
      p <- 'n.s.'
    } else {
      if(p < 0.000001) {
        p <- 'p < 0.000001'
      } else {
        p <- paste('p =',p)
      }
    }
    paste('vs Control:',p)
  },USE.NAMES=TRUE,simplify=FALSE)
  stats2 <- sapply(plot_order,function(i) {
    mrk <- FindMarkers(ast,features='FABP7',ident.1=i,ident.2='Periplaque',min.pct=0,logfc.threshold=0)
    p <- mrk$p_val_adj
    p <- signif(p,digits=2)
    if(p > 0.05) {
      p <- 'n.s.'
    } else {
      if(p < 0.000001) {
        #p <- 'p < 0.000001'
        p <- paste('p = ',p)
      } else {
        p <- paste('p = ',p)
      }
    }
    paste('vs Periplaque:',p)
  },USE.NAMES=TRUE,simplify=FALSE)
  stats1$Control <- ''
  stats2$Periplaque <- ''
  # UMAP + expression
  DefaultAssay(ast) <- 'SCT'
  reduc <- Reductions(ast,'umap')
  rangeX <- range(reduc[[]][,1])
  rangeY <- range(reduc[[]][,2])
  expr_upper_lim <- 0.9*max(FetchData(ast,'FABP7')$FABP7)
  l <- SplitObject(ast,split.by='Lesion_Compartment')
  l <- l[plot_order]
  p <- lapply(seq_along(l),function(i) {
    FeaturePlot(l[[i]],features='FABP7',order=TRUE,pt.size=0.5) +
      ggtitle(label=names(l)[i],subtitle=paste0(stats1[[i]],'\n',stats2[[i]])) +
      theme(plot.title=element_text(size=12,face='plain'),
            plot.subtitle=element_text(size=8)) + 
      xlim(rangeX) +
      ylim(rangeY) +
      NoAxes() + 
      NoLegend() +
      scale_color_distiller(palette='Blues',direction=1,limits=c(0,expr_upper_lim))
  })
  wrap_plots(p,ncol=2)
  ### get color legend
  legend <- as_ggplot(get_legend(p[[1]] +
                                   RestoreLegend() +
                                   labs(color='FABP7 Expression    ') +
                                   theme(
                                     legend.text=element_text(size=10),
                                     legend.title=element_text(size=12,face='plain')
                                     ),
                                 position='top')
                      )
  legend
