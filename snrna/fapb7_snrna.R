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
  library(plyr)
  

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



  ### Immune cells
  imm <- readRDS('imm_final_032823.RDS')



  ### Immune cell overview (Fig. 6f)
  # Determine optimal clustering resolution
  ElbowPlot(imm,ndims=30,reduction='pca') # choose 13
  DefaultAssay(imm) <- 'integrated'
  imm <- RunUMAP(imm,reduction='pca',dims=1:13)
  imm <- FindNeighbors(imm,reduction='pca',dims=1:13)
  imm <- FindClusters(imm,resolution=seq(from=0.2,to=0.8,by=0.2))
  DefaultAssay(imm) <- 'SCT'

  # Cluster colors and order
  imm$Lesion_Compartment <- factor(imm$Lesion_Compartment,levels=plot_order)
  cluster_order <- c('8','11','9','2','10','5','3','4','7','0','1','6')
  colors <- setNames(brewer.pal(length(cluster_order),name='Set3'),cluster_order)
  imm$seurat_clusters <- factor(imm$seurat_clusters,levels=cluster_order)
  DimPlot(imm,label=T) + coord_fixed() + NoLegend() +
    scale_color_manual(values=colors) +
    annotate('text',x=Inf,y=Inf,label=paste(ncol(imm),'cells'),vjust=10,hjust=1) +
    ggtitle(label='Immune Cells') +
    theme(text = element_text(size = 12)) +
    theme(plot.title=element_text(size=14,face='plain',hjust=0.5))
  double_pos <- WhichCells(imm,expression = CD14 > 0 & FCGR3A > 0)
  DimPlot(imm,label=F,cells.highlight=double_pos,split.by='Lesion_Compartment',cols.highlight='dodgerblue',pt.size=0.2,sizes.highlight=0.4) + coord_fixed() + NoLegend() 



  ### Cell Types (Fig. 6g)
  T_cell <- c('CD2','CD3E','CD4','CD8A','IL7R')
  B_cell <- c('IGHM','IGHG1','IGHA1','IL2RA','CD79A','CD38','SLAMF7') # https://www.abcam.com/primary-antibodies/b-cells-basic-immunophenotyping
  Mono <- c('CD14','FCGR3A','CD36','CD68','CD80','CD83','CD86','CD163','MRC1','CD209','SPI1','MYB')
  homeo <- c('PTPRC','ITGAM','C1QA','P2RY12','CSF1R','FOS','GPR34','MERTK','PROS1','TGFBR1')
  activ <- c('IL1B','IL1RAP','IL17RA','AIF1','APOE','AXL','SPP1','ITGAX','TREM2','CHIT1','GPNMB')
  Iron <- c('FTL','FTH1')
  # CD206 (MRC1) is also highly expressed by the perivascular macrophages that reside at the CNS border, enabling them to be discriminated from microglia in the CNS parenchyma.   
  imm_genes <- c(T_cell,B_cell,Mono,homeo,activ,Iron)
  imm_idents <- c('8' ='8-T_cells',
                  '11'='11-B_cells',
                  '0' ='0-Homeostatic_microglia',
                  '1' ='1-Homeostatic_microglia',
                  '3' ='3-Homeostatic_microglia',
                  '7' ='7-CD11c_microglia',
                  '9' ='9-Infiltrating/perivascular_macrophage',
                  '6' ='6-Early_activated_microglia',
                  '5' ='5-Activated_microglia_IL1B',
                  '4' ='4-Activated_microglia_CD83',
                  '2' ='2-Activated_microglia_MIMS',
                  '10'='10-Activated_microglia_MIMS-iron')
  imm$seurat_clusters <- revalue(imm$seurat_clusters,imm_idents)
  Idents(imm) <- 'seurat_clusters'

  DotPlot(imm,assay='SCT',features=imm_genes,scale.min = 0,scale.max=40,dot.scale=4,group.by='seurat_clusters',cluster.idents=TRUE) +
    coord_flip() +
    geom_point(aes(size=pct.exp),shape=21,color='gray30',stroke=0.5) +
    scale_color_viridis_b(option='B') + # try B or D
    guides(size=guide_legend(title = 'Perc.\nexpr.',override.aes=list(shape=21,color='gray30',fill='white')),
           colour=guide_colorbar(title='Avg.\nexpr.')) +
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=9),
          axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
          legend.text=element_text(size=9),
          legend.title=element_text(size=10)) + # Export 300x500 (or 300x800)
    coord_fixed()

  Idents(imm) <- 'seurat_clusters'
  imm <- PrepSCTFindMarkers(imm,assay='SCT')
  # Clustree
  DefaultAssay(imm) <- 'integrated'
  clustree(imm)
  clustree_overlay(imm,'pca1','pca2',red_dim='pca',assay='integrated')
  clustree_overlay(imm,'umap1','umap2',red_dim='umap',assay='integrated')



  ### Relative proportion of immune cell types (Fig. 6h)
  # first, get cd14 cd16 double pos cells by cluster
  DefaultAssay(imm) <- 'SCT'
  ic <- SplitObject(imm,split.by='seurat_clusters')
  perc_doublepos <- sapply(ic,function(i) {
    expr <- FetchData(i,vars=c('CD14','FCGR3A'))
    (sum(apply(expr>0,1,all)) / nrow(expr))*100
  })
  mrk_expr <- sapply(ic,function(i) {
    colSums(FetchData(i,vars=c('CD14','FCGR3A','CD163','MRC1'))) / nrow(i)
  })
  mrk_expr <- as.data.frame(t(mrk_expr))
  perc_doublepos <- as.data.frame(perc_doublepos)
  colnames(perc_doublepos) <- 'CD14/CD16'
  dat <- SplitObject(imm,split.by='Lesion_Compartment')
  dat <- sapply(dat,function(i) {
    tab <- table(i$seurat_clusters) + 1 #pseudocount
    tab / sum(tab)
  })
  dat <- as.data.frame(dat)
  pc <- apply(dat,2,function(col) {
    old <- dat[,'Control']
    ((col-old)/old)*100
  })  
  COL <- rev(brewer.pal(9,name='RdBu'))
  cluster_order2 <- c('9-Infiltrating/perivascular_macrophage','5-Activated_microglia_IL1B','2-Activated_microglia_MIMS','10-Activated_microglia_MIMS-iron',
                      '6-Early_activated_microglia','0-Homeostatic_microglia','1-Homeostatic_microglia','3-Homeostatic_microglia',
                      '4-Activated_microglia_CD83','7-CD11c_microglia','8-T_cells','11-B_cells')
  pheatmap(pc[cluster_order2,],annotation_row=perc_doublepos,cluster_rows=FALSE,annotation_legend=FALSE,breaks=seq(from=-150,to=150,length.out=length(COL)+1),color=COL,main='Immune cell\ncomposition')



  ### Tissue microenvironments (Fig. 6i)
  DefaultAssay(ast) <- 'SCT'
  ae <- AverageExpression(ast,features='FABP7',group.by='GEO_ID',assays='SCT',slot='data')$SCT
  ae <- t(ae)
  colnames(ae) <- 'FABP7'
  ae <- as.data.frame(ae)
  ig <- SplitObject(imm,split.by='GEO_ID')
  ig <- sapply(ig,function(i) {
    tab <- table(i$seurat_clusters) # + 1 #pseudocount
    tab / sum(tab)
  })
  ig <- as.data.frame(t(ig))
  shared <- intersect(row.names(ae),row.names(ig))
  ae <- ae[row.names(ae) %in% shared,,drop=FALSE]
  ig <- ig[row.names(ig) %in% shared,,drop=FALSE]
  mrg <- merge(ae,ig,by=0)
  mrg$Row.names <- NULL
  cm <- cor(mrg,method='spearman')
  cv <- cm[,'FABP7']
  cv <- cv[order(cv,decreasing=TRUE)]
  pheatmap(t(cv[-1]),cluster_rows=FALSE,cluster_cols=FALSE,breaks=seq(from=-0.5,to=0.5,length.out=10),color=rev(brewer.pal(9,'RdYlBu')),main='Astrocyte FABP7 expression\nassociation with immune cell subsets')



  ### Association between FABP7 and cell types (Fig. 6j)
  mrg$MIMS <- mrg$`2-Activated_microglia_MIMS` + mrg$`10-Activated_microglia_MIMS-iron`
  mrg$Homeostatic_microglia <- mrg$`0-Homeostatic_microglia` + mrg$`1-Homeostatic_microglia` + mrg$`3-Homeostatic_microglia`
  mrg$Activated_microglia <- mrg$`2-Activated_microglia_MIMS` + mrg$`4-Activated_microglia_CD83` + mrg$`5-Activated_microglia_IL1B` + mrg$`10-Activated_microglia_MIMS-iron`
  mrg$Inflamed_microglia <- mrg$`2-Activated_microglia_MIMS` + mrg$`10-Activated_microglia_MIMS-iron`
  
  p1_stat <- cor.test(mrg$FABP7,mrg$Homeostatic_microglia,method='spearman',exact = FALSE)
  p1 <- ggplot(mrg,aes(x=FABP7,y=Homeostatic_microglia)) + geom_point() + geom_smooth(method='lm') + theme_bw() + ylab('') + xlab('Astrocyte\nFABP7 expr.') +
    ggtitle('Homeostatic\nmicroglia',subtitle=paste('rho =',signif(p1_stat$estimate,2),'\np =',signif(p1_stat$p.value,2)))
  
  p2_stat <- cor.test(mrg$FABP7,mrg$Inflamed_microglia,method='spearman',exact = FALSE)
  p2 <- ggplot(mrg,aes(x=FABP7,y=Inflamed_microglia)) + geom_point() + geom_smooth(method='lm') + theme_bw() + ylab('') + xlab('Astrocyte\nFABP7 expr.') +
    ggtitle('Inflamed\nmicroglia',subtitle=paste('rho =',signif(p2_stat$estimate,2),'\np =',signif(p2_stat$p.value,2)))
  
  p1+p2

