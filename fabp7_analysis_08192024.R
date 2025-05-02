
  setwd('I:/Work (E)/FABP7/fabp7/')
    
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
  
  # Install seurat version 4.4.0
  # remotes::install_version('Seurat',version='4.4.0',lib='C:/Users/Devin/Desktop/FABP7/lib')
  
  
  # Astrocytes
  ast <- readRDS('ast_final2_040323.RDS')
  plot_order <- c('Control','Periplaque','Lesion-core','Rim')
  
  ### Astrocyte overview# Export 300x300
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
  ### FABP7 expression # Export 300x300 (or 300x250)
  # # Statistics
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
  wrap_plots(p,ncol=2) # Export 300x400
  ### get color legend # Export 300x100
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
  # Immune cell overview
  imm <- readRDS('imm_final_032823.RDS')
  
  # Determine optimal clustering resolution using scSHC
  #library(scSHC)
  tmp <- imm
  DefaultAssay(tmp) <- 'RNA'
  # x <- GetAssayData(tmp,assay='RNA',layer='counts')
  # scs <- scSHC(x,num_features = 2000)
  # clus.scs <- scs[[1]]
  # 
  # 
  # # 
  # data(p3cl)
  # seu <- p3cl
  # multik <- MultiK(seu, reps=100)
  # DiagMultiKPlot(multik$k, multik$consensus)
  # 
  # 
  # imm.norm <- imm
  # DefaultAssay(imm.norm) <- 'RNA'
  # mk <- MultiK(imm.norm,resolution = seq(0.1,2,0.2),nPC=30,reps=200,pSample=0.2,seed=17)
  # 
  # DiagMultiKPlot(mk$k,mk$consensus)
  # 
  # clusters <- getClusters(imm.norm, 14)
  # 
  # pval <- CalcSigClust(imm.norm, clusters$clusters)
  # 
  # PlotSigClust(imm.norm, clusters$clusters, pval)
  # 
  # tog <- as.data.frame(table(mk$k)[table(mk$k) > 1])
  # colnames(tog)[1] <- "ks"
  # pacobj <- CalcPAC(x1=0.1, x2=0.9, xvec = tog$ks, ml = multik$consensus)
  # tog$rpac <- pacobj$rPAC
  # tog$one_minus_rpac  <- 1-tog$rpac
  # optK <- findOptK(tog)
  # 
  # 
  # 
  # 
  # 
  # tog.f <- tog[tog$Freq > 100 | tog$Freq == 100, ]
  # hpts <- chull(tog.f[, c("one_minus_rpac", "Freq")])
  # hpts <- c(hpts, hpts[1])
  # ch.df <- tog.f[hpts, ]
  # df <- ch.df[, c("ks", "one_minus_rpac", "Freq")]
  # colnames(df) <- c("k", "x", "y")
  # b <- c()
  # end_points <- c()
  # for (i in 1:(nrow(df) - 1)) {
  #   end_points[i] <- paste(as.character(df[i, ]$k), as.character(df[(i + 
  #                                                                      1), ]$k), sep = "-")
  #   b[i] <- (df[(i + 1), ]$y - df[i, ]$y)/(df[(i + 1), ]$x - 
  #                                            df[i, ]$x)
  # }
  # lineseg.df <- data.frame(end_points = end_points, slope = b)
  # lineseg.df$p1 <- do.call("rbind", strsplit(lineseg.df$end_points, 
  #                                            "-"))[, 1]
  # lineseg.df$p2 <- do.call("rbind", strsplit(lineseg.df$end_points, 
  #                                            "-"))[, 2]
  # which.k <- as.character(ch.df[which.max(ch.df$Freq), ]$ks)
  # if (all(lineseg.df[lineseg.df$p1 == which.k | lineseg.df$p2 == 
  #                    which.k, ]$slope > 0)) {
  #   optK <- which.k
  # }
  # else {
  #   tmp <- which(lineseg.df[lineseg.df$p1 == which.k | lineseg.df$p2 == 
  #                             which.k, ]$slope < 0)
  #   tmp <- lineseg.df[lineseg.df$p1 == which.k | lineseg.df$p2 == 
  #                       which.k, ][tmp, ]
  #   which.k2 <- as.character(c(tmp$p1, tmp$p2)[which(c(tmp$p1, 
  #                                                      tmp$p2) != which.k)])
  #   lineseg.df.sub <- lineseg.df[lineseg.df$p1 != which.k & 
  #                                  lineseg.df$p2 != which.k, ]
  #   if (lineseg.df.sub[lineseg.df.sub$p1 == which.k2 | lineseg.df.sub$p2 == 
  #                      which.k2, ]$slope > tmp$slope) {
  #     optK <- c(which.k, which.k2)
  #   }
  #   else {
  #     tmp <- which(lineseg.df.sub[lineseg.df.sub$p1 == 
  #                                   which.k2 | lineseg.df.sub$p2 == which.k2, ]$slope < 
  #                    0)
  #     tmp <- lineseg.df.sub[lineseg.df.sub$p1 == which.k2 | 
  #                             lineseg.df.sub$p2 == which.k2, ][tmp, ]
  #     which.k3 <- as.character(c(tmp$p1, tmp$p2)[which(c(tmp$p1, 
  #                                                        tmp$p2) != which.k & c(tmp$p1, tmp$p2) != which.k2)])
  #     lineseg.df.sub <- lineseg.df[lineseg.df$p1 != which.k & 
  #                                    lineseg.df$p2 != which.k & lineseg.df$p1 != which.k2 & 
  #                                    lineseg.df$p2 != which.k2, ]
  #     if (lineseg.df.sub[lineseg.df.sub$p1 == which.k3 | 
  #                        lineseg.df.sub$p2 == which.k3, ]$slope > tmp$slope) {
  #       optK <- c(which.k, which.k2, which.k3)
  #     }
  #     else {
  #       optK <- c(which.k, which.k2, which.k3)
  #     }
  #   }
  # }
  # 
  # 
  # 
  # 
  # 
  # mk.old <- mk
  # 
  
  
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
  ### Overview UMAP
  DimPlot(imm,label=T) + coord_fixed() + NoLegend() +
    scale_color_manual(values=colors) +
    annotate('text',x=Inf,y=Inf,label=paste(ncol(imm),'cells'),vjust=10,hjust=1) +
    ggtitle(label='Immune Cells') +
    theme(text = element_text(size = 12)) +
    theme(plot.title=element_text(size=14,face='plain',hjust=0.5)) 
  # DimPlot(imm,label=T,split.by='MS_Control') + coord_fixed() + NoLegend() + scale_color_manual(values=colors)
  # DimPlot(imm,label=T,split.by='Lesion_Compartment') + coord_fixed() + NoLegend() + scale_color_manual(values=colors)
  double_pos <- WhichCells(imm,expression = CD14 > 0 & FCGR3A > 0)
  DimPlot(imm,label=F,cells.highlight=double_pos,split.by='Lesion_Compartment',cols.highlight='dodgerblue',pt.size=0.2,sizes.highlight=0.4) + coord_fixed() + NoLegend() 

  ### Immune marker dot plot
  
  T_cell <- c('CD2','CD3E','CD4','CD8A','IL7R')
  B_cell <- c('IGHM','IGHG1','IGHA1','IL2RA','CD79A','CD38','SLAMF7') # https://www.abcam.com/primary-antibodies/b-cells-basic-immunophenotyping
  Mono <- c('CD14','FCGR3A','CD36','CD68','CD80','CD83','CD86','CD163','MRC1','CD209','SPI1','MYB')
  homeo <- c('PTPRC','ITGAM','C1QA','P2RY12','CSF1R','FOS','GPR34','MERTK','PROS1','TGFBR1')
  activ <- c('IL1B','IL1RAP','IL17RA','AIF1','APOE','AXL','SPP1','ITGAX','TREM2','CHIT1','GPNMB')
  Iron <- c('FTL','FTH1')
  
  # CD206 (MRC1) is also highly expressed by the perivascular macrophages that reside at the CNS border, enabling them to be discriminated from microglia in the CNS parenchyma. 
  
  imm_genes <- c(T_cell,B_cell,Mono,homeo,activ,Iron)
  
  # https://www.biocompare.com/Editorial-Articles/586779-A-Guide-to-Microglial-Markers/#:~:text=In%20DAM%2C%20homeostatic%20microglial%20markers,in%20other%20neurodegenerative%20contexts%2C%20too.
  
  
  # imm_genes <- c('AIF1','LGALS3','PTPRC','P2RY12','SALL1','TMEM119','CD209','MRC1','CD163','CD86',
  #                'CD83','CD80','CD68','CD36','FCGR3A','CD14','IL1B','C1QA','ITGAM',
  #                'CD79A','IL7R','CD8A','CD4','CD3E','NTM','IL1RAPL1',
  #                'LPL','NUPR1','MERTK','NLRP3','TREM2','APOE','FTL','FTH1','S100A9','IGHG1','IGHM')
  # see this reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8557350/


  library(plyr)
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
  
  # Absinta labels
  Idents(imm) <- 'seurat_clusters'
  imm <- PrepSCTFindMarkers(imm,assay='SCT')
  mrk <- FindAllMarkers(imm)
  library(openxlsx)
  abs_markers <- read.xlsx('absinta_supplementary_tables_S1-S13.xlsx',sheet=5,startRow=2,rowNames=TRUE)
  # scenrich_v0.1.R
  tmp <- scenrich(imm,seurat_marker_list=mrk,reference_marker_list=abs_markers) # export 900x200
  # new cluster 8 (abs7)
  # new cluster 11 (abs10)
  # new cluster 9 (unknown)
  # new cluster 2 (abs8 or abs1)
  # new cluster 10 (abs8)
  # new cluster 5 (abs0)
  # new cluster 3 (abs1)
  # new cluster 4 (abs6 or abs2)
  # new cluster 7 (unknown)
  # new cluster 0 (abs0)
  # new cluster 1 (unknown)
  # new cluster 6 (abs0)
  # MIMS-iron (cluster 10, and to a lesser extent 2)
  # MIMS-foamy (cluster 3, and to a lesser extent 2)
  
  ### Clustree
  DefaultAssay(imm) <- 'integrated'
  clustree(imm)
  clustree_overlay(imm,'pca1','pca2',red_dim='pca',assay='integrated')
  clustree_overlay(imm,'umap1','umap2',red_dim='umap',assay='integrated')
  
  ### Relative proportion of immune cell types
  # first, get cd14/cd16 double pos cells by cluster
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
  # 
  dat <- SplitObject(imm,split.by='Lesion_Compartment')
  dat <- sapply(dat,function(i) {
    tab <- table(i$seurat_clusters) + 1 #pseudocount
    tab / sum(tab)
  })
  dat <- as.data.frame(dat)
  # lfc <- apply(dat,2,function(col) {
  #   log2((col)/(dat[,'Control'])) # log2-FC
  # })
  pc <- apply(dat,2,function(col) {
    old <- dat[,'Control']
    ((col-old)/old)*100
  })
  #perc_doublepos.order <- row.names(perc_doublepos[order(perc_doublepos[,1],decreasing=TRUE),,drop=FALSE])
  #pheatmap(lfc[perc_doublepos.order,],annotation_row=perc_doublepos,cluster_rows=FALSE)
  
  COL <- rev(brewer.pal(9,name='RdBu'))
  #pheatmap(pc,annotation_row=mrk_expr,annotation_legend=FALSE,breaks=seq(from=-150,to=150,length.out=length(COL)+1),color=COL)
  cluster_order2 <- c('9-Infiltrating/perivascular_macrophage','5-Activated_microglia_IL1B','2-Activated_microglia_MIMS','10-Activated_microglia_MIMS-iron',
                      '6-Early_activated_microglia','0-Homeostatic_microglia','1-Homeostatic_microglia','3-Homeostatic_microglia',
                      '4-Activated_microglia_CD83','7-CD11c_microglia','8-T_cells','11-B_cells')
  pheatmap(pc[cluster_order2,],annotation_row=perc_doublepos,cluster_rows=FALSE,annotation_legend=FALSE,breaks=seq(from=-150,to=150,length.out=length(COL)+1),color=COL,main='Immune cell\ncomposition') # export 400x350
  #pheatmap(pc,annotation_row=perc_doublepos,cluster_rows=TRUE,annotation_legend=FALSE,breaks=seq(from=-150,to=150,length.out=length(COL)+1),color=COL,main='Immune cell abundance (percent change)')
  
  
  #pheatmap(lfc,annotation_row=mrk_expr,annotation_legend=FALSE,breaks=seq(from=-2,to=2,length.out=length(COL)+1),color=COL)
  # pheatmap(lfc,annotation_row=perc_doublepos,annotation_legend=FALSE,breaks=seq(from=-2,to=2,length.out=length(COL)+1),color=COL)
  
  
  
  # ggplot(dat2,aes(x=Lesion_Compartment,y=Cluster,color=double_pos,size=Relative_Proportion)) + geom_point() + scale_color_distiller(palette='Reds',direction=1) + theme_bw()
  
  
  ### Tissue microenvironments
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
  # Export 360x320
  
  ggplot(mrg,aes(x=FABP7,y=`2-Activated_microglia_MIMS`)) + geom_point() + geom_smooth(method='lm') + theme_bw()
  # Export 300x300
  
  mrg$MIMS <- mrg$`2-Activated_microglia_MIMS` + mrg$`10-Activated_microglia_MIMS-iron`
  mrg$Homeostatic_microglia <- mrg$`0-Homeostatic_microglia` + mrg$`1-Homeostatic_microglia` + mrg$`3-Homeostatic_microglia`
  mrg$Activated_microglia <- mrg$`2-Activated_microglia_MIMS` + mrg$`4-Activated_microglia_CD83` + mrg$`5-Activated_microglia_IL1B` + mrg$`10-Activated_microglia_MIMS-iron`
  mrg$Inflamed_microglia <- mrg$`2-Activated_microglia_MIMS` + mrg$`10-Activated_microglia_MIMS-iron`
  
  # ggplot(mrg,aes(x=FABP7,y=MIMS)) + geom_point() + geom_smooth(method='lm') + theme_bw()
  
  p1_stat <- cor.test(mrg$FABP7,mrg$Homeostatic_microglia,method='spearman',exact = FALSE)
  p1 <- ggplot(mrg,aes(x=FABP7,y=Homeostatic_microglia)) + geom_point() + geom_smooth(method='lm') + theme_bw() + ylab('') + xlab('Astrocyte\nFABP7 expr.') +
    ggtitle('Homeostatic\nmicroglia',subtitle=paste('rho =',signif(p1_stat$estimate,2),'\np =',signif(p1_stat$p.value,2)))
  
  #p2_stat <- cor.test(mrg$FABP7,mrg$Activated_microglia,method='spearman',exact = FALSE)
  #p2 <- ggplot(mrg,aes(x=FABP7,y=Activated_microglia)) + geom_point() + geom_smooth(method='lm') + theme_bw() + ylab('') + xlab('Astrocyte\nFABP7 expr.') +
  #  ggtitle('Activated microglia',subtitle=paste('rho =',signif(p2_stat$estimate,2),'\np =',signif(p2_stat$p.value,2)))
  p2_stat <- cor.test(mrg$FABP7,mrg$Inflamed_microglia,method='spearman',exact = FALSE)
  p2 <- ggplot(mrg,aes(x=FABP7,y=Inflamed_microglia)) + geom_point() + geom_smooth(method='lm') + theme_bw() + ylab('') + xlab('Astrocyte\nFABP7 expr.') +
    ggtitle('Inflamed\nmicroglia',subtitle=paste('rho =',signif(p2_stat$estimate,2),'\np =',signif(p2_stat$p.value,2)))
  
  p1+p2 # Export 300 x 200

  ### Pseudobulk
  library(pheatmap)
  DefaultAssay(imm) <- 'SCT'
  Idents(imm) <- 'seurat_clusters'
  ae <- as.data.frame(AverageExpression(imm,assays='SCT',slot='data',group.by='seurat_clusters')$SCT)
  bulk <- read.csv('compare_control_vs_FABP.csv')
  
  # baseMean
  dat <- merge(bulk[,c('symbol','baseMean','log2FoldChange')],ae,by.x='symbol',by.y='row.names')
  cm <- cor(dat[,c(-1,-3)],method='spearman')[,'baseMean']
  cm <- as.data.frame(cm)[-1,,drop=FALSE]
  pheatmap(cm,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=TRUE)
  
  pheatmap(cm,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=TRUE,breaks=seq(from=0,to=0.5,length.out=length(COL)+1),color=COL)
  
  # log2FC
  cm <- cor(dat[,c(-1,-2)],method='spearman')[,'log2FoldChange']
  cm <- as.data.frame(cm)[-1,,drop=FALSE]
  pheatmap(cm,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=TRUE)
  
  # differential expression between cluster 10 and 0
  DefaultAssay(imm) <- 'SCT'
  Idents(imm) <- 'seurat_clusters'
  mrk <- FindMarkers(imm,ident.1='10',ident.2='0',method='DESeq2',min.pct=0,logfc.threshold=0)
  #dat2 <- merge(bulk[,c('symbol','log2FoldChange')],mrk[,'avg_log2FC',drop=FALSE],by.x='symbol',by.y='row.names')
  dat2 <- merge(bulk[,c('symbol','padj')],mrk[,'p_val_adj',drop=FALSE],by.x='symbol',by.y='row.names')
  
  
  ast <- readRDS('ast_final2_040323.RDS')
  Idents(ast) <- 'MS_Control'
  FindMarkers(ast,features='FABP7',ident.1='MS',ident.2='Control',min.pct=0,logfc.threshold=0)
  Idents(ast) <- 'Lesion_Compartment'
  FindMarkers(ast,features='FABP7',ident.1='Periplaque',ident.2='Control',min.pct=0,logfc.threshold=0)
  FindMarkers(ast,features='FABP7',ident.1='Rim',ident.2='Control',min.pct=0,logfc.threshold=0)
  FindMarkers(ast,features='FABP7',ident.1='Lesion-core',ident.2='Control',min.pct=0,logfc.threshold=0)
  FindMarkers(ast,features='FABP7',ident.1='Rim',ident.2='Periplaque',min.pct=0,logfc.threshold=0)
  FindMarkers(ast,features='FABP7',ident.1='Rim',ident.2='Lesion-core',min.pct=0,logfc.threshold=0)
  
  
  ast <- readRDS('~/Desktop/astros.remerged.reassigned.mvimolg.finalrecluster.allenSingleR.subtypegroups.Absfixed.sobj.RDS')