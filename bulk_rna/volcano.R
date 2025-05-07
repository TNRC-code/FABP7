## Load Libraries
library(EnhancedVolcano)
library(tidyverse)
library(ggplot2)

## --- Read Data 
new_data <- read_csv('compare_control_vs_FABP.csv')

## -- Set volcano plot Object 
volcano_plot <- EnhancedVolcano(new_data,
                lab = new_data$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-12,
                FCcutoff = 1.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                title = 'Differential Expression Analysis',
                subtitle = 'Control vs FABP7',
                pointSize = 2.1,
                labSize = 2.5,
                colAlpha = 0.5) + coord_flip()

## --Save file
png('volcano_plot_control_vs_fabp7.png', height = 8, width=10,units="in",res=500)
print(volcano_plot)
dev.off()

