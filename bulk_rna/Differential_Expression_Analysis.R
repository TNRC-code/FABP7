## Load Libraries
library(tidyverse) ## Data Manipulation
library(janitor) ## Cleaning column names
library(org.Hs.eg.db) ## Annotations
library(DESeq2) ## Differential Gene Expression 


## -- Phenotype loading --- 
pheno <- read.table('lps_fabp_sample_info.txt', header=T) 

## -- Count loading ---
count_data <- read.csv('count_matrix.csv', header=T,row.names = 1)

## -- Setting colnames --
colnames(count_data) <- gsub(".bam",'', colnames(count_data))
colnames(count_data) <- gsub('-','_', colnames(count_data))


## -- Updating pheno file to match count data and filtering --

pheno_fabp <- pheno %>%
  dplyr::mutate(Barcode = gsub('-','_',Barcode)) %>%
  dplyr::mutate(Group = as.factor(Group)) %>%
  dplyr::mutate(Sex = as.factor(Sex)) %>%
  dplyr::filter(Group %in% c('Control','FABP7')) %>%
  dplyr::mutate(Subject_id = as.factor(Subject_id)) %>%
  droplevels()


## -- selecting columns matching control and FABP7 --

count_fabp <- count_data %>%
  dplyr::select(pheno_fabp$Barcode) 

#all(colnames(count_data) == pheno$Barcode)


## -- Setting up DESeq2 Object --
dds <- DESeqDataSetFromMatrix(countData = count_fabp, 
                              colData = pheno_fabp,
                              design = ~0 + Group)

# --- Filter out low-count genes ---
# Keep genes that have at least 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## -- Run DESeq2 -- 

dds <- DESeq(dds, parallel = TRUE)

## -- Processing results and setting contrast to capture control vs fabp7
results_control_vs_fabp7 <- results(dds, contrast = c('Group','Control','FABP7')) %>%
  data.frame() %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(pvalue)


## -- Sex Specific Analysis --

## Male Control vs Male FABP7 

dds_men <- subset(dds, select = colData(dds)$Sex == 'Male') ## Subset the object ## 

dds_men <- DESeq(dds_men, parallel = TRUE) ## Run DESeq2 

## Extract results ## 
male_compare_control_vs_FABP <- results(dds_men,contrast = c('Group','Control','FABP7')) %>%
  data.frame() %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj)



## Female Control vs Female FABP7 

dds_female <- subset(dds, select = colData(dds)$Sex == 'Female') ## Subset the object ## 
dds_female$Sex <- droplevels(dds_female$Sex)

dds_female <- DESeq(dds_female, parallel = TRUE) ## Run DESeq2 


female_compare_control_vs_FABP <- results(dds_female,contrast = c('Group','Control','FABP7')) %>%
  data.frame() %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj)


