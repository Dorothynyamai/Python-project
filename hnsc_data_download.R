
BiocManager::install("maftools") 
BiocManager::install("TCGAbiolinks") 
BiocManager::install("SummarizedExperiement") 
BiocManager::install("GDCquery")

BiocManager::install("DESeq2")

library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)


# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-HNSC')

# building a query
query_TCGA <- GDCquery(project = 'TCGA-HNSC',
                       data.category = 'Transcriptome Profiling',
                       data.type = 'miRNA Expression Quantification',
                       experimental.strategy = 'miRNA-Seq',
                       workflow.type = 'BCGSC miRNA Profiling',
                       access = 'open')


output_query_TCGA <- getResults(query_TCGA)

#query_TCGA <- GDCquery(project = 'TCGA-HNSC',
#                       data.category = 'Transcriptome Profiling',
 #                      data.type = 'miRNA Expression Quantification',
 #                      experimental.strategy = 'miRNA-Seq',
  #                     workflow.type = 'BCGSC miRNA Profiling',
  #                     access = 'open')
                       #barcode = c('TCGA-LL-A73Y-01A-11R-A33J-07', 'TCGA-E2-A1IU-01A-11R-A14D-07','TCGA-AO-A03U-01B-21R-A10J-07'))


#getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)

# prepare data
tcga_hnsc_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)

#tcga_hnsc_data  <- na.omit(tcga_hnsc_data)
library(Biobase)
myexpr  <- ExpressionSet(assayData = list(counts = as.matrix(tcga_hnsc_data )))

assay(myexpr)
#library(Biobase) 
#tcga_hnsc_datab  <- read.csv(tcga_hnsc_data)

#write.csv(tcga_hnsc_data, "hnsc_data.csv")
hnsc_matrix <- assay(tcga_hnsc_data, 'read_count')






