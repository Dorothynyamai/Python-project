#script for differential analysis of mouse data using DESeq2
#setwd("/Users/dorothynyamai/Documents/Prof_Tao/Ctnnb1_analysis/TPM") 

# load libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(org.Mm.eg.db)


# read in counts data
counts <- read.csv('gene_tpm.csv', header = TRUE, row.names = 1)

#filter out counts with rowsums less than 10

#counts <-counts[which(rowSums(counts)>10),]


#define our groups
condition <- factor(c("ctnnb","ctnnb","ctnnb","Control","Control","Control"))

#makedataframe for the conditions and samples
coldata <- data.frame(row.names=colnames(counts), condition)

##convert RPM values to integers so that it doesn't throw an error when creating the dds object
##removes decimal points
counts <- round(counts)

#make an deseq object

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)



#set factor level to indicate reference level
#dds$condition <- relevel(dds$condition, ref="Control")

dds <- DESeq(dds)


#create a new variable vsdata to perform VST of the data
#replace `vst(dds, blind = FALSE)` with `varianceStabilizingTransformation(dds)` in the R code if there not enough 
#samples in the dataset to perform a variance stabilizing transformation (VST) using the DESeq2 package 

vsdata <- varianceStabilizingTransformation(dds, blind = FALSE)

#do QC of the data

#generate a PCA plot

png("PCAplot_ctnnbTPM.png", res = 300, width = 1500, height =1500)
plotPCA(vsdata, intgroup = "condition")

dev.off()


#look at dispersion of the data .shows variability between replicates as a function of normalized counts

png("dispersion_ctnnbTPM.png", res = 300, width = 1500, height =1500)
plotDispEsts(dds)

dev.off()

#create a table of differentially expressed genes
#Use contrast to indicate which conditions you are comparing

res <- results(dds, contrast = c("condition", "ctnnb", "Control"))

res

res <- na.omit(res)

write.csv(res, file = 'ctnnb_TPMdeseq_results.csv')

#padjusted allows for correction for multiple tests

#resb <- na.omit(res)

#resb.df <-  as.data.frame(res)
#create a volcano plot

#create a dataframe of res

res.df <- as.data.frame(res)
#map ensembl ids to symbol ids

res.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res.df), keytype = "ENSEMBL", column = "SYMBOL")

#resb.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(resb.df), keytype = "ENSEMBL", column = "SYMBOL")

#resb.df
res.df

#use Enhnancedvolcano to generate volcano plot
png("volcano_ctnnbTPM.png", res = 300, width = 2000, height =3500)
#EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol)
EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol)
#pCutoff = 1e-4, FCcutoff = 1)
dev.off()

####################################
#make a list to label a few select genes
selected = c("ND5", "ND4L", "ND6", "CTYB", "COX2", "Polr2a")


#png("volcano_cul3.png", res = 300, width = 2000, height =1500)

pdf("volcano_cul3b.pdf", width = 8, height = 10)
#EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol)
EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 1e-4, FCcutoff = 1, selectLab = selected)
dev.off()
###############################


#filter genes with NA entries

sigs <- na.omit(res)

sigs
#filter genes with a p adjusted value less than 0.05
#very high adjusted P values >0.99

sigs <- sigs[sigs$padj<0.05,]

#write the sig genes to a csv file 

write.csv(sigs, file = "ctnnbTPM_sigs_genes.csv")



#create a data frame of significant genes 

sigs.df<- as.data.frame(sigs)

#convert Ensembl to symbol and save in a new coulmn

sigs.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")


#write sigs.df with the the mapped symbol ids to a file
write.csv(sigs.df, file = "ctnnb_ensembl_to_symbol.csv")

#sym_id <- read.csv('ensembl_to_symbol.csv', header = TRUE, row.names = 1)

sigs.df <- na.omit(sigs.df)
sigs.df

write.csv(sigs.df, file = "ctnnb_sigs_genes_symbol.csv")

########################################
#you can also use a dictionary to map symbol IDS

ensembl_map <- read.csv('ensembl_key_mapper.csv', header = FALSE)

keys <- ensembl_map$V1
values <- ensembl_map$V2

l <- list()
for(i in 1:length(keys)){
  l[keys[i]] <- values[i]
}

#for non-mapped labels #in this case have already removed them

no_values <- setdiff(rownames(sigs.df), keys)
for (i in 1:length(no_values)){
  l[no_values[i]] <- "NA"
}


#add a symbol column to the dataframe 
sigs.df <- unlist(l[rownames(df)], use.names= FALSE)
######################################################################

#select the top differentially expressed genes to reduce the genes 

df.top <- sigs.df[(sigs.df$padj < 0.05) & (abs(sigs.df$log2FoldChange)>0.5),] 

dim(df.top)
#order them in the decreasing order based on log2foldchange

df.top <- df.top[order(df.top$log2FoldChange,decreasing = TRUE),]

df.top
#get normalized count data from dds object using rlog

rlog_out <- rlog(dds, blind = FALSE)

#create a matrix of significant genes
mat <- assay(rlog_out)[rownames(df.top), rownames(coldata)]

colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)


#get center and scale each column (z-score) then transpose
#get the Z score of each row, scale it and transpose it because the scale applies to the row 
mat.scaled <- t(apply(mat, 1, scale))

#rename the columns by getting the name from the mat object
colnames(mat.scaled) <- colnames(mat)


#because we ordered the miRNA based on the log2foldchange we are going to keep the top 25 and bottom 25 miRNAs
num_keep <- 30

#1 to num-keep (top) and len-num_keep to length of the dataframe (bottom)
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))

rows_keep


##get the values for the two extra columns in the data frame
#get log2 values for each gene we are keeping and set it as a matrix

l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange)
colnames(l2_val) <-"logFC"

#get the mean for each miRNA we are keeping 
mean <- as.matrix(df.top[rows_keep,]$baseMean)
colnames(mean) <- "AveExpr"

library(RColorBrewer)
library(circlize)


#make color map that maps values between blue, white and red for min and max l2 values

col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red"))

#maps between 0% quantile and 75 quantile of mean values ---- 0, 25, 50, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

#make the heatmap



ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill=2),
                                               height = unit(2, "cm")))

##do not cluster genes but you can cluster samples
h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F,
              column_labels = colnames(mat.scaled), name = "Z-score",
              cluster_columns = T)

h2 <- Heatmap(l2_val, row_labels = df.top$symbol[rows_keep],
              cluster_rows = F, name = "logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid based on l2_val
                grid.text(round(l2_val[i, j],2), x, y)
              })

h3 <- Heatmap(mean, row_labels = df.top$symbol[rows_keep],
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid based on mean
                grid.text(round(mean[i, j],2), x, y)
              })

#h <- h1+h2+h3
h <- h1+h2
png("geneTPMheatmap.png",  res = 300, width = 3500, height =4500)
print(h)
dev.off()


#################################
#simple heatmap
#sigtop.df <- sigs.df[(sigs.df$padj < 0.05) & (sigs.df$baseMean>2500),] 

sigtop.df <- sigs.df[(sigs.df$padj < 0.05) & (abs(sigs.df$log2FoldChange)>1),] 

dim(sigtop.df)

write.csv(sigtop.df, file = "top_sig_genes.csv")

mats <- counts(dds, normalized = T)[rownames(sigtop.df),]

mats.z <- t(apply(mats, 1, scale))
colnames(mats.z) <- rownames(coldata)
#mats.z

#make the heatmap
png("simple_heatmap.png",  res = 300, width = 3500, height =4500)
Heatmap(mats.z, cluster_rows = F, cluster_columns = T, column_labels = colnames(mats.z),
        name="z-score", row_labels = sigtop.df[rownames(mats.z),]$symbol)
dev.off()

#############################################


#KEGG analysis

#load libraries
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(dplyr)
library(AnnotationDbi)



##convert Ensembl to Entrez ID and save in a new column

kegg.df <- sigs.df

kegg.df$ENTREZID <- mapIds(org.Mm.eg.db, keys = rownames(kegg.df), keytype = "ENSEMBL", column = "ENTREZID")

kegg.df

topkegg.df <- kegg.df[(kegg.df$padj < 0.05) & (abs(kegg.df$log2FoldChange)>0.5),] 

#write the genes to use for KEGG and GO analysis to a file

write.csv(topkegg.df, file = "kegg_Go_genes_symbol.csv")

entrez_genes <- as.vector(topkegg.df$ENTREZID)
entrez_genes <- as.data.frame(topkegg.df$ENTREZID)



enrichment <- function(x,y){
  plot = enrichGO(
    x,
    org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = y,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    readable = FALSE,
    pool = FALSE
  )
  dotplot(plot)
}

data("geneList")


gene = names(geneList)[geneList>1]

#pdf("ctnnb_MF.pdf", width = 6, height = 2)
png("ctnnb_MF.png", res = 300, width =2500, height = 1800)
enrichment(entrez_genes, "MF")
dev.off()

png("ctnnbTPM_CC.png", res = 300, width =2200, height = 1500)
enrichment(entrez_genes, "CC")
dev.off()

png("ctnnbTPM_BP.png", res = 300, width =2200, height = 1800)
enrichment(entrez_genes, "BP")
dev.off()



#convert entrez ids to kegg ids

library(KEGGREST)

#entrez_ids <- as.vector(topkegg.df$ENTREZID)

#entrez_ids <- na.omit(entrez_ids)
##Get the Entrez gene IDs associated with those symbols

#EG_IDs <- as.character(topkegg.df$ENTREZID)

#EG_IDs = mget(sym, revmap(org.Mm.egSYMBOL),ifnotfound=NA)

##Then get the KEGG IDs associated with those entrez genes.
#KEGG_IDs <- mget(as.character(EG_IDs), org.Mm.egPATH,ifnotfound=NA)

#KEGG_IDs <- na.omit(KEGG_IDs)
#kegg_ids <- keggConv("ncbi-geneid:entrez", "kegg", entrez_ids$entrez_ids)

#write.table(KEGG_IDs, "entrez_keggids.csv", row.names = TRUE)

#print(KEGG_IDs)


#use select function from the "org.Mm.eg.db" to retrieve gene annotations for the entrez IDs

Mm <- org.Mm.eg.db
gene_anno <- AnnotationDbi::select(org.Mm.eg.db, keys = kegg.df$ENTREZID, keytype = "ENTREZID", columns= c("ENTREZID", "PATH"))

gene_annoGO <- AnnotationDbi::select(org.Mm.eg.db, keys = kegg.df$ENTREZID, keytype = "ENTREZID", columns= c("ENTREZID", "GO"))
gene_annoSY <- AnnotationDbi::select(org.Mm.eg.db, keys = kegg.df$ENTREZID, keytype = "ENTREZID", columns= c("ENTREZID", "SYMBOL"))

#perform KEGG pathway enrichment analysis using enrichKEGG

#kegg_enrichm <- enrichKEGG(gene = gene_anno$SYMBOL, )

#kegg_pw <- enrichKEGG(gene = gene_anno$PATH, organism = 'mmu', keyType = "PATH", pvalueCutoff = 1)

genelist <- as.vector(gene_anno$PATH)

genelist <-na.omit(genelist)

png("Ctnnb_KEGG.png", res = 300, width =2500, height = 2000)
Kegg <- enrichKEGG(gene_anno$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
Kegg@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", Kegg@result$Description, fixed = T)
dotplot(Kegg, showCategory = 20)
dev.off()

#png("ctnnb_KEGG_boxplot.png", res = 300, width =2000, height = 1500)
#Kegg <- enrichKEGG(as.vector(gene_anno$ENTREZID), pvalueCutoff = 0.05)
#boxplot(Kegg)
#plot(barplot(Kegg, showCategory = 20))
#dev.off()
