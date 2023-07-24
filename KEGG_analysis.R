#a script for performing functional enrichment analysis of differentially expresssed miRNA gene targets
# setwd("~/Documents/MRNA_MiRNA/HNSC/GSE49190_WCGA/KEGG_Analysis")

#load libraries
library(DOSE)
library(enrichplot)
library(clusterProfiler)
#library(AnnotationDbi)
library(org.Hs.eg.db)

enrichment <- function(x,y){
  plot = enrichGO(
    x,
    org.Hs.eg.db,
    keyType = "ENSEMBL",
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
eid <- c("ENSG00000167601", "ENSG00000134852", "ENSG00000101266", "ENSG00000100697", "ENSG00000133216", "ENSG00000151491", "ENSG00000068024",
         "ENSG00000017427", "ENSG00000073792", "ENSG00000101384", "ENSG00000134333", "ENSG00000288299", "ENSG00000157227", "ENSG00000057935", "ENSG00000150593",
         "ENSG00000145675", "ENSG00000181690", "ENSG00000069667", "ENSG00000187764", "ENSG00000141646", "ENSG00000072274")
gene = names(geneList)[geneList>1]
#pdf("GSE49190_MF.pdf", width = 6, height = 2)
png("GSE49190_MFb.png", res = 300, width =2000, height = 800)
enrichment(eid, "MF")
dev.off()

png("GSE49190_CCb.png", res = 300, width =2200, height = 1200)
enrichment(eid, "CC")
dev.off()

png("GSE49190_BP.png", res = 300, width =2200, height = 1500)
enrichment(eid, "BP")
dev.off()

Ent_ids <- c("23405", "2048","9575", "2059","5324","6095","27250", "558","3939", "182", "4089", "9759", "10644", "7037", "1457", "10507","4323","57504","5295","3479")

png("GSE49190_KEGGb.png", res = 300, width =2000, height = 1500)
Kegg <- enrichKEGG(Ent_ids, pvalueCutoff = 0.05)
dotplot(Kegg, showCategory = 20)
dev.off()

png("GSE49190_KEGG_boxplot.png", res = 300, width =2000, height = 1500)
Kegg <- enrichKEGG(Ent_ids, pvalueCutoff = 0.05)
#boxplot(Kegg)
plot(barplot(Kegg, showCategory = 20))
dev.off()
