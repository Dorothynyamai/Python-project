
##script for converting gene symbol to ensembl ID
BiocManager::install("biomaRT") 

#library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")

#Read the CSV file containing gene symbols
gene_symbols <- read.csv(file = "GSE49190_gene_list.csv", header = TRUE)
#gene_symbols <- as.data.frame(gene_symbols)

#Create a vector of gene symbols
symbols <- gene_symbols$SYMBOL

#Use the `select` function from `org.Hs.eg.db` to convert gene symbols to Ensembl IDs:
ensembl_ids <- select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", column = "ENSEMBL")

#Merge the original data frame with the Ensembl IDs:

merged_data <- merge(gene_symbols, ensembl_ids, by.x = "SYMBOL", by.y = "SYMBOL")


#Save the updated data frame to a new CSV file

write.csv(merged_data, "gene_symbols_ensembl_ids.csv", row.names = FALSE)

#get KEGG IDS 
#ens_IDS <- c("ENSG00000167601", "ENSG00000134852", "ENSG00000101266", "ENSG00000100697", "ENSG00000133216", "ENSG00000151491", "ENSG00000068024",
#                 "ENSG00000017427", "ENSG00000073792", "ENSG00000101384", "ENSG00000134333", "ENSG00000288299", "ENSG00000157227", "ENSG00000057935", "ENSG00000150593",
 #                "ENSG00000145675", "ENSG00000181690", "ENSG00000069667", "ENSG00000187764", "ENSG00000141646", "ENSG00000072274")

#KEGG_IDs = mget(as.character(ens_IDS), org.Hs.egPATH,ifnotfound=NA)


 ##A list of gene symbols:
sym = c("DICER1",  "EPHB2",   "CLOCK",   "EPS8",    "PLAG1",   "RORA",    "PDCD4",   "AXL",     "LDHA",    "JAG1",    "SMAD4",  
          "HDAC4",   "IGF2BP2", "TFRC",    "CSNK2A1", "SEMA4D",  "MMP14",   "MTA3",    "PIK3R1",  "IGF1" )

##Get the Entrez gene IDs associated with those symbols
EG_IDs = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)

##Then get the KEGG IDs associated with those entrez genes.
KEGG_IDs = mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA)

#"23405", "2048","9575", "2059","5324","6095","27250", "558","3939", "182", "4089", "9759", "10644", "7037", "1457", "10507","4323","57504","5295","3479"
