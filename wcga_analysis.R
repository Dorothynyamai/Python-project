#manipulation miRNA expression data retrieved from GEO (GSE49190)

#setwd("~/Documents/MRNA_MiRNA/HNSC/GSE49190_WCGA/WCGNA")

##Load the necessary libraries

library(dplyr)
library(tidyverse)
library(GEOquery)


dat <- read.csv(file = "GSE49190_non_normalized_data.csv")


# get metadata --------
gse <- getGEO(GEO = 'GSE49190', GSEMatrix = TRUE)

gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)


metadata.modified <- metadata %>%
  select(2,10,12) %>%
  rename(condition = characteristics_ch1) %>%
  rename(tissue = characteristics_ch1.2) %>%
  mutate(condition = gsub("disease state: ", "", condition)) %>%
  mutate(condition = gsub("cell line: ", "", condition)) %>%
  mutate(tissue = gsub("tissue: ", "", tissue))


write.csv(metadata.modified, file = 'metadata.modified.csv')

metadata.modified <- read.csv(file = "metadata.modified.csv")
# looking at gene expression data ---------
head(dat)

# reshaping count data - from wide to long--------
dat.long <- dat %>%
  rename(ILMNID = SYMBOL) %>%
  gather(key = 'samples', value = 'counts', -ILMNID)


# join dataframes = dat.long + metadata.modified

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "geo_accession")) 

# explore data ------
# write the data into a csv file

write.csv(dat.long, file = 'gse49190_counts_data.csv')

##Prepare data for WCGNA
#First select only ilumna ID, counts and sample columns
wcgna.data <- dat.long %>% 
  select(1,2,3) %>% 
  spread(key = 'samples', value = 'counts') %>% #change to the wide format. column will be samples
  column_to_rownames(var = 'ILMNID') #set first column to be the rownames

###############################################
#Now perform wcgna analysis 

#load the necessary libraries
library(WGCNA)
library(DESeq2)
#library(GEOquery) #already loaded
#library(tidyverse) #already loaded
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()          # allow multi-threading (optional)

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes using goodSamplesGenes from wcgna which requires the rows
#to be the samples and the coulmns to be the genes so we transpose the data

gsg <- goodSamplesGenes(t(wcgna.data))
summary(gsg)
gsg$allOK ##if its TRUE it means none of the genes and samples are outliers

#in this case all sampls and genes have passed the test
#to count genes that are outliers
#table(gsg$goodGenes)
#to count samples that are outliers
#table(gsg$goodSamples)

# remove genes that are detectd as outliers
#wcgna.data <- wcgna.data[gsg$goodGenes == TRUE,]

# You can also detect outlier samples - using hierarchical clustering method 
pdf("hiecluster_gse49190.pdf", width = 10, height = 12)
htree <- hclust(dist(t(wcgna.data)), method = "average")
plot(htree)

dev.off()


# checking for outliers using pca method 

pca <- prcomp(t(wcgna.data))
#the information about principal components is stored n a slot called X
pca.dat <- pca$x

##calculate variance explained by each principal component (square of std deviation)
pca.var <- pca$sdev^2
#get percentage of each variance and round it to 2 digits
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

#first convert pca data into a dataframe
pca.dat <- as.data.frame(pca.dat)

#plot first two principal components
pdf("PCA_gse49190.pdf", width = 10, height = 12)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()

### NOTE: If there are batch effects observed, correct for them before moving ahead

# exclude outlier samples
#Example
#samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
#data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset
#phenotypic data
colData <- metadata.modified

#remove the two technical replicates in coldata
write.csv(colData, file = 'coldata.csv')


##read in the data 

colData <- read.csv(file = "coldata_modified.csv", header = TRUE, row.names = 1,)
# incase of outliers then exclude outlier samples
#colData <- metadata %>% 
  #filter(!row.names(.) %in% samples.to.be.excluded)

# fixing column names in colData
names(colData)

#incase of special characters in the heading fix them as follows
#names(colData) <- gsub(':ch1', '', names(colData))
#names(colData) <- gsub('\\s', '_', names(colData)) #to remove space

# making the rownames and column names identical
#if you get an error make sure you specify the first row is a header
all(rownames(colData) %in% colnames(wcgna.data))
#ensure they are in the same order
all(rownames(colData) == colnames(wcgna.data))

wcgna.data <- round(wcgna.data) ##convert data to integers

# create a dds object
dds <- DESeqDataSetFromMatrix(countData = wcgna.data,
                              colData = colData,
                              design = ~ 1) # not specifying model

head(dds)

dds <- DESeq(dds)

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

#dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
#nrow(dds75) # 13284 genes

# perform variance stabilization

#dds_norm <- vst(dds)
#create a new variable vsdata to perform VST of the data
#replace `vst(dds, blind = FALSE)` with `varianceStabilizingTransformation(dds)` in the R code if there not enough 
#samples in the dataset to perform a variance stabilizing transformation (VST) using the DESeq2 package 
#dds_norm <- vst(dds, blind = FALSE)

dds_norm <- varianceStabilizingTransformation(dds, blind = FALSE)

#assay(dds_norm) %>% 
  #head()

# get normalized counts (and transpose the data)
norm.counts <- assay(dds_norm) %>% 
  t()
#head(norm.counts)

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
#create a vector of powers for which the  networks will be calculated
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function #verbose shows output once it is running
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

#use a matrix caclulated by this function to decide which power gives a scale free topology
sft.data <- sft$fitIndices


# visualization to pick power
#plot of rsquare values (a dot scatter plot)

pdf("rsquareplot_gse49190.pdf", width = 10, height = 12)
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +      ###label slightly above the points
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()
#print(a1)
#dev.off()


#Plot for mean connectivity
#pdf("mean_connect_gse49190.pdf", width = 10, height = 12)
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

#visualize as rows
grid.arrange(a1, a2, nrow = 2)

dev.off()

#inorder to run the blockwise expression model functions we need the data to be in numeric form
# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

#choose a power above the red line
soft_power <- 12
#assign cor function to a temp_cor variable
temp_cor <- cor
#assign the cor function to the WCGNA correlation function
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000, ##how many genes should be assigned to one block
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234, #for reproducibility
                          verbose = 3)

#assign the temp_cor function to the original cor function after running the above code
cor <- temp_cor


# 5. Get Module Eigengenes information (stored as MEs) from bwnet object---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview #since we set the names of the numeric names as false the 
#names of  module Eigengenes will be the name of the colors
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

#visualizations
# Plot the dendrogram and the module colors before and after merging underneath

pdf("module_dendogram_gse49190.pdf", width = 10, height = 12)
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

dev.off()


# grey module = all genes that doesn't fall into other modules were assigned to the grey module
#fewer colors in the merged than unmerged indicates that there were several similar modules 
#in the unmerged that were merged

#use the merged modules for further analysis

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

#colData <- read.csv(file = "coldata_modified.csv", header = TRUE, row.names = 1,)

# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(condition_bin = ifelse(grepl('tumor', condition), 1, 0)) %>% 
  select(3) #to select the condition bin


# binarize all other categorical variables as follows

#colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

#severity.out <- binarizeCategoricalColumns(colData$severity,
 #                                          includePairwise = FALSE,
  #                                         includeLevelVsAll = TRUE,
  #                                         minCount = 1)

#then join them
#traits <- cbind(traits, severity.out)


# Define numbers of  samples and genes
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

##calculate the correlation between moduleeigengenes and traits
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
#calculate the p values for the correlation to see the modules that are significantly correlated with disease state
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

#covert first column into row names
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[7],
             y = names(heatmap.data)[1:6],
             col = c("blue1", "skyblue", "white", "pink", "red"))


#identify which miRNAs that are associated with module blue since it has highest correlation
#create a dataframe for the miRNAs and associated modules
module.gene.mapping <- as.data.frame(bwnet$colors)

#select only miRNAs associated with module blue
select_miRNA <-module.gene.mapping %>% 
  filter(`bwnet$colors` == 'blue') %>% 
  rownames()

#write these miRNAs to a csv file
write.csv(select_miRNA, file = 'tumor_assoc_miRNAs.csv')

# 6B. Intramodular analysis: Identifying driver genes ---------------



# Calculate the module membership and the associated p-values
#module eigengene are representative miRNA expression profiles of a cluster

# The module membership/intramodular connectivity is calculated as the correlation of the module eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')

#calculate p values for these measures
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:6,1:10]


# Calculate the miRNA significance and associated p-values
#use Pearsons correlation
gene.signf.corr <- cor(norm.counts, traits$data.tumor.vs.control, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

#top 25 miRNAs
top25_miRNAs <-gene.signf.corr.pvals #%>% 
top25_miRNAs  <- as.data.frame(top25_miRNAs) #%>% 
  head()
  arrange(V1) %>% 
  head(25)


# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

