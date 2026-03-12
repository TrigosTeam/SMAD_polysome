library(anota2seq)
library(biomaRt)
library(limma)

data_raw <- read.delim("salmon.merged.gene_counts.tsv", sep="\t")
rownames(data_raw) <- data_raw$gene_id
data_raw$gene_id <- NULL
data_raw$gene_name <- NULL

##Get list of coding genes

#mart <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl' )
#genes <- biomaRt::getBM(attributes = c("external_gene_name", "transcript_biotype"),  mart = mart)
#coding_genes <- genes[genes$transcript_biotype == "protein_coding",]
#save(coding_genes, file="coding_genes.Rdata")
load("coding_genes.Rdata")

##Select coding genes from the data
data <- data_raw[rownames(data_raw) %in% coding_genes$external_gene_name,]
colnames(data) <- substr(colnames(data), 10, nchar(colnames(data)))
colnames(data) <- gsub("^\\.", "", colnames(data))

data <- data[, order(colnames(data))]

data[,c("KO2.POLY.N1")] <- NULL
data[,c("KO1.POLY.N1")] <- NULL
data[,c("KO2.INPUT.N1")] <- NULL
data[,c("KO1.INPUT.N1")] <- NULL

data_I <- data[, grep("INPUT", colnames(data))]
data_P <- data[, grep("POLY", colnames(data))]

pheno_vec <- colnames(data_I)
pheno_vec <- substr(pheno_vec, 1, 3)
pheno_vec <- gsub("\\.", "", pheno_vec)


ads <- anota2seqDataSetFromMatrix(
  dataP = data_P,
  dataT = data_I,
  phenoVec = pheno_vec,
  dataType = "RNAseq",
  filterZeroGenes = FALSE,
  normalize = TRUE,
  transformation = "TMM-log2",
  varCutOff = NULL)

ads <- anota2seqPerformQC(Anota2seqDataSet = ads,
                          generateSingleGenePlots = TRUE)

data_norm <- data.frame(ads@dataT, ads@dataP)
boxplot(data_norm)
plotMDS(data_norm)

ads <- anota2seqResidOutlierTest(ads, residFitPlot = FALSE,
                                 generateSingleGenePlots = TRUE, nGraphs = 12)

phenoLev <- levels(as.factor(pheno_vec))
myContrast <- matrix(nrow =length(phenoLev),ncol=length(phenoLev)-1)
rownames(myContrast) <- phenoLev
colnames(myContrast) <- c("TDvsKO", "TDvsCAS", "TDvsKO1")

myContrast[,1] <- c(0, -1, -1, 2)
myContrast[,2] <- c(-1, 0, 0, 1)
myContrast[,3] <- c(0, -1, 0, 1)

ads <- anota2seqAnalyze(Anota2seqDataSet = ads,
                        analysis = c("total mRNA", "translated mRNA",
                                     "translation", "buffering"),
                        contrasts = myContrast)
ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            analysis = c("total mRNA", "translated mRNA",
                                         "translation", "buffering"),
                            selContrast = 1,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.05)

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            analysis = c("total mRNA", "translated mRNA",
                                         "translation", "buffering"),
                            selContrast = 2,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.05)

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            analysis = c("total mRNA", "translated mRNA",
                                         "translation", "buffering"),
                            selContrast = 3,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.05)


ads <- anota2seqRegModes(ads)
head(anota2seqGetOutput(object = ads, output="regModes",
                        selContrast = 1, analysis="buffering",
                        getRVM = TRUE))[, c("apvSlope", "apvEff", "apvRvmP",
                                            "apvRvmPAdj", "singleRegMode")]

anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE)

output_contrast_TDvsKO <-  anota2seqGetOutput(ads,output="singleDf",selContrast=1)
output_contrast_TDvsCAS <-  anota2seqGetOutput(ads,output="singleDf",selContrast=2)
output_contrast_TDvsKO1 <-  anota2seqGetOutput(ads,output="singleDf",selContrast=3)


write.table(output_contrast_TDvsKO, file="TDvsKO_no_outliers.txt",
            sep='\t', quote=FALSE, row.names=FALSE)
write.table(output_contrast_TDvsCAS, file="TDvsCAS_no_outliers.txt",
            sep='\t', quote=FALSE, row.names=FALSE)
write.table(output_contrast_TDvsKO1, file="TDvsKO1_no_outliers.txt",
            sep='\t', quote=FALSE, row.names=FALSE)




###Excluding TD

library(anota2seq)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(limma)

data_raw <- read.delim("salmon.merged.gene_counts.tsv", sep="\t")
rownames(data_raw) <- data_raw$gene_id
data_raw$gene_id <- NULL
data_raw$gene_name <- NULL

##Get list of coding genes
load("coding_genes.Rdata")

##Select coding genes from the data
data <- data_raw[rownames(data_raw) %in% coding_genes$external_gene_name,]
colnames(data) <- substr(colnames(data), 10, nchar(colnames(data)))
colnames(data) <- gsub("^\\.", "", colnames(data))

data <- data[, order(colnames(data))]

data <- data[,-grep("TD", colnames(data))]

data[,c("KO2.POLY.N1")] <- NULL
data[,c("KO1.POLY.N1")] <- NULL
data[,c("KO2.INPUT.N1")] <- NULL
data[,c("KO1.INPUT.N1")] <- NULL

data_I <- data[, grep("INPUT", colnames(data))]
data_P <- data[, grep("POLY", colnames(data))]

pheno_vec <- colnames(data_I)
pheno_vec <- substr(pheno_vec, 1, 3)
pheno_vec <- gsub("\\.", "", pheno_vec)


ads <- anota2seqDataSetFromMatrix(
  dataP = data_P,
  dataT = data_I,
  phenoVec = pheno_vec,
  dataType = "RNAseq",
  filterZeroGenes = FALSE,
  normalize = TRUE,
  transformation = "TMM-log2",
  varCutOff = NULL)

ads <- anota2seqPerformQC(Anota2seqDataSet = ads,
                          generateSingleGenePlots = TRUE)

data_norm <- data.frame(ads@dataT, ads@dataP)
boxplot(data_norm)

plotMDS(data_norm)

ads <- anota2seqResidOutlierTest(ads, residFitPlot = FALSE,
                                 generateSingleGenePlots = TRUE, nGraphs = 12)

phenoLev <- levels(as.factor(pheno_vec))
myContrast <- matrix(nrow =length(phenoLev),ncol=length(phenoLev)-1)
rownames(myContrast) <- phenoLev
colnames(myContrast) <- c("KO1vsCAS", "KO2vsCAS")

myContrast[,1] <- c(-1, 1, 0)
myContrast[,2] <- c(-1, 0, 1)

ads <- anota2seqAnalyze(Anota2seqDataSet = ads,
                        analysis = c("total mRNA", "translated mRNA",
                                     "translation", "buffering"),
                        contrasts = myContrast)
ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            analysis = c("total mRNA", "translated mRNA",
                                         "translation", "buffering"),
                            selContrast = 1,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.05)

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            analysis = c("total mRNA", "translated mRNA",
                                         "translation", "buffering"),
                            selContrast = 2,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.05)

ads <- anota2seqRegModes(ads)
head(anota2seqGetOutput(object = ads, output="regModes",
                        selContrast = 1, analysis="buffering",
                        getRVM = TRUE))[, c("apvSlope", "apvEff", "apvRvmP",
                                            "apvRvmPAdj", "singleRegMode")]

anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE)

output_contrast_KO1vsCAS <-  anota2seqGetOutput(ads,output="singleDf",selContrast=1)
output_contrast_KO2vsCAS <-  anota2seqGetOutput(ads,output="singleDf",selContrast=2)
write.table(output_contrast_KO1vsCAS, file="KO1vsCAS_no_outliers_no_TD.txt",
            sep='\t', quote=FALSE, row.names=FALSE)
write.table(output_contrast_KO2vsCAS, file="KO2vsCAS_no_outliers_no_TD.txt",
            sep='\t', quote=FALSE, row.names=FALSE)
