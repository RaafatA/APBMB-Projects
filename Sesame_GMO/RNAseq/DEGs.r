# #############################################################
#                                                             #
# project:      Sesamum transcriptome                         #
# Copywrites:   Raafat A. Fahmy                               #
# Title:        Finding DEGs from RNA-seq Data using DESeq2   #
# Date:         24, May 2023                                  #
#                                                             #
###############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
# Load packages
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
install.packages("tidyverse")
library(EnhancedVolcano)
library(dplyr)
library(DESeq2)
library(tidyverse)

# 1. Import and prepare the count data --------------------------------------------------------------------------------------------------
setwd("/media/raafat/927d53ff-9823-4e2d-9c98-cee61adf83f7/home/rafat/APBMB/sesamum/DESeq2")
### 1.1. Read the conut Data 
count_data = as.matrix(read.csv("S_Indicum_read_count_output.csv", row.names=1))
colnames(count_data)

#count_data = count_data[which(rowSums(count_data)>50),]
head(count_data)

### 1.2. Read the experimental meta-data
sample_info <- read.csv("exp_data.csv", header = TRUE, row.names = 1)

### Convert the factors to appropriate levels
sample_info$condition <- factor(sample_info$condition)
sample_info$timepoint <- factor(sample_info$timepoint)
sample_info$condition2 <- factor(sample_info$condition2)
### 

#count_data=apply(count_data,2,as.integer)
#row.names(count_data)=genes
###
#hist(count_data, col = "blue", main="Histogram", breaks = 100)
#hist(log2(count_data+1), col = "blue", main="Histogram")
# 2. Construct a DESeqDataset object ----------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = sample_info, 
                              design = ~ condition)

#dds$timepoint <- factor(dds$timepoint, levels = c("lateT","earlyT"))
dds <- dds[!is.na(dds$timepoint), ]

### 2.1 cleaning the data 
keep = rowSums(counts(dds)) >= 30
dds = dds[keep,]
#### 2.1.1 Removing the null and low gene counts

#### 2.1.2 Setting the factor level: the factor level is the level of the factor that is used as the reference level for all comparisons. 

# 3 Running DESeq -----------------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds)
?results
# contrast = c("condition", "earlyT", "lateT") meaning that genes with logFC > 0 are overexpressed in earlyT
comp1 <- results(dds,contrast = c("condition", "trans", "control"))
write.csv(as.data.frame(comp1),file="comp1_unfiltered.csv", quote=F,row.names=T)

earlyc_earlyt_deseq_result <- results(dds,contrast = c("timepoint", "earlyT", "lateT"))
earlyc_earlyt_deseq_result2 <- results(dds,contrast = c("timepoint", "earlyC", "lateC"))
c_t_deseq_result
earlyc_earlyt_deseq_result
comp1 <- as.data.frame(comp1)
class(c_t_deseq_result)
# order the genes based on the pvalue
earlyc_earlyt_deseq_result <- earlyc_earlyt_deseq_result[order(earlyc_earlyt_deseq_result$pvalue),]
#earlyc_earlyt_deseq_result2 <- earlyc_earlyt_deseq_result2[order(earlyc_earlyt_deseq_result2$pvalue),]

# Filter the most differenanially expressed genes 
# Select Genes With a significant change in gene expression (adjusted pvalue below 0.05)
filterd <- comp1 %>% filter(comp1$padj<0.05)

# Filter based on the fold changes threshold of 2 
filterd <- filterd %>% filter(abs(filterd$log2FoldChange)>2)
dim(filterd)

write_csv(as.matrix(filterd),file="comp2.csv", quote=F,row.names=T)
write.csv(as.matrix(filterd),file="comp7.csv", quote=F,row.names=T)  

EnhancedVolcano(comp1, x= "log2FoldChange", y="padj" ,
                lab=rownames(comp1), 
                title="Control Vs. Trans ",
                subtitle = bquote(italic(padj<0.05 | log2FoldChange > 2 )),
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 6.0,)

?EnhancedVolcano



tophits <- earlyc_earlyt_deseq_result[order(earlyc_earlyt_deseq_result$padj),][1:200,]
tophits <- row.names(tophits)

tophits
library(pheatmap)
rld <- rlog(dds,blind = FALSE)
rld
pheatmap(assay(rld[,c("FS11A_TS1","FS11A_TS2","FS21A_CS1","FS22A_CS2")])[tophits,], )
pheatmap(t(assay(rld[,c("FS21A_CS1","FS22A_CS2","FS23A_CS3","FS24A_CS4")])[tophits,]))



## Data preparation
pval_threshold <- 0.05
t.degs <- row.names(earlyc_earlyt_deseq_result[which(earlyc_earlyt_deseq_result$padj <= pval_threshold), ])
c.degs <- row.names(earlyc_earlyt_deseq_result2[which(earlyc_earlyt_deseq_result2$padj <= pval_threshold), ])

## Venn-diagram using the `VennDiagram` library (see below for alternative method)
library(VennDiagram)
# Arguments for a pairwise (two-sets) venn-diagram are sizes for set1, set2 and overlap (intersect)
# Many more functions are available for triple, quad and quantuple diagrams (starting with 'draw.***')
venn.plot <- draw.pairwise.venn(length(c.degs),
                                length(t.degs),
                                # Calculate the intersection of the two sets
                                length( intersect(c.degs, t.degs) ),
                                category = c("Early_Late Control vs Control 4", "Early_Late Trans "), scaled = F,
                                fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))

# Actually plot the plot
grid.draw(venn.plot)








### 3.1 Create a DESeq2 object
dds.run = DESeq(dds)
vstdata = vst(dds.run ,blind = FALSE)
plotPCA(vstdata, intgroup = "timepoint")
plotDispEsts(dds)
### 3.2 Perform differential expression analysis
res <- results(dds.run)
res <- results(dds.run, contrast=c("condition", "trans", "control"))
res_time <- results(dds.run, contrast = c("timepoint", "earlyT", "earlyC"))
# 4 Exploring the results ---------------------------------------------------------------------------------------------------------------
plotMA(res_time)
plotMA(results)


res_time=as.data.frame(res_time[complete.cases(res), ])
deseq.degtime=res_time[res$padj < 0.05 & abs(res$log2FoldChange)>2,]
write.csv(as.matrix(deseq.degtime),file="deseq.deg.csv", quote=F,row.names=T)  
res.df=as.data.frame(res_time)
EnhancedVolcano(res.df, x= "log2FoldChange", y="padj" ,
                lab=rownames(res.df), 
                title="Early Control Vs. Early Trans",pCutoff = 1e-4,
                FCcutoff = 5,
                pointSize = 3.0,
                labSize = 6.0,)
