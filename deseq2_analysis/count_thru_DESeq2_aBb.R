#########################################################################################################################
# purpose: generate DESeq2 results objects of DE genes
# input: path to directory which contains (at any level) counts files from the same experiment
# output: results tables, MA plots
# written by: anne sapiro, annesapiro@gmail.com with help from calla martyn and chase mateusiak

#########################################################################################################################
#
# experiment: aBb over tick bloodmeal days 1-4
#
################################################## user input #############################################################

# output directory
mappingDir <- '~/Desktop/aBb_tick_feeding_Bb_time_course/deseq2_analysis'
analysisDir <- '~/Desktop/aBb_tick_feeding_Bb_time_course/deseq2_analysis/results'

# path to quant dir (below will work if the (singular) quants dir is in your workingDir. Won't work if there are many quants dirs)
quantDir <- file.path(mappingDir, 'deseq2_files')

############################################### end user input ##############################################################
#############################################################################################################################
###################################### install packages and functions #######################################################

# description of packages
# bioconductor -- repo, biocmanager is used to download
# DESeq2 -- differential gene expression, used throughough
# pheatmap and RColorBrewer -- used to create heatmaps
# ggplot2 -- used specifically for PCA plot, but generally has nice plotting functionality
# tidyverse -- 

# load packages

librarian::shelf(ggplot2, tidyverse, DESeq2, pheatmap, RColorBrewer)

# print session info. Suggestion - save the output of this in your workingDir.
sessionInfo()
################################################### start main script #######################################################

# set working directory
setwd(mappingDir)

sampleFiles <- grep('counts.txt', list.files(quantDir), value = TRUE)
sampleFiles
sampleCondition<-c("day1", "day1", "day1", "day1", "day2", 'day2', 'day2', 'day2',
                   "day3", "day3", "day3", "day3", "day4", "day4", "day4", "day4")

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
sampleTable
# names of sample_replicate -- must be in the same order as in the data directory
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=quantDir, design=~condition)

#### look at dds object
colSums(assay(ddsHTSeq))
colData(ddsHTSeq)
rowData(ddsHTSeq)

############################### use heatmap, pca and dispersion to visualize data ########################################
setwd(analysisDir)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
#### check the percent of genes without zero counts.
# count number of genes
geneCounts <- counts(ddsHTSeq)
# make logical maxtrix of genes with (true) and without (false) zero count
indexNotZero <- apply(geneCounts, 1, function(x) {all (x > 0)})
# display percent of genes with zero counts
sum(indexNotZero == TRUE)/nrow(assay(ddsHTSeq))

#### create heatmap to examine similiarities btwn samples
rld <- rlogTransformation(ddsHTSeq, blind=TRUE)

distsRL <- dist(t(assay(rld)))
distsRL
mat <- as.matrix(distsRL)
rownames(mat) <- colData(rld)$condition
colnames(mat) <- colData(rld)$condition

hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
pheatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
heatmap <- pheatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
pdf(file="heatmap.pdf")
heatmap
dev.off() 

# PCA
# from the vignette re: vst : DESeq2 offers transformations for count data that stabilize the variance across the mean: the regularize logarithm (rlog) and the variance stabilizing transformation (VST). 
# These have slightly different implementations, discussed a bit in the DESeq2 paper and in the vignette, but a similar goal of stablizing the variance across the range of values. 
# Both produce log2-like values for high counts
#vsd <- vst(ddsHTSeq)
vsd <- varianceStabilizingTransformation(ddsHTSeq)
plotPCA(vsd, 'condition')
PCAplot_result <- plotPCA(vsd, 'condition')
pdf(file="PCAplot.pdf",useDingbats=FALSE)
PCAplot_result
dev.off() 


#rld for alt heatmaps
rld <- rlogTransformation(ddsHTSeq, blind=FALSE)
deseq2rld <- assay(rld)
deseq2rld <- as.data.frame(deseq2rld)
deseq2rld$Gene <- rownames(deseq2rld)
head(deseq2rld)
write.csv(as.data.frame(deseq2rld),file="timecourse_deseq2rld.csv", quote=F)

### plot dispersion
ddsHTSeq_disp <- estimateSizeFactors(ddsHTSeq)
ddsHTSeq_disp <- estimateDispersions(ddsHTSeq_disp)
plotDispEsts(ddsHTSeq_disp, main = "Dispersion")

############################################# begin DE analysis #####################################################

# create DESeq object
dds <-  DESeq(ddsHTSeq)
#ddsbp <- DESeq(ddsHTSeq, betaPrior = TRUE)

### output and plot PC1 and PC2 loadings for genes that drive the most variation along those axes
md <- read.csv("~/Desktop/aBb_tick_feeding_Bb_time_course/deseq2_analysis/metadata.csv")
x <-  SummarizedExperiment::assay(dds) # get gene expression values
pca_res <- prcomp(t(vst(x))) #take log, transform and compute principle components
rownames(pca_res$x)=colData(dds)$sample # switch to whatever your "sample column is called"
pca_coord <- as_tibble(pca_res$x, rownames='sample') # extract PC coordinates as tibble
pca_coord <- full_join(md, pca_coord, by='sample') # add back in metadata
summary(pca_res)

data <- as_tibble(pca_res$rotation,rownames='gene')%>%
  select(gene,PC1)%>% #select the PC you are interested in
  arrange(desc(abs(PC1)))# arrange by abolute value of the PC loading

df <- data.frame(data)
p<-ggplot(data=df, aes(x=reorder(gene, -abs(PC1)), y=PC1)) +
  geom_bar(stat = "identity")
p

write.csv(df,file="timecourse_PC1_loadings.csv", quote=F)


df2 <- head(df,100)
p<-ggplot(data=df2, aes(x=reorder(gene, -abs(PC1)), y=PC1)) +
  geom_bar(stat = "identity")
p

pdf(file="PC1_loadings.pdf",useDingbats=FALSE)
p
dev.off() 


### PC2
data2 <- as_tibble(pca_res$rotation,rownames='gene')%>%
  select(gene,PC2)%>% #select the PC you are interested in
  arrange(desc(abs(PC2)))# arrange by abolute value of the PC loading

df3 <- data.frame(data2)
p2<-ggplot(data=df3, aes(x=reorder(gene, -abs(PC2)), y=PC2)) +
  geom_bar(stat = "identity")
p2

write.csv(df,file="timecourse_PC2_loadings.csv", quote=F)


df4 <- head(df3,100)
p<-ggplot(data=df4, aes(x=reorder(gene, -abs(PC2)), y=PC2)) +
  geom_bar(stat = "identity")
p

pdf(file="PC2_loadings.pdf",useDingbats=FALSE)
p
dev.off() 
########################### create results objects ##################################
# 
# # compare day 1 to day 2
# with lfcShrink
res_12_lfc <- lfcShrink(dds, coef = "condition_day2_vs_day1", type="apeglm")
res_12_lfc <- res_12_lfc[order(res_12_lfc$padj),]
head(res_12_lfc)
# make an MA plot to visualize data
plotMA(res_12_lfc, ylim=c(-6,6), main="day2 v day1")

pdf(file="day1_v_day2_apeglm.pdf")
plotMA(res_12_lfc, ylim=c(-4,4), main="day2 v day1")
dev.off() 
# print out the results file as a csv
write.csv(as.data.frame(res_12_lfc),file="day2_v_day1.csv", quote=F)

# # compare day 1 to day 3
# with lfcShrink
res_13_lfc <- lfcShrink(dds, coef = "condition_day3_vs_day1", type="apeglm")
res_13_lfc <- res_13_lfc[order(res_13_lfc$padj),]
head(res_13_lfc)
# make an MA plot to visualize data
plotMA(res_13_lfc, ylim=c(-6,6), main="day3 v day1")

pdf(file="day1_v_day3_apeglm.pdf")
plotMA(res_13_lfc, ylim=c(-6,6), main="day3 v day1")
dev.off() 
# print out the results file as a csv
write.csv(as.data.frame(res_13_lfc),file="day3_v_day1.csv", quote=F)

# # compare day 1 to day4
# with lfcShrink
res_1e_lfc <- lfcShrink(dds, coef = "condition_day4_vs_day1", type="apeglm")
res_1e_lfc <- res_1e_lfc[order(res_1e_lfc$padj),]
head(res_1e_lfc)
# make an MA plot to visualize data
plotMA(res_1e_lfc, ylim=c(-6,6), main="day4 v day1")

pdf(file="day1_v_day4_apeglm.pdf")
plotMA(res_1e_lfc, ylim=c(-6,6), main="day4 v day1")
dev.off() 
# print out the results file as a csv
write.csv(as.data.frame(res_1e_lfc),file="day4_v_day1.csv", quote=F)


######################################################3
# # compare day 2 to day 3 
# first make day2 the reference day
dds$condition <- relevel(dds$condition, ref = "day2")
dds <- nbinomWaldTest(dds)
resultsNames(dds)

# with lfcShrink
res_23_lfc <- lfcShrink(dds, coef = "condition_day3_vs_day2", type="apeglm")
res_23_lfc <- res_23_lfc[order(res_23_lfc$padj),]
head(res_23_lfc)
# make an MA plot to visualize data
plotMA(res_23_lfc, ylim=c(-6,6), main="day2 v day3")

pdf(file="day2_v_day3_apeglm.pdf")
plotMA(res_23_lfc, ylim=c(-6,6), main="day2 v day3")
dev.off() 
# print out the results file as a csv
write.csv(as.data.frame(res_23_lfc),file="day2_v_day3.csv", quote=F)

# # compare day 3 to day4 
# first make day3 the reference day
dds$condition <- relevel(dds$condition, ref = "day3")
dds <- nbinomWaldTest(dds)
resultsNames(dds)

# with lfcShrink
res_3e_lfc <- lfcShrink(dds, coef = "condition_day4_vs_day3", type="apeglm")
res_3e_lfc <- res_3e_lfc[order(res_3e_lfc$padj),]
head(res_3e_lfc)
# make an MA plot to visualize data
plotMA(res_3e_lfc, ylim=c(-6,6), main="day3 v day4")

pdf(file="day3_v_day4_apeglm.pdf")
plotMA(res_3e_lfc, ylim=c(-6,6), main="day3 v day4")
dev.off() 
# print out the results file as a csv
write.csv(as.data.frame(res_3e_lfc),file="day3_v_day4.csv", quote=F)
####################################################################

# make a file of normalized counts, which can help look at genes over time course
dds$condition <- relevel(dds$condition, ref = "day1")
dds <- nbinomWaldTest(dds)
normCounts <- counts(dds, normalized=TRUE)
write.csv(normCounts, file='day1_thru_day4_normCounts.csv', quote=F)

# thess snippets make graphs of 1 gene over all samples, pick your favorites
plotCounts(dds, gene='BB_0001')
plotCounts(dds, gene='BB_A15')
plotCounts(dds, gene='BB_B19')


