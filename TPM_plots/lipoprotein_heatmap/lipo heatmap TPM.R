# this script plots a heatmap from TPMs of Bb lipoproteins

#load packages for plot
librarian::shelf(pheatmap, RColorBrewer, viridis)

#import avg norm count data
timecourse_TPMs_avg <- read.delim("~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/lipoprotein_heatmap/avg_TPMs_annotes.txt")

#pull out outer surface lipoproteins
oslipo <- subset(timecourse_TPMs_avg, timecourse_TPMs_avg$lipo_annote == "S")

#remove extraneous info from matrix and set row names to genes
row.names(oslipo) <- oslipo$symbol
oslipo_avgs <- oslipo[,4:7]
oslipo_avgs <- as.data.frame(oslipo_avgs)
head(oslipo_avgs)
df <- data.frame(t(oslipo_avgs))
dim(df)
dataForHeat<-t(df)

#heatmap of log(TPMs)
oslipo_heat_log <- pheatmap(log(dataForHeat), color = viridis(171, direction = -1), show_rownames = T, cluster_rows = TRUE, cluster_cols = FALSE, border_color = NA)
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/lipoprotein_heatmap/oslipo_TPM_heatmap_log.pdf", width= 10, height = 14, useDingbats=FALSE)
oslipo_heat_log
dev.off() 

#normalized heatmap of TPMs to show how expression changes
oslipo_heat_norm <- pheatmap(dataForHeat, color = viridis(171, direction = -1), show_rownames = T, cluster_rows = TRUE, cluster_cols = FALSE, border_color = NA, scale = 'row')
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/lipoprotein_heatmap/oslipo_TPM_heatmap_normalized.pdf", width= 10, height = 14, useDingbats=FALSE)
oslipo_heat_norm
dev.off() 

