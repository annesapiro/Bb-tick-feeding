# this R script plots volcano plots from DESeq2 data across feeding timepoints

# load ggplot2 and ggrepel for labeling
librarian::shelf(ggplot2, ggrepel)

#load DE data to plot
timecourse_DE <- read.delim("~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/timecourse_DE_volcano.regulon.txt")

### day 1 vs subsequent days with all colors in all graphs

#day 2 v day 1
#identify significant genes for each day
sig4 <- subset(timecourse_DE, timecourse_DE$padj_14 < 0.05 & abs(timecourse_DE$log2fc_14) >= 1 & abs(timecourse_DE$log2fc_13) < 1 & abs(timecourse_DE$log2fc_12) < 1)
sig3 <- subset(timecourse_DE, timecourse_DE$padj_13 < 0.05 & abs(timecourse_DE$log2fc_13) >= 1 & abs(timecourse_DE$log2fc_12) < 1)
sig2 <- subset(timecourse_DE, timecourse_DE$padj_12< 0.05 & abs(timecourse_DE$log2fc_12) >= 1)
nonsig <- subset(timecourse_DE, abs(timecourse_DE$log2fc_14) < 1 & abs(timecourse_DE$log2fc_13) < 1 & abs(timecourse_DE$log2fc_12) < 1)

plot12 <- ggplot(data = timecourse_DE, aes(x=log2fc_12, y=padj_12), axis.ticks=10) + 
  xlab("Log 2 Fold Change") + ylab("-log10(padj)") +
  xlim(-4, 4) + ylim(0, 60) + geom_hline(yintercept=1.3) + geom_vline(xintercept=1) + geom_vline(xintercept=-1) +
  geom_point(data=nonsig, aes(x=log2fc_12, y=-log10(padj_12)), shape=20, colour='grey60', size=1.5) +
  geom_point(data=sig4, aes(x=log2fc_12, y=-log10(padj_12)), shape=20, colour='purple4', size=1.5) +
  geom_point(data=sig3, aes(x=log2fc_12, y=-log10(padj_12)), shape=20, colour='firebrick2', size=1.5) +
  geom_point(data=sig2, aes(x=log2fc_12, y=-log10(padj_12)), shape=20, colour='goldenrod1', size=1.5) +
  #label day 2 sig genes
  geom_text_repel(data=sig2, label=sig2$gene, y=-log10(sig2$padj_12), aes(size = "8")) + 
  theme_bw() + 
  theme(axis.text = element_text(size = '8', face = 'bold'),
        panel.border = element_rect(colour='black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.background = element_blank(),
        legend.position = "none")
plot12
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/plots/day2_v_day1_all_colors.pdf", width= 6, height = 6,useDingbats=FALSE)
plot12
dev.off() 

# day 3 v day 1

sig4 <- subset(timecourse_DE, timecourse_DE$padj_14 < 0.05 & abs(timecourse_DE$log2fc_14) >= 1 & abs(timecourse_DE$log2fc_13) < 1 & abs(timecourse_DE$log2fc_12) < 1)
sig2 <- subset(timecourse_DE, timecourse_DE$padj_12< 0.05 & abs(timecourse_DE$log2fc_12) >= 1)
sig3 <- subset(timecourse_DE, timecourse_DE$padj_13 < 0.05 & abs(timecourse_DE$log2fc_13) >= 1 & abs(timecourse_DE$log2fc_12) < 1)
nonsig <- subset(timecourse_DE, abs(timecourse_DE$log2fc_14) < 1 & abs(timecourse_DE$log2fc_13) < 1 & abs(timecourse_DE$log2fc_12) < 1)

plot12 <- ggplot(data = timecourse_DE, aes(x=log2fc_13, y=padj_13), axis.ticks=10) + 
  xlab("Log 2 Fold Change") + ylab("-log10(padj)") +
  xlim(-4, 4) + ylim(0, 60) + geom_hline(yintercept=1.3) + geom_vline(xintercept=1) + geom_vline(xintercept=-1) +
  geom_point(data=nonsig, aes(x=log2fc_13, y=-log10(padj_13)), shape=20, colour='grey60', size=1.5) +
  geom_point(data=sig4, aes(x=log2fc_13, y=-log10(padj_13)), shape=20, colour='purple4', size=1.5) +
  geom_point(data=sig3, aes(x=log2fc_13, y=-log10(padj_13)), shape=20, colour='firebrick2', size=1.5) +
  geom_point(data=sig2, aes(x=log2fc_13, y=-log10(padj_13)), shape=20, colour='goldenrod1', size=1.5) +
  #label day 3 sig genes
  geom_text_repel(data=sig3, label=sig3$gene, y=-log10(sig3$padj_13), aes(size = "8")) + 
  theme_bw() + 
  theme(axis.text = element_text(size = '8', face = 'bold'),
        panel.border = element_rect(colour='black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.background = element_blank(),
        legend.position = "none")
plot12
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/plots/day3_v_day1_all_colors.pdf", width= 6, height = 6,useDingbats=FALSE)
plot12
dev.off() 


# day 1 v day 4 colors

sig4 <- subset(timecourse_DE, timecourse_DE$padj_14 < 0.05 & abs(timecourse_DE$log2fc_14) >= 1 & abs(timecourse_DE$log2fc_13) < 1 & abs(timecourse_DE$log2fc_12) < 1)
#cap -log10(padj) at 60
highUP <- subset(timecourse_DE, -log10(timecourse_DE$padj_14) > 60 & timecourse_DE$log2fc_14 > 0)
highDOWN <- subset(timecourse_DE, -log10(timecourse_DE$padj_14) > 60 & timecourse_DE$log2fc_14 < 0)
sig3 <- subset(timecourse_DE, timecourse_DE$padj_13 < 0.05 & abs(timecourse_DE$log2fc_13) >= 1 & abs(timecourse_DE$log2fc_12) < 1)
sig2 <- subset(timecourse_DE, timecourse_DE$padj_12< 0.05 & abs(timecourse_DE$log2fc_12) >= 1)
nonsig <- subset(timecourse_DE, abs(timecourse_DE$log2fc_14) < 1 & abs(timecourse_DE$log2fc_13) < 1 & abs(timecourse_DE$log2fc_12) < 1)
#cap log2fc at 4 (two points)
right <- subset(timecourse_DE, timecourse_DE$log2fc_14 > 4)

plot12 <- ggplot(data = timecourse_DE, aes(x=log2fc_14, y=padj_14), axis.ticks=10) + 
  xlab("Log 2 Fold Change") + ylab("-log10(padj)") +
  xlim(-4, 4) + ylim(0, 60) + geom_hline(yintercept=1.3) + geom_vline(xintercept=1) + geom_vline(xintercept=-1) +
  geom_point(data=nonsig, aes(x=log2fc_14, y=-log10(padj_14)), shape=20, colour='grey60', size=1.5) +
  geom_point(data=sig4, aes(x=log2fc_14, y=-log10(padj_14)), shape=20, colour='purple4', size=1.5) +
  geom_point(data=sig3, aes(x=log2fc_14, y=-log10(padj_14)), shape=20, colour='firebrick2', size=1.5) +
  geom_point(data=sig2, aes(x=log2fc_14, y=-log10(padj_14)), shape=20, colour='goldenrod1', size=1.5) +
  geom_point(data=highUP, aes(x=log2fc_14, y=60), shape=20, colour='goldenrod1', size=1.5) + 
  geom_point(data=highDOWN, aes(x=log2fc_14, y=60), shape=20, colour='purple4', size=1.5) + 
  geom_point(data=right, aes(x=4, y=-log10(padj_14)), shape=20, colour='purple4', size=1.5) +
  #label day 4 sig genes
  geom_text_repel(data=sig4, label=sig4$gene, y=-log10(sig4$padj_14), aes(size = "8")) + 
  theme_bw() + 
  theme(axis.text = element_text(size = '8', face = 'bold'),
        panel.border = element_rect(colour='black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.background = element_blank(),
        legend.position = "none")
plot12
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/plots/day4_v_day1_all_colors.pdf", width= 6, height = 6,useDingbats=FALSE)
plot12
dev.off() 


#############################################################################################
############################################################################################
# day 4 vs day 1 RpoS-regulated (Caimano et al. 2019) highlight

rposUP <- subset(timecourse_regulon, timecourse_regulon$rpos_annote == 'up')
rposDOWN <- subset(timecourse_regulon, timecourse_regulon$rpos_annote == 'down')
high <- subset(timecourse_regulon, -log10(timecourse_regulon$padj_1e) >= 60)
highDOWN <- subset(rposDOWN, -log10(rposDOWN$padj_1e) >= 60)
highUP <- subset(rposUP, -log10(rposUP$padj_1e) >= 60)
right <- subset(timecourse_regulon, timecourse_regulon$log2fc_1e > 6)

plot12 <- ggplot(data = timecourse_regulon, aes(x=log2fc_1e, y=padj_1e), axis.ticks=10) + 
  xlab("Log 2 Fold Change") + ylab("-log10(padj)") +
  xlim(-6, 6) + ylim(0, 60) + geom_hline(yintercept=1.3) + 
  geom_point(data=timecourse_regulon, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='grey60', size=2) +
  geom_point(data=rposUP, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='dodgerblue', size=2) +
  geom_point(data=rposDOWN, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='deeppink', size=2) +
  geom_point(data=high, aes(x=log2fc_1e, y=60), shape=20, colour='grey60', size=2) +
  geom_point(data=highDOWN, aes(x=log2fc_1e, y=60), shape=20, colour='deeppink', size=2) +
  geom_point(data=highUP, aes(x=log2fc_1e, y=60), shape=20, colour='dodgerblue', size=2) +
  geom_point(data=right, aes(x=6, y=-log10(padj_1e)), shape=20, colour='grey60', size=2) +
  theme_bw() +
  theme(axis.text = element_text(size = '8', face = 'bold'),
        panel.border = element_rect(colour='black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.background = element_blank(),
        legend.position = "none")
plot12
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/plots/day4_v_day1_rpos.pdf", width= 6, height = 6,useDingbats=FALSE)
plot12
dev.off() 


# day 4 vs day 1 Rrp1 regulated (Caimano et al. 2015) highlight

rrp1UP <- subset(timecourse_regulon, timecourse_regulon$rrp1_annote == 'up')
rrp1DOWN <- subset(timecourse_regulon, timecourse_regulon$rrp1_annote == 'down')
high <- subset(timecourse_regulon, -log10(timecourse_regulon$padj_1e) >= 60)
highDOWN <- subset(rrp1DOWN, -log10(rrp1DOWN$padj_1e) >= 60)
highUP <- subset(rrp1UP, -log10(rrp1UP$padj_1e) >= 60)
right <- subset(timecourse_regulon, timecourse_regulon$log2fc_1e > 6)


plot12 <- ggplot(data = timecourse_regulon, aes(x=log2fc_1e, y=padj_1e), axis.ticks=10) + 
  xlab("Log 2 Fold Change") + ylab("-log10(padj)") +
  xlim(-6, 6) + ylim(0, 60) + geom_hline(yintercept=1.3) + 
  geom_point(data=timecourse_regulon, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='grey60', size=2) +
  geom_point(data=rrp1UP, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='dodgerblue', size=2) +
  geom_point(data=rrp1DOWN, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='deeppink', size=2) +
  geom_point(data=high, aes(x=log2fc_1e, y=60), shape=20, colour='grey60', size=2) +
  geom_point(data=highDOWN, aes(x=log2fc_1e, y=60), shape=20, colour='deeppink', size=2) +
  geom_point(data=highUP, aes(x=log2fc_1e, y=60), shape=20, colour='dodgerblue', size=2) +
  geom_point(data=right, aes(x=6, y=-log10(padj_1e)), shape=20, colour='grey60', size=2) +
  theme_bw() + 
  theme(axis.text = element_text(size = '8', face = 'bold'),
        panel.border = element_rect(colour='black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.background = element_blank(),
        legend.position = "none")
plot12
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/plots/day4_v_day1_rrp1.pdf", width= 6, height = 6,useDingbats=FALSE)
plot12
dev.off() 

# day 4 v day 1 RelBbu regulated (Drecktrah et al. 2015) highlight

relUP <- subset(timecourse_regulon, timecourse_regulon$rel_recov_annote == 'up' | timecourse_regulon$rel_starve_annote == 'up' | timecourse_regulon$rel_station_annote == 'up')
relDOWN <- subset(timecourse_regulon, timecourse_regulon$rel_recov_annote == 'down' | timecourse_regulon$rel_starve_annote == 'down' | timecourse_regulon$rel_station_annote == 'down')
high <- subset(timecourse_regulon, -log10(timecourse_regulon$padj_1e) >= 60)
highDOWN <- subset(relDOWN, -log10(relDOWN$padj_1e) >= 60)
highUP <- subset(relUP, -log10(relUP$padj_1e) >= 60)
right <- subset(timecourse_regulon, timecourse_regulon$log2fc_1e > 6)

#plotting DOWN first because there is 1 gene that is both up and down and it actually goes up
plot12 <- ggplot(data = timecourse_regulon, aes(x=log2fc_1e, y=padj_1e), axis.ticks=10) + 
  xlab("Log 2 Fold Change") + ylab("-log10(padj)") +
  xlim(-6, 6) + ylim(0, 60) + geom_hline(yintercept=1.3) + 
  geom_point(data=timecourse_regulon, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='grey60', size=2) +
  geom_point(data=relDOWN, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='deeppink', size=2) +
  geom_point(data=relUP, aes(x=log2fc_1e, y=-log10(padj_1e)), shape=20, colour='dodgerblue', size=2) +
  geom_point(data=high, aes(x=log2fc_1e, y=60), shape=20, colour='grey60', size=2) +
  geom_point(data=highDOWN, aes(x=log2fc_1e, y=60), shape=20, colour='deeppink', size=2) +
  geom_point(data=right, aes(x=6, y=-log10(padj_1e)), shape=20, colour='grey60', size=2) +
  theme_bw() + 
  theme(axis.text = element_text(size = '8', face = 'bold'),
        panel.border = element_rect(colour='black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.background = element_blank(),
        legend.position = "none")
plot12
pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/volcano_plots/plots/day4_v_day1_relbbu.pdf", width= 6, height = 6,useDingbats=FALSE)
plot12
dev.off() 
