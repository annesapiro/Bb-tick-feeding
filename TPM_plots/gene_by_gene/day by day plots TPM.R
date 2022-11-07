#########################################################################################
## this is for making TPM plots of any random gene
#########################################################################################

#load packages
librarian::shelf(ggplot2, ggrepel)

##
timecourse_TPMs_all_clean <- read.delim("~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/gene_by_gene/TPMs_all_clean.txt")

# make a plot with Tukey boxplot of ospA TPMs
bba15 <- subset(timecourse_TPMs_all_clean, timecourse_TPMs_all_clean$gene == "BB_A15")

pbba15 <- ggplot(bba15, aes(x=day, y=TPM, fill=gene)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scale="free", nrow=1) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(legend.position="none") +
  theme(aspect.ratio = .86)
pbba15   

pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/gene_by_gene/ospA.pdf", width= 3, height = 3,useDingbats=FALSE)
pbba15
dev.off() 

# rpoS
bb0771 <- subset(timecourse_TPMs_all_clean, timecourse_TPMs_all_clean$gene == "BB_0771")

pbb0771 <- ggplot(bb0771, aes(x=day, y=TPM, fill=gene)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scale="free", nrow=1) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #labs(title="sig genes on cp26") +
  theme(legend.position="none") +
  theme(aspect.ratio = .86)
pbb0771   

pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/gene_by_gene/rpos.pdf", width= 3, height = 3,useDingbats=FALSE)
pbb0771
dev.off() 

# ospC
bbb19 <- subset(timecourse_TPMs_all_clean, timecourse_TPMs_all_clean$gene == "BB_B19")

pbbb19 <- ggplot(bbb19, aes(x=day, y=TPM, fill=gene)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scale="free", nrow=1) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #labs(title="sig genes on cp26") +
  theme(legend.position="none") +
  theme(aspect.ratio = .86)
pbbb19  

pdf(file="~/Desktop/aBb_tick_feeding_Bb_time_course/TPM_plots/gene_by_gene/ospC.pdf", width= 3, height = 3,useDingbats=FALSE)
pbbb19
dev.off() 






