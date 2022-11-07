# Bb-tick-feeding
This repo contains the intermediate files and code used to generate figures for Longitudinal map of transcriptome changes in the Lyme pathogen Borrelia burgdorferi during tick-borne transmission by Anne L. Sapiro, Beth M. Hayes, Regan F. Volk, Jenny Y. Zhang, Diane M. Brooks, Calla Martyn, Atanas Radkov, Ziyi Zhao, Margie Kinnersley, Patrick R. Secor, Balyn W. Zaro, and Seemay Chou.

Differential expression analysis was performed using DESeq2 with R script `count_thru_DESeq2_aBb.R`

Volcano plots of DESeq output were created with R script `volcanoPlotsColors.R`

Boxplots of TPM data across days for specific genes were created with R script `day by day plots TPM.R`

A heatmap of TPMs of all outer surface lipoproteins was created with R script `lipo heatmap TPM.R`
