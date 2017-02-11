# For linux
library(multicore,  lib.loc="/home/dennis/scratch/R/linux")
library(mgcv,  lib.loc="/home/dennis/scratch/R/linux")
library(doMC,  lib.loc="/home/dennis/scratch/R/linux")
library(lars)
library(e1071)
library(vegan)


load("/scratch/dennis/projects/TFBS_cell_specificity_hg18/K562_GM12878_whole_genome_500bp_hg18.RData")

source("/home/dennis/scratch/projects/expression_prediction/predict_expression_funcs.R")
source("/scratch/dennis/projects/insulators/distr_correlation.R")
source("/scratch/dennis/projects/MK_ChIPseq_analysis/cooccurrence_funcs.R")


gm.chromatin <- as.matrix(read.table("/scratch/dennis/projects/chromatin_states/wgEncodeBroadHmmGm12878HMM.bed", sep="\t", header=FALSE))
colnames(gm.chromatin)[1:4] <- c("chr", "start","end", "state")

gm.crm.chrom <- overlap.crms( gm.combined.mat, gm.chromatin)

chrom.states <- unique(gm.crm.chrom[,"state"])
chrom.count <- c()
for(i in 1:length(chrom.states)){
  chrom.count <- c(chrom.count, length(which(gm.crm.chrom[,"state"] == chrom.states[i])))
}
names(chrom.count) <- chrom.states

png("/scratch/dennis/projects/chromatin_states/output/states_count.png")
par(mar=c(4,9,1,1))
par(las=2)
barplot(sort(chrom.count), horiz=TRUE, xlab="Number of CRMs")
dev.off()


load("/scratch/dennis/projects/chromatin_states/GM12878_CRM_chromatinStates_hg18.RData")
