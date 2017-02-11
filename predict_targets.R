setwd("C:")
# Load CRMs
load("/scratch/dennis/projects/target_gene_prediction/GM12878_whole_genome_500bp_hg18.RData") # 29 TFs in GM12878 - CRMs with 2 or more TFs

load("/scratch/dennis/projects/target_gene_prediction/GM12878_whole_genome_incl_singles_500bp_hg18.RData") # 29 TFs - CRMs with 1 or more TFs

load("/scratch/dennis/projects/target_gene_prediction/GM12878_whole_genome_incl_singles_500bp_hg18_predict_Fibroblast_eQTLs.RData") # Fibroblast eQTLs and CRM co-localized with eQTLs 
load("/scratch/dennis/projects/target_gene_prediction/GM12878_whole_genome_incl_singles_500bp_hg18_predict_Tcells_eQTLs.RData") # Tcell eQTLs and CRM co-localized with eQTLs 
load("/scratch/dennis/projects/target_gene_prediction/TF_cooccurrence_matrix.RData") # TF co-occurrence matrix

library(multicore,  lib.loc="/home/dennis/scratch/R/linux")
library(mgcv,  lib.loc="/home/dennis/scratch/R/linux")
library(doMC,  lib.loc="/home/dennis/scratch/R/linux")
library(lars)
library(e1071)
library(vegan)
library(randomForest)
library(ROCR)
library(party)
source("/home/dennis/scratch/projects/expression_prediction/predict_expression_funcs.R")
source("/scratch/dennis/projects/insulators/distr_correlation.R")
source("/scratch/dennis/projects/MK_ChIPseq_analysis/cooccurrence_funcs.R")
source("/scratch/dennis/projects/target_gene_prediction/classification_performance.R")

gene.pos <- read.table("/home/dennis/scratch/projects/expression_prediction/gene_positions_hg18.txt", header=TRUE, sep="\t")
#gene.pos <- read.table("/scratch/dennis/projects/target_gene_prediction/gene_positions_hg18.txt", header=TRUE, sep="\t")


# Filter for distal CRMs
#gm.combined.mat.promoter <- map.gene.crm(gm.combined.mat, gene.pos, 5000)
#gm.combined.mat.enhancer <- gm.combined.mat[-as.numeric(gm.combined.mat.promoter[,"crm.ind"]),]


# Load FAIRE
gm.faire <- as.matrix(read.table("/scratch/dennis/projects/FAIRE/GM12878_FAIREseq_hg18.pk", sep="\t", header=TRUE))

# Load CTCF
gm.ctcf <-  as.matrix(read.table("/scratch/dennis/projects/chromatin_states/GM12878_CTCF.bed", sep="\t", header=TRUE))

# Load chromatin states
gm.chromatin <- as.matrix(read.table("/scratch/dennis/projects/chromatin_states/wgEncodeBroadHmmGm12878HMM.bed", sep="\t", header=FALSE))
colnames(gm.chromatin)[1:4] <- c("chr", "start","end", "state")

# load SNPs
dimas.eqtl <- as.matrix(read.table("/scratch/dennis/projects/TFBS_cell_specificity_hg18/eQTL/dimas_eQTLs_data.txt", sep="\t", header=TRUE))
LCL.eqtl <- dimas.eqtl[dimas.eqtl[,"Tissue"] == "L",]
dimas.proxy <- as.matrix(read.table("/scratch/dennis/projects/TFBS_cell_specificity_hg18/eQTL/dimas_eQTL_proxies.txt", sep="\t", header=TRUE))
dimas.proxy <- dimas.proxy[-which(dimas.proxy[,"Coordinate_HG18"] == "N/A"),]
dimas.proxy <- cbind(dimas.proxy, Gene=dimas.eqtl[match( dimas.proxy[,"SNP"], dimas.eqtl[,"SNP_Label"]), c("Associated.Gene.Name", "pvalue", "Tissue")])
LCL.proxy <- dimas.proxy[dimas.proxy[,"Tissue"]=="L",]

fibro.proxy <- dimas.proxy[dimas.proxy[,"Tissue"]=="F",]
tcell.proxy <- dimas.proxy[dimas.proxy[,"Tissue"]=="T",]

# load ChIA-PET data for CD4 Tcells
tcell.interact <- as.matrix(read.table("/scratch/dennis/projects/target_gene_prediction/tcell_interact.txt", sep="\t", header=TRUE, skip=1))
colnames(tcell.interact)[1:6] <- c("chr", "start", "end", "chr", "start", "end")
interact.1 <- cbind(region.id= 1:nrow(tcell.interact), tcell.interact[,1:3]) 
interact.2 <-  cbind(region.id= 1:nrow(tcell.interact), tcell.interact[,4:6])



#####################
# load Broad data#
#####################
haem.express <- as.matrix(read.table("/scratch/dennis/projects/expression_prediction/Broad/DMap_data.gct", sep="\t", header=TRUE))
colnames(haem.express)[2] <- "Gene"

all.samples <- 3:dim(haem.express)[2]
#colnames(haem.express) <- unlist(lapply(strsplit(colnames(haem.express), split="_"), function(x){x[1]}))
training.express <- haem.express


#####################
# load LCL expression data for 85 individuals   
#####################
lcl.express <- as.matrix(read.table("/scratch/dennis/projects/eQTL_CRM/LCL_expression/Dimas_LCL_85_individuals.txt", sep="\t", header=TRUE))
all.samples <- 3:dim(lcl.express)[2]
training.express <- lcl.express

# mean centre
training.express[, all.samples] <- t(apply(training.express[, all.samples], 1, function(v){return(  (as.numeric(v) - mean(as.numeric(v))) )}))


                                        # distribution of distance between SNPs and target genes
snp.gene.dist <- c()
for(i in 1:dim(crm.snps)[1]){
  snp.gene.dist <- c(snp.gene.dist, as.numeric(crm.snps[i, "Coordinate_HG18"]) - gene.pos[which(as.character(gene.pos[,"Associated.Gene.Name"]) == crm.snps[i, "Associated.Gene.Name"]), "Transcript.Start..bp."])
}

boxplot(abs(snp.gene.dist) , outline=FALSE, ylab="Distance to target gene (bp)", main="Distances of associated SNP to target gene")

###############
# Find CRMs   #
###############

gm.combined.mat.enhancer <- overlap.crms(gm.combined.mat, gm.chromatin)

crm.snps <- overlap.eqtl(gm.combined.mat.enhancer, LCL.proxy) # only enhancer regions
crm.snps <- overlap.eqtl(gm.combined.mat.enhancer, fibro.proxy)
crm.snps <- overlap.eqtl(gm.combined.mat.enhancer, tcell.proxy)


# remove duplicates
keep.crms <- which(!duplicated(crm.snps[,c("start", "end", "pattern")]))
crm.snps <- crm.snps[keep.crms,]
crm.mat <- crm.snps
crm.gene.preds <- crm.gene.preds[keep.crms]
gene.sim <- gene.sim[keep.crms]
insulators <- insulators[keep.crms]


# print CRM data
crm.dists <- abs(as.numeric(crm.mat[,"Coordinate_HG18"]) - gene.pos[match(crm.mat[,"Associated.Gene.Name"], gene.pos[,"Associated.Gene.Name"]), "Transcript.Start..bp."])
crm.out <- cbind(crm.mat[,c("chr", "start","end", "state", "Proxy", "Associated.Gene.Name")], crm.dists)
crm.out <- cbind(unlist(apply(crm.mat[,1:length(gm.tfs)], 1, function(x){paste(names(x)[which(x==1)], collapse=", ")})) , crm.out)
colnames(crm.out) <- c("TFs", "chromosome", "start", "end", "chromatin state", "SNP", "target_gene", "distance")
crm.out <- crm.out[which(!duplicated(crm.out[,c("TFs","start", "end", "target_gene")])),]
write.table(crm.out, file="/scratch/dennis/projects/target_gene_prediction/output/CREs_table.txt", sep="\t", row.names=FALSE, quote=FALSE)

########################################################################################################
# Analysis
#########################################################################################################

# plot co-localization of TF binding regions, chromatin signatures, open chromatin, and eQTLs
tf.dists <- c()
enhancer.dists <- c()
faire.dists <- c()
snp.dists <- c()
for(chr in unique(gm.combined.mat[,"chr"])){
 tf.dists <- c(tf.dists, sapply(as.numeric(LCL.eqtl[which(LCL.eqtl[,"SNP_chr"] == chr),"SNP_loc"]) , function(x){ as.numeric(gm.combined.mat[which(gm.combined.mat[,"chr"] == chr), "start"]) - x}))
 enhancer.dists <- c(enhancer.dists,  sapply(as.numeric(LCL.eqtl[which(LCL.eqtl[,"SNP_chr"] == chr),"SNP_loc"]) , function(x){ as.numeric(gm.chromatin[which(gm.chromatin[,"chr"] == chr), "start"]) - x}))
 faire.dists <- c(faire.dists,  sapply(as.numeric(LCL.eqtl[which(LCL.eqtl[,"SNP_chr"] == chr),"SNP_loc"]) , function(x){ as.numeric(gm.faire[which(gm.faire[,"chr"] == chr), "start"]) - x}))
 snp.dists <- c(snp.dists,  sapply(as.numeric(LCL.eqtl[which(LCL.eqtl[,"SNP_chr"] == chr),"SNP_loc"]) , function(x){ as.numeric(LCL.proxy[which(LCL.proxy[,"Chromosome"] == chr), "Coordinate_HG18"]) - x}))
}

png("/scratch/dennis/projects/target_gene_prediction/output/distances_from_eQTLs.png")
par(mar=c(5.1,5.1,4.1,2.1))
plot(density(unlist(enhancer.dists),bw=10^5, from=-1*10^7, to=1*10^7), xlab="distance from SNP (kb)", cex.axis=1.5, cex.lab=1.5, main="", col="blue", xaxt='n')
lines(density(unlist(tf.dists),bw=10^5, from=-1*10^7, to=1*10^7), col="red")
lines(density(unlist(faire.dists) ,bw=10^5, from=-1*10^7, to=1*10^7), col="green")
legend("topright", legend=c("chromatin states", "TF binding sites", "open chromatin"), fill=c( "blue", "red", "green"), cex=1.2)
axis(1, at=seq(-1*(10^7), 1*(10^7), by=5*10^6), seq(-1000, 1000, by=500), cex.axis=1.5)
dev.off()


# plot distribution of crms from target gene
crm.centres <- (as.numeric(crm.mat[,"end"]) - as.numeric(crm.mat[,"start"]))/2 + as.numeric(crm.mat[,"start"])
gene.match <- match(crm.mat[,"Associated.Gene.Name"], gene.pos[,"Associated.Gene.Name"])
dist.tmp <- abs(crm.centres[!is.na(gene.match)] - as.numeric(gene.pos[na.omit(gene.match), "Transcript.Start..bp."]))
png("/scratch/dennis/projects/target_gene_prediction/output/crm_target_distances.png")
par(mar=c(5.1,5.1,4.1,2.1))
hist(dist.tmp, xlim=c(0,500000), breaks=10000, xaxt='n', ylab="number of target genes", xlab="distance (kb)", cex.axis=1.5, cex.lab=1.5, main="")
axis(1, at=seq(0, 500000, 50000), seq(0,500, 50), cex.axis=1.5)
dev.off()

# check significance of co-localization between CRMs and eQTL SNPs

a.100 <- c(nrow(overlap.eqtl(gm.combined.mat.enhancer, LCL.proxy, 100)), nrow(overlap.eqtl(gm.combined.mat.enhancer, fibro.proxy, 100)), nrow(overlap.eqtl(gm.combined.mat.enhancer, tcell.proxy, 100)))
a.500 <- c(nrow(overlap.eqtl(gm.combined.mat.enhancer, LCL.proxy, 500)), nrow(overlap.eqtl(gm.combined.mat.enhancer, fibro.proxy, 500)), nrow(overlap.eqtl(gm.combined.mat.enhancer, tcell.proxy, 500)))
a.1000 <- c(nrow(overlap.eqtl(gm.combined.mat.enhancer, LCL.proxy, 1000)), nrow(overlap.eqtl(gm.combined.mat.enhancer, fibro.proxy, 1000)), nrow(overlap.eqtl(gm.combined.mat.enhancer, tcell.proxy, 1000)))

crm.num <- rbind(a.100, a.500, a.1000)
colnames(crm.num) <- c("LCL", "Fibroblast", "T-cell")
rownames(crm.num) <- c("100", "500", "1000")


png("/scratch/dennis/projects/target_gene_prediction/output/eQTL_colocalization.png")
par(mar=c(5,5,2.1,2.1), mgp=c(3,1,0))
barplot2(t(crm.num), beside=TRUE, ylab="Number of elements", xlab="distance from CRE (bases)", cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
legend("topleft", c("LCL", "Fibroblast", "T-cell"), fill=c("red", "orange", "yellow"), cex=1.5, bty='n')
dev.off()



# TF-gene expression regression over distance
#############################################
load("/scratch/dennis/projects/target_gene_prediction/GM12878_whole_genome_predict_expression.RData")

gene.match <- match(training.express[,"Gene"], gene.pos[,"Associated.Gene.Name"])
gene.dists <- rep(NA, dim(training.express)[1])
gene.dists[which(!is.na(gene.match))] <- gene.pos[na.omit(gene.match), "Transcript.Start..bp."]

for(i in 1:length(crm.gene.preds)){
  crm.gene.preds[[i]][,"CRM_ind"] <- rep(i, dim(crm.gene.preds[[i]])[1])
}

tf.gene.dists <- unlist(lapply(crm.gene.preds, function(x){mean(c(as.numeric(crm.mat[as.numeric(x[,"CRM_ind"]), "start"]), as.numeric(crm.mat[as.numeric(x[,"CRM_ind"]), "end"]))) -  gene.dists[as.numeric(x[,"target_expression_ind"])]}))
mean.r2.fits <- unlist(lapply(crm.gene.preds, function(x){x[,"mean.r2.fit"]}))
mean.r2.preds <- unlist(lapply(crm.gene.preds, function(x){x[,"mean.r2.pred"]}))

results <- cbind(tf.gene.dists, mean.r2.fits, mean.r2.preds)
class(results) <- "numeric"
results <- results[order(results[,"tf.gene.dists"]),]
results <- results[which(results[,"tf.gene.dists"] < 2000000 & results[,"tf.gene.dists"] > -2000000),]

a <- split(as.data.frame(results), rep(1:308, each=dim(results)[1]/308))

a.2 <- cbind(unlist(lapply(a,function(x){mean(x[,1])})), unlist(lapply(a,function(x){mean(x[,2])})), unlist(lapply(a,function(x){mean(x[,3])})))

pdf("/scratch/dennis/projects/target_gene_prediction/output/expression_regression_pred.pdf")
par(mar=c(5,5.2,2.1,2.1))
smoothScatter(a.2[,1], a.2[,3], xaxt='n', xlab="distance from TSS (Mb)", ylab=expression(paste(r^2)), cex.lab=1.5, cex.axis=1.5)
lines(lowess(a.2[,1], a.2[,3]), col="red",)
axis(1,seq(-2*10^6, 2*10^6, 10^6), -2:2, cex.axis=1.5)
dev.off()

results <- matrix(NA, nrow=dim(crm.mat)[1], ncol=3)
for(i in 1:dim(crm.mat)[1]){
  if(length(crm.gene.preds[[i]]) > 1){
  pred.acc <- crm.gene.preds[[i]][which(crm.gene.preds[[i]][,"target_gene"] == crm.mat[i,"Associated.Gene.Name"]), c("mean.r2.fit", "mean.r2.pred"), drop=FALSE]

  if(length(pred.acc) > 0){
    results[i,2:3] <- pred.acc[1,]
    results[i, 1] <-  mean(c(as.numeric(crm.mat[i, "start"]), as.numeric(crm.mat[i, "end"]))) - gene.pos[gene.pos[,"Associated.Gene.Name"] == crm.mat[i,"Associated.Gene.Name"], "Transcript.Start..bp."][1]
  }
}
}

png("/scratch/dennis/projects/target_gene_prediction/output/expression_regression_fit_onlyeQTLs.png")
par(mar=c(5,5.2,2.1,2.1))
plot(as.numeric(results[,1]), as.numeric(results[,2]), pch="+", xaxt='n', xlab="distance from TSS (kb)", ylab=expression(paste(R^2)), cex.lab=1.5, cex.axis=1.5)
axis(1,seq(-5*10^5, 5*10^5, 5*10^5), c(-500,0,500), cex.axis=1.5)
lines(lowess(na.omit(as.numeric(results[,1])), na.omit(as.numeric(results[,2]))), col="red")
dev.off()


###################
# Analyze Features#
###################
train.set <- feature.select(crm.gene.preds, crm.mat, 1:length(crm.gene.preds), faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat, dist.thres=0)
class(train.set) <- "numeric"
colnames(train.set) <- c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
train.set[train.set[,"distance"] > 10^6,"distance"] <- 10^6 

keep.crms <- which(!duplicated(crm.mat[,c("start", "end", "pattern", "Associated.Gene.Name")]))
crms.tmp <- crm.mat[keep.crms,]
target.count <- c()
crms.id <- c()
for(i in 1:nrow(crms.tmp)){
 target.count <- length(which( crms.tmp[,"start"] == crms.tmp[i, "start"] &  crms.tmp[,"end"] == crms.tmp[i, "end"] & crms.tmp[,"Associated.Gene.Name"] != crms.tmp[i,"Associated.Gene.Name"]))
 if(target.count > 0){crms.id <- c(crms.id, crms.tmp[i,"start"])}
}

keep.crms <- which(!duplicated(crm.mat[,c("start", "end", "pattern")]))
nontarget.num <- unlist(lapply(crm.gene.preds[keep.crms], function(x){nrow(x)}))
nontarget.num <- nontarget.num - 1
nontarget.num[match(crms.id, crm.mat[keep.crms, "start"])] <- nontarget.num[match(crms.id, crm.mat[keep.crms, "start"])] - 1
names(nontarget.num) <- 1:length(nontarget.num)

png("/scratch/dennis/projects/target_gene_prediction/output/nontargets_distribution.png")
barplot(sort(nontarget.num), ylim=c(0,35), xlab="CRE", ylab="number of non-targets", cex.lab=1.5)
dev.off()
    
# plot distribution of features for target and non-targets
png("/scratch/dennis/projects/target_gene_prediction/output/features_boxplots.png", width=1000, height=500)
par(mar=c(9,6.2,2.1,2.1), mfrow=c(1,6), mgp=c(4.2,1,0))
boxplot(list(train.set[train.set[,"target"] == 1, "distance"], train.set[train.set[,"target"] == 0, "distance"]), ylim=c(0, 10^6), names=c("target", "non-target"),  cex.axis=2, ylab="distance (kb)", cex.lab=2, las=2, yaxt='n', col=c("yellow", "red"))
axis(2, at=seq(0,10^6, 200000), seq(0,1000, 200), cex.axis=2)
boxplot(list(train.set[train.set[,"target"] == 1, "mean.r2.fit"], train.set[train.set[,"target"] == 0, "mean.r2.fit"]),  names=c("target", "non-target"), cex.axis=2, ylab=expression(paste(R^2)), cex.lab=2, las=2, col=c("yellow", "red"))
boxplot(list(train.set[train.set[,"target"] == 1, "crm.pattern"], train.set[train.set[,"target"] == 0, "crm.pattern"]),  names=c("target", "non-target"), cex.axis=2, ylab="co-occurrence score", cex.lab=2, las=2, col=c("yellow", "red"))
boxplot(list(train.set[train.set[,"target"] == 1, "faire.signal"] , train.set[train.set[,"target"] == 0, "faire.signal"] ),  names=c("target", "non-target"), cex.axis=2, ylab="FAIRE signal", cex.lab=2, las=2, col=c("yellow", "red"))
boxplot(list(train.set[train.set[,"target"] == 1, "GO.similarity"], train.set[train.set[,"target"] == 0, "GO.similarity"]),  names=c("target", "non-target"),  cex.axis=2, ylab="GO similarity", cex.lab=2, las=2, col=c("yellow", "red"))
boxplot(list(train.set[train.set[,"target"] == 1, "insulator"], train.set[train.set[,"target"] == 0, "insulator"]),  names=c("target", "non-target"),  cex.axis=2, ylab="CTCF signal", cex.lab=2, las=2, col=c("yellow", "red"))
dev.off()

t.test(train.set[train.set[,"target"] == 1, "insulator"] , train.set[train.set[,"target"] == 0, "insulator"], alternative="two.sided")

##########################
# get variable importance#
##########################

cv.fold <- 10
crm.index <-  sample(1:nrow(crm.mat), nrow(crm.mat))
registerDoMC()
var.acc <- foreach(cv.trial = 1:cv.fold, .combine=rbind)%dopar%{    
  test.genes <- crm.index[ceiling(length(crm.index)*((cv.trial-1)/cv.fold)+1):ceiling(length(crm.index)*(cv.trial/cv.fold))]  # k-fold cv
  train.genes <- crm.index[-match(test.genes, crm.index)]
  
  return(var.importance(gene.preds=crm.gene.preds, crm.mat=crm.mat, train.samps=train.genes, faire.gene=faire.gene, gene.crms=gene.crms, gene.sim=gene.sim, insulators=insulators, score.mat=score.mat, conditional=TRUE))
}

means <-  apply(var.acc ,2, mean)
stdev <- sqrt(apply(var.acc,2, var))
ciw   <- qt(0.975, 10) * stdev / sqrt(nrow(var.acc))
names(means) <- c("gene expression", "distance", "TF co-occurrence", "open chromatin", "GO similarity", "insulator")
png("/scratch/dennis/projects/target_gene_prediction/output/variable_importance.png")
par(mar=c(12.1,7.1,4.1,2.1), mgp=c(5,1,0))
barplot2(means, plot.ci=TRUE, ci.l=means-ciw, ci.u=means+ciw, xpd=FALSE, las=2, ylab="importance score", cex.lab=2, cex.axis=1.5, cex.names=1.5)
dev.off()

#################################
# predict target genes of eQTLs #
#################################

crm.mat <- crm.snps

# predict target gene expression
crm.gene.preds <- predict.targets(crm.mat, training.express, training.express, all.samples, gene.pos=gene.pos, regress.method="RF")

##
# Get GO similarity for TF-gene pairs
##
library(GOSim)
tf.entrez <- as.numeric(training.express[match(names(gm.tfs), training.express[,"Gene"]),"NAME"])
names(tf.entrez) <- names(gm.tfs)
registerDoMC()
gene.sim <- foreach(i = 1:length(crm.gene.preds))%dopar%{
  if(!is.na(crm.gene.preds[[i]])){
    matched.genes <- as.numeric(training.express[as.numeric(crm.gene.preds[[i]][,"target_expression_ind"]), "NAME"])
    names(matched.genes) <- crm.gene.preds[[i]][,"target_gene"]
    crm.tfs <- tf.entrez[ which(crm.mat[i,1:length(gm.tfs)]==1)]
    
    crm.gene.sim <- tryCatch(getGeneSim(as.character(crm.tfs), as.character(matched.genes)), error = function(e){NA})
    if(length(crm.gene.sim) > 1){
      rownames(crm.gene.sim) <- names(crm.tfs)[match(rownames(crm.gene.sim), as.character(crm.tfs))]  
      colnames(crm.gene.sim) <- names(matched.genes[match(colnames(crm.gene.sim), matched.genes)])
    }
    return(crm.gene.sim)
  }
  else{return(NA)}
}

# Find insulators#
registerDoMC()
insulators <- foreach(i = 1:length(crm.gene.preds))%dopar%{
  if(!is.na(crm.gene.preds[[i]])){
    matched.genes <- gene.pos[match(crm.gene.preds[[i]][,"target_gene"], gene.pos[,"Associated.Gene.Name"]),]
    ctcf.chr <- gm.ctcf[gm.ctcf[,"chr"] == crm.mat[i,"chr"],]
    ctcf.upstream <- ctcf.chr[which(as.numeric(ctcf.chr[,"start"]) > as.numeric(crm.mat[i, "Coordinate_HG18"])),]
    ctcf.downstream <- ctcf.chr[which(as.numeric(ctcf.chr[,"start"]) < as.numeric(crm.mat[i, "Coordinate_HG18"])),]
    gene.ctcf <- rep(0, dim(crm.gene.preds[[i]])[1])
    names(gene.ctcf) <- crm.gene.preds[[i]][,"target_gene"]
    for(j in 1:dim(matched.genes)[1]){
      if(matched.genes[j,"Transcript.Start..bp."] > crm.mat[i, "Coordinate_HG18"]){
        gene.ctcf[names(gene.ctcf) == matched.genes[j, "Associated.Gene.Name"]] <- mean(as.numeric(ctcf.upstream[which(as.numeric(ctcf.upstream[,"start"]) < matched.genes[j,"Transcript.Start..bp."]), "signalValue"]))
      }
      if(matched.genes[j,"Transcript.Start..bp."] < crm.mat[i, "Coordinate_HG18"]){
        gene.ctcf[names(gene.ctcf) == matched.genes[j, "Associated.Gene.Name"]] <- mean(as.numeric(ctcf.downstream[which(as.numeric(ctcf.downstream[,"start"]) > matched.genes[j,"Transcript.Start..bp."]), "signalValue"]))
      }
    }
    gene.ctcf[is.na(gene.ctcf)] <- 0       
    return(gene.ctcf)
  }
  else{return(NA)}
}

# TF co-occurrence matrix
score.mat <- cooccur.score(gm.combined.mat[,1:length(gm.tfs)])



# cross-validation of target predictions
crm.index <-  sample(1:length(crm.gene.preds), length(crm.gene.preds))
classification <- FALSE
registerDoMC()
result <- foreach(min.dist = seq(0,150000, 5000), .combine=rbind)%dopar%{
  perf.vars <- c()
  perf.prec <- c()
  perf.rec <- c()
  pred.lst <- list()
  obs.lst <- list()
  samps.lst <- list()
 i.var <-  c("distance",  "mean.r2.fit", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
 # for(trial in 1:10){
 for(var.ind in 1:length(i.var)){
   
  #for(i.var in c("distance",  "crm.pattern", "mean.r2.fit", "faire.signal", "GO.similarity", "insulator")){
    cv.preds <- list()
    cv.obs <- list()
    cv.acc <- list()
    cv.samps <- list()
    cv.fold <- 10
    for (cv.trial in 1:cv.fold){    
                                        # test.genes <- sample(1:length(crm.gene.preds), length(crm.gene.preds)*0.2) # random resampling cv
      
      test.genes <- crm.index[ceiling(length(crm.index)*((cv.trial-1)/cv.fold)+1):ceiling(length(crm.index)*(cv.trial/cv.fold))]  # k-fold cv
      train.genes <- crm.index[-match(test.genes, crm.index)]
  
      target.preds <- classify.targets(gene.preds=crm.gene.preds, crm.mat=crm.mat, train.samps=train.genes, test.samps=test.genes, faire.gene=faire.gene, gene.crms=gene.crms, gene.sim=gene.sim, insulators=insulators, method="RF-regression", model.var=c("target", i.var[var.ind:length(i.var)]), score.mat=score.mat, dist.thres=min.dist)
      
      cv.preds[[cv.trial]] <- as.numeric(target.preds[,"pred"])
      cv.obs[[cv.trial]] <- as.numeric(target.preds[,"obs"])
      cv.samps[[cv.trial]] <- target.preds[,"crm.genes"]
      if(classification){cv.acc[[cv.trial]] <- prec.rec(target.preds[,"pred"], target.preds[,"obs"])}
    }
    
    if(classification){
      perf.vars <- cbind(perf.vars, mean(unlist(lapply(cv.acc, mean))))
      perf.prec <- cbind(perf.prec, mean(unlist(lapply(cv.acc, function(x){x[1]}))))
      perf.rec <- cbind(perf.rec, mean(unlist(lapply(cv.acc, function(x){x[2]}))))
    }else{
      #perf.vars <- cbind(perf.vars, mean(unlist( performance(prediction(cv.preds, cv.obs) , "auc")@y.values)))
      pred.lst[[var.ind]] <- cv.preds
      obs.lst[[var.ind]] <- cv.obs
      samps.lst[[var.ind]] <- cv.samps
    }
  }
  
  return(perf.vars)
}
colnames(result) <- c("distance", "gene expression", "TF co-occurrence", "open chromatin", "GO similarity", "insulator")

# plot AUCs or Precision & Recall
png("/scratch/dennis/projects/target_gene_prediction/output/auc_genes_by_distance.png")
par(mar=c(5.1,5.1,4.1,2.1))
matplot( seq(0,150, 10), result[seq(1, 31, 2),], type="l", col=c(1, rainbow(6)), ylim=c(0.2,1), lwd=3, xlab="minimum distance from CRE (kb)", ylab="average AUC", xaxt='n', cex.lab=1.5, cex.axis=1.5)
matpoints(seq(0,150, 10), result[seq(1, 31, 2),], pch="+", col=c(1, rainbow(6)))
axis(1, seq(0,150, 10), cex.axis=1.5)
legend("bottomright", legend=colnames(result), fill=c(1, rainbow(6)), cex=1.3, bty="n")
dev.off()

##################
# plot ROC curves#
##################
perf <- performance(prediction(cv.preds,cv.obs) , "tpr", "fpr")
perf.dist <- performance(prediction(cv.preds, cv.obs) , "tpr", "fpr")

png("/scratch/dennis/projects/target_gene_prediction/output/ROC_full_model.png")
par(mar=c(5.1,5.1,4.1,2.1))
plot(perf, col="grey82", lyt=3, ylab="Sensitivity", xlab="1-Specificity", cex.lab=1.5)
plot(perf, lwd=3, avg="vertical", col=1, add=TRUE)
plot(perf.dist, col=colors()[86], lyt=3, add=TRUE,  cex.lab=1.5)
plot(perf.dist, lwd=3, avg="vertical", col=colors()[258], add=TRUE)
abline(0:1, col="red")
legend("bottomright", legend=c("all six features", "only genomic distance"), lty=1, lwd=3, col=c("black", colors()[258]), cex=1.3, bty="n")
dev.off()

classify.eval(unlist(cv.preds), unlist(cv.obs))

png("/scratch/dennis/projects/target_gene_prediction/output/eval_full_model_150kbp_limit.png")
binary_eval(unlist(cv.preds), unlist(cv.obs), cutoff=0.4664)
dev.off()

png("/scratch/dennis/projects/target_gene_prediction/output/ROC_step_models.png")
perf <- performance(prediction(pred.lst[[1]], obs.lst[[1]]) , "tpr", "fpr")
auc <- c(round(mean(unlist( performance(prediction(pred.lst[[1]], obs.lst[[1]]) , "auc")@y.values)), 3))
par(mar=c(5.1,5.1,4.1,2.1))
plot(perf, lwd=3, avg="vertical", ylab="Sensitivity", xlab="1-Specificity", cex.lab=1.5)
for(i in 2:length(pred.lst)){
  perf <- performance(prediction(pred.lst[[i]], obs.lst[[i]]) , "tpr", "fpr")
  plot(perf, lwd=3, avg="vertical", add=TRUE, col=i)
  auc <- c(auc, c(round(mean(unlist( performance(prediction(pred.lst[[i]], obs.lst[[i]]) , "auc")@y.values)), 3)))
}
abline(0:1, col="red", lty=2)
legend("bottomright", c(paste("Model6 = ", auc[1]), paste("Model5 = ", auc[2]), paste("Model4 = ", auc[3]), paste("Model3 = ", auc[4]), paste("Model2 = ", auc[5]), paste("Model1 = ", auc[6])), title="AUC:", lty=1, lwd=3, col=1:6, cex=1.3)
dev.off()


##
# Predict eQTLs in Fibroblasts and Tcells
##
perf <- performance(prediction(unlist(cv.preds), unlist(cv.obs)), "tpr", "fpr")

load("/scratch/dennis/projects/target_gene_prediction/perf_tcell.RData") # Performance data for predicting T-cell eQTLs: perf.tcell
load("/scratch/dennis/projects/target_gene_prediction/perf_fibroblast.RData") # predicting fibroblast eQTLs: perf.fibro

png("/scratch/dennis/projects/target_gene_prediction/output/ROC_different_cells2.png")
par(mar=c(5.1,5.1,4.1,2.1))
plot(perf, lwd=3,, ylab="Sensitivity", xlab="1-Specificity", cex.lab=1.5)
plot(perf.tcell, lwd=3, avg="vertical", add=TRUE, col="blue")
plot(perf.tcell.dist, lwd=3, avg="vertical", add=TRUE, col="blue", lty=3)
plot(perf.fibro, lwd=3, avg="vertical", add=TRUE, col="green")
plot(perf.fibro.dist, lwd=3, avg="vertical", add=TRUE, col="green", lty=3)
abline(0:1, col="red")
legend("bottomright", c("targets of LCL eQTLs", "targets of T-cell eQTLs", "targets of fibroblast eQTLs"), fill=c("black", "blue", "green"), cex=1.5, bty='n')
dev.off()

              

# map target genes
crm.index <-  sample(1:length(crm.gene.preds), length(crm.gene.preds))
test.genes <- crm.index[ceiling(length(crm.index)*((1-1)/10)+1):ceiling(length(crm.index)*(1/10))]
train.genes <-   crm.index[-match(test.genes, crm.index)]
target.preds <- classify.targets(gene.preds=crm.gene.preds, crm.mat=crm.mat, train.samps=train.genes, test.samps=test.genes, faire.gene=faire.gene, gene.crms=gene.crms, gene.sim=gene.sim, insulators=insulators, score.mat=score.mat, method="RF-regression", model.var=c("target", "crm.pattern"),  dist.thres=0)
prec.rec(target.preds[,"pred"], target.preds[,"obs"])
unlist(performance(prediction(as.numeric(target.preds[,"pred"]), as.numeric(target.preds[,"obs"])) , "auc")@y.values)


class(target.preds) <- "numeric"
plot(performance(prediction(target.preds[,"pred"], target.preds[,"obs"]) , "tpr", "fpr"))
abline(0:1, col="red")

# plot prediction accuracy by chromatin state element type
trial.lst <- list()
for(trial in 1:10){
  cv.preds <- pred.lst[[trial]]
  cv.obs <- obs.lst[[trial]]
  cv.samps <- samps.lst[[trial]]
  
  crm.samps <- unlist(cv.samps)
  crm.ind <- as.numeric(unlist(lapply(sapply(crm.samps, strsplit, ":"), function(x){x[1]})))
  pred.perf <- list()
  
for(i in 1:length(unique(crm.mat[,"state"]))){
  crm.preds <- unlist(cv.preds)[ which(crm.mat[crm.ind,"state"] == unique(crm.mat[,"state"])[i])]
  crm.obs <- unlist(cv.obs)[ which(crm.mat[crm.ind,"state"] == unique(crm.mat[,"state"])[i])]
   
  if(length(which(crm.obs == 1)) > 0 & length(crm.obs) > 30){ 
    pred.perf[[i]] <- binary_eval(crm.preds, crm.obs)
  } else{ pred.perf[[i]] <- NA}

}
  names(pred.perf) <- unique(crm.mat[,"state"])
  pred.perf <- pred.perf[!is.na(pred.perf)]
  perf.mat <- cbind(unlist(lapply(pred.perf, function(x){x["Sensitivity"]})), unlist(lapply(pred.perf, function(x){x["Specificity"]})),  unlist(lapply(pred.perf, function(x){x["Precision"]})))
  rownames(perf.mat) <- names(pred.perf)
  colnames(perf.mat) <- c("Sensitivity", "Specificity", "Precision")
  trial.lst[[trial]] <- perf.mat
}

crm.sens <- c()
crm.spec <- c()
crm.prec <- c()
for(i in 1:length(trial.lst)){
  crm.sens <- cbind(crm.sens, trial.lst[[i]][,"Sensitivity"])
  crm.spec <- cbind(crm.spec, trial.lst[[i]][,"Specificity"])
  crm.prec <- cbind(crm.prec, trial.lst[[i]][,"Precision"])
}


means <- rbind(apply(crm.sens,1, mean), apply(crm.spec,1, mean), apply(crm.prec,1, mean))
stdev <- rbind(sqrt(apply(crm.sens,1, var)), sqrt(apply(crm.spec,1, var)), sqrt(apply(crm.prec,1, var)))
ciw   <- rbind(qt(0.975, 10) * stdev[1,] / sqrt(10), qt(0.975, 10) * stdev[2,] / sqrt(10), qt(0.975, 10) * stdev[3,] / sqrt(10))

pdf("/scratch/dennis/projects/target_gene_prediction/output/performance_chromatin_states.pdf", width=1000, height=400)
par(mar=c(12.1,5.1,4.1,2.1))
par(mfrow=c(1,3))
barplot2(means[1,], plot.ci=TRUE, ci.l=means[1,]-ciw[1,], ci.u=means[1,]+ciw[1,], ylim=c(0,1.05), xpd=FALSE, las=2, ylab="Sensitivity", cex.lab=2, cex.axis=1.5, cex.names=1.5)
barplot2(means[2,], plot.ci=TRUE, ci.l=means[2,]-ciw[2,], ci.u=means[2,]+ciw[2,], ylim=c(0,1.05), xpd=FALSE, las=2, ylab="Specificity", cex.lab=2, cex.axis=1.5, cex.names=1.5)
barplot2(means[3,], plot.ci=TRUE, ci.l=means[3,]-ciw[3,], ci.u=means[3,]+ciw[3,], ylim=c(0,1.05), xpd=FALSE, las=2, ylab="Precision", cex.lab=2, cex.axis=1.5, cex.names=1.5)
dev.off()

#########################
# Validate with ChIA-PET#
#########################

# eqtl matches
crm.interact.1 <- overlap.crms(crm.mat, interact.1, dist.thres=100)
crm.interact.2 <-  overlap.crms(crm.mat, interact.2, dist.thres=100)

crm.interact <- rbind(crm.interact.1, crm.interact.2)
crm.interact <- crm.interact[which(!duplicated(crm.interact[,c("start", "end", "pattern", "chr", "Associated.Gene.Name")])),]

tcell.eqtl <- dimas.eqtl[dimas.eqtl[,"Tissue"] == "T",]
target.pos <- tcell.eqtl[ match(crm.interact[,"Associated.Gene.Name"], tcell.eqtl[, "Associated.Gene.Name"]), "Gene_loc"]

matched.regions <- find.interact(target.pos, crm.interact)
crm.interact[matched.regions,]

# print chromatin interactions
colnames(crm.interact)[55:57] <- c("element_chr", "element_start", "element_end")
write.table(crm.interact[matched.regions,c("Associated.Gene.Name","Proxy", "chr", "start", "end", "element_chr", "element_start", "element_end")], file="/scratch/dennis/projects/target_gene_prediction/output/ChIA-pet_mapped.txt", row.names=FALSE, quote=FALSE)


# CRM matches
target.preds <- cbind(pred=unlist(cv.preds), obs=unlist(cv.obs), crm.genes=unlist(cv.samps))
target.preds <- target.preds[as.numeric(target.preds[,"pred"]) > 0,]

crm.mat.2 <- cbind(crm.mat[as.numeric(unlist(lapply(strsplit(target.preds[, "crm.genes"], split=":"), function(x){x[1]}))),, drop=FALSE], target.preds)
crm.interact.1 <- overlap.crms(crm.mat.2, interact.1, dist.thres=100)
crm.interact.2 <-  overlap.crms(crm.mat.2, interact.2, dist.thres=100)
crm.interact <- rbind(crm.interact.1, crm.interact.2)

target.pos <- c()
keep.crms <- c()
for(i in 1:nrow(target.preds)){
  crm.ind <- as.numeric(unlist(lapply(strsplit(target.preds[i, "crm.genes"], split=":"), function(x){x[1]})))
  gene.ind <- as.numeric(unlist(lapply(strsplit(target.preds[i, "crm.genes"], split=":"), function(x){x[2]})))
 
  crm.match <- which(as.numeric( crm.interact[,"start"]) == as.numeric(crm.mat[crm.ind, "start"] ) &  as.numeric(crm.interact[,"end"]) == as.numeric(crm.mat[crm.ind, "end"] ) &  crm.interact[,"chr"] == crm.mat[crm.ind, "chr"] )
  gene.match <- match( crm.gene.preds[[crm.ind]][gene.ind,"target_gene"] ,gene.pos[,"Associated.Gene.Name"])
  if(length(crm.match) > 0 & length(gene.match) > 0){
    keep.crms <- c(keep.crms, crm.match[1])
    target.pos <- c(target.pos, gene.pos[gene.match, "Transcript.Start..bp."])
  }
}
crm.interact <- crm.interact[keep.crms,]
target.pos <- target.pos[which(!duplicated(crm.interact[,c("start", "end", "chr")]))]
crm.interact <- crm.interact[which(!duplicated(crm.interact[,c("start", "end", "chr")])),]

matched.regions <- find.interact(target.pos, crm.interact)
crm.interact[matched.regions,]

##############################
# Validate with 5C - in hg19 #
##############################

chrom.inter <- read.table("/scratch/dennis/projects/chromatin_states/5C/GM12878_hg19_matrix.matrix", header=TRUE, row.names=1)
crms1 <- lapply(strsplit(colnames(chrom.inter), split='[.]'), function(x){x[3:5]})
crms1 <- cbind(1:length(crms1), unlist(lapply(crms1, function(x){x[1]})), unlist(lapply(crms1, function(x){x[2]})), unlist(lapply(crms1, function(x){x[3]})))
colnames(crms1) <- c("col_index", "chr", "start", "end")

crms2 <- unlist(lapply(strsplit(rownames(chrom.inter), split='[|]'), function(x){x[3]}))
crms2 <- cbind(1:length(crms2), unlist(lapply(strsplit(crms2, split=':'), function(x){x[1]})), unlist(lapply(strsplit(unlist(lapply(strsplit(crms2, split=':'), function(x){x[2]})), split='-'), function(x){x[1]})), unlist(lapply(strsplit(unlist(lapply(strsplit(crms2, split=':'), function(x){x[2]})), split='-'), function(x){x[2]})))
colnames(crms2) <- c("row_index", "chr", "start", "end")

# get predicted target and non-targets
samps.tested <- c()
for(i in 1:length(unlist(cv.samps))){
  crm.ind <- as.numeric(unlist(strsplit(unlist(cv.samps)[i], ":"))[1])
 test.crm <- crm.mat[crm.ind, c("chr", "start", "end")]
  test.gene <- crm.gene.preds[[crm.ind]][as.numeric(unlist(strsplit(unlist(cv.samps)[i], ":"))[2]) , "target_gene"]
  samps.tested <- rbind(samps.tested, c(test.crm, test.gene, unlist(cv.preds)[i], unlist(cv.obs)[i]))
}
colnames(samps.tested) <- c("chr", "start", "end", "target_gene", "pred", "obs")

# liftover from hg18 to hg19
samps.tested[,"start"] <- as.integer(samps.tested[,"start"])
samps.tested[,"end"] <- as.integer(samps.tested[,"end"]) +1
write.table(samps.tested, file="/scratch/dennis/projects/target_gene_prediction/predictions_hg18.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
samps.tested <- as.matrix(read.table("/scratch/dennis/projects/target_gene_prediction/predictions_hg18.bed"))
colnames(samps.tested) <- c("chr", "start", "end", "target_gene", "pred", "obs")

# add gene positions
gene.pos.hg19 <- as.matrix(read.table("/scratch/dennis/projects/target_gene_prediction/gene_positions_hg19.txt", header=TRUE, sep="\t"))
samps.tested <- cbind(samps.tested, TSS=gene.pos.hg19[match(samps.tested[,"target_gene"], gene.pos.hg19[,"Associated.Gene.Name"]),"Transcript.Start..bp."])

crm.inter <- overlap.crms(crms2, samps.tested, dist.thres=1000)



find.interact <- function(target.pos, crm.interact){
  crm.target <- c()
  for(i in 1:nrow(crm.interact)){
    if(as.numeric(interact.1[as.numeric(crm.interact[i, "region.id"]), "end"]) == as.numeric(crm.interact[i, ncol(crm.interact)])){
      crm.target <- rbind(crm.target, interact.2[as.numeric(crm.interact[i, "region.id"]),])
    }else{
      crm.target <- rbind(crm.target, interact.1[as.numeric(crm.interact[i, "region.id"]),])
    }
  }
  matched.regions <- c()
  for(i in 1:length(target.pos)){
    if(target.pos[i] > (as.numeric(crm.target[i, "start"]) - 1000) &&  target.pos[i] < (as.numeric(crm.target[i, "end"]) + 1000)){
      matched.regions <- c(matched.regions, i)
    }
  }
  return(matched.regions)
}

predict.targets <- function(crm.mat, training.express, test.express, all.samples, max.distance=1000000, regress.method="RF", gene.pos){
registerDoMC()
crm.gene.preds <- foreach(crm.ind = 1:dim(crm.mat)[1])%dopar%{
  gene.preds <- c()
 gene.chr<- gene.pos[gene.pos[,"Chromosome.Name"] == crm.mat[crm.ind, "chr"],]
 matched.genes <- as.character(gene.chr[which(abs(gene.chr[,"Transcript.Start..bp."] - mean(c(as.numeric(crm.mat[crm.ind,"start"]), as.numeric(crm.mat[crm.ind,"end"])))) <= max.distance), "Associated.Gene.Name"])
 if(length(matched.genes) > 0){
   matched.genes <- matched.genes[!duplicated(matched.genes)]
   combined.mat <- crm.mat[crm.ind,, drop=FALSE]
   combined.mat <- cbind(combined.mat[rep(nrow(combined.mat), each=length(matched.genes)),, drop=FALSE], matched.genes)
   colnames(combined.mat)[dim(combined.mat)[2]] <- "gene"
 
   tfs <- gm.tfs
  
   patterns <- unique(as.numeric(combined.mat[,"pattern"]))
   names(patterns) <- patterns
   pred.results.lst <- list()
   for(k in 0:4){
     cv.samps <- ceiling(length(all.samples)/5)
     test.express.cols <- c(all.samples[max(cv.samps*k, 1) : min(cv.samps*(k+1), length(all.samples))])
     training.cols <- all.samples[-match(test.express.cols, all.samples)]
     
     pred.results <- predict.gene.express.2(patterns, combined.mat, tfs, training.express, test.express, training.cols, test.express.cols, method=regress.method, step.wise=FALSE, bootstrap=FALSE, null.tfs=all.tfs)
     if(length(pred.results) < 1){break}
     obs.express <- test.express[as.numeric(pred.results[,"target_expression_ind"]), test.express.cols, drop=FALSE]
     r.squareds <- c()
     for(i in 1:dim(pred.results)[1]){
       r.squareds <- c(r.squareds, summary(lm(as.numeric(pred.results[i, 4:(dim(pred.results)[2] -1)]) ~ as.numeric(obs.express[i,])))$r.squared)
     }
     pred.results.lst[[k+1]] <- cbind(pred.results, r.squareds)
   }
   
   if(length(pred.results.lst) < 1){return(NA)}
   
   gene.preds <- pred.results.lst[[1]][,c("CRM_ind", "target_expression_ind", "target_gene"), drop=FALSE]
   mean.r2.fit <- c()
   mean.r2.pred <- c()
   for(i in 1:dim(gene.preds)[1]){
     mean.r2.fit <- c(mean.r2.fit, mean(as.numeric(unlist(lapply(pred.results.lst, function(x){x[i, "R.squared_fit"]})))))
     mean.r2.pred <- c(mean.r2.pred, mean(as.numeric(unlist(lapply(pred.results.lst, function(x){x[i, "r.squareds"]})))))
   }
   gene.preds <- cbind(gene.preds, mean.r2.fit, mean.r2.pred)
 }
  return(gene.preds)
}
return(crm.gene.preds)
}


classify.targets <- function(gene.preds, crm.mat, train.samps, test.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat, method="RF-classification", model.var=c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator"), dist.thres=0){
  trainset <- feature.select(gene.preds, crm.mat, train.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat)
  class(trainset) <- "numeric"
  colnames(trainset) <- c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
  rownames(trainset) <- 1:dim(trainset)[1]
  trainset <- trainset[,match(model.var, colnames(trainset))]
  model.obj <- NULL
  if(method=="svm"){
    model.obj <- try(svm( target ~ . , data=trainset))
  }
  if(method=="gam"){
    model.obj <- gam(target ~ s(mean.r2.fit) + enhancer.type + crm.pattern + s(faire.signal) + s(GO.similarity) + s(insulator), data=as.data.frame(trainset))
  }
  if(method =="RF-classification"){
    model.obj <- try(randomForest( trainset[,-1, drop=FALSE] , as.factor(trainset[,1])))
  }
  if(method=="RF-regression"){
    model.obj <- try(randomForest( target ~ . , data=trainset))
  }

  #load("/scratch/dennis/projects/target_gene_prediction/LCL_model_100kb.RData")
  #gene.preds <- crm.gene.preds; test.samps <- 1:length(gene.preds);  dist.thres <- 0; model.var=c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator"); method <- "RF-regression";
  
  
  testset <- feature.select(gene.preds, crm.mat, test.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat, dist.thres)
  test.response <- as.numeric(testset[,1])
  test.params <- testset[,-1]
  
  colnames(test.params) <- c("mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
  class(test.params) <- "numeric"
  test.params <- test.params[,match(model.var[-1], colnames(test.params)), drop=FALSE]
  rownames(test.params) <- 1:dim(test.params)[1]
  target.pred <- NULL
  if(method == "svm"){
    target.pred <- predict(model.obj,test.params)
  }
  if(method == "gam"){
    target.pred <- predict(model.obj, as.data.frame(test.params) )
  }
  if(method == "RF-classification"){
    target.pred <- predict(model.obj,test.params, type="response")
    target.pred <- as.numeric(levels(target.pred)[target.pred])
  }
   if(method == "RF-regression"){
    target.pred <- predict(model.obj,test.params)
  }
  
  return(cbind(pred=target.pred, obs=test.response, crm.genes=rownames(testset)))
}

feature.select <- function( gene.preds, crm.mat, samps, faire.gene, gene.crms, gene.sim, insulators, score.mat, dist.thres=0){
train.params <- c()
train.params2 <- c()
train.params3 <- c()
train.params4 <- c()
train.params5 <- c()
train.params6 <- c()
train.response <- c()
gene.names <- c()
samps.tested <- c()
for(i in 1:length(samps)){
  crm.centre <- mean( as.numeric(crm.mat[samps[i], "end"]), as.numeric(crm.mat[samps[i], "start"]))
  if(length(gene.preds[[samps[i]]]) > 1){
    if( length(which(gene.preds[[samps[i]]][,"target_gene"] == crm.mat[samps[i],"Associated.Gene.Name"])) > 0){
      
      crm.gene.dist <- abs(gene.pos[match(gene.preds[[samps[i]]][,"target_gene"], gene.pos[,"Associated.Gene.Name"]),"Transcript.Start..bp."] - crm.centre )
      gene.filter <- which(crm.gene.dist >= dist.thres)
      
      train.params <- c(train.params, gene.preds[[samps[i]]][gene.filter,"mean.r2.fit"])
      
      train.params2 <- c(train.params2, crm.gene.dist[gene.filter])
      
      gene.match <- gene.crms[match( gene.preds[[samps[i]]][,"target_gene"], gene.crms[,"gene"]), 1:length(gm.tfs),drop=FALSE]
      gene.match[is.na(gene.match)] <- 0

      crm.dis <- apply(gene.match, 1, function(x){ crm.cooccur( as.numeric(x), score.mat)})
      train.params3 <- c(train.params3, crm.dis[gene.filter])
      
      gene.match <- match(gene.preds[[samps[i]]][,"target_gene"], faire.gene[,"gene"])
      gene.match[is.na(gene.match)] <- 0
      gene.match[which(gene.match != 0)] <- faire.gene[gene.match, "signalValue"]
      train.params4 <- c(train.params4, gene.match[gene.filter])

      if(length(gene.sim[[samps[i]]]) == 1 & is.na(gene.sim[[samps[i]]]))
        { gene.match <- rep(0, dim(gene.preds[[samps[i]]])[1])}
      else{
        sim.score <- apply(gene.sim[[samps[i]]], 2, mean)
        sim.score[is.na(sim.score)] <- 0
        gene.match <- rep(0, dim(gene.preds[[samps[i]]])[1])
        gene.match[match(  names(sim.score) , gene.preds[[samps[i]]][,"target_gene"])] <- sim.score
      }
      train.params5 <- c(train.params5, gene.match[gene.filter])
      
      
      gene.match <- insulators[[samps[i]]][match( gene.preds[[samps[i]]][,"target_gene"], names(insulators[[samps[i]]]) )]
      train.params6 <- c(train.params6, gene.match[gene.filter])
      
      response.tmp <- rep(0, dim(gene.preds[[samps[i]]])[1])
      response.tmp[which(gene.preds[[samps[i]]][,"target_gene"] == crm.mat[samps[i],"Associated.Gene.Name"])] <-  1
      train.response <- c(train.response, response.tmp[gene.filter])
      gene.names <- c(gene.names, gene.preds[[samps[i]]][,"target_gene"][gene.filter])

      if(length(gene.filter) > 0){ samps.tested <- c(samps.tested, paste(samps[i],":", gene.filter, sep=""))}
    }
  }
}
features.mat <- cbind(train.response, train.params, train.params2, train.params3, train.params4, train.params5, train.params6)
rownames(features.mat) <- samps.tested
return(features.mat)
}

# compute precision and recall
classify.eval <- function(preds, obs){
 tp <- length(which(!is.na(match(which(preds == 1), which(obs == 1)))))
 fp <- length(which(is.na(match(which(preds == 1), which(obs == 1)))))
 fn <- length(which(!is.na(match(which(preds == 0), which(obs == 1)))))
 tn <- length(which(!is.na(match(which(preds == 0), which(obs == 0)))))
 
 precision <- tp/(tp+fp)
 recall <- tp/(tp+fn)
 specificity <- tn/(fp+tn)
 if(!is.finite(precision)){ precision <- 1}
 if(!is.finite(recall)){ recall <- 1}
 if(!is.finite(specificity)){ specificity <- 1}
 
 return(c(precision = precision, recall = recall, specificity=specificity))
}

var.importance <- function(gene.preds, crm.mat, train.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat, model.var=c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator"), conditional=TRUE){
  trainset <- feature.select(gene.preds, crm.mat, train.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat)
  class(trainset) <- "numeric"
  colnames(trainset) <- c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
  rownames(trainset) <- 1:dim(trainset)[1]
  trainset <- trainset[,match(model.var, colnames(trainset))]
  model.obj <- try(cforest( target ~ . , data=as.data.frame(trainset)))
  return(varimp(model.obj, conditional=conditional))
}

# match crms
crm1 <- crm.mat.fibro
crm2 <- crm.match
crm.match <- c()
for(i in 1:nrow(crm1)){
  crm.match <- rbind(crm.match, crm2[which(crm2[,"start"] == crm1[i,"start"] & crm2[,"end"] == crm1[i,"end"]),])
}
crm.match <- crm.match[!duplicated(crm.match[,c("start", "end", "pattern")]),]
