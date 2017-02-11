# For linux
library(multicore,  lib.loc="/home/dennis/scratch/R/linux")
library(mgcv,  lib.loc="/home/dennis/scratch/R/linux")
library(doMC,  lib.loc="/home/dennis/scratch/R/linux")
library(lars)
library(e1071)
library(vegan)
library(gplots)
library(lattice)
library(randomForest)
library(tools)
library(ROCR)

load("/scratch/dennis/projects/TFBS_cell_specificity_hg18/K562_GM12878_whole_genome_500bp_hg18.RData")

source("/home/dennis/scratch/projects/expression_prediction/predict_expression_funcs.R")
source("/scratch/dennis/projects/insulators/distr_correlation.R")
source("/scratch/dennis/projects/MK_ChIPseq_analysis/cooccurrence_funcs.R")

# load chromosome conformation file
load.chr <- function(chr, file.type="obs"){
gm.chromatin.tmp <- as.matrix(read.table(paste("/scratch/dennis/projects/chromatin_states/HiC/GM06690_heatmaps/HIC_gm06690_", chr,"_", chr, "_100000_",file.type, ".txt", sep=""), sep="\t", header=TRUE, skip=1, fill=TRUE))
chr.pos <- gm.chromatin.tmp[,1]
gm.chromatin <- matrix(NA, nrow=length(chr.pos), ncol=length(chr.pos))
gm.chromatin <- gm.chromatin.tmp[1:length(chr.pos), 2:(length(chr.pos) + 1)]
chr.pos <- unlist(lapply(strsplit(chr.pos, split=":"), function(x){x[2]}))
colnames(gm.chromatin) <- chr.pos
rownames(gm.chromatin) <- chr.pos
class(gm.chromatin) <- "numeric" 
return(gm.chromatin)
}

# load all chromosomes
load.genome <- function(chromosomes, file.type="obs"){
registerDoMC()
chrom.assoc <- foreach(chr = chromosomes)%dopar%{
return(load.chr(chr, file.type))
}
names(chrom.assoc) <- chromosomes
return(chrom.assoc)
}

# filer for all interactions with CRMs
find.crm.inter <- function(crms, genome.pos){
crm.genome.pos <- matrix(NA, nrow=dim(crms)[1], ncol=dim(genome.pos[["chr1"]])[2])
for(i in 1:dim(crms)[1]){
 chrom.chr <- genome.pos[[crms[i,"chr"]]] 
 crm.centre <- mean(as.numeric(crms[i,"end"]) , as.numeric(crms[i,"start"]))
 crm.genome.pos[i,1:dim(chrom.chr)[2]] <-  chrom.chr[ceiling(crm.centre/100000),]
}
colnames(crm.genome.pos) <- colnames(genome.pos[["chr1"]])
return(crm.genome.pos)
}

# modify matrix for levelplot
prep.levelplot <- function(input.mat)
{
colnames(input.mat) <- c()
rownames(input.mat) <- c()
input.mat <- input.mat[,-which(apply(input.mat, 2, function(x){sum(x)}) == 0)]
input.mat <- input.mat[-which(apply(input.mat, 1, function(x){sum(x)}) == 0),]
return(input.mat)
}
  

#gm.combined.mat <- gm.combined.mat[-which(gm.combined.mat[,"chr"] == "chrY"),]


genome.pos <- load.genome(unique(gm.combined.mat[,"chr"]), file.type="obs")
crm.chrom.pos <- find.crm.inter(gm.combined.mat, genome.pos)
corr.mat <- crm.chrom.pos[which(gm.combined.mat[,"chr"] == "chr14"), 1:dim(genome.pos[["chr14"]])[2]]
# plot observed CRM reads along chromosome 14
corr.mat.2 <- corr.mat[,-which(apply(corr.mat, 2, function(x){sum(x)}) == 0)]
colnames(corr.mat.2) <- c()
rownames(corr.mat.2) <- c()
rgb.palette <- colorRampPalette(c("white", "red"), space="rgb")

png("/scratch/dennis/projects/chromatin_states/output/CRM_positioning_obs.png")
levelplot(corr.mat.2, col.regions=rgb.palette(120), cuts=10, at=c(0,1,500), xlab="CRMs", ylab="Nucleotide position (every 100kb)", colorkey=FALSE)
dev.off()

corr.mat.2 <- prep.levelplot(genome.pos[["chr14"]])
png("/scratch/dennis/projects/chromatin_states/output/chr14_positioning_obs.png")
levelplot(corr.mat.2, col.regions=rgb.palette(120), cuts=10, at=c(0,1,500), xlab="Nucleotide position", ylab="Nucleotide position (every 100kb)", colorkey=FALSE)
dev.off()

# plot observed/expected CRM reads along chromosome
genome.pos.exp <- load.genome(unique(gm.combined.mat[,"chr"]), file.type="exp")
crm.chrom.pos.exp <- find.crm.inter(gm.combined.mat, genome.pos.exp)
exp.mat <- crm.chrom.pos.exp[which(gm.combined.mat[,"chr"] == "chr14"), 1:dim(genome.pos.exp[["chr14"]])[2]]

enrich.mat <- corr.mat/exp.mat
ennrich.mat.2 <- prep.levelplot(enrich.mat)
# plot enrichment for CRMs
png("/scratch/dennis/projects/chromatin_states/output/CRM_positioning_enrich.png")
levelplot(enrich.mat.2, col.regions=rgb.palette(120), pretty=TRUE,at=c(0,2,50), xlab="CRMs", ylab="Nucleotide position (every 100kb)", colorkey=FALSE)
dev.off()

enrich.mat.2 <- prep.levelplot(genome.pos[["chr14"]] / genome.pos.exp[["chr14"]])
# plot enrichment for whole chromosome
png("/scratch/dennis/projects/chromatin_states/output/chr14_positioning_enrich.png")
levelplot(enrich.mat.2, col.regions=rgb.palette(120), cuts=10, at=c(0,2,50), xlab="Nucleotide position", ylab="Nucleotide position (every 100kb)", colorkey=FALSE)
dev.off()

load("/scratch/dennis/projects/chromatin_states/LCL_whole_genome_500bp_hg18.RData")



# find chromatin associations for target genes
genome.pos <- load.genome(unique(gm.combined.mat[,"chr"]), file.type="obs")
crm.chrom.pos <- find.crm.inter(crm.mat[,1:15], genome.pos)

genome.pos.exp <- load.genome(unique(gm.combined.mat[,"chr"]), file.type="exp")
crm.chrom.pos.exp <- find.crm.inter(crm.mat[,1:15], genome.pos.exp)

enrich.mat <- crm.chrom.pos/crm.chrom.pos.exp

registerDoMC()
pred.enrich <- foreach(i = 1:dim(crm.mat)[1], .combine=rbind)%dopar%{
  pred.enrich.tmp <- c()
crm.target <- crm.gene.preds[[1]]
for(j in 1:dim(crm.target)[1])
  {
    pred.enrich.tmp <- rbind(pred.enrich.tmp, c(crm.target[j,"mean.r2.fit"], crm.target[j,"mean.r2.pred"], enrich.mat[i, ceiling(gene.pos[gene.pos[,"Associated.Gene.Name"] == crm.target[j, "target_gene"], "Transcript.Start..bp."][1] / 100000)]))
  }
  return(pred.enrich.tmp)
}

class(pred.enrich) <- "numeric"

pred.enrich.1 <- pred.enrich
plot(pred.enrich.1[,2], pred.enrich.1[,3])
plot(loess.smooth( pred.enrich.1[,1], pred.enrich.1[,3], span=8.5/10, degree=1), type="l")

pred.enrich.2 <- pred.enrich
plot(pred.enrich.2[,2], pred.enrich.2[,3])
plot(loess.smooth( pred.enrich.2[,1], pred.enrich.2[,3], span=7/10, degree=1), type="l")

load("/scratch/dennis/projects/chromatin_states/LCL_whole_genome_500bp_crm_predictions.RData")

#########################################
# chromosome interactions for CRM pairs #
#########################################

# find genes with TF binding sites
crm.genes <- map.gene.crm(gm.combined.mat, gene.pos)

crm.crm <- list()
for(i in 1:nrow(crm.mat)){
 crm.chr <- crm.genes[which(crm.genes[,"chr"] == crm.mat[i,"chr"]),]
 crm.crm[[i]] <- crm.chr[which(abs(as.numeric(crm.chr[,"start"]) - as.numeric(crm.mat[i,"start"])) > 5000000),]
}

########################
# CRM-gene interactions#
########################

crm.mat <- samps.tested[as.numeric(samps.tested[,"pred"]) == 1 & as.numeric(samps.tested[,"obs"]) == 0, ]
genome.pos <- load.genome(unique(crm.mat[,"chr"]), file.type="obs")
genome.pos.exp <- load.genome(unique(crm.mat[,"chr"]), file.type="exp")
crm.inter.obs <- find.crm.inter(crm.mat, genome.pos)
crm.inter.exp <- find.crm.inter(crm.mat, genome.pos.exp)
enrich.mat <- crm.inter.obs

crm.gene.inter <- c()
for(i in 1:nrow(crm.mat)){
  crm.gene.inter <- c(crm.gene.inter, enrich.mat[i,gene.pos[gene.pos[,"Associated.Gene.Name"] == crm.mat[i, "target_gene"], "Transcript.Start..bp."]/100000])
}
mean(crm.gene.inter)


# extract features for each pairwise CRM
########################################
genome.pos <- load.genome(unique(crm.mat[,"chr"]), file.type="obs")
genome.pos.exp <- load.genome(unique(crm.mat[,"chr"]), file.type="exp")

crm.inter.obs <- find.crm.inter(crm.mat, genome.pos)
crm.inter.exp <- find.crm.inter(crm.mat, genome.pos.exp)

enrich.mat <- crm.inter.obs

# gene co-expression
crm.gene.preds <- predict.targets(crm.mat, training.express, training.express, all.samples, gene.pos=gene.pos, max.distance=300000000, regress.method="RF")

# filter for genes with chromatin interactions
gene.preds <- crm.gene.preds
for(i in 1:nrow(crm.mat)){
  if(!is.na(gene.preds[[i]])){
    tmp <- merge(gene.preds[[i]], crm.crm[[i]], by.x="target_gene", by.y="gene")
    tmp <- tmp[!duplicated(tmp[,c("target_gene", "target_expression_ind")]),1:ncol(gene.preds[[i]])]

    sample.size <- 20
    if(nrow(tmp) < sample.size){sample.size <- nrow(tmp)}
    gene.preds[[i]] <- tmp[sample(nrow(tmp), sample.size),]
  }
}


# Get GO similarity for TF-gene pairs
library(SparseM, lib.loc="/home/dennis/scratch/R/linux")
library(GOSim)
tf.entrez <- as.numeric(training.express[match(names(gm.tfs), training.express[,"Gene"]),"NAME"])
names(tf.entrez) <- names(gm.tfs)
registerDoMC()
psnice(value=1)
gene.sim <- foreach(i = 1:length(gene.preds))%dopar%{
  if(!is.na(gene.preds[[i]])){
    matched.genes <- as.numeric(training.express[as.numeric(gene.preds[[i]][,"target_expression_ind"]), "NAME"])
    names(matched.genes) <- gene.preds[[i]][,"target_gene"]
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
insulators <- foreach(i = 1:length(gene.preds))%dopar%{
  if(!is.na(gene.preds[[i]])){
    matched.genes <- gene.pos[match(gene.preds[[i]][,"target_gene"], gene.pos[,"Associated.Gene.Name"]),, drop=FALSE]
    ctcf.chr <- gm.ctcf[gm.ctcf[,"chr"] == crm.mat[i,"chr"],]
    ctcf.upstream <- ctcf.chr[which(as.numeric(ctcf.chr[,"start"]) > as.numeric(crm.mat[i, "Coordinate_HG18"])),, drop=FALSE]
    ctcf.downstream <- ctcf.chr[which(as.numeric(ctcf.chr[,"start"]) < as.numeric(crm.mat[i, "Coordinate_HG18"])),, drop=FALSE]
    gene.ctcf <- rep(0, dim(gene.preds[[i]])[1])
    names(gene.ctcf) <- gene.preds[[i]][,"target_gene"]
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


# extract features
for(i in 1:nrow(crm.mat)){
  gene.preds[[i]] <- as.matrix(gene.preds[[i]])
  crm.mat[i,"Associated.Gene.Name"] <- gene.preds[[i]][1,"target_gene"]
}

train.set <- feature.select(gene.preds, crm.mat, 1:length(gene.preds), faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat, dist.thres=0)

load("/scratch/dennis/projects/chromatin_states/chromatin_interaction_prediction.RData")


##############
# Predict HiC#
##############

# set response as chromatin interaction enrichment
for(i in 1:nrow(train.set)){
 crm.ind <- as.numeric(unlist(strsplit(rownames(train.set)[i], split=":"))[1])
 train.set[i,"train.response"] <- enrich.mat[crm.ind, ceiling(gene.pos[gene.pos[,"Associated.Gene.Name"] == train.set[i,"gene.names"],"Transcript.Start..bp."]/ 100000)[1]]
# if(log(as.numeric(train.set[i, "train.response"])) > mean(na.omit(enrich.mat[crm.ind, ]))){
#   train.set[i, "train.response"] <- 1
# }else{
#    train.set[i, "train.response"] <- 0
# }
}
colnames(train.set) <- c("gene.names", "read enrichment", "gene co-expression", "distance", "crm pattern", "faire signal", "GO similarity", "insulator")
train.set <- train.set[!is.na(train.set[,"read enrichment"]),]

# check features
plot( train.set[,"faire signal"], train.set[,"read enrichment"])

png("/scratch/dennis/projects/target_gene_prediction/output/features_HiC.png" )
par(mfrow=c(3,2), mar=c(5.1,5.1,4.1,2.1))
plot(smooth.spline(as.numeric(train.set[,"distance"]), as.numeric(train.set[,"read enrichment"]), df=3), type="l", ylab="interaction frequency", xlab="distance (kb)", cex.lab=1.5, cex.axis=1.5)
plot(smooth.spline(as.numeric(train.set[,"gene co-expression"]), as.numeric(train.set[,"read enrichment"]), df=3), type="l", ylab="interaction frequency", xlab=expression(R^2), cex.lab=1.5, cex.axis=1.5)
plot(smooth.spline(as.numeric(train.set[,"crm pattern"]), as.numeric(train.set[,"read enrichment"]), df=3), type="l", ylab="interaction frequency", xlab="co-occurrence score", cex.lab=1.5, cex.axis=1.5)
plot(smooth.spline(as.numeric(train.set[,"faire signal"]), as.numeric(train.set[,"read enrichment"]), df=3), type="l", ylab="interaction frequency", xlab="FAIRE signal", cex.lab=1.5, cex.axis=1.5)
plot(smooth.spline(as.numeric(train.set[,"GO similarity"]), as.numeric(train.set[,"read enrichment"]), df=3), type="l", ylab="interaction frequency", xlab="GO similarity", cex.lab=1.5, cex.axis=1.5)
plot(smooth.spline(as.numeric(train.set[,"insulator"]), as.numeric(train.set[,"read enrichment"]), df=3), type="l", ylab="interaction frequency", xlab="CTCF signal", cex.lab=1.5, cex.axis=1.5)
dev.off()


#obs.tmp <- genome.pos[[crm.mat[1, "chr"]]] [ ceiling(as.numeric(crm.mat[1, "start"])/100000),    ceiling(gene.pos[gene.pos[,"Associated.Gene.Name"] == "C12orf32","Transcript.Start..bp."]/100000) ]
#exp.tmp <- genome.pos.exp[[crm.mat[1, "chr"]]] [ ceiling(as.numeric(crm.mat[1, "start"])/100000),    ceiling(gene.pos[gene.pos[,"Associated.Gene.Name"] == "C12orf32","Transcript.Start..bp."]/100000) ]

# predict chromatin interactions
predict.chrom.interactions <- function(response.var, feature.set, train.samps, model.var=c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator"), classification=FALSE){
  trainset <- cbind(response.var[train.samps], feature.set[train.samps,])
  class(trainset) <- "numeric"
  colnames(trainset) <- c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
  rownames(trainset) <- 1:dim(trainset)[1]
  trainset <- trainset[,match(model.var, colnames(trainset))]

  model.obj <- NULL
  if(classification == TRUE){
    model.obj <- try(randomForest( trainset[,-1, drop=FALSE] , as.factor(trainset[,1])))
  }else{
    model.obj <- try(randomForest( target ~ . , data=trainset, ntree=100))
  }
  
  test.params <-  feature.set[-train.samps,]
  test.response <- as.numeric(response.var[-train.samps])
  colnames(test.params) <- c("mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
  class(test.params) <- "numeric"
  gene.distance <- test.params[,"distance"]
  test.params <- test.params[,match(model.var[-1], colnames(test.params)), drop=FALSE]
  rownames(test.params) <- 1:dim(test.params)[1]

  target.pred <- NULL
  if(classification == TRUE){
    target.pred <- predict(model.obj,test.params, type="response")
    target.pred <- as.numeric(levels(target.pred)[target.pred])
  }else{
     target.pred <- predict(model.obj,test.params)
  }
  
  return(cbind(pred=target.pred, obs=test.response, dist=gene.distance))
}

predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=1:5000, model.var=c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator"), classification=FALSE)
 cor(predict.results[,"pred"], predict.results[,"obs"])^2

crm.index <- 1:nrow(train.set)
cv.fold <- 10
registerDoMC()
cv.acc <- foreach (cv.trial = 1:cv.fold) %dopar%{    
  
  test.genes <- crm.index[ceiling(length(crm.index)*((cv.trial-1)/cv.fold)+1):ceiling(length(crm.index)*(cv.trial/cv.fold))]  # k-fold cv
  train.genes <- crm.index[-match(test.genes, crm.index)]

  r.sqreds <- c()
  
  predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target", "mean.r2.fit", "crm.pattern", "faire.signal", "GO.similarity", "insulator"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)

  return(cbind(predict.results[,"pred"], predict.results[,"obs"], predict.results[,"dist"]))
}
  
  print(paste(cv.trial, "1"))
  
  predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target",  "distance"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)

 print(paste(cv.trial, "2"))

  predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target", "mean.r2.fit"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)
print(paste(cv.trial, "3"))
  
  predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target",  "crm.pattern"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)

  print(paste(cv.trial, "4"))
  
  predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target", "faire.signal"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)

print(paste(cv.trial, "5"))
  
  predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target", "GO.similarity"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)

  print(paste(cv.trial, "6"))
  
   predict.results <- predict.chrom.interactions(train.set[,"read enrichment"], train.set[,3:8], train.samps=train.genes, model.var=c("target", "insulator"), classification=FALSE)
  r.sqreds <- c(r.sqreds, cor(predict.results[,"pred"], predict.results[,"obs"])^2)
  

  print(paste(cv.trial, "7"))
  return(r.sqreds)
}


png("/scratch/dennis/projects/target_gene_prediction/output/HiC_pred_vs_obs.png" )
par( mar=c(5.1,5.1,4.1,2.1))
smoothScatter(preds ~ obs, nrpoints=Inf, bandwidth=1, xlab="observed interaction frequency", ylab="predicted interaction frequency", cex.lab=1.5, cex.axis=1.5)
abline(lm(preds ~ obs), col="red")
dev.off()


target.enrich <- test.response[which(target.pred > 0)]
target.dist <- gene.distance[which(target.pred > 0)]

nontarget.enrich <- test.response[which(target.pred == 0)]
nontarget.dist <- gene.distance[which(target.pred == 0)]
nontarget.enrich <- nontarget.enrich[sample(which(nontarget.dist < max(target.dist)), length(target.dist))]

png("/scratch/dennis/projects/target_gene_prediction/output/HiC_eqtl_targets.png" )
par( mar=c(5.1,5.1,4.1,2.1))
t.test(target.enrich, nontarget.enrich, alternative="greater")
boxplot(list(Targets = target.enrich, Non.targets = nontarget.enrich), ylab="interaction frequency", cex.lab=1.5, cex.axis=1.5)
dev.off()

r.sqreds <- c()
for(i in 1:length(cv.acc)){
preds <- cv.acc[[i]][,1]
obs <- cv.acc[[i]][,2]
dists <- cv.acc[[i]][,3]

r.sqreds <- c(r.sqreds, cor(preds, obs)^2)

}


