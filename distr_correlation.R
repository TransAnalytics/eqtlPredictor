library(IRanges)
library(flexmix)
library(fields)

procPeakFile <- function(raw.data){
for (i in 1:length(raw.data)) {
  raw.data[[i]]$chr <- raw.data[[i]][,1]
  raw.data[[i]]$start <- raw.data[[i]][,2]
  raw.data[[i]]$end <- raw.data[[i]][,3]
  raw.data[[i]]$center <- (raw.data[[i]]$start + raw.data[[i]]$end)/2
}
return(raw.data)
}

#
# Plot distributions
#

#
# function to plot distribution of binding sites
# inputs: chips - chipseq data
#         ch - chromosome
#
plotBSdistr <- function(chips, ch, intensity=FALSE){
pos.lst <- list()
all.pos <- c()

for (chip in chips) {
  chip <- chip[chip$chr == ch,]
  center <- chip[,"center"]
  if(intensity==TRUE){
   center <- cbind(center, chip[,5])
   colnames(center)[2] <- "intensity"
  }
  pos.lst <- c(pos.lst,list(center))
 
  all.pos <- c(all.pos,chip[,"center"])
}
#png("gata_nfe2_myc_profiles.png")
#colors <- rainbow(length(pos.lst))
#plot(density(pos.lst[[1]],from=0,to=max(all.pos),bw=100, n=100),  ylim=c(0, 4.0*10^-5), main="Density distributions of binding sites", col=colors[1]);# rug(all.pos)

#for (i in (2:(length(pos.lst)))) {
#  pos <- pos.lst[[i]]
#  lines(density(pos,from=0,to=max(pos), bw=100, n=10000),col=colors[i])
#  rug(pos,col=colors[i])
# }
#legend("topleft", names(chips), cex=0.8, col=colors[1:length(pos.lst)], lty=1)
#dev.off()
 pos.lst <- c(pos.lst, list(sort(all.pos)))
 names(pos.lst) <- c(names(chips), "all")
 return(pos.lst)
}


#
# find correlation
#

# function to compute KL divergence between two distributions
# input: list of histone mods sites
KL.dist <- function(h1, bw=10){
  distr.mat <- c()
  for (i in 1:length(h1)){
    d1 <- density(h1[[i]], from=0, to=max(h1[[i]]), bw=bw)
    distr.mat <- cbind(distr.mat, d1$y)
  }
  distr.dist <- KLdiv(distr.mat, method="discrete")
  colnames(distr.dist) <- names(h1)

  rownames(distr.dist) <- names(h1)
  return(distr.dist)
}


#
#  mean of distances from TF1 to nearest TF2 binding site 
nearest.dist <- function(bs1, bs2){
  dists <- c()
  if(length(bs1) > length(bs2)){
    dists <- sapply(bs2, function(x, y)min(abs(y - x)), bs1)
  } else{
    dists <- sapply(bs1, function(x, y)min(abs(y - x)), bs2)
  }
#  return(mean(dists))
#  return(median(dists))
  return(dists)
}

#
# generates distance matrix based on given distance measure
#
get.dist.mat <- function(h1, method){
  dist.mat <- matrix(NA, nrow=length(h1), ncol=length(h1))
  for(i in 1:length(h1)){
    for(j in i:length(h1)){
      if(method=="nearest"){
        dist.mat[i, j] <- nearest.dist(h1[[i]], h1[[j]])
      }
      if(method=="KS"){        
       dist.mat[i, j] <- ks.test(h1[[i]], h1[[j]], alternative="two.sided", exact=TRUE)$statistic 
      }
    }
  }
  if(method == "KL"){
   dist.mat <- KL.dist(h1)
  }
  colnames(dist.mat) <- names(h1)
  rownames(dist.mat) <- names(h1)
  return(dist.mat)
}
  
# function to resampling (1000x) from pooled sites to generate KL distances between two random histone mods
resample <- function(ha, h1, h2, method, bw=10){

  distr.dists <- c()
  for (i in 1:1000){
    h1.rand <- ha[sample(length(ha), length(h1), replace=F)]
    h2.rand <- ha[sample(length(ha), length(h2), replace=F)]
    if(method=="KL"){
      distr.dists <- c(distr.dists, KL.dist(list(h1.rand, h2.rand), bw=bw)[1,2])
    }
    if(method=="KS"){
      distr.dists <- c(distr.dists, get.dist.mat(list(h1.rand, h2.rand), method="KS")[1,2])
    }
    if(method=="nearest"){
      distr.dists <- c(distr.dists, nearest.dist(h1.rand, h2.rand))
    }
  }
  return(distr.dists)
}


# test KL distance measures
testKL.dist <- function(pos.lst, hist.dists=NULL, bw=10){
  if(is.null(hist.dists)){
    hist.dists <- KL.dist(pos.lst[-length(pos.lst)], bw=bw)
  }
#samp.dists.KL <- c() # for testing
corr.diff <- matrix(NA, nrow=dim(hist.dists)[1], ncol=dim(hist.dists)[2])
colnames(corr.diff) <- colnames(hist.dists)
rownames(corr.diff) <- rownames(hist.dists)
for(i in 1:dim(hist.dists)[1]){
 for(j in 1:dim(hist.dists)[1])
   {
     h1 <- pos.lst[[rownames(hist.dists)[i]]]
     h2 <- pos.lst[[rownames(hist.dists)[j]]]

     samp.dists <- resample(pos.lst$all, h1, h2, method = "KL", bw=bw)
#     samp.dists.KL <- c(samp.dists.KL, samp.dists) # for testing
#     cdf <- ecdf(samp.dists)
  #  corr.diff[i,j] <- cdf(hist.dists[i,j])
     corr.diff[i,j] <- length(which(samp.dists < hist.dists[i,j])) / length(samp.dists)
   }
}
return(list(corr.diff, hist.dists))
}


#
# test KS distance measure
# DEPRECATED
testKS.dist <- function(pos.lst, hist.dists=NULL){
  if(is.null(hist.dists)){
    hist.dists <- KS.dist(pos.lst[-length(pos.lst)])
  }
samp.dists.KS <- c() # for testing
corr.diff <- matrix(NA, nrow=dim(hist.dists)[1], ncol=dim(hist.dists)[2])
colnames(corr.diff) <- colnames(hist.dists)
rownames(corr.diff) <- rownames(hist.dists)
for(i in 1:(dim(hist.dists)[1] -1)){
 for(j in i:dim(hist.dists)[1])
   {
     h1 <- pos.lst[[rownames(hist.dists)[i]]]
     h2 <- pos.lst[[rownames(hist.dists)[j]]]

     samp.dists <- resample(pos.lst$all, h1, h2, method="KS")
     samp.dists.KS <- c(samp.dists.KS, samp.dists) # for testing
     cdf <- ecdf(samp.dists)
     corr.diff[i,j] <- cdf(hist.dists[i,j])
   }
}
return(corr.diff, hist.dists)
}

test.dist <- function(pos.lst, method, hist.dists=NULL){

if(is.null(hist.dists)){hist.dists <- get.dist.mat(pos.lst[-length(pos.lst)], method)}

samp.dists.lst <- c() # for testing
corr.diff <- matrix(NA, nrow=dim(hist.dists)[1], ncol=dim(hist.dists)[2])
colnames(corr.diff) <- colnames(hist.dists)
rownames(corr.diff) <- rownames(hist.dists)

comp.time <- system.time(for(i in 1:(dim(hist.dists)[1] -1)){
 for(j in i:dim(hist.dists)[1])
   {
     h1 <- pos.lst[[rownames(hist.dists)[i]]]
     h2 <- pos.lst[[rownames(hist.dists)[j]]]

     samp.dists <- resample(pos.lst$all, h1, h2, method=method)
     samp.dists.lst <- c(samp.dists.lst, samp.dists) # for testing
     
     corr.diff[i,j] <- length(which(samp.dists < hist.dists[i,j])) / length(samp.dists)
   }
})

return(corr.diff, hist.dists, samp.dists.lst)
}

#
# sum KL distances
#
add.KL <- function(tf1, tf2, KL.dists){
 
}

          
#
# Profile of KL distance as a function of position
#
profile.KL <- function(sites.lst, window=1000000){
  tf1 <- sites.lst[[1]]
  tf2 <- sites.lst[[2]]
  all.pos <- sites.lst$all
  KL.dists <- c()
  for(i in seq(1, max(all.pos), by=100000)){
    if(i > max(tf1) - window  || i > max(tf2) - window ){ break}
    tf1.win <- tf1[which(tf1 > i & tf1 < i+window)]
    tf2.win <- tf2[which(tf2 > i & tf2 < i+window)]
    all.win <- all.pos[which(all.pos > i & all.pos < i+window)]
 #   d1 <- density(tf1.win, from=0, to=max(tf1.win), bw=10)
 #   d2 <- density(tf2.win, from=0, to=max(tf2.win), bw=10)
 #   distr.mat <- cbind(d1$y, d2$y)
 #   dist.mat <- KLdiv(distr.mat)
    dist.score = 1
    dist.sig = 1
    if(length(tf1.win) != 0 & length(tf2.win) != 0){
    pos.win <- list(tf1.win, tf2.win, all.win)
    names(pos.win) <- names(sites.lst)
  #  distr.stat <- test.dist(pos.win, method="KS") # try KS distance
    distr.stat <- testKL.dist(pos.win)
    dist.score <- distr.stat$hist.dists[1,2]
    dist.sig <-  distr.stat$corr.diff[1,2]
  }
    KL.dists <- rbind(KL.dists, c(i+(window/2),dist.score , dist.sig))
  }
  colnames(KL.dists) <- c("position", "KL_dist", "significance")
  return(KL.dists)
}


#
# combines chromosomes
#
combine.chr <- function(chips, buffer.size=10000){

chromosomes <- c("chr1", "chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
buffer.lst <- c(0)
for(chr in chromosomes){
    buffer.lst <-  c(buffer.lst, max(plotBSdistr(chips , chr)$all) + buffer.size + buffer.lst[length(buffer.lst)])
  }
  
for(i in 1:length(chips)){
  chip <- chips[[i]]
  chip$chr <-  substr(chip$chr, 4, 100)
  chip$chr[chip$chr=="X"] <- 23 #recode X as 23 for sorting
  chip$chr[chip$chr=="Y"] <- 24 
  chip = chip[order(as.numeric(chip$chr),as.numeric(chip$start)),]
  chip$chr <- paste("chr", chip$chr, sep="")
 
  for(j in 1:length(unique(chip$chr))){
    chr <- unique(chip$chr)[j]
    chip[chip$chr == chr,]$start <-  chip[chip$chr == chr,]$start + buffer.lst[j]
    chip[chip$chr == chr,]$end <-  chip[chip$chr == chr,]$end  + buffer.lst[j]
    chip[chip$chr == chr,]$center <-  chip[chip$chr == chr,]$center + buffer.lst[j]
  }

  chip$chr <- "combined"
  chips[[i]] <- chip
}
return(chips)
}


next.tss.chr <- function(centers,tss) {
  pts <- c(tss,centers)
  is.center <- c(rep(FALSE,length(tss)),rep(TRUE,length(centers)))
  inds <- c(1:length(tss),1:length(centers))
  o <- order(pts)
  pts.o <- pts[o]
  inds.o <- inds[o]
  pos.c <- which(is.center[o])
  inds.c <- inds[o][pos.c]
  dist <- rep(NA,length(centers))
  for (i in inds.c) {
    peak <- centers[i]
    pos <- pos.c[i]
    j <- pos - 1
    while (j >= 1 && is.center[j]) j <- j - 1
    if (j >= 1)
      low.gene <- pts.o[j]
    else
      low.gene <- -Inf
    j <- pos + 1
    while (j <= length(pts.o) && is.center[j]) j <- j + 1
    if (j <= length(pts.o))
      high.gene <- pts.o[j]
    else
      high.gene <- Inf
    dist[i] <- min(peak - low.gene, high.gene - peak)
  }
  dist
}


getBlockSigs <- function(blk, mod.lst){
 for(i in 1:length(mod.lst)){
   mod.lst[[i]] <- mod.lst[[i]][which(mod.lst[[i]] >= (blk[1]-10000) & mod.lst[[i]] <=( blk[2] + 10000))]
 }
 return(mod.lst)
}


getCTCF.blocks <- function(enh, chr, random=FALSE){
 ctcf.lst.chr <- plotBSdistr(enh, chr)
   
########################
# Testing random blocks#  
########################
 if(random == TRUE){
 ctcf.lst.chr$enh <- sample(min(ctcf.lst.chr$enh):max(ctcf.lst.chr$enh), size=length(ctcf.lst.chr$enh) )
 ctcf.lst.chr$blks <- sample(min(ctcf.lst.chr$blks):max(ctcf.lst.chr$blks), size=length(ctcf.lst.chr$blks) )
 ctcf.lst.chr$all <- c(ctcf.lst.chr$enh, ctcf.lst.chr$blks)
 }
#######################
 
 ctcf.sites <- sort(ctcf.lst.chr$blks)
 ctcf.lst.chr$enh <- sort(ctcf.lst.chr$enh)
 blk.lst <- c()

 for(i in 1:length(ctcf.lst.chr$enh)){
  e1 <- ctcf.lst.chr$enh[i]
  if(e1 > max(ctcf.sites)){ break}
  block.end <- which(ctcf.sites >= e1)[1]
  if(block.end > 1){
  block.start <- block.end - 1
  block <- c(ctcf.sites[block.start], ctcf.sites[block.end])
  blk.lst <- rbind(blk.lst, block)
  ctcf.sites <- ctcf.sites[block.end:length(ctcf.sites)]
  }
}
blk.lst <- unique(blk.lst)
blk.lst <- IRanges(start=blk.lst[,1], end=blk.lst[,2])
 return(blk.lst)
}
  
# gets all ChIP positions in CTCF blocks around predicted enhancers
# input: chips - ChIP profiles
#        blk.lst - IRange of ctcf blocks
#        chr - chromosome of interest
filter.profiles.blocks <- function(chips, blk.lst, chr){

# get ChIP profiles only at CTCF blocks  
blk.chips <- NULL

for(chip in chips){
 
  inds <- chip$chr == chr
  chip.lst <- IRanges(chip[inds, "start"], chip[inds,"end"] )
  blk.profile <- intersect(chip.lst, blk.lst)
  blk.chips <- c(blk.chips, list((end(blk.profile) + start(blk.profile))/2))  
}

names(blk.chips) <- names(chips)
 
chips <- plotBSdistr(chips, chr) 
blk.chips <- c(blk.chips, list(chips$all))
names(blk.chips)[length(blk.chips)] <-  "all"

return(blk.chips)
}


#####################################
#Associate with Protein Interactions#
##################################### 


inter.dist.table <- function(KL_dists, KS_dists, KL_sig, KS_sig, inter){
# make symmetric
KS_dists[is.na(KS_dists)] <- 0
KS_dists<- KS_dists + t(KS_dists)
KL_sig[is.na(KL_sig)] <- 0
KL_sig<- KL_sig+ t(KL_sig)
KS_sig[is.na(KS_sig)] <- 0
KS_sig<- KS_sig+ t(KS_sig) 
  
inter.dists <- cbind(inter, matrix(NA, nrow=dim(inter)[1], ncol=4))
colnames(inter.dists) <- c(colnames(inter), "KL_dists", "KL_sig", "KS_dists", "KS_sig")
tfs <- colnames(KS_dists)
for(i in 1:dim(inter.dists)[1]){
 tf1 <- tolower(as.character(inter.dists[i,"node1"]))
 tf2 <- tolower(as.character(inter.dists[i,"node2"]))
  if(length(c(which(tfs == tf1), which(tfs == tf2))) > 1){ 
    inter.dists[i,"KL_dists"] <- KL_dists[tf1,tf2]
    inter.dists[i,"KL_sig"] <- KL_sig[tf1,tf2]
    inter.dists[i,"KS_dists"] <- KS_dists[tf1,tf2]
    inter.dists[i,"KS_sig"] <- KS_sig[tf1,tf2]
  }
}
return(inter.dists)
}

# generate distance matrix from interaction summary data
getInterMat <- function(inter, evidence){
   proteins <- unique(c(as.character(unique(inter$node1)),  as.character(unique(inter$node2))))
   inter.mat <- matrix(0, ncol= length(proteins), nrow=length(proteins), dimnames=list(proteins, proteins))
   for(i in 1:dim(inter)[1]){
     p1 <- as.character(inter$node1[i])
     p2 <- as.character(inter$node2[i])
     inter.mat[ p1, p2] <- inter[i, evidence]
     inter.mat[ p2, p1] <- inter[i, evidence]
   }
   diag(inter.mat) <- 1
   return(inter.mat)
}


# finds indirect relationships between two proteins
# input: matrix of protein-protein interactions
# returns matrix of indirect relationships and combined scores 
find.indirect <- function(inter){
  relations.mat <- c()
  for(i in 1:dim(inter)[2]){
    for(j in i:dim(inter)[1]){
      if(inter[j,i] != 0 & i != j){
        for(k in 1:i){
          if(inter[j, k] != 0 & j != k & i != k){
            relations.mat <- rbind(relations.mat,  c(colnames(inter)[i], colnames(inter)[j], colnames(inter)[k], inter[j, i]*inter[j, k]))
          }
        }

      }
    }
  }
colnames(relations.mat) <- c("node1", "intermediate", "node2", "combined_score")
return(relations.mat)
}

# tf2 is usually the intermediate
sig.relationships <- function(tf1, tf2, tf3, pos.lst, chr){ 
tf1.tf2.dist <- profile.KL(c(pos.lst[tf1], pos.lst[tf2], pos.lst["all"]))
tf2.tf3.dist <- profile.KL(c(pos.lst[tf2], pos.lst[tf3], pos.lst["all"]))
tf3.tf1.dist <- profile.KL(c(pos.lst[tf1], pos.lst[tf3], pos.lst["all"]))

write.table(tf1.tf2.dist, file=paste(tf1, "_", tf2, "_", chr, ".txt", sep=""), sep="\t") 
write.table(tf2.tf3.dist, file=paste(tf2, "_", tf3, "_", chr, ".txt", sep=""), sep="\t")
write.table(tf3.tf1.dist, file=paste(tf1, "_", tf3, "_", chr, ".txt", sep=""), sep="\t") 

sig.regions <- tf3.tf1.dist[which(tf3.tf1.dist[, "significance"] < 0.000036),]
sig.overlap.lst <- rep(0, dim(sig.regions)[1])
names(sig.overlap.lst) <- sig.regions[,"position" ]

for(i in 1:length(sig.regions[,"position"])){
  sig1 <- tf2.tf3.dist[which((tf2.tf3.dist[,"position"] < sig.regions[i,"position"]+50000) & ( tf2.tf3.dist[,"position"] > sig.regions[i,"position"]-50000)), "significance"]
  sig2 <- tf1.tf2.dist[which((tf1.tf2.dist[,"position"] < sig.regions[i,"position"]+50000) & ( tf1.tf2.dist[,"position"] > sig.regions[i,"position"]-50000)), "significance"]
  sig.overlap.lst[i] <- length(which(sig1 < 0.000036)) + length(which(sig2 < 0.000036))
}
return(sig.overlap.lst)
}


#inter.combo.score <- cbind(as.character(inter$tfs$node1), as.character(inter$tfs$node2), inter$tfs$combined_score)
  
##########################################
#Integrate gene expression in CTCF blocks#
##########################################

library(Biobase)
library(GEOquery)
library(hgu95av2.db)
library(AnnotationDbi)

load.geo <- function(geofile){
geo.data <- getGEO(filename=geofile, destdir="/scratch/dennis/projects/insulators/Gene_expression/", AnnotGPL=TRUE)
geo.plat <- getGEO(Meta(geo.data)$platform, destdir="/scratch/dennis/projects/insulators/Gene_expression/")
MA <- GDS2MA(geo.data, GPL = geo.plat)
expr.set <- MA$M
expr.set <- cbind( as.character(MA$genes$Gene.Symbol), as.character(MA$genes$Gene.Title), expr.set)
colnames(expr.set)[1:2] <- c('Gene.Symbol', 'Gene.Title')
rownames(expr.set) <- as.character(MA$genes$ID)
return(expr.set)
}

load.anno <- function(chr){
probes.chr <- mget(chr, revmap(hgu95av2CHR))
probes.loc <- abs(unlist(mget(unlist(probes.chr), env=hgu95av2CHRLOC)))
probes.loc <- sort(probes.loc[!is.na(probes.loc)])
return(probes.loc)
}


# gets TSS range  
get.TSS.range <- function(chr){
require(BSgenome.Hsapiens.UCSC.hg19)
hs.chr <- Hsapiens[[chr]]

}

# filters 
filter.pos <- function(pos.lst, filter.ranges){


}

