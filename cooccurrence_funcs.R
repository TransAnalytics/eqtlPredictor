procPeakFile <- function(raw.data, type="bed"){
  if(type == "bed"){
    for (i in 1:length(raw.data)) {
      raw.data[[i]]$chr <- raw.data[[i]][,1]
      raw.data[[i]]$start <- raw.data[[i]][,2]
      raw.data[[i]]$end <- raw.data[[i]][,3]
      raw.data[[i]]$center <- (raw.data[[i]]$start + raw.data[[i]]$end)/2
    }
  }
  if(type == "ucsc"){
     for (i in 1:length(raw.data)) {
       raw.data[[i]]$chr <- raw.data[[i]][,"chrom"]
       raw.data[[i]]$start <- raw.data[[i]][,"chromStart"]
       raw.data[[i]]$end <- raw.data[[i]][,"chromEnd"]
       raw.data[[i]]$center <- (raw.data[[i]]$start + raw.data[[i]]$end)/2
     }
   }
return(raw.data)
}

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
  center <-as.numeric(chip[,"center"])
  if(intensity==TRUE){
   center <- cbind(center, chip[,5])
   colnames(center)[2] <- "intensity"
  }
  pos.lst <- c(pos.lst,list(center))
 
  all.pos <- c(all.pos, as.numeric((chip[,"center"])))
}

 pos.lst <- c(pos.lst, list(sort(all.pos)))
 names(pos.lst) <- c(names(chips), "all")
 return(pos.lst)
}

shuffle.patterns <- function(combined.mat, tf.num){
  samp.mat <- matrix(0, nrow=dim(combined.mat)[1], ncol=tf.num)
  colnames(samp.mat)=colnames(combined.mat)[1:tf.num]

  for(i in 1:tf.num){

    samp.mat[sample(dim(combined.mat)[1], length(which(combined.mat[,i] > 0))), i] <- 1
  }

  samp.mat <- getPattern(samp.mat, tf.num)
  samp.patterns <- sort(table(samp.mat[,"pattern"]), TRUE)

  return(samp.patterns)
}

pattern2tfs <- function(patterns, tf.names){
  tf.patterns <- list()
  for(i in 1:length(patterns)){
    bin.pat <-  integer.base.b(as.integer(names(patterns)[i]))

    bin.pat <- c(rep(0, length(tf.names) - length(bin.pat)), bin.pat)

    tf.patterns[[i]] <- tf.names[which(bin.pat == 1)]
  }
  names(tf.patterns) <- patterns
  return(tf.patterns)
}

pattern2crm <- function(pattern, tf.names){
  bin.pat <-  integer.base.b(as.integer(pattern))

  bin.pat <- c(rep(0, length(tf.names) - length(bin.pat)), bin.pat)

  return(tf.names[which(bin.pat == 1)])


}


integer.base.b <- function(x, b=2){
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base.b <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base.b[, ndigits-i+1] <- (x %% b)
                x <- (x %/% b)
        }
        if(N ==1) Base.b[1, ] else Base.b
}

bin2dec <- function(x) sum(2^(which(rev(x) == 1)-1))

getPattern <- function(occur.mat, tf.num){
  pattern <- c()
  occur.mat.tmp <- occur.mat[,1:tf.num]
  pattern <- apply(occur.mat.tmp, 1, bin2dec)
  return(cbind(occur.mat, pattern))
}

maskControl<- function(tf.pos, control.pos){

# Mask with IgG
  for(peak in control.pos$igG){
    for(j in 1:(length(tf.pos)-1)){
      diff <- abs(tf.pos[[j]] - peak)
      if(length(which(diff < 200)) > 0){
        tf.pos[[j]] <- tf.pos[[j]][-which(diff < 200)]
      }
    }
  }
  all.pos <- c()
  for(j in 1:(length(tf.pos)-1)){
    all.pos <- c(all.pos, tf.pos[[j]])
  }
  tf.pos$all <- all.pos
  return(tf.pos)
}

# analyzes whole chromosome (genic and intergenic regions)
analyzeChr <- function(tf.pos, tf.dist=400){
  tf.pos$all <- sort(tf.pos$all)
  regions <- c()
  start <- tf.pos$all[1]
  for(i in 1:length(tf.pos$all)){
    if( i == length(tf.pos$all) | (tf.pos$all[i+1] - tf.pos$all[i]) > tf.dist ){

      tf.region <- tf.pos[-length(tf.pos)]
      for(j in 1:(length(tf.region) )){
        tf.region[[j]] <- tf.pos[[j]][which(tf.pos[[j]] > start-1 & tf.pos[[j]] < tf.pos$all[i]+1)]
      }
      
      regions <- rbind(regions, c(countOccur(tf.region), start, tf.pos$all[i]))
      start <- tf.pos$all[i+1]
    }
  }
  return(regions)
}

analyzeGenes <- function(genes.chr, tf.pos, tf.dist=100, win.size=1000,  peak.score=TRUE){

 gene.tab <- c()

 for(gene in 1:dim(genes.chr)[1]){
   # define promoter region
   start <- 0
   end <- 0
   if( genes.chr[gene,"strand"] == "-"){
     start <- genes.chr[gene,"txEnd"] - win.size
     end <- start + (win.size * 2)
     }
   if(genes.chr[gene,"strand"] == "+"){
     start <- genes.chr[gene,"txStart"] - 1000
     end <- start + (win.size * 2)
   }
   tf.pos$all <- sort(tf.pos$all)
   # examine binding sites in region
         
   # extend region if there's nearby binding sites within tf.dist
   region.ind <- which(tf.pos$all >= start & tf.pos$all <= end)
   if(length(region.ind) > 0){
  
   while(region.ind[1] > 1){  
     if((tf.pos$all[region.ind[1]] - tf.pos$all[region.ind[1] - 1]) > tf.dist){
       break
     }
     region.ind <- c(region.ind[1] - 1, region.ind)
   }
   while(region.ind[length(region.ind)] < length(tf.pos$all)){
     if( (tf.pos$all[region.ind[length(region.ind)] + 1] - tf.pos$all[region.ind[length(region.ind)]]) > tf.dist)
       {
         break
       }
     region.ind <- c(region.ind, region.ind[length(region.ind)] + 1)
   }
   }
   
   start <- tf.pos$all[region.ind[1]]
   
   
   if(length(region.ind) > 0){
     for(i in 1:length(region.ind)){
       if(  i == length(region.ind) | (tf.pos$all[region.ind[i+1]] - tf.pos$all[region.ind[i]]) > tf.dist ){
         tf.region <- tf.pos[-length(tf.pos)]
   
         for(j in 1:(length(tf.region) )){
           tf.region[[j]] <- tf.pos[[j]][which(tf.pos[[j]] > (start-1) & tf.pos[[j]] < (tf.pos$all[region.ind[i]]+1))]
       
         }
         gene.tab <- rbind(gene.tab, c(countOccur(tf.region), start, tf.pos$all[region.ind[i]], as.character(genes.chr$name2[gene])))
         start <- tf.pos$all[region.ind[i+1]]
       }
     }
   }
 }
 if(length(gene.tab) != 0){
   colnames(gene.tab) <- c(names(tf.pos)[-length(tf.pos)], "start", "end", "gene")
 }
 return(gene.tab)
}



countOccur <- function(tfs.win){
  tf.count <-  rep(NA, length(tfs.win))
  names(tf.count) <- names(tfs.win)

  for(i in 1:length(tfs.win)){
    if(length(tfs.win[[i]]) > 0){
      tf.count[i] <- 1
    }
    else{
      tf.count[i] <- 0
    }
  }
    return(tf.count)
}

sampleBS <- function(samp.chip, test.chip, chr){
  pseudo.pos <- list()
  samp.pos <- plotBSdistr(c(samp.chip, test.chip), chr)
  test.pos <- plotBSdistr(test.chip, chr)

  for(i in test.pos){
    pseudo.pos <- c(pseudo.pos, list(sample(samp.pos$all, length(i))))
  }
  names(pseudo.pos) <- names(test.pos)
  return(pseudo.pos)
}

# find expression values for genes in MK
# input: file - gene expression data
#        genes - Hs genes with Ensembl ids
# output: matrix of Hs genes with probe annotations and expression values
load.express <- function(file="/scratch/dennis/projects/Haematlas/geneExpress.txt", genes){
  express <- read.csv(file, sep="\t")
  express.val <- rep(NA, dim(genes)[1])
  for(i in 1:dim(genes)[1]){
    gene.express <- express[which( express[,"Ensembl_gene_id"] == rownames(genes)[i] & express[,"MK"] == TRUE),]
    if(dim(gene.express)[1] > 0){
      express.val[i] <- gene.express[1,"MK.1"]
    }
  }
  genes <- cbind(genes, express.val)
  return(genes)
}

get.sig.occur <- function(obs.patterns, boot.patterns){
sig.patterns <- c()
for( i in 1:length(obs.patterns)){
 pattern.code <- names(obs.patterns)[i]
 occur.distr <- boot.patterns[which(names(boot.patterns) == pattern.code)]
 sig.patterns <- c(sig.patterns, length(which(occur.distr > obs.patterns[i]))/ length(occur.distr))
}
names(sig.patterns) <- names(obs.patterns)
#sig.patterns <- sort(sig.patterns)
return(sig.patterns)
}

# detects possible CRMs
get.CRMs <- function(type="promoter", tfs, controls=NULL, genes, tf.dist=100){
 genome.mat <- list()
  for(i in c(1:22, "X", "Y")){
    chr <- paste("chr", i, sep="")
    tf.pos <- plotBSdistr(tfs, chr)

    if(length(controls) > 0){
      control.pos <- plotBSdistr(controls, chr)
      
      tf.pos <- maskControl(tf.pos, control.pos)
    }
    if(type=="promoter"){
      genes.chr <- genes[genes[,"chrom"] == chr,]
      genes.chr <- subset(genes.chr, !duplicated(genes.chr[,"name2"]))
      genome.mat[[chr]] <- analyzeGenes(genes.chr, tf.pos, tf.dist=tf.dist)
    }
 
    if(type=="genome"){
      genome.mat[[chr]] <- analyzeChr(tf.pos, tf.dist=tf.dist)
    }
  }
 
 #combine occurrence matrix for all chromosomes
  combined.mat <- c()
  for(i in 1:length(genome.mat)){
    combined.mat <- rbind(combined.mat, cbind(genome.mat[[i]], names(genome.mat)[i]))
  }
 colnames(combined.mat)[dim(combined.mat)[2]] <- "chr"
 # get patterns of binding sites
 combined.mat <- getPattern(combined.mat, length(tfs))

# remove singleton pattern
 singletons <- c()
 for (i in 1:length(tfs)){
   pattern.tmp <- rep(0, length(tfs))
   pattern.tmp[i] <- 1
   singletons <- c(singletons, which(combined.mat[,"pattern"] == bin2dec(pattern.tmp)))
 }
 combined.mat <- combined.mat[-singletons,]
 obs.patterns <- sort(table(combined.mat[,"pattern"]), TRUE)

return(obs.patterns)
}


# detect CRMs and find significance
run.detection <- function(type="promoter", tfs, controls=NULL, genes, singles=FALSE, tf.dist=100, win.size=1000,  parallel.comp=FALSE, bootstrap=TRUE){

# get occurrences for the whole genome
 genome.mat <- list()
  for(chr.num in c(1:22, "X", "Y")){
    chr <- paste("chr", chr.num, sep="")
    tf.pos <- plotBSdistr(tfs, chr, intensity=FALSE)
    if(length(tf.pos$all) > 0){
      if(length(controls) > 0){
        control.pos <- plotBSdistr(controls, chr)
        
        tf.pos <- maskControl(tf.pos, control.pos)
      }
      if(type=="promoter"){
        genes.chr <- genes[genes[,"chrom"] == chr,]
        genes.chr <- subset(genes.chr, !duplicated(genes.chr[,"name2"]))
        if(dim(genes.chr)[1] > 0){
        genome.mat[[chr]] <- analyzeGenes(genes.chr, tf.pos, tf.dist=tf.dist, win.size=win.size)
      }
      }
      if(type=="genome"){
        genome.mat[[chr]] <- analyzeChr(tf.pos, tf.dist=tf.dist)
      }
    }
  }
 
 #combine occurrence matrix for all chromosomes
  combined.mat <- c()
  for(i in 1:length(genome.mat)){
    combined.mat <- rbind(combined.mat, cbind(genome.mat[[i]], names(genome.mat)[i]))
  }
 colnames(combined.mat)[dim(combined.mat)[2]] <- "chr"
 # get patterns of binding sites
  combined.mat <- getPattern(combined.mat, length(tfs))

                                        # remove singleton pattern
 if(singles == FALSE){
   singletons <- c()
   for (i in 1:length(tfs)){
     pattern.tmp <- rep(0, length(tfs))
     pattern.tmp[i] <- 1
     singletons <- c(singletons, which(combined.mat[,"pattern"] == bin2dec(pattern.tmp) | combined.mat[,"pattern"] == 0))
   }
   combined.mat <- combined.mat[-singletons,]
 }
 
obs.patterns <- sort(table(combined.mat[,"pattern"]), TRUE)
sig.patterns <- c()
boot.patterns <- c()

# get occurrences for sampled binding sites
 if(bootstrap==TRUE){
   boot.patterns <- c()
   if (parallel.comp == FALSE){
     for(boot.trial in 1:1000){
       boot.patterns <- c(boot.patterns, shuffle.patterns(combined.mat, length(tfs)))
     }
   }
   else{
                                        # parallel method
     registerDoMC()
     boot.patterns <- foreach(boot.trial = 1:1000, .combine=c) %dopar% {
       shuffle.patterns(combined.mat, length(tfs))
     }
   }
 }
 
# get significance of observed occurrences
  sig.patterns <- get.sig.occur(obs.patterns, c(boot.patterns, obs.patterns))

 colnames(combined.mat)[(1:2)+length(tfs)] <- c("start","end")
 return(list(sig.patterns=sig.patterns, obs.patterns=obs.patterns, boot.patterns=boot.patterns, combined.mat=combined.mat))
}

# number of each TF binding site in co-localization region
detect.binding.count <- function(type="promoter", tfs, controls=NULL, genes, tf.dist=100, win.size=1000){

# get occurrences for the whole genome
 genome.mat <- list()
  for(i in c(1:22, "X", "Y")){
    chr <- paste("chr", i, sep="")
    tf.pos <- plotBSdistr(tfs, chr, intensity=FALSE)
    if(length(tf.pos$all) > 0){
      if(length(controls) > 0){
        control.pos <- plotBSdistr(controls, chr)
        
        tf.pos <- maskControl(tf.pos, control.pos)
      }
      if(type=="promoter"){
        genes.chr <- genes[genes[,"chrom"] == chr,]
        genes.chr <- subset(genes.chr, !duplicated(genes.chr[,"name2"]))
        if(dim(genes.chr)[1] > 0){
        genome.mat[[chr]] <- analyzeGenes(genes.chr, tf.pos, tf.dist=tf.dist, win.size=win.size)
      }
      }
      if(type=="genome"){
        genome.mat[[chr]] <- analyzeChr(tf.pos, tf.dist=tf.dist)
      }
    }
  }
 
 #combine occurrence matrix for all chromosomes
  combined.mat <- c()
  for(i in 1:length(genome.mat)){
    combined.mat <- rbind(combined.mat, cbind(genome.mat[[i]], names(genome.mat)[i]))
  }
 colnames(combined.mat)[dim(combined.mat)[2]] <- "chr"


 colnames(combined.mat)[(1:2)+length(tfs)] <- c("start","end")
 return(combined.mat)
}

                                  

examinePattern <- function(pattern, tfs, genes,  hists, controls, cell2.tfs=NULL, cell2.combined.mat=NULL, cell2.sig.patterns=NULL){
         
crm.tfs <- pattern2tfs(pattern, names(tfs))
crm.mat <- combined.mat[which(combined.mat[,"pattern"] == names(pattern)),]

# filter for only CRMs with high peaks
tfs.tmp <- tfs
for(j in 1:length(tfs)){
  peak.thres <- summary(tfs[[j]][,5])["Median"]
  tfs.tmp[[j]] <- tfs.tmp[[j]][which(tfs.tmp[[j]][,5] > peak.thres),]
}

detection.results <- run.detection(type="promoter", tfs.tmp, controls, genes, tf.dist=500)
sig.patterns <- detection.results[["sig.patterns"]]
obs.patterns <- detection.results[["obs.patterns"]]
combined.mat <- detection.results[["combined.mat"]]

crm.hpeaks <- combined.mat[which(combined.mat[,"pattern"] == names(pattern)),]

high.peaks <- rep(0, dim(crm.mat)[1])
high.peaks[match( crm.hpeaks[,"gene"], crm.mat[,"gene"], nomatch=0)] <- 1
crm.mat <- cbind(crm.mat, high.peaks)

# CRMs with histone markers
detection.results <- run.detection(type="promoter", hists, controls, genes, tf.dist=2000, singles=TRUE)
sig.patterns <- detection.results[["sig.patterns"]]
obs.patterns <- detection.results[["obs.patterns"]]
combined.mat <- detection.results[["combined.mat"]]

hist.mat <- matrix(0, nrow=dim(crm.mat)[1], ncol=length(hists), dimnames=list(c(),names(hists)))

hist.matches <- match(crm.mat[,"gene"], combined.mat[,"gene"], nomatch=0)
for(i in (which(hist.matches > 0))){
  hist.mat[i,] <- combined.mat[hist.matches[i], 1:length(hists)]
}
crm.mat <- cbind(crm.mat,hist.mat)

if(length(cell2.tfs) > 0){
crm.match <- match(crm.mat[,"gene"], cell2.combined.mat[,"gene"])
K562.crm <-  rep(NA, dim(crm.mat)[1])
K562.pvalue <- rep(NA, dim(crm.mat)[1])
for(i in 1:length(crm.match)){
  if(!is.na(crm.match[i])){
    pvalue.match <- cell2.sig.patterns[which(names(cell2.sig.patterns) == cell2.combined.mat[crm.match[i], "pattern"])]
  
    K562.crm[i] <- paste(pattern2tfs(pvalue.match, names(cell2.tfs))[[1]], collapse="-")
    K562.pvalue[i] <- pvalue.match
  }
}
crm.mat <- cbind(crm.mat, K562.crm, K562.pvalue)
}

write.table(crm.mat, file=paste("./output/CRMs_of_interest/", paste(crm.tfs[[1]], collapse="-"), ".txt", sep=""), sep="\t")
return(crm.mat)
}


# gets width (relative to each other) of CRM near a gene
get.crm.size <- function(tfs.names, gene, crm.data){
  gene.match <- crm.data[which(crm.data[,"gene"] == gene),, drop =FALSE]
  crm.width <- NA
  if(length(which(tfs.names == "TAL1")) > 0){
    tfs.names[which(tfs.names == "TAL1")] <- "scl"
  }
     
  for(i in 1:dim(gene.match)[1]){
  
    if( sum(as.numeric(gene.match[i,match(tfs.names, colnames(gene.match)) ])) == length(tfs.names)){
      crm.width <- as.numeric( gene.match[i, "end"]) - as.numeric(gene.match[i,"start"]) 
    }
  }
  return(crm.width)
}

# returns
analyze.CRM <- function(crm.mat, tfs, genes){
crm.info <- list()
for(i in 1:dim(crm.mat)[1]){
  tf.names <- colnames(crm.mat)[which(crm.mat[i,1:length(tfs)] == 1)]
  tf.dists <- c()
  tf.peaks <- c()
  tf.names.tmp <- c()
  for(j in 1:length(tf.names)){
    matched.site <-  tfs[[tf.names[j]]][ tfs[[tf.names[j]]][,"chr"] == crm.mat[i,"chr"] &  tfs[[tf.names[j]]][,"center"] >= as.numeric(crm.mat[i,"start"]) &  tfs[[tf.names[j]]][,"center"] <= as.numeric(crm.mat[i,"end"]),]
    tf.dists <- c(tf.dists, matched.site[,"center"])
    tf.peaks <- c(tf.peaks, matched.site[,7])
    tf.names.tmp <- c(tf.names.tmp, rep(tf.names[j], dim(matched.site)[1]))
  }
  info.tmp <- rbind(tf.dists, tf.peaks)
  colnames(info.tmp) <- tf.names.tmp
  rownames(info.tmp) <- c("center.pos", "peak.height")
  crm.info[[i]] <- info.tmp
}
return(crm.info)
}


# create list of genes and their nearby crms
# input:
# crm matrix
# matrix of gene name and TSS position
# window threshold

dist.gene.crm <- function(combined.mat, gene.dat, win.thres=500000){

 combined.mat <- cbind(combined.mat, 1:dim(combined.mat)[1])
 
 registerDoMC(20)
 gene.lst <- foreach(i = 1:dim(gene.dat)[1])%dopar%{
   nearby.crms <- c()
   crms <- combined.mat[which(combined.mat[,"chr"] == gene.dat[i, "Chromosome.Name"]),]
   if(length(crms) > 0){
     crms.centre <- ((as.numeric(crms[,"end"]) -as.numeric( crms[,"start"]))/2) +as.numeric( crms[,"start"])
     crms <- crms[which(crms.centre >  (as.numeric(gene.dat[i, "Transcript.Start..bp."]) - win.thres) & crms.centre <  (as.numeric(gene.dat[i, "Transcript.Start..bp."]) + win.thres)),, drop=FALSE]
     if(length(crms) > 0){
       nearby.crms <- rbind(crms[,dim(crms)[2]], crms[,"pattern"])
       rownames(nearby.crms) <- c("crm.index", "pattern")
     }
   }
  return(nearby.crms) 
 }
 return(gene.lst)
}

# map crms to promoters of genes
map.gene.crm <- function(combined.mat, gene.dat, win.thres=1000){
 registerDoMC()
 gene.crms <- foreach(i = 1:dim(gene.dat)[1], .combine=rbind)%dopar%{
   crm.ind <- which(combined.mat[,"chr"] == gene.dat[i, "Chromosome.Name"])
    crms <- combined.mat[crm.ind,]
     if(length(crms) > 0){
       crms.centre <- ((as.numeric(crms[,"end"]) -as.numeric( crms[,"start"]))/2) +as.numeric( crms[,"start"])
       crm.ind2 <- which(crms.centre >  (as.numeric(gene.dat[i, "Transcript.Start..bp."]) - win.thres) & crms.centre <  (as.numeric(gene.dat[i, "Transcript.Start..bp."]) + win.thres))
       crms <- crms[crm.ind2,, drop=FALSE]
       crm.ind <- crm.ind[crm.ind2]
     }
     if(length(crms) > 0){
       return(cbind(crms, as.character(gene.dat[i, "Associated.Gene.Name"]), crm.ind))
     }else{return(c())}
 }
 colnames(gene.crms)[c(dim(gene.crms)[2]-1,dim(gene.crms)[2])] <- c("gene", "crm.ind")
 gene.crms <- gene.crms[!duplicated(gene.crms[,c("gene","start")]),]
 return(gene.crms)
}


# find regions where CRMs overlap between two datasets
overlap.crms <- function(crm1, crm2, indexes=FALSE, dist.thres=0){
  crm.centre2 <- (as.numeric(crm2[, "end"]) - as.numeric(crm2[, "start"]))/2 + as.numeric(crm2[, "start"])
  
  registerDoMC()
  matched.crms <- foreach(i = 1:dim(crm1)[1], .combine=rbind)%dopar%{
    crm.centre1 <- (as.numeric(crm1[i, "end"]) - as.numeric(crm1[i, "start"]))/2 + as.numeric(crm1[i, "start"])
    matched.crm <- which(crm2[,"chr"] == crm1[i, "chr"] & (as.numeric(crm2[,"start"]) - dist.thres) < crm.centre1  & (as.numeric(crm2[,"end"]) + dist.thres) > crm.centre1 )
    if(length(matched.crm) > 0){
      return(c(i, matched.crm[1]))
    }else{
      matched.crm <- which(crm2[,"chr"] == crm1[i, "chr"] & (crm.centre2 + dist.thres) > as.numeric(crm1[i,"start"]) & (crm.centre2 - dist.thres) < as.numeric(crm1[i,"end"]))
      if(length(matched.crm) > 0){
        return(c(i, matched.crm[1]))
      }
      return(c())
    }
    
  }
  results <- NULL
  if(indexes == TRUE){
    results <- matched.crms # return only indexes of crm matrix
  }else{
    results <- cbind(crm1[matched.crms[,1],], crm2[matched.crms[,2],]) # return combined crm matrices
  }
 return(results)
}

# functions locates all CRMs with at least one eQTL SNP
overlap.eqtl <- function(combined.mat, matched.snps, dist.thres=100){
registerDoMC()
snp.crms <- foreach(i = 1:dim(matched.snps)[1], .combine=rbind)%dopar%{

  snp.crm <- combined.mat[which( combined.mat[,"chr"] == matched.snps[i,"Chromosome"] &   as.numeric(combined.mat[,"start"]) - dist.thres < as.numeric(matched.snps[i,"Coordinate_HG18"]) & as.numeric(combined.mat[,"end"]) + dist.thres  >  as.numeric(matched.snps[i,"Coordinate_HG18"]) ),, drop=FALSE]
  if(length(snp.crm) > 0){
    return( c(snp.crm[1,], matched.snps[i,]) )
  }
  return(c())
}
keep.crms <- which(!duplicated(snp.crms[,c( "start", "end", "pattern")]))
snp.crms <- snp.crms[keep.crms,]
return(snp.crms)
}

# calculate dissimilarity between a two lists of overlapping CRMs
crm.dissimilarity <- function(crm1, crm2, method="hamming"){ 
  crm.dists <- c()
  for(i in 1:dim(crm1)[1]){
    gm.bin.pattern <- as.numeric(crm1[i,])
    k562.bin.pattern <-  as.numeric(crm2[i,])
    if(method == "jaccard"){
      crm.dists <- c(crm.dists,  vegdist(rbind(k562.bin.pattern, gm.bin.pattern), method="jaccard", binary=TRUE))
    }
    if(method == "hamming"){
      crm.dists <- c(crm.dists,  hamming.distance(k562.bin.pattern, gm.bin.pattern))
    } 
  }
  return(crm.dists)
}
       
# removes TFs from CRM matrix
reduce.crm <- function(combined.mat, orig.tfs, new.tfs){
 matched.crms <- combined.mat[,match( new.tfs, colnames(crm.mat))]
 class(matched.crms) <- "numeric"

 new.mat <- combined.mat[which(apply(matched.crms, 1, sum) == 1), -which(is.na(match(orig.tfs, new.tfs)))]
 
 new.mat[,"pattern"] <- apply(new.mat[,1:length(new.tfs)], 1,   bin2dec)
 return(new.mat)
}


# generates a matrix of TF co-occurrence scores
cooccur.score <- function(combined.mat){
  expect.mat <- matrix(0, nrow=ncol(combined.mat), ncol=ncol(combined.mat), dimnames=list(colnames(combined.mat), colnames(combined.mat)))
  tf.num <- ncol(combined.mat)
  trial.num <- 100
  for(trial in 1:trial.num){  
    samp.mat <- matrix(0, nrow=nrow(combined.mat), ncol=tf.num)
    for(i in 1:tf.num){
      samp.mat[sample(nrow(combined.mat), length(which(combined.mat[,i] > 0))), i] <- 1
    }
    for(i in 1:tf.num){
      expect.mat[,i] <- expect.mat[,i] + apply(samp.mat[which(samp.mat[,i] > 0), ], 2, function(x){length(which(x > 0))})
    }
  }
  expect.mat <- expect.mat/trial.num
  
  observe.mat <- matrix(0, nrow=ncol(combined.mat), ncol=ncol(combined.mat), dimnames=list(colnames(combined.mat), colnames(combined.mat)))
  for(i in 1:tf.num){
    observe.mat[,i] <-  apply(combined.mat[which(combined.mat[,i] > 0), ], 2, function(x){length(which(x > 0))})
  }
  score.mat <- log(observe.mat/expect.mat)
  score.mat[row(score.mat) == col(score.mat)] <- 1
  return(score.mat)
}

# calculate co-occurrence score for two CRMs
crm.cooccur <- function(crm2, score.mat){
  names(crm2) <- colnames(score.mat)
  if(sum(as.numeric(crm2)) > 0){
    mat.brief <- score.mat[names(crm2)[crm2 > 0], names(crm2)[crm2 > 0], drop=FALSE]
    return(sum(apply(mat.brief, 2, max)))
  }else{return(-2)}
}
