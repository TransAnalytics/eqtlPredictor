 getPattern <- function(occur.mat, tf.num){
  pattern <- c()
  occur.mat.tmp <- occur.mat[,1:tf.num]
  pattern <- apply(occur.mat.tmp, 1, bin2dec)
  return(cbind(occur.mat, pattern))
}

# functions locates all CRMs with at least one eQTL SNP
overlap.eqtl <- function(combined.mat, matched.snps, dist.thres=100){
  
  snp.crms = c()  
  chrs = unique(matched.snps[,"Chromosome"])
  
  for(i in chrs){
    
    snp.chrs = matched.snps[which(matched.snps[,"Chromosome"] == i),]
    crm.chrs = combined.mat[which(combined.mat[,"chr"] == i),]
    
    matched = apply(crm.chrs, 1, function(x){which(findInterval(as.numeric(snp.chrs[,"Coordinate_HG18"]), c(as.numeric(x["start"]) - dist.thres, as.numeric(x["end"]) + dist.thres) ,rightmost.closed= TRUE)==1)[1] })
    snp.crm = cbind(crm.chrs[which(!is.na(matched)), ], snp.chrs[na.omit(matched),]) 
    
    if(nrow(snp.crm) > 0){snp.crms = rbind(snp.crms, snp.crm)}
  }
  keep.crms <- which(!duplicated(snp.crms[,c( "start", "end", "pattern")]))
  snp.crms <- snp.crms[keep.crms,]
  return(snp.crms)
}

# deprecated
regress.expression <- function (crm.mat, training.express, all.samples, gene.pos, max.distance = 1e+06) 
{
  test.express = training.express
  crm.gene.preds <- list()
  for(crm.ind in 1:nrow(crm.mat)) 
  {
    gene.preds <- c()
    gene.chr <- gene.pos[gene.pos[, "Chromosome.Name"] == crm.mat[crm.ind, "chr"], ]
    matched.genes <- as.character(gene.chr[which(abs(gene.chr[, "Transcript.Start..bp."] - mean(c(as.numeric(crm.mat[crm.ind, "start"]), as.numeric(crm.mat[crm.ind, "end"])))) <= max.distance), "Associated.Gene.Name"])
  
    if (length(matched.genes) > 0) {
    matched.genes <- matched.genes[!duplicated(matched.genes)]
    combined.mat <- crm.mat[crm.ind, , drop = FALSE]
    combined.mat <- cbind(combined.mat[rep(nrow(combined.mat), each = length(matched.genes)), , drop = FALSE], matched.genes)
    
    colnames(combined.mat)[dim(combined.mat)[2]] <- "gene"
    tfs <- gm.tfs
    patterns <- unique(as.numeric(combined.mat[, 
                                               "pattern"]))
    names(patterns) <- patterns
    pred.results.lst <- list()
    for (k in 0:4) {
      cv.samps <- ceiling(length(all.samples)/5)
      test.express.cols <- c(all.samples[max(cv.samps * 
                                               k, 1):min(cv.samps * (k + 1), length(all.samples))])
      training.cols <- all.samples[-match(test.express.cols, 
                                          all.samples)]
      pred.results <- predict.gene.express.2(patterns, 
                                             combined.mat, tfs, training.express, test.express, 
                                             training.cols, test.express.cols, method = "RF", 
                                             step.wise = FALSE, bootstrap = FALSE, null.tfs = all.tfs)
      if (length(pred.results) < 1) {
        break
      }
      obs.express <- test.express[as.numeric(pred.results[, "target_expression_ind"]), test.express.cols, drop = FALSE]
      r.squareds <- c()
      for (i in 1:dim(pred.results)[1]) {
        r.squareds <- c(r.squareds, summary(lm(as.numeric(pred.results[i, 4:(dim(pred.results)[2] - 1)]) ~ as.numeric(obs.express[i, 
                                                                                                                               ])))$r.squared)
      }
      pred.results.lst[[k + 1]] <- cbind(pred.results, r.squareds)
    }
    if (length(pred.results.lst) < 1) {
      crm.gene.preds[[crm.ind]] <- NA
    }else{
    gene.preds <- pred.results.lst[[1]][, c("CRM_ind", 
                                            "target_expression_ind", "target_gene"), drop = FALSE]
    mean.r2.fit <- c()
    mean.r2.pred <- c()
    for (i in 1:dim(gene.preds)[1]) {
      mean.r2.fit <- c(mean.r2.fit, mean(as.numeric(unlist(lapply(pred.results.lst, 
                                                                  function(x) {
                                                                    x[i, "R.squared_fit"]
                                                                  })))))
      mean.r2.pred <- c(mean.r2.pred, mean(as.numeric(unlist(lapply(pred.results.lst, 
                                                                    function(x) {
                                                                      x[i, "r.squareds"]
                                                                    })))))
      }
      gene.preds <- cbind(gene.preds, mean.r2.fit, mean.r2.pred)
      crm.gene.preds[[crm.ind]] <- gene.preds
    }
    }else{ crm.gene.preds[[crm.ind]] <- NA}
}
  return(crm.gene.preds)
}

# function to find CRMs located within 1kb of every gene's TSS
genes.to.crms <- function(gm.combined.mat, gene.pos, dist.thres=500){
  gene.crms = c()  
  chrs = unique(gm.combined.mat[,"chr"])
  
  for(i in chrs){
    
    
    crm.chrs = gm.combined.mat[which(gm.combined.mat[,"chr"] == i),, drop=FALSE]
    gene.chrs = gene.pos[which(gene.pos[,"Chromosome.Name"] == i),, drop=FALSE]
    
    matched = apply(crm.chrs, 1, function(x){which(findInterval(as.numeric(gene.chrs[,"Transcript.Start..bp."]), c(as.numeric(x["start"]) - dist.thres, as.numeric(x["end"]) + dist.thres) ,rightmost.closed= TRUE)==1)[1] })
    gene.crm = cbind(crm.chrs[which(!is.na(matched)), ], gene.chrs[na.omit(matched),]) 
    
    if(nrow(gene.crm) > 0){gene.crms = rbind(gene.crms, gene.crm)}
  }
  keep.crms = which(!duplicated(gene.crms[,c( "start", "end", "chr")]))
  gene.crms = gene.crms[keep.crms,]
  return(gene.crms)
}

  
##
# Get GO similarity for TF-gene pairs
##
GO.similarity <- function(crm.snps, gene.pos, tf.names, dist.thres=1000000){
  tf.entrez = as.numeric(gene.pos[match(tf.names, gene.pos[,"Associated.Gene.Name"]), "EntrezGene.ID"])
  names(tf.entrez) = tf.names
  
  gene.sim <- list()
  for(i in 1:nrow(crm.snps)){
   gene.chr = gene.pos[which( gene.pos[,"Chromosome.Name"] == crm.snps[i,"Chromosome"]),]
   matched.ids = gene.chr[which( abs(as.numeric(gene.chr[,"Transcript.Start..bp."]) - as.numeric(crm.snps[i, "Coordinate_HG18"])) < dist.thres),]
   matched.genes = as.numeric(matched.ids[,"EntrezGene.ID"])
   names(matched.genes) = matched.ids[,"Associated.Gene.Name"]
   matched.genes <- matched.genes[-which(is.na(matched.genes))]
   matched.genes <- matched.genes[!duplicated(matched.genes)]
   
   matched.tfs <- tf.entrez[which(crm.snps[i,tf.names] == 1)]
 
   
  crm.gene.sim <- tryCatch(getGeneSim(as.character(matched.tfs), as.character(matched.genes), similarity="OA", method="Lin"), error = function(e){NA})
  if(length(crm.gene.sim) > 1){
     rownames(crm.gene.sim) <- names(matched.tfs)[match(rownames(crm.gene.sim), as.character(matched.tfs))]  
     colnames(crm.gene.sim) <- names(matched.genes[match(colnames(crm.gene.sim), matched.genes)])
   }
   gene.sim[[i]] <- crm.gene.sim
  }
  return(gene.sim)
}
 
map.ctcf <- function(crm.snps, gene.pos, gm.ctcf){
  gene.ctcfs <- list()
  for(i in 1:nrow(crm.snps)){
    gene.chr = gene.pos[which( gene.pos[,"Chromosome.Name"] == crm.snps[i,"Chromosome"]),]
    ctcf.chr = gm.ctcf[gm.ctcf[,"chr"] == crm.snps[i,"chr"],]
    
    intervals <- cbind(crm.snps[i,"Coordinate_HG18"], gene.chr[,"Transcript.Start..bp."])
    
    gene.ctcf = apply(intervals, 1, function(x){mean(as.numeric(ctcf.chr[which(findInterval(as.numeric(ctcf.chr[,"start"]), sort(as.numeric(x)) ,rightmost.closed= TRUE)==1), "signalValue"]))})
    gene.ctcf[is.na(gene.ctcf)] <- 0   
    gene.ctcfs[[i]] <- gene.ctcf
  }
  return(gene.ctcfs)
}

# generates a matrix of TF co-occurrence scores
 cooccur.score <- function(combined.mat){
   expect.mat <- matrix(0, nrow=ncol(combined.mat), ncol=ncol(combined.mat), dimnames=list(colnames(combined.mat), colnames(combined.mat)))
   tf.num <- ncol(combined.mat)
   trial.num <- 1000
   for(trial in 1:trial.num){  
     
     samp.mat <- matrix(0, nrow=nrow(combined.mat), ncol=tf.num)
     for(i in 1:tf.num){
       samp.mat[sample(nrow(combined.mat), length(which(combined.mat[,i] > 0))), i] <- 1
     }
     for(i in 1:tf.num){
       expect.mat[,i] <- expect.mat[,i] + apply(samp.mat[which(samp.mat[,i] > 0), ], 2, function(x){length(which(x > 0))})
     }
     
     #progress
     if(trial%%100 == 0){print(paste((trial/trial.num)*100, "% complete", sep=""))}
     
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
  
 
 feature.select <- function(gene.lst, gene.preds, crm.mat, samps, faire.gene, gene.crms, gene.sim, insulators, score.mat, dist.thres=0){
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
     #progress
     if(i%%10 == 0){print(paste(round((i/length(samps))*100), "% complete", sep=""))}
    
   
     crm.centre <- mean( as.numeric(crm.mat[samps[i], "end"]), as.numeric(crm.mat[samps[i], "start"]))

     crm.gene.dist <- abs(gene.pos[match(gene.preds[[samps[i]]][,"target_gene"], gene.pos[,"Associated.Gene.Name"]),"Transcript.Start..bp."] - crm.centre )
     gene.filter <- which(crm.gene.dist >= dist.thres)
         
       
        train.params <- c(train.params, gene.preds[[samps[i]]][gene.filter,"mean.r2.fit"])
         
         
         train.params2 <- c(train.params2, crm.gene.dist[gene.filter])
         
         gene.match <- gene.crms[match( gene.preds[[samps[i]]][,"target_gene"], gene.crms[,"gene"]), 1:length(gm.tfs),drop=FALSE]
         gene.match[is.na(gene.match)] <- 0
         
         crm.dis <- apply(gene.match, 1, function(x){ crm.cooccur( as.numeric(x), score.mat)})
         train.params3 <- c(train.params3, crm.dis[gene.filter])
         
         gene.match <- match(gene.preds[[samps[i]]][,"target_gene"], faire.gene[,"Associated.Gene.Name"])
         gene.match[is.na(gene.match)] <- 0
         gene.match[which(gene.match != 0)] <- faire.gene[gene.match, "signalValue"]
         train.params4 <- c(train.params4, gene.match[gene.filter])
         
         if(length(gene.sim[[samps[i]]]) == 1 & length(which(is.na(gene.sim[[samps[i]]]))) > 0 )
         { gene.match <- rep(0, nrow(gene.preds[[samps[i]]]))}else{
           sim.score <- apply(gene.sim[[samps[i]]], 2, mean)
           sim.score[is.na(sim.score)] <- 0
           gene.match <- rep(0, dim(gene.preds[[samps[i]]])[1])
           match.ind <- match(  names(sim.score) , gene.preds[[samps[i]]][,"target_gene"])
           gene.match[na.omit(match.ind)] <- sim.score[which(!is.na(match.ind))]
         }
         train.params5 <- c(train.params5, gene.match[gene.filter])
         
         
         gene.match <- insulators[[samps[i]]][match( gene.preds[[samps[i]]][,"target_gene"], names(insulators[[samps[i]]]) )]
         gene.match[is.na(gene.match)] <- 0
         train.params6 <- c(train.params6, gene.match[gene.filter])
         
         response.tmp <- rep(0, dim(gene.preds[[samps[i]]])[1])
         response.tmp[which(gene.preds[[samps[i]]][,"target_gene"] == crm.mat[samps[i],"Associated.Gene.Name"])] <-  1
         train.response <- c(train.response, response.tmp[gene.filter])
         gene.names <- c(gene.names, gene.preds[[samps[i]]][,"target_gene"][gene.filter])
         
         if(length(gene.filter) > 0){ samps.tested <- c(samps.tested, paste(samps[i],":", gene.filter, sep=""))}
       }
     
   
   features.mat <- cbind(train.response, train.params, train.params2, train.params3, train.params4, train.params5, train.params6)
   rownames(features.mat) <- samps.tested
   return(features.mat)
 }

 
 
 classify.targets <- function(gene.preds, crm.mat, train.samps, test.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat, model.var=c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator"), dist.thres=0){
   trainset <- feature.select(gene.preds, crm.mat, train.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat)
   class(trainset) <- "numeric"
   colnames(trainset) <- c("target", "mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
   rownames(trainset) <- 1:dim(trainset)[1]
   trainset <- trainset[,match(model.var, colnames(trainset))]
   model.obj <- NULL
 

  model.obj <- try(randomForest( trainset[,-1, drop=FALSE] , as.factor(trainset[,1])))
   
   
   
   testset <- feature.select(gene.preds, crm.mat, test.samps, faire.gene, gene.crms, gene.sim, insulators, score.mat=score.mat, dist.thres)
   test.response <- as.numeric(testset[,1])
   test.params <- testset[,-1]
   
   colnames(test.params) <- c("mean.r2.fit", "distance", "crm.pattern", "faire.signal", "GO.similarity", "insulator")
   class(test.params) <- "numeric"
   test.params <- test.params[,match(model.var[-1], colnames(test.params)), drop=FALSE]
   rownames(test.params) <- 1:dim(test.params)[1]
   target.pred <- NULL

     target.pred <- predict(model.obj,test.params, type="response")
     target.pred <- as.numeric(levels(target.pred)[target.pred])
   
 
   return(cbind(pred=target.pred, obs=test.response, crm.genes=rownames(testset)))
 }

  