# removes duplicate gene entries and keeps gene with lowest associated value
lowest.gene <- function(gene.pvals){
  filtered.genes <- c()
  for(i in 1:length(gene.pvals)){
   gene.tmp <- gene.pvals[which(names(gene.pvals) == names(gene.pvals)[i])]
   filtered.genes <- c(filtered.genes, gene.tmp[which.min(gene.tmp)])
  }
  filtered.genes <- filtered.genes[!duplicated(names(filtered.genes))]
  return(filtered.genes)
}


train.gam <- function(response.express, predictor.express, stepwise=TRUE, smoothing="s", k=-1, inter.effects=FALSE){
  

    gam.dat <- c()
    gam.obj <- NULL

    for(h in 1:dim(predictor.express)[1]){
      gam.dat <- cbind(gam.dat, as.numeric(predictor.express[h,]))
    }
  
    gam.dat <- cbind(as.numeric(response.express), gam.dat)
    colnames(gam.dat) <- c("target", rownames(predictor.express))
    if(stepwise==TRUE){
          
    gam.obj <- gam(as.formula(gam.formula(rownames(predictor.express), smoothing=FALSE, inter.effects=FALSE, k=k)), data=as.data.frame(gam.dat))
    gam.obj <- step.gam(gam.obj, scope=step.list( terms(as.formula(gam.formula(rownames(predictor.express), smoothing=FALSE, inter.effects=FALSE, k=k))), crm.tfs=rownames(predictor.express), target.intercept=-1, inter.effects=inter.effects), trace=FALSE, data=as.data.frame(gam.dat))
    
  }else{ 
    gam.obj <- gam(as.formula(gam.formula(rownames(predictor.express), smoothing=smoothing, inter.effects=inter.effects, k=k))    , data = as.data.frame(gam.dat))
  }
    
   
  return(gam.obj)
}


predict.express <- function(gam.obj,  predictor.express){
  
     haem.params <- c()
     for(k in 1:dim(predictor.express)[1]){
       haem.params <- cbind(haem.params, as.numeric(predictor.express[k,]))
     }
     colnames(haem.params) <- rownames(predictor.express )   
   
   return( predict.gam(gam.obj, as.data.frame(haem.params)))
     
 }

step.list <- function(formula.terms, crm.tfs, target.intercept, inter.effects=FALSE){
 term.lst <- attr(formula.terms, "term.labels")
 form.lst <- list()
  inter.iter <-  1 +length(term.lst)
 for(i in 1:length(term.lst)){
   if(length(term.lst) > 1){    
     if(inter.effects == TRUE & i < length(term.lst)){
       
       form.lst <- c(form.lst, as.formula(paste("~", target.intercept, " + ", term.lst[i], " + ", attr(terms(as.formula(gam.formula(crm.tfs, smoothing="s", inter.effects=TRUE))), "term.labels")[i] , "+", paste(attr(terms(as.formula(gam.formula(crm.tfs, smoothing="s", inter.effects=TRUE))), "term.labels")[inter.iter : (inter.iter + length(term.lst) -1 - i)], collapse="+")   )))
       inter.iter <-  (inter.iter + length(term.lst) - i) 
     }
     else{
       form.lst <- c(form.lst, as.formula(paste("~", target.intercept, " + ", term.lst[i], " + ", attr(terms(as.formula(gam.formula(crm.tfs, smoothing="s", inter.effects=TRUE))), "term.labels")[i] )))
     }
   }else{
     form.lst <- c(form.lst, as.formula(paste("~", target.intercept, " + ", term.lst[i], " + ", attr(terms(as.formula(gam.formula(crm.tfs, smoothing="s", inter.effects=FALSE))), "term.labels")[i] )))
   }
 }
 names(form.lst) <- crm.tfs
 return(form.lst)
}


                         

gam.formula <- function(tf.names, smoothing=FALSE, inter.effects=FALSE, k=-1){
  model.formula <- ""
  if(smoothing == FALSE){
    model.formula <- paste("target", "~", paste( tf.names,  collapse=" + "))
  }else{
    model.formula <- paste("target", "~", paste(paste(smoothing, "(", tf.names, ",k=", k, ")", sep=""),  collapse=" + "))
  }
  if(inter.effects == TRUE){
    inter.effects <- ""
    for(i in 1:(length(tf.names)-1)){
      for(j in (i+1):length(tf.names)){
        if(smoothing == FALSE){
          inter.effects <- paste(inter.effects, " + ", tf.names[i], ":", tf.names[j], sep="")
        }
        if(smoothing == "s")
        {
          inter.effects <- paste(inter.effects, " + ", smoothing, "(", tf.names[i], "*", tf.names[j], ",k=", k,")", sep="")
        }
        if(smoothing == "lo")
        {
          inter.effects <- paste(inter.effects, " + ", smoothing, "(", tf.names[i], ",", tf.names[j], ")", sep="")
        }
      }
    }
    model.formula <- paste(model.formula, inter.effects, sep="")
  }
  return(model.formula)
}


tissue.predictions <- function(crm.mat, training.express, tissue.ind, all.cols=2:159, step.wise=FALSE, inter.effects=FALSE){

  preds <- c()
  observs <- c()
  training.cols <- all.cols[-which(all.cols == tissue.ind)]
  test.express.cols <- tissue.ind
  predictions <- c()
  observations <- c()
  pred.genes <- c()

    for(i in 1:nrow(crm.mat)){
    tf.names <- names(tfs)[which(crm.mat[i, 1:length(tfs)] == 1)]
    
    target.gene <- crm.mat[i, "gene"]    
 
    tf.express <- as.matrix(training.express[match( toupper(tf.names), training.express[,"Gene"]),, drop=FALSE])
    
    rownames(tf.express) <- tf.names
    target.express <- training.express[which( training.express[,"Gene"] == target.gene),,drop=FALSE]
  
    if(dim(target.express)[1] > 1){
      target.express <- target.express[which.max(apply(target.express[, all.cols], 1, var)),, drop=FALSE] 
    }
 
    if(dim(target.express)[1] > 0){
      
      gam.obj <- NULL
      gam.obj <- train.gam(target.express[, training.cols, drop=FALSE], tf.express[,training.cols, drop=FALSE], stepwise=step.wise, k=-1)
      
      if(length( attr(terms(as.formula(gam.obj)), "term.labels")) > 0){
        if(length(gam.obj) > 0 ){
          haem.params <- c()
          for(k in 1:dim(tf.express)[1]){
            haem.params <- cbind(haem.params, as.numeric(tf.express[k,test.express.cols]))
          }
          colnames(haem.params) <- tf.names
        
          tar.gene.pred <- predict.gam(gam.obj, as.data.frame(haem.params))
          
          observations <- c(observations, as.numeric(target.express[, test.express.cols]))
          predictions <- c(predictions, tar.gene.pred)
          pred.genes <- c(pred.genes, crm.mat[i,"gene"])
        }
      }
    }
  }
  names(observations) <- pred.genes
  names(predictions) <- pred.genes
  return(rbind(observations, predictions))
}


step.gam <- function (object, scope, scale.aic,   direction = c("both", "backward", 
    "forward"), trace = TRUE, keep = NULL, steps = 1000, ...) 
{
    scope.char <- function(formula) {
        formula = update(formula, ~-1 + .)
        tt <- terms(formula)
        tl <- attr(tt, "term.labels")
        if (attr(tt, "intercept")) 
            c("1", tl)
        else tl
    }
    re.arrange <- function(keep) {
        namr <- names(k1 <- keep[[1]])
        namc <- names(keep)
        nc <- length(keep)
        nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
            namc))
    }
    untangle.scope <- function(terms, regimens) {
        a <- attributes(terms)
        response <- deparse(a$variables[[2]])
        term.labels <- a$term.labels
        if (!is.null(a$offset)) {
            off1 <- deparse(a$variables[[a$offset]])
        }
        nt <- length(regimens)
        select <- integer(nt)
        for (i in seq(nt)) {
            j <- match(regimens[[i]], term.labels, 0)
            if (any(j)) {
                if (sum(j > 0) > 1) 
                  stop(paste("The elements of a regimen", i, 
                    "appear more than once in the initial model", 
                    sep = " "))
                select[i] <- seq(j)[j > 0]
                term.labels <- term.labels[-sum(j)]
            }
            else {
                if (!(j <- match("1", regimens[[i]], 0))) 
                  stop(paste("regimen", i, "does not appear in the initial model", 
                    sep = " "))
                select[i] <- j
            }
        }
        if (length(term.labels)) 
            term.labels <- paste(term.labels, "+")
        if (!is.null(a$offset)) 
            term.labels <- paste(off1, term.labels, sep = " + ")
        return(list(response = paste(response, term.labels, sep = " ~ "), 
            select = select))
    }

   deviance.lm <- function (object, ...) {
if (is.null(w <- object$weights)) sum(object$residuals^2) else sum(w * 
    object$residuals^2)
}

   as.anova <- function (df, heading) 
{
    if (!inherits(df, "data.frame")) 
        stop("df must be a data frame")
    attr(df, "heading") <- heading
    if (inherits(df, "anova")) {
        dfClasses <- attr(df, "class")
        if (dfClasses[1] == "anova") 
            return(df)
    }
    class(df) <- unique(c("anova", class(df)))
    df
}

    make.step <- function(models, fit, scale, object) {
        chfrom <- sapply(models, "[[", "from")
        chfrom[chfrom == "1"] <- ""
        chto <- sapply(models, "[[", "to")
        chto[chto == "1"] <- ""
        dev <- sapply(models, "[[", "deviance")
        df <- sapply(models, "[[", "df.resid")
        ddev <- c(NA, diff(dev))
        ddf <- c(NA, diff(df))
        AIC <- sapply(models, "[[", "AIC")
        heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
            "\nInitial Model:", deparse(as.vector(formula(object))), 
            "\nFinal Model:", deparse(as.vector(formula(fit))), 
            paste("\nScale: ", format(scale), "\n", sep = ""))
        aod <- data.frame(From = chfrom, To = chto, Df = ddf, 
            Deviance = ddev, `Resid. Df` = df, `Resid. Dev` = dev, 
            AIC = AIC, check.names = FALSE)
        fit$anova <- as.anova(aod, heading)
        fit
      }

    
    direction <- match.arg(direction)
    if (missing(scope)) 
        stop("you must supply a scope argument to step.gam(); the gam.scope() function might be useful")
    if (!is.character(scope[[1]])) 
        scope <- lapply(scope, scope.char)
    response <- untangle.scope(object$terms, scope)
    form.y <- response$response
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    items <- response$select
    family <- family(object)
    Call <- object$call
    term.lengths <- sapply(scope, length)
    n.items <- length(items)
    visited <- array(FALSE, term.lengths)
    visited[array(items, c(1, n.items))] <- TRUE
    if (!is.null(keep)) {
        keep.list <- vector("list", length(visited))
        nv <- 1
    }
    models <- vector("list", length(visited))
    nm <- 2
    form.vector <- character(n.items)
    for (i in seq(n.items)) form.vector[i] <- scope[[i]][items[i]]
    form <- deparse(object$formula)
    if (trace) 
        cat("Start: ", form)
    fit <- object
    n <- length(fit$fitted)
    if (missing(scale.aic)) {
        famname <- family$family["name"]
        scale.aic <- switch(famname, Poisson = 1, Binomial = 1, deviance.lm(fit)/fit$df.resid)
    }
    else if (scale.aic == 0) 
        scale.aic <- deviance.lm(fit)/fit$df.resid
    bAIC <- fit$aic
    if (trace) 
        cat("; AIC=", format(round(bAIC, 4)), "\n")
    models[[1]] <- list(deviance = deviance(fit), df.resid = fit$df.resid, 
        AIC = bAIC, from = "", to = "")
    if (!is.null(keep)) {
        keep.list[[nv]] <- keep(fit, bAIC)
        nv <- nv + 1
    }
    AIC <- bAIC + 1
    while (bAIC < AIC & steps > 0) {
        steps <- steps - 1
        AIC <- bAIC
        bitems <- items
        bfit <- fit
        for (i in seq(n.items)) {
            if (backward) {
                trial <- items
                trial[i] <- trial[i] - 1
                if (trial[i] > 0 && !visited[array(trial, c(1, 
                  n.items))]) {
                  visited[array(trial, c(1, n.items))] <- TRUE
                  tform.vector <- form.vector
                  tform.vector[i] <- scope[[i]][trial[i]]
                  form <- paste(form.y, paste(tform.vector, collapse = " + "))
                  if (trace) 
                    cat("Trial: ", form)
                  tfit <- update(object, eval(parse(text = form)), 
                    trace = FALSE, ...)
                  tAIC <- tfit$aic
                  if (!is.null(keep)) {
                    keep.list[[nv]] <- keep(tfit, tAIC)
                    nv <- nv + 1
                  }
                  if (tAIC < bAIC) {
                    bAIC <- tAIC
                    bitems <- trial
                    bfit <- tfit
                    bform.vector <- tform.vector
                    bfrom <- form.vector[i]
                    bto <- tform.vector[i]
                  }
                  if (trace) 
                    cat("; AIC=", format(round(tAIC, 4)), "\n")
                }
            }
            if (forward) {
                trial <- items
                trial[i] <- trial[i] + 1
                if (trial[i] <= term.lengths[i] && !visited[array(trial, 
                  c(1, n.items))]) {
                  visited[array(trial, c(1, n.items))] <- TRUE
                  tform.vector <- form.vector
                  tform.vector[i] <- scope[[i]][trial[i]]
                  form <- paste(form.y, paste(tform.vector, collapse = " + "))
                  if (trace) 
                    cat("Trial: ", form)
                  tfit <- update(object, eval(parse(text = form)), 
                    trace = FALSE, ...)
                  tAIC <- tfit$aic
                  if (!is.null(keep)) {
                    keep.list[[nv]] <- keep(tfit, tAIC)
                    nv <- nv + 1
                  }
                  if (tAIC < bAIC) {
                    bAIC <- tAIC
                    bitems <- trial
                    bfit <- tfit
                    bform.vector <- tform.vector
                    bfrom <- form.vector[i]
                    bto <- tform.vector[i]
                  }
                  if (trace) 
                    cat("; AIC=", format(round(tAIC, 4)), "\n")
                }
              }
          }
        if (bAIC >= AIC | steps == 0) {
            if (!is.null(keep)) 
                fit$keep <- re.arrange(keep.list[seq(nv - 1)])
            return(make.step(models[seq(nm - 1)], fit, scale.aic, 
                object))
        }
        else {
            if (trace) 
                cat("Step : ", deparse(bfit$formula), "; AIC=", 
                  format(round(bAIC, 4)), "\n\n")
            items <- bitems
            models[[nm]] <- list(deviance = deviance(bfit), df.resid = bfit$df.resid, 
                AIC = bAIC, from = bfrom, to = bto)
            nm <- nm + 1
            fit <- bfit
            form.vector <- bform.vector
        }
    }
}


parse.sample.names <- function(ts.names){
for(i in 1:length(ts.names)){
 ts.names[i] <- paste(unlist(strsplit(unlist(strsplit(ts.names[i], "*_")), ".CEL")[-1]), collapse="_")
}
return(ts.names)
}

# predict a gene's expression profile using a CRM
predict.gene.express <- function(patterns, combined.mat, tfs, training.express, test.express, training.cols, test.express.cols, step.wise=FALSE, inter.effects=FALSE, bootstrap=FALSE, smoothing.k=-1){

# set parameters
training.params <- c()
test.params <- c()

for(i in 1:length(names(tfs))){

  param.express <- training.express[which(training.express[,"Gene"] == names(tfs)[i]),, drop=FALSE]
  if(dim(param.express)[1] > 1){
    param.express <- param.express[which.max(apply(param.express[, training.cols, drop=FALSE], 1, var)),]
  }
  training.params <- rbind(training.params, param.express)
 
   
  
  param.express <- test.express[which(test.express[,"Gene"] == names(tfs)[i]),, drop=FALSE]
  if(dim(param.express)[1] > 1){
    test.params <- rbind(test.params, param.express[which.max(apply(param.express[, training.cols, drop=FALSE], 1, var)),])
  }else{
    test.params <- rbind(test.params, param.express)
  }
}


 predictions <- c()

# include TF binding distance
#crm.dat <- analyze.CRM(combined.mat, tfs, genes)


    for(j in 1:dim(combined.mat)[1]){
      if(length(which(names(patterns) == combined.mat[j, "pattern"])) > 0) {
        tf.names <- pattern2tfs(patterns[which(names(patterns) == combined.mat[j, "pattern"])], names(tfs))[[1]]  
        target.gene <- combined.mat[j, "gene"]
      if(length(which(training.express[,"Gene"] == target.gene)) > 0 & length(which(test.express[,"Gene"] == target.gene)) > 0 ){

        if(bootstrap == TRUE){
          tf.names <- names(tfs)[sample(length(tfs), length(tf.names))]
     
        }

              
        tf.express <- training.params[match( toupper(tf.names), training.params[,"Gene"]),, drop=FALSE]
        rownames(tf.express) <- tf.names

   # Include TF binding distance
 # dist.mat <- as.matrix(dist(crm.dat[[j]]["center.pos",]))
 # diag(dist.mat) <- Inf
 # tf.dists <- apply(dist.mat, 1, min)
 # tf.dists <- tf.dists[!duplicated(names(tf.dists))]
 # tf.express[,training.cols] <- apply(tf.express[,training.cols], 2, function(x){ as.numeric(x)*exp(-tf.dists/100)})  
                        
        target.express <- training.express[which( training.express[,"Gene"] == target.gene),, drop=FALSE]
        target.express <- target.express[!duplicated(target.express[,"Transcript.Start..bp."]),, drop=FALSE]
        if(dim(target.express)[1] > 0){
           for(tar.iter in 1:dim(target.express)[1]){
             
             gam.obj <- NULL
             gam.obj <- train.gam(target.express[tar.iter, training.cols, drop=FALSE], tf.express[,training.cols, drop=FALSE], stepwise=step.wise, k=smoothing.k)

             test.tfs <-  test.params[match( toupper(tf.names), test.params[,"Gene"]),, drop=FALSE]
             if(length(gam.obj)>0){
               if(length( attr(terms(as.formula(gam.obj)), "term.labels")) > 0){
                 haem.params <- test.params[match( toupper(tf.names), test.params[,"Gene"]), test.express.cols, drop=FALSE]
                 haem.params <- t(apply(haem.params, 2, as.numeric, drop=FALSE))
                 colnames(haem.params) <- tf.names

               #  haem.params <- t(apply(haem.params, 1, function(x){ as.numeric(x)*exp(-tf.dists/100)}))
                 
                 tar.gene.pred <- predict.gam(gam.obj, as.data.frame(haem.params))

                 test.gene.ind <- which( test.express[,"Gene"] == target.gene)
                 test.gene.ind <- na.omit(test.gene.ind[!duplicated(test.express[test.gene.ind,"Transcript.Start..bp."])])
              
                 test.target <- test.express[test.gene.ind,,drop=FALSE]
                 for(k in 1:length(test.gene.ind)){
                   predictions <- rbind(predictions, c(test.gene.ind[k], test.target[k,"Gene"],test.target[k, "Transcript.Start..bp."], tar.gene.pred, try(summary(gam.obj)$r.sq, silent=TRUE)) )
                 }
               }
             }
           
           }
         }
      }
    }
  }
return(predictions)
}


predict.sample.express <- function(combined.mat, param.dat, express.mat, express.cols, sample.cols, training.genes, test.genes, step.wise=FALSE, inter.effects=FALSE, smoothing.k=60){

response.var <- c()
matched.genes <- c()
observed.express <- c()
for(j in 1:dim(combined.mat)[1]){
  
  tf.names <- colnames(combined.mat)[ which(combined.mat[j,1:length(tfs)] == 1)]
  target.gene <- combined.mat[j, "gene"]
  if(length(which(express.mat[,"Gene"] == target.gene)) > 0){
   
    target.express <- express.mat[which( express.mat[,"Gene"] == target.gene),, drop=FALSE]

    
    target.express <- target.express[which.max(apply(target.express[,express.cols, drop=FALSE], 1, var)),, drop=FALSE]
        
    target.express <- target.express[!duplicated(target.express[,"Transcript.Start..bp."]),, drop=FALSE]

    response.var <- c(response.var, as.numeric(target.express[,sample.cols] ))
    matched.genes <- c(matched.genes, j)
    observed.express <- rbind(observed.express, target.express)
  }
}
training.params <- param.dat[,intersect(matched.genes, training.genes), drop=FALSE]


gam.obj <- NULL
gam.obj <- train.gam(response.var, training.params, stepwise=step.wise, k=smoothing.k, inter.effects=inter.effects)

test.params <- t(param.dat[,intersect(matched.genes, test.genes), drop=FALSE])
predictions <- predict.gam(gam.obj, as.data.frame(test.params))
observed.express <- observed.express[na.omit(match(test.genes, matched.genes)),]

return(list(predictions=predictions, observations=observed.express, genes=intersect(matched.genes, test.genes)))
}


predict.gene.express.2 <- function(patterns, combined.mat, tfs, training.express, test.express, training.cols, test.express.cols, method="gam", step.wise=FALSE, inter.effects=FALSE, bootstrap=FALSE, smoothing.k=-1, null.tfs=NULL){

# set parameters
training.params <- c()
test.params <- c()

crm.tfs <- names(tfs)

if(bootstrap ==TRUE){
  if(!is.null(null.tfs)){
    crm.tfs <- null.tfs # for bootstrapping with Broad TFs
  }
}

for(i in 1:length(crm.tfs)){

  param.express <- training.express[which(training.express[,"Gene"] == crm.tfs[i]),, drop=FALSE]
  if(dim(param.express)[1] > 1){
    param.express <- param.express[which.max(apply(param.express[, training.cols, drop=FALSE], 1, var)),]
  }
  training.params <- rbind(training.params, param.express)
 
   
  
  param.express <- test.express[which(test.express[,"Gene"] == crm.tfs[i]),, drop=FALSE]
  if(dim(param.express)[1] > 1){
    if(dim(param.express[,test.express.cols, drop=FALSE])[2] < 2){
      test.params <- rbind(test.params, param.express[which.max(param.express[,test.express.cols]),])
    }
    else{ test.params <- rbind(test.params, param.express[which.max(apply(param.express[, test.express.cols, drop=FALSE], 1, var)),])}
  }else{
    test.params <- rbind(test.params, param.express)
  }
}

 predictions <- c()

    for(j in 1:dim(combined.mat)[1]){
      if(length(which(names(patterns) == combined.mat[j, "pattern"])) > 0) {
        tf.names <- pattern2tfs(patterns[which(names(patterns) == combined.mat[j, "pattern"])], names(tfs))[[1]]  
        target.gene <- combined.mat[j, "gene"]
      if(length(which(training.express[,"Gene"] == target.gene)) > 0 & length(which(test.express[,"Gene"] == target.gene)) > 0 ){

        if(bootstrap == TRUE){
          tf.names <- names(tfs)[sample(length(tfs), length(tf.names))]
          if(!is.null(null.tfs)){
            tf.names <- null.tfs[sample(length(null.tfs), length(tf.names))]
          }
        }

              
        tf.express <- training.params[match( tf.names, training.params[,"Gene"]),, drop=FALSE]
        rownames(tf.express) <- tf.names

        target.express <- training.express[which( training.express[,"Gene"] == target.gene),, drop=FALSE]

        if(dim(target.express)[1] > 0){
           for(tar.iter in 1:dim(target.express)[1]){
             
             model.obj <- NULL
             training.response <- as.numeric(target.express[tar.iter, training.cols, drop=FALSE])
             covariate.params <- t(tf.express[,training.cols, drop=FALSE])
             class(covariate.params) <- "numeric"
             trainset <- cbind(training.response, covariate.params)
             
             if(method == "gam"){
               model.obj <- train.gam(target.express[tar.iter, training.cols, drop=FALSE], tf.express[,training.cols, drop=FALSE], stepwise=step.wise, k=smoothing.k)
             }
             if(method == "svm"){
               model.obj <- try(svm( training.response ~ . , data=trainset, cost=1, gamma=1e-04))
             }
             if(method == "lasso"){
               model.obj <- try(lars( covariate.params , training.response, type="lasso", use.Gram=FALSE ), silent=T)
             }
             if(method == "RF"){
               model.obj <- try(randomForest( training.response ~ . , data=trainset))
             }
             
             
             test.tfs <-  test.params[match( tf.names, test.params[,"Gene"]),, drop=FALSE]
             if(length(model.obj)>0){
              
                 haem.params <- t(test.params[match( tf.names, test.params[,"Gene"]), test.express.cols, drop=FALSE])
                 class(haem.params) <- "numeric"
                 colnames(haem.params) <- tf.names
                 
                 r.sqr.fit <- NULL
                 tar.gene.pred <- NULL
                 if(method == "gam"){
                   tar.gene.pred <- predict.gam(model.obj, as.data.frame(haem.params))
                    r.sqr.fit <- try(summary(model.obj)$r.sq, silent=TRUE)
                 }
                 if(method == "svm"){
                   tar.gene.pred <- predict(model.obj, haem.params)
                 }
                 if(method == "lasso"){
                   tar.gene.pred <- predict.lars(model.obj, haem.params, type="fit")$fit[,as.numeric(names(which.max(model.obj$R2)))]
                   r.sqr.fit <- model.obj$R2[length(model.obj$R2)]
                 }
                 if(method == "RF"){
                   tar.gene.pred <- predict(model.obj, haem.params)
                   r.sqr.fit <- model.obj$rsq[length(model.obj$rsq)]
                 }
                                
                 if(length(r.sqr.fit) == 0){r.sqr.fit <- NA}

                 test.gene.ind <- which( test.express[,"Gene"] == target.gene)[tar.iter]
                 
                 test.target <- test.express[test.gene.ind,,drop=FALSE]
                 
                 for(k in 1:length(test.gene.ind)){
                   predictions <- rbind(predictions, c(j, test.gene.ind[k], test.target[k,"Gene"], tar.gene.pred, r.sqr.fit))
                 }
               
             }
           
           }
         }
      }
    }
  }
if(length(predictions) > 0){
colnames(predictions)[c(1:3, dim(predictions)[2])] <- c("CRM_ind", "target_expression_ind", "target_gene", "R.squared_fit")
}
return(predictions)
}
