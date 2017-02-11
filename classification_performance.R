library(pastecs)

binary_eval <- function(pred,labels, cutoff='naive', repar=TRUE, ...) {
  # Various binary classification evaluation plots and metrics
  library(ROCR)
  # plot(performance(prediction(pred,y),'acc'))
  rocr_pred = prediction(pred,labels)
  acc = performance(rocr_pred,'acc')
  f1 = performance(rocr_pred,'f')
  auc = performance(rocr_pred,'auc')@y.values[[1]]
  roc = performance(rocr_pred,'rec','spec')
  bac = if (rocr_pred@n.pos[[1]] != rocr_pred@n.neg[[1]])
      sapply(1:length(roc@x.values[[1]]), function(i)
        mean(c(roc@x.values[[1]][i], roc@y.values[[1]][i])))
    else
      rep(-1,length(pred))
  # sensspec = performance(rocr_pred,'rec','spec')
  pr_curve = performance(rocr_pred,'prec','rec')
  rp_curve = performance(rocr_pred,'rec','prec')

  printf("AUC = %.3f\n", auc)

  if (cutoff=='naive') {
    if (all(pred>=0) & all(pred<=1)) {
      printf("Predictions seem to be probabilities, so ")
      cutoff = 0.5
    } else if (any(pred<0) & any(pred>0)) {
      printf("Predictions seem to be real-valued scores, so ")
      cutoff = 0
    } else {
      warning("cant tell what naive cutoff should be")
      cutoff = NULL
    }
    printf("using naive cutoff %s:\n", cutoff)
  } else if (class(cutoff)=='character') {
    printf("Using %s-best cutoff ", cutoff)
    if (cutoff=='bac') {
      perf = NULL
      perf_y = bac
    } else {
      perf = performance(rocr_pred, cutoff, ...)
      perf_y = perf@y.values[[1]]
    }
    cutoff_ind = which.max(perf_y)
    cutoff = if (cutoff=='prbe') perf@x.values[[1]][1] else rocr_pred@cutoffs[[1]][cutoff_ind]
    printf("%f\n", cutoff)
  } else {
    printf("For cutoff %s:\n", cutoff)
  }
  cutoff_ind = last(which(rocr_pred@cutoffs[[1]] >= cutoff))

  if (repar) par(mfrow=c(2,2))

  pp = function(perf) {
    if (length(cutoff_ind)>0 && is.finite(cutoff_ind)) {
      x=perf@x.values[[1]][cutoff_ind]
      y=perf@y.values[[1]][cutoff_ind]
      points(x,y, col='blue')
      linelight(x,y, col='lightblue')
    }
  }
  plot(acc); pp(acc)
  plot(f1); pp(f1)
  plot(roc); pp(roc)
  abline(a=1,b=-1,lty='dashed',col='gray')
  legend('bottomleft',legend=sprintf("AUC = %.3f",auc))
  plot(rp_curve); pp(rp_curve)
  pp = function(ind,...) points(rp_curve@x.values[[1]][ind], rp_curve@y.values[[1]][ind], ...)
  best_f1 = which.max(f1@y.values[[1]])
  pp(best_f1, pch=2,col='green')
  f05 = performance(rocr_pred,'f',beta=0.5)
  best_f05 = which.max(f05@y.values[[1]])
  pp(best_f05,pch=2,col='green')
  f2 = performance(rocr_pred,'f',beta=2)
  best_f2 = which.max(f2@y.values[[1]])
  pp(best_f2,pch=2,col='green')

  prbe = performance(rocr_pred,'prbe')@y.values[[1]]
  linelight(prbe,prbe,col='lightgray')

  # printf("Acc = %.3f\n", mean((pred >= cutoff) == (labels > 0)))
  printf("Acc %.3f, ", acc@y.values[[1]][cutoff_ind])

  printf(" F %.3f, Prec %.3f, Rec %.3f, Spec %.3f",
    f1@y.values[[1]][cutoff_ind],
    pr_curve@y.values[[1]][cutoff_ind],
    pr_curve@x.values[[1]][cutoff_ind],
    roc@x.values[[1]][cutoff_ind])
  # printf(" Prec = %.3f\n", pr_curve@y.values[[1]][cutoff_ind])
  # printf(" Rec = %.3f\n", pr_curve@x.values[[1]][cutoff_ind])
  # printf(" Spec = %.3f\n", roc@x.values[[1]][cutoff_ind])

  if (bac[1] != -1)
    printf(", BalAcc %.3f", mean(bac))
  printf("\n")


  invisible(rocr_pred)
  return(c( "Accuracy" = acc@y.values[[1]][cutoff_ind], "Precision/Sensitivity"= f1@y.values[[1]][cutoff_ind],
    "Precision" = pr_curve@y.values[[1]][cutoff_ind],
    "Sensitivity" = pr_curve@x.values[[1]][cutoff_ind],
    "Specificity" = roc@x.values[[1]][cutoff_ind], "Balanced Accuracy"= mean(bac)))
}

printf <- function(...) cat(sprintf(...))

linelight <- function(x,y, lty='dashed', col='lightgray', ...) {
  # highlight a point with lines running to the axes.
  left = par('usr')[1]
  bot = par('usr')[3]
  segments(left,y, x,y, lty=lty, col=col, ...)
  segments(x,bot, x,y, lty=lty, col=col, ...)
}
