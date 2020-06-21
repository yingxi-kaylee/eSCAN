#' @references Gogarten, S.M., Sofer, T., Chen, H., Yu, C., Brody, J.A., Thornton, 
#' T.A., Rice, K.M., and Conomos, M.P. (2019). Genetic association testing using 
#' the GENESIS R/Bioconductor package. Bioinformatics.
testVariantSetSKAT <- function(nullmod, G, weights, neig = 200, ntrace = 500, verbose = FALSE){

  G <- genoAsMatrix(nullmod, G)

  # multiply G by weights
  if(is(G, "Matrix")){
    G <- G %*% Diagonal(x = weights)
  }else{
    G <- t(t(G) * weights)
  }

  # scores
  U <- as.vector(crossprod(G, nullmod$resid)) # WGPY
  # SKAT test statistic
  Q <- sum(U^2)

  # adjust G for covariates and random effects
  G <- calcGtilde(nullmod, G) # P^{1/2}GW

  # compute the p-value
  out <- calcPvalVCTest(Q = Q, G = G, neig = neig, ntrace = ntrace, verbose = verbose)

  return(list(Q = Q, pval = out$pval, err = out$err, pval.method = out$pval.method))
}


calcGtilde <- function(nullmod, G){
  C <- nullmod$cholSigmaInv

  if(length(C) > 1){ # n by n cholSigmaInv (may be Diagonal)
    CG <- crossprod(C, G)
  }else{ # cholSigmaInv is a scalar
    CG <- C*G
  }

  # calculate Gtilde
  nrowG <- as.numeric(nrow(CG))
  ncolG <- as.numeric(ncol(CG))
  if(length(C) == 1 || nrowG*ncolG <= 2^31){
    Gtilde <- CG - tcrossprod(nullmod$CXCXI, crossprod(CG, nullmod$CX))
    # base::qr.resid(nullmod$qr, CG) # QR seems to be slower unexpectedly

  }else{
    # too large when G sparse; break into multiple blocks
    nblock <- ceiling(nrowG*ncolG/2^31)
    blocks <- unname(split(1:ncolG, cut(1:ncolG, nblock)))
    Gtilde <- list()
    for(i in 1:length(blocks)){
      Gtilde[[i]] <- as.matrix(CG[,blocks[[i]]] - tcrossprod(nullmod$CXCXI, crossprod(CG[,blocks[[i]]], nullmod$CX)))
    }
    Gtilde <- do.call(cbind, Gtilde)
  }

  return(Gtilde)
}


calcPvalVCTest <- function(Q, G, neig = 200, ntrace = 500, verbose = FALSE){
  if(!requireNamespace("survey")) stop("package 'survey' must be installed to calculate p-values for SKAT or SMMAT")
  if(!requireNamespace("CompQuadForm")) stop("package 'CompQuadForm' must be installed to calculate p-values for SKAT or SMMAT")

  ncolG <- ncol(G) # number of snps
  nrowG <- nrow(G) # number of samples
  if (verbose) message('nsamp = ', nrowG, '; nsnp = ', ncolG)

  if(min(ncolG, nrowG) < 6000 + 20*neig){
    if(ncolG <= nrowG){
      V <- crossprod(G) # WGPGW
    }else{
      V <- tcrossprod(G) # same eigenspace but smaller matrix
    }
    if(mean(abs(V)) < sqrt(.Machine$double.eps)){
      return(list(pval = NA_real_, pval.method = NA_character_, err = 1))
    }

    if(min(ncolG, nrowG) < 2*neig){
      # use "regular" method
      pv <- regular(Q, V, ncolG)

    }else{
      # use "fast H" method
      if (verbose) message("using method fast_H")
      pv <- fastH(Q, V, neig)
    }

  }else{
    # use "fast G" method
    if (verbose) message("using method fast_G")
    pv <- fastG(Q, G, neig, ntrace)
  }

  return(list(pval = pv$pval, pval.method = pv$method, err = pv$err))
}

regular <- function(Q, V, ncolG) {
  if(ncolG == 1){
    pv <- list(pval = pchisq(as.numeric(Q/V), df=1, lower.tail=FALSE), method = "integration")

  }else{
    lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
    pv <- pchisqsum(x = Q, df = rep(1, length(lambda)), a = lambda)
  }
  pv$err <- ifelse(is.na(pv$pval), 1, 0)
  return(pv)
}


fastH <- function(Q, V, neig) {
  pv <- list(pval = NA_real_, method = NA_character_)
  pval.try = 0
  while(is.na(pv$pval) & pval.try < 10){
    pv <- tryCatch( pchisqsum_ssvd(x = Q, M = as.matrix(V), n = neig, p = 10, q = 1),
                    error = function(e){ list(pval = NA_real_, method = "error") } )
    pval.try <- pval.try + 1
  }
  pv$method <- paste0('ssvd_', pv$method)
  if(is.na(pv$pval)){
    err <- 1
  }else if(pval.try > 1){
    err <- 2
  }else{
    err <- 0
  }
  pv[["err"]] <- err
  return(pv)
}



fastG <- function(Q, G, neig, ntrace) {
  pv <- list(pval = NA_real_, method = NA_character_)
  pval.try = 0
  while(is.na(pv$pval) & pval.try < 10){
    pv <- tryCatch( pchisqsum_rsvd(x = Q, M = as.matrix(G), n = neig, p = 10, q = 3, tr2.sample.size = ntrace),
                    error = function(e){ list(pval = NA_real_, method = "error") }   )
    pval.try <- pval.try + 1
  }
  pv$method <- paste0('rsvd_', pv$method)
  if(is.na(pv$pval)){
    err <- 1
  }else if(pval.try > 1){
    err <- 2
  }else{
    err <- 0
  }
  pv[["err"]] <- err
  return(pv)
}


pchisqsum_rsvd <- function(x,M,n=100,p=10,q=2, tr2.sample.size=100){
  Q<-srfht2(M,n+p,q=q)
  B<-M%*%Q
  ## svd decomposition
  ee<-svd(B,nu=0,nv=0)$d[1:n]^2
  diags <- colSums(M^2)
  tr<-sum(diags)
  if (tr2.sample.size>0){
    tr2<-tracefht(M,k=tr2.sample.size,trace.full=tr)
  } else {
    Ms<-crossprod(M)
    tr2<- sum(Ms^2)
  }
  tr.small<-tr-sum(ee)
  tr2.small<-tr2-sum(ee^2)
  scale<-tr2.small/tr.small
  nu<-(tr.small^2)/tr2.small
  pchisqsum(x, c(rep(1,n), nu), c(ee, scale))
}


pchisqsum_ssvd <- function(x,M,n=100,p=10,q=0){
  Q<-srfht(M,n+p,q=q)
  B<-t(Q)%*%M%*%Q
  ## eigen decomposition
  ee<-eigen(B,symmetric=TRUE,only.values=TRUE)$values[1:n]
  tr<-sum(diag(M))
  tr2<-sum(M^2)
  tr.small<-tr-sum(ee)
  tr2.small<-tr2-sum(ee^2)
  scale<-tr2.small/tr.small
  nu<-(tr.small^2)/tr2.small
  pchisqsum(x, c(rep(1,n),ceiling(nu)), c(ee, scale))
}

srfht <- function(A,k,q=0){
  m<-NROW(A)
  mbig<-2^ceiling(log2(m))
  R<-sample(c(-1,1),m,replace=TRUE)
  Astar<-A*R
  idx<-sample(mbig,k)
  AOmega<-fht(Astar)[idx,]/sqrt(k)
  Q<-qr.Q(qr(t(AOmega)))
  for(i in seq_len(q)){
    tildeQ<-qr.Q(qr(crossprod(A,Q)))
    Q<-qr.Q(qr(A%*%tildeQ))
  }
  Q
}


srfht2 <-function(A,k,q=0){
  m<-NROW(A)
  mbig<-2^ceiling(log2(m))
  R<-sample(c(-1,1),m,replace=TRUE)
  Astar<-A*R
  idx<-sample(mbig,k)
  AOmega<-fht(Astar)[idx,]/sqrt(k)
  Q<-qr.Q(qr(t(AOmega)))
  for(i in seq_len(q)){
    tildeQ<-qr.Q(qr(A%*%Q))
    Q<-qr.Q(qr(crossprod(A,tildeQ)))
  }
  Q
}

tracefht <- function(A,k,trace.full=NULL){
  m<-NROW(A)
  mbig<-2^ceiling(log2(m))
  R<-sample(c(-1,1),m,replace=TRUE)
  Astar<-A*R
  idx<-sample(mbig,k)
  AOmega<-fht(Astar)[idx,]
  AAOmega<-tcrossprod(AOmega,A)
  if (!is.null(trace.full))
    tr<-sum(rowSums(AOmega*AOmega))/k
  trsquared<-sum(rowSums(AAOmega*AAOmega))/k
  if (is.null(trace.full))
    trsquared
  else
    trsquared*(trace.full/tr)^2
}

fht<-function (A,big=TRUE) {
  m <- NROW(A)
  mbig <- 2^ceiling(log2(m))
  Abig <- matrix(0, nrow = mbig, ncol = NCOL(A))
  Abig[1:m, ] <- A
  if (big){  return(big_mfwht(Abig)) }
  else { return(mfwht(Abig, as.integer(mbig), as.integer(NCOL(Abig)))) }
}


genoAsMatrix <- function(nullmod, G) {
  if (is(nullmod$cholSigmaInv, "Matrix") & is.matrix(G)) {
    G <- Matrix(G)
  }
  G
}



pchisqsum <- function(x, df, a){

  ## check for bad.df
  ## can happen with randomised trace estimator if most remaining singular values are very small
  ## leads to unreliable p-values
  if(any(df < 1)){
    stop("Negative/fractional df")
  }

  df<-round(df)

  ## try integration
  f <- suppressWarnings(CompQuadForm::davies(x, a, df, acc = 1e-9))
  if((f$ifault > 0) | (f$Qq < 1e3*.Machine$double.eps) | (f$Qq > 1)){
    ## try saddlepoint
    pval <- survey:::saddle(x, rep(a, df))
    #pval <- Saddle(x, rep(a, df))
    method <- "saddlepoint"

    #if(pval==2){
    #  f <- suppressWarnings(CompQuadForm::davies(x, a, df, acc = 1e-4))
    #  pval <- f$Qq
    #  method <- "integration"
    #}

  }else{
    pval <- f$Qq
    method <- "integration"
  }

  return(list(pval = pval, method = method))
}
