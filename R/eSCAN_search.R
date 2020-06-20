#' @title eSCAN_SEARCH
#' @description compute p-value for each enhancer
#' @param genotype an n*p genotype matrix (dosage matrix)
#' @param nullmod an object from fitNullModel() function
#' @param new_enloc an m*3 matrix of enhancer information
#' @param fam a logical value to indicate phenotype family: 1 indicates binomial, 0 indicates gaussian
#' @param weights vector of weights for variants of interest
#' @return a vector with len=n_en, listing p-value for each enhancer
#' @import Rcpp
#' @import RcppArmadillo
#' @import Matrix
eSCAN_SEARCH <- function(genotype, nullmod, new_enloc, fam=0, weights){

  G <- genotype
  ## variants number
  p <-  ncol(G)
  ## sample size
  n <-  nrow(G)

  ## X is design matrix
  X <- nullmod$model.matrix

  ## number of weights used
  w_num <- ncol(weights)

  n0 <- nrow(new_enloc) # enhancer number

  ## create a vector to record test statistics Q
  Q <- rep(0, n0)

  ## create a vector to record the p-value
  CCT_p0 <- matrix(1, nrow = w_num, ncol = n0)

  for (k in 1:w_num) {
    for (i in 1:n0) {
      begin <- new_enloc[i, 4]
      end <- new_enloc[i, 5]

      G_sub <- G[, begin:end] # sub genotype matrix
      weight_sub <- weights[begin:end, k] # sub weights

      # multiply G_sub by weight_sub
      if(is(G_sub, "Matrix")){
        G_sub <- G_sub %*% Diagonal(x = weight_sub)
      }else{
        G_sub <- t(t(G_sub) * weight_sub)
      }

      # scores
      U <- as.vector(crossprod(G_sub, nullmod$resid)) # WGPY
      # SKAT test statistic
      Q[i] <- sum(U^2)

      # adjust G for covariates and random effects
      G_sub <- calcGtilde(nullmod, G_sub) # P^{1/2}GW

      if(ncol(G_sub) <= nrow(G_sub)){
        V <- crossprod(G_sub) # WGPGW
      }else{
        V <- tcrossprod(G_sub) # same eigenspace but smaller matrix
      }

      set.seed(123)
      if(min(ncol(G_sub), nrow(G_sub)) < 200){
        pv <- regular(Q[i], V, ncol(G_sub))
      } else{
        pv <- fastH(Q[i], V, 120)
      }
      CCT_p0[k, i] <- pv$pval
    }
  }

  ## cauchy method
  pval0 <- rep(1, n0)
  CCT_weights <- rep(1, w_num)
  for (j in 1:n0) {
    pval0[j] <- CCT_pval(CCT_p0[,j], CCT_weights)
  }

  return(pval0)
}
