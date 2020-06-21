#' @references Zilin Li, Xihao Li, Yaowu Liu, Jincheng Shen, Han Chen, Hufeng Zhou, 
#' Alanna C. Morrison, Eric Boerwinkle, and Xihong Lin (2019) "Dynamic Scan Procedure
#' for Detecting Rare-Variant Association Regions in Whole-Genome Sequencing Studies". 
#' The American Journal of Human Genetics, 104(5), 802-814.
#' @references Gogarten, S.M., Sofer, T., Chen, H., Yu, C., Brody, J.A., Thornton, 
#' T.A., Rice, K.M., and Conomos, M.P. (2019). Genetic association testing using 
#' the GENESIS R/Bioconductor package. Bioinformatics.
eSCAN_thres <- function(genotype, nullmod, new_enloc, fam=0, times, weights){

  seed <- 123
  G <- genotype
  ## variants number
  p <-  ncol(G)
  ## sample size
  n <-  nrow(G)

  ## number of enhancers
  n_en <- nrow(new_enloc)

  ## X is design matrix
  X <- nullmod$model.matrix
  ## covariates number
  q <-  ncol(X)

  ## number of weights used
  w_num <- ncol(weights)


  working <- nullmod$working
  sigma <- sqrt(nullmod$varComp)


  ## create a vector to record the p-value
  CCT_p <- matrix(1, nrow = n_en, ncol = times)

  for (i in 1:n_en) {
        begin <- new_enloc[i, 4]
        end <- new_enloc[i, 5]

        G_sub <- G[, begin:end] # sub genotype matrix
        weight_sub <- weights[begin:end, ] # sub weights matrix
        set.seed(seed)
        x_sub <- compx(as.matrix(G_sub), X, working, sigma, fam, times) # sub pseudo-score

        Cov_sub <- compCov(as.matrix(G_sub), X, working, sigma, fam)

        p_sub <- ncol(G_sub)
        Covw_sub <- compCovw(p_sub, w_num, Cov_sub, weight_sub)

        for (j in 1:times) {
          Q_sub <- rep(0, w_num)
          pval_temp <- rep(1, w_num)
          for (k in 1:w_num) {
            Q_sub[k] <- sum((weight_sub[,k]*x_sub[j,])^2)
            if(p_sub<200){
              pv <- regular(Q_sub[k], Covw_sub[((k-1)*p_sub+1):(k*p_sub), ], p_sub)
            } else{
              pv <- fastH(Q_sub[k], Covw_sub[((k-1)*p_sub+1):(k*p_sub), ], 150)
            }
            pval_temp[k] <- pv$pval
          }

          # cauchy method
          CCT_weights <- rep(1, w_num)
          CCT_p[i,j] <- CCT_pval(pval_temp, CCT_weights)
       }

  }
  return(CCT_p)
}
