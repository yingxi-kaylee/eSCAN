#' @title eSCAN
#' @description scan the enhancer regions
#' @param genotype n*p genotype matrix
#' @param nullmod fitted null model from GENESIS package
#' @param new_enloc NULL or matrix whose number of rows is enhancer number, and has 6 columns, chr No., start_pos, end_pos, start_adj, end_adj, length
#' @param gap threshold to split independent loci, with default value 10^4
#' @param times simulation times for MC, with default value=1000
#' @param alpha significance level
#' @param analy TRUE indicates analytical threshold, FALSE indicates empirical threshold
#' @return a list containing 'res', 'res0', 'thres_mat', 'thres'
#' @export
#' @import Rcpp
#' @import RcppArmadillo
#' @import Matrix
eSCAN <- function(genotype, nullmod, new_enloc=NULL, gap=10^4, times = 1000, alpha=0.05, analy=FALSE)
{
  if(nullmod$family[1] == "binomial")
  {
    fam <- 1
  }
  if(nullmod$family[1] == "gaussian")
  {
    fam <- 0
  }

  new_list <- preprocess(genotype, new_enloc, gap)
  genotype <- new_list[[1]]
  MAF <- new_list[[2]]
  new_enloc <- new_list[[3]]

  # calculate weights
  weights <- cbind(dbeta(MAF,1,1),dbeta(MAF,1,25))

  # creat a n_en*times matrix to keep pesudo-pvalue
  thres_mat <- matrix(0,nrow=nrow(new_enloc),ncol=times)

  # get the pseudo p-value
  thres_mat <- eSCAN_THRES(genotype, nullmod, new_enloc, fam, times, weights)

  thres_temp <- apply(thres_mat, 2, min)
  thres <- quantile(thres_temp,alpha,na.rm = TRUE)

  if(analy){
    Get.Gumbel.M<-function(re.Q){
      re.mean<-mean(re.Q,na.rm = TRUE)
      re.variance<-var(re.Q,na.rm = TRUE)
      beta<-sqrt(re.variance*6/pi^2)
      mu<-re.mean-beta*0.5772156649 #Euler Mascheroni constant

      threshold.Q<-mu-log(-log(1-0.05))*beta
      M<-0.05/exp(-threshold.Q)
      return(M)
    }
    M <- Get.Gumbel.M(-log(thres_temp))
    thres <- alpha/M
  }
  # calculate the pvalue using real phenotype
  pvalue <- eSCAN_SEARCH(genotype, nullmod, new_enloc, fam, weights)

  # keep track of p-value for each enhancer
  res0 <- cbind(pvalue, new_enloc[, 2:3])

  res_en <- res0[res0[,1] < thres, ]
  if(length(res_en)==0)
  {
    res_en <- c(0,0,0,0)
  }

  Lst <- list(res=res_en, res0=res0, thres=thres)
  return(Lst)
}
