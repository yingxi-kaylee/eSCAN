#' @title Scan the enhancers
#' @description The \code{eSCAN} function is the main function in the package. It takes in 
#' genotype matrix, null model object and annotation information which could be specified 
#' by users, and then detect rare variant association between a quantitative/dichotomous 
#' phenotype and a regulatory region in a sequence by using eSCAN procedure.
#' @param genotype an n*p genotype matrix, where n is the sample size and p is the number
#' of rare varints included.
#' @param nullmod a null model object returned from \code{\link{fitNull}}. Note that 
#' \code{\link{fitNull}} is a wrapper of \code{\link[GENESIS]{fitNullModel}} function from the 
#' \code{\link[GENESIS]{GENESIS}} package.
#' @param new_enloc a data frame of annotation information with dimension q*6, where
#' q is the number of the target regulatory regions. The six columns indicates no., physical
#' start position, physical end position, start position in the \code{genotype} matrix, 
#' end position in the \code{genotype} matrix and region length, respectively. 
#' If annotation information is not specified by users (default=NULL), a data frame of 
#' enhancer location will be automatically created by \code{\link{eGenerator}}.
#' @param gap if \code{new_enloc} is not specified by users, this parameter will be used 
#' to generate enhancer locations in functoin \code{\link{eGenerator}}, where \code{gap}
#' is the threshold to split independent loci (default=10^4)
#' @param times simulation times for MC (default=1000)
#' @param alpha significance level (default=0.05)
#' @param analy TRUE indicates analytical threshold, FALSE indicates empirical threshold
#'  (default=FALSE)
#' @return The function returns a list containing the following elements:
#' @return \code{res}: A matrix of significant regions detected by eSCAN.
#' The first column is the p-value of the detected region.
#' The next two columns are the location of the detected rgeion (physical position on
#'  chromosome).
#' @return \code{res0}: A matrix to summarise all the regions included in the analysis.
#' The first column is the p-value of the detected region.
#' The next two columns are the location of the detected rgeion (physical position on
#' chromosome).
#' @return \code{thres}: threshold of eSCAN to control the family-wise/genome-wide error
#' at \code{alpha} level.
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

  if(times>0){
    # creat a n_en*times matrix to keep pesudo-pvalue
    thres_mat <- matrix(0,nrow=nrow(new_enloc),ncol=times)
    
    # get the pseudo p-value
    thres_mat <- eSCAN_thres(genotype, nullmod, new_enloc, fam, times, weights)
    
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
    pvalue <- eSCAN_search(genotype, nullmod, new_enloc, fam, weights)
    
    # keep track of p-value for each enhancer
    res0 <- cbind(pvalue, new_enloc[, 2:3])
    
    res_en <- res0[res0[,1] < thres, ]
    if(length(res_en)==0)
    {
      res_en <- c(1,0,0)
    }
    
  } else{
    # for test purposes only
    # calculate p-values
    pvalue <- eSCAN_SEARCH(genotype, nullmod, new_enloc, fam, weights)
    
    res0 <- cbind(pvalue, new_enloc[, 2:3])
    res_en <- NULL
    thres <- NULL
  }
 

  Lst <- list(res=res_en, res0=res0, thres=thres)
  return(Lst)
}
