#' @title Fit generalized linear model under the null hypothesis for unrelated samples.
#' @description The \code{fitNull} function is a wrapper of \code{\link[GENESIS]{fitNullModel}} from the 
#' \code{\link[GENESIS]{GENESIS}} package. It fits a regression model under the null hypothesis for 
#' unrelated samples which is the preliminary step for subsequent analysis.
#' @param x a data frame containing outcome variable and covariates.
#' @param outcome a character string specifying the name of outcome variable in \code{x}.
#' @param covars a vector of character strings specifying the name of covariates in \code{x}.
#' @param fam Can be either "gaussian" for continuous phenotype or "binomial" for binary phenotype.
#' @return The function returns an object of model fit from \code{\link[GENESIS]{fitNullModel}}. See
#' \code{\link[GENESIS]{fitNullModel}} in the \code{\link[GENESIS]{GENESIS}} package for more details.
#' @references Gogarten, S.M., Sofer, T., Chen, H., Yu, C., Brody, J.A., Thornton, 
#' T.A., Rice, K.M., and Conomos, M.P. (2019). Genetic association testing using 
#' the GENESIS R/Bioconductor package. Bioinformatics.
#' @export
#' @importFrom GENESIS fitNullModel
fitNull <- function(x, outcome=NULL, covars=NULL, fam = "gaussian"){
  if((fam != "gaussian") & (fam != "binomial")){
    stop("'fam' can only be gaussian or binomial!")
  }

  if(!is.data.frame(x)){
    stop("'x' must be a data frame!")
  }

  if(is.null(outcome)){
    outcome <- colnames(x)[1]
  }

  if(is.null(covars)){
    covars <- colnames(x)[2:ncol(x)]
  }

  col_name <- colnames(x)

  if(!(outcome %in% col_name)){
    stop("Please give a valid name of outcome variable!")
  }

  if(!(outcome %in% col_name)){
    stop("Undefined covariates columns selected!")
  }

  nullmod <- fitNullModel(x, outcome = outcome, covars=covars, family=fam)

  if(fam=="binomial"){
    X <- nullmod$model.matrix
    logist <- glm(x[,outcome]~-1 + X, family = binomial(link = logit))
    nullmod$working <- logist$weights
    nullmod$varComp <- 1
  } else{
    samplesize <- nrow(x)
    nullmod$working <- rep(1, samplesize)
  }

  return(nullmod)
}
