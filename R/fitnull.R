#' @title fitNull
#' @description Fit generalized linear model under the null hypothesis. It's a wrapper of fitNullModel from GENESIS package.
#' @param x a data frame containing outcome variable and covariates
#' @param outcome a character string specifying the name of outcome variable in x
#' @param covars a vector of character strings specifying the name of covariates in x
#' @param fam a description of the error distribution and link function to be used in the model. Can only be either "gaussian" for continuous phenotype or "binomial" for binary phenotype
#' @return The function returns an object of model fit from \code{fitNullModel}. See \code{fitNullModel} in GENESIS package for more details.
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
