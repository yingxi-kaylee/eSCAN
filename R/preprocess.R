#' @title Data preprocessing
#' @description This function is a preliminary step for subsequent analysis using eSCAN.
#' @param geno an n*p genotype matrix, where n is the sample size and p is the number
#' of rare varints included.
#' @param enhancer a data frame of annotation information with dimension q*6, where
#' q is the number of the target regulatory regions. The six columns indicates no., physical
#' start position, physical end position, start position in the \code{geno} matrix, 
#' end position in the \code{geno} matrix and region length, respectively. 
#' If annotation information is not specified by users (default=NULL), a data frame of 
#' enhancer location will be automatically created by \code{\link{eGenerator}} function.
#' @param gap threshold to split independent loci (default=10^4)
#' @return The functoin returns a list with the following elements:
#' \code{genotype}: An n*p genotype matrix by increasing order of variant position.
#' \code{MAF}: A vector (length=p) of the minor allele frequency.
#' \code{new_enloc}: A data frame of enhancer location of dimension q*6, where
#' q is the number of the target regulatory regions. The six columns indicates No., physical
#' start position, physical end position, start position in the \code{genotype} matrix, 
#' end position in the \code{genotype} matrix and region length, respectively.
#' @export
preprocess <- function(geno, enhancer, gap){
  # re-arrange the genotype matrix by the SNP position
  numSNP <- dim(geno)[2]
  posiSNP <- as.numeric(colnames(geno))

  orde <- order(posiSNP)
  posiSNP <- posiSNP[orde]
  geno <- geno[,orde]
  colnames(geno) <- as.character(posiSNP)

  # compute MAF
  MAF <- colMeans(geno)/2

  #generate enhancer list
  if(is.null(enhancer)){
    enhancer <- as.matrix(eGenerator(geno, gap))
  }

  return(list(genotype=geno, MAF=MAF, new_enloc=enhancer))

}
