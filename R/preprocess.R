#' @title Data preprocessing
#' @description This function is data preparation for subsequent analysis using eSCAN.
#' @param geno an n*p genotype matrix, where n is the sample size and p is the number
#' of rare variants included.
#' @param enhancer a data frame of annotation information with dimension q*6, where
#' q is the number of candidate regulatory regions. The six columns indicate index, 
#' start position, end position, start index, end index
#' the \code{geno} matrix sorted by genomic positions and length of the enhancers, respectively. 
#' If annotation information is not specified by users (default=NULL), a data frame of 
#' locations of the enhancers will be automatically created by \code{\link{eGenerator}}.
#' @param gap threshold to split independent loci (default=10^4).
#' @return The function returns a list with the following elements:
#' @return \code{genotype}: an n*p genotype matrix sorted by genomic positions.
#' @return \code{MAF}: a vector (length=p) of minor allele frequencies.
#' @return \code{new_enloc}: a data frame of locations of the enhancers of dimension q*6, where
#' q is the number of candidate regulatory regions. The six columns indicate index, 
#' start position, end position, start index, end index 
#' in the \code{genotype} matrix sorted by genomic postions and length of the enhancers, respectively.
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
