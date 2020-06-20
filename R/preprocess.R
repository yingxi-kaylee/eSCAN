#' @title preprocess
#' @description pre-process data for eSCAN
#' @param geno n*p genotype matrix
#' @param enhancer NULL or matrix whose number of rows is enhancer number, and has 6 columns, chr No., start_pos, end_pos, start_adj, end_adj, length
#' @param gap threshold to split independent loci, with default value 10^4
#' @return a list containing 'genotype', 'MAF' and 'new_enloc'
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
    enhancer <- as.matrix(generate_enhancer(geno, gap))
  }

  return(list(genotype=geno, MAF=MAF, new_enloc=enhancer))

}
