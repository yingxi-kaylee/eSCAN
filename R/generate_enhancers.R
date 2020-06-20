#' @title generate enhancer
#' @description generate enhancer list for further use
#' @param geno genotype matrix
#' @param maxgap threshold to split independent loci, with default value 10^4
#' @return a data frame involving ennhancer information
#' @export
generate_enhancer <- function(geno, maxgap=10^4){
  numSNP <- dim(geno)[2]
  posiSNP <- as.numeric(colnames(geno))

  loci <- data.frame(CHR=rep(0,numSNP), POS=posiSNP, gap=rep(0,numSNP), numLoci=rep(1,numSNP))
  loci$gap[-1] <- loci$POS[-1]-loci$POS[-numSNP]
  loci$numLoci[-1] <- loci$gap[-1] > maxgap
  loci$numLoci <- cumsum(loci$numLoci)

  num_enhancer <- max(loci$numLoci)
  enhancer_list <- as.data.frame(matrix(0, nrow = num_enhancer, ncol = 6))
  for (i in 1:num_enhancer) {
    start <- min(loci[loci$numLoci==i,2])
    end <- max(loci[loci$numLoci==i,2])
    enhancer_list[i,2] <- start
    enhancer_list[i,3] <- end
    enhancer_list[i,4] <- which(posiSNP==start)[1]
    enhancer_list[i,5] <- which(posiSNP==end)[1]
  }
  enhancer_list[,6] <- enhancer_list[,5] - enhancer_list[,4] + 1
  colnames(enhancer_list) <- c("No.", "Start_pos", "End_pos", "Start_adj", "End_adj", "Length")

  short <- enhancer_list$Length > 40
  enhancer_list <- enhancer_list[short,]

  super <- enhancer_list$Length > 5000
  if(sum(super)>0){
    super_enhancer <- as.matrix(enhancer_list[super,])
    enhancer_list <- enhancer_list[!super,]
    #split super enhancers
    for(k in 1:sum(super)){
      num_split <- as.integer(super_enhancer[k,6]/5000) + 1
      add_pos <- 0:num_split
      new_start <- super_enhancer[k,4] + add_pos[-(num_split+1)]*5000
      new_end <- new_start + 5000 - 1
      new_end[num_split] <- min(new_end[num_split], super_enhancer[k,5])
      new_len <- new_end - new_start + 1

      new_mat <- cbind(0,posiSNP[new_start], posiSNP[new_end],new_start, new_end, new_len)
      super_enhancer <- rbind(super_enhancer,new_mat)
    }

    super_enhancer <- super_enhancer[-(1:sum(super)),]
    enhancer_list <- rbind(enhancer_list, super_enhancer)
    enhancer_list <- enhancer_list[order(enhancer_list[,4]),]
  }
  num_enhancer <- dim(enhancer_list)[1]
  enhancer_list[,1] <- 1:num_enhancer
  return(enhancer_list)
}


